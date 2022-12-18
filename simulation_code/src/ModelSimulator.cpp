#include <execution> // execution policy

#include "ModelSimulator.h"
#include "Cell.h"
#include "Utils.h"

using Point = Geometry::Point;
using Line = Geometry::Line;
using Polar = Material::Polar;

namespace {
    constexpr double SCALING_FACTOR{1e9}; // Factor to scale the scattering time (ns per second)
    constexpr std::size_t PHONON_CUTOFF{5'000'000}; // Switch from individual phonons to builders at this point
    constexpr std::size_t BUILDER_MAX_PHONONS{100'000};
    // Consider phonon velocity to be 0 in this direction if it is less than this value. Prevents some FP issues
    // Phonon would move at most 5 nm in this direction over a 50ns simulation, so this is safe to do in most cases
    // Should perhaps be scaled based on simulation settings (step_time_ and cell dimensions)
    constexpr double VELOCITY_EPS{0.01};
    // Prevent phonons from endlessly bouncing in tight corners. Consider scaling this based on step_time_?
    constexpr std::size_t MAX_COLLISIONS{100};
}

ModelSimulator::ModelSimulator(std::size_t measurement_steps, double simulation_time, bool phasor_sim)
    : step_time_{simulation_time / measurement_steps},
      phasor_sim_{phasor_sim}
{
    // Set up timing vector - each entry is the time at which a measurement will take place
    step_times_.resize(measurement_steps);
    std::generate(std::begin(step_times_), std::end(step_times_), [&, n=1] () mutable {
        return n++ * simulation_time / measurement_steps;
    });
}

void ModelSimulator::runSimulation(double t_eq) {
    (total_phonons_ < PHONON_CUTOFF) ? runPhononByPhonon(t_eq) : runUsingBuilders(t_eq);
}

void ModelSimulator::initPhononBuilders(std::vector<Cell>& cells, double t_eq, double eff_energy) noexcept {
    auto getPhonons = [&eff_energy](const double fractional_energy) {
        double temp_phonons = 0;
        const auto frac_phonons = std::modf(fractional_energy / eff_energy, &temp_phonons);
        auto num_phonons = static_cast<std::size_t>(temp_phonons); // rounds down
        return (Utils::urand() < frac_phonons) ? ++num_phonons : num_phonons; // 'round' up or stay rounded down
    };
    CellOriginBuilder cb{};
    for (auto& cell : cells) {
        const auto mat = cell.getMaterial();
        // Add initial phonon builders
        const auto init_energy = cell.getInitEnergy(t_eq);
        if (const auto init_phonons = getPhonons(init_energy); init_phonons > 0) {
            total_phonons_ += init_phonons;
            // Using max_phonons/2 since init_phonons generally have longer simulation times that emitted phonons
            if (const auto phonons = cb.totalPhonons(); phonons + init_phonons > BUILDER_MAX_PHONONS / 2 && phonons != 0) {
                phonon_builders_.emplace_back(std::move(cb));
                cb = CellOriginBuilder{};
            }
            cb.addCellPhonons(&cell, init_phonons);
        }
        // Add emit surface phonon builders
        for (const auto& boundary : cell.getBoundaries()) {
            for (const auto& es : boundary.getEmitSurfaces()) {
                const auto temp = es.getTemp();
                const auto energy_factor = mat.emitEnergy(temp) * es.getEmitDuration() * es.getLength() / 4.;
                const auto emit_energy = (t_eq == 0.) ? energy_factor : energy_factor * std::fabs(t_eq - temp);
                auto emit_phonons = getPhonons(emit_energy);
                total_phonons_ += emit_phonons;
                while (emit_phonons > BUILDER_MAX_PHONONS) {
                    (phasor_sim_) ? phonon_builders_.emplace_back(PhasorBuilder{cell, es, BUILDER_MAX_PHONONS})
                                  : phonon_builders_.emplace_back(SurfaceOriginBuilder{cell, es, BUILDER_MAX_PHONONS});
                    emit_phonons -= BUILDER_MAX_PHONONS;
                }
                (phasor_sim_) ? phonon_builders_.emplace_back(PhasorBuilder{cell, es, emit_phonons})
                              : phonon_builders_.emplace_back(SurfaceOriginBuilder{cell, es, emit_phonons});
            }
        }
    } if (cb.hasPhonons()) {
        phonon_builders_.emplace_back(std::move(cb));
    }
}

std::optional<double> ModelSimulator::nextImpact(Phonon& p, double time) const noexcept {
    // Get some necessary information
    const auto& [px, py] = p.getPosition();
    const auto& [vx, vy] = p.getVelVector();
    const Point start_point{px, py};
    const Point end_point{px+time*vx, py+time*vy};
    const auto boundaryLines = p.getCellBoundaryLines();
    if (start_point == end_point) {
        return std::nullopt;
    }
    const Line phonon_path{start_point, end_point};

    auto getTime = [](double start_coord, double end_coord, double velocity, double max_time) {
        return (velocity > VELOCITY_EPS || velocity < -VELOCITY_EPS) ? (end_coord - start_coord) / velocity : max_time;
    };

    // Find the nearest impact time and corresponding impact point
    std::optional<Point> impact_point = std::nullopt;
    for (const auto& line : boundaryLines) {
        // If there is a point of intersection that is not the start point
        if (const auto poi = line.getIntersection(phonon_path); poi && (*poi != start_point)) {
            // If the time taken to hit this POI is <= previous shortest time -> store POI and time
            const auto impact_time_x = getTime(start_point.x, (*poi).x, vx, time);
            const auto impact_time_y = getTime(start_point.y, (*poi).y, vy, time);
            const auto impact_time = (impact_time_x <= impact_time_y) ? impact_time_x : impact_time_y;
            if (impact_time <= time) {
                time = impact_time;
                impact_point = poi;
            }
        }
    }
    if (impact_point) {
        p.setPosition((*impact_point).x, (*impact_point).y);
        p.handleSurfaceCollision(*impact_point, step_time_);
        return std::make_optional(time);
    }
    return std::nullopt;
}

void ModelSimulator::scatter(Phonon& p, const Phonon::RelaxRates& relax_rates) noexcept {
    const auto [tau_N_inv, tau_U_inv, tau_I_inv] = relax_rates;
    const double tau_inv = std::accumulate(std::cbegin(relax_rates), std::cend(relax_rates), 0.);
    const double rand = Utils::urand();
    if (rand <= (tau_N_inv + tau_U_inv) / tau_inv) { // Not an impurity scatter
        // Resample the new phonon (freq, vel & polarization)
        p.scatterUpdate();
        if (rand > tau_N_inv / tau_inv) { // Umklapp scatter -> change direction vector
            p.setRandDirection(Utils::urand(), Utils::urand());
        }
    } else if (tau_I_inv > 0.) { // Impurity scatter
        p.setRandDirection(Utils::urand(), Utils::urand());
    }
}

void ModelSimulator::simulatePhonon(Phonon&& p, std::size_t measurement_steps) const {
    bool phonon_alive = true;
    double phonon_age = p.getLifetime();
    auto step = static_cast<std::size_t>(phonon_age / step_time_);
    p.setLifeStep(step); // For transient simulations
    Phonon::RelaxRates relax_rates{};
    double time_to_scatter = 0.;
    double time_to_measurement = 0.;

    auto get_scatter_info = [&step](const Phonon& p) {
        const auto relax_rates = p.getRelaxRates(step);
        return std::make_pair(relax_rates,
                              SCALING_FACTOR * -log(Utils::urand()) / std::accumulate(std::cbegin(relax_rates),
                                                                                      std::cend(relax_rates), 0.));
    };

    while (phonon_alive) {
        // If the phonon has scattered on the previous iteration recalculate new scattering rates and
        // find the time to the next scattering event
        if (time_to_scatter <= 0.) {
            std::tie(relax_rates, time_to_scatter) = get_scatter_info(p);
        }
        // If a measurement event occurred -> find the time to the next measurement event
        if (time_to_measurement <= 0.) {
            time_to_measurement = step_times_[step] - phonon_age;
        }
        auto drift_time = std::min(time_to_scatter, time_to_measurement); // Drift time until next non-impact event
        const auto sensor_id = p.getCellSensorID();
        // drifted_time is how long the phonon drifts before an impact event
        // Will be equal to drift_time if there is no impact
        // Will be false/null if the phonon impacts an emitting surface - signals it should be removed from system
        const std::optional<double> drifted_time = handleImpacts(p, drift_time, sensor_id);

        if (drifted_time) { // If the phonon had a transition/boundary surface collision
            // If the phonon has transitioned to a new sensor area (impact with transition surface)
            // Adjust drift_time to reflect there may be additional impacts but, first we need to find
            // a new scattering time before continuing
            if (p.getCellSensorID() != sensor_id) {
                // reduce drift_time to the amount of time the phonon has drifted, so we can start
                // the process over with a fresh scattering time as we have entered a different sensor area
                // i.e. old time_to_scatter is no longer valid
                drift_time = *drifted_time;
            }
            p.drift(drift_time-*drifted_time);
            phonon_age += drift_time;
            time_to_measurement -= drift_time;
            time_to_scatter -= drift_time;
            if (time_to_measurement == 0.) { // Take a measurement
                if (++step < measurement_steps) { // Simulation time has not been exceeded
                    p.setLifeStep(step);
                    if (step >= step_adjustment_) {
                        p.updateCellHeatParams(step - step_adjustment_);
                    }
                } else { // Exceeds simulation time
                    phonon_alive = false;
                }
            } else if (!phasor_sim_ && time_to_scatter == 0.) { //
                scatter(p, relax_rates);
            } else { // This is the condition when the phonon transitions to a new sensor area
                time_to_scatter = 0.; // Reset scattering time based on new sensor area properties
            }
        } else { // Phonon made impact with an emitting surface (left system)
            phonon_alive = false;
        }
    }
}

// returning 0 means the calling function will drift the phonon for drift_time (it does not drift here)
// returning a number will reduce the amount of time the calling function drifts the phonon by that amount
// The next impact method places the phonon on the poi of the impacted surface effectively drifting it for impact_time
// The higher drifted time is here, the less amount of time the calling function drifts the phonon
// returning std::nullopt will kill the phonon
std::optional<double> ModelSimulator::handleImpacts(Phonon& p, double drift_time, std::size_t sensor_id) const {
    auto impact_time = nextImpact(p, drift_time);
    double drifted_time = 0.;
    std::size_t collision_counter = 0;
    // Impact with an emitting surface will set the phonon cell to nullptr
    while (impact_time) {
        if (p.outsideCell()) { return std::nullopt; }
        drifted_time += *impact_time;
        // If phonon is stuck, move it to a random location in the cell - primarily used to handle FP issues
        if (++collision_counter > MAX_COLLISIONS) {
            p.setRandPoint(Utils::urand(), Utils::urand());
            return std::make_optional<double>(drift_time); // calling function will not further drift the phonon
        }
        // If the phonon has changed sensor areas - return immediately as scatter time must be reset
        const auto cur_sensor_id = p.getCellSensorID();
        if (sensor_id != cur_sensor_id ) {
            return std::make_optional<double>(drifted_time);
        }
        impact_time = nextImpact(p, drift_time - drifted_time);
    }
    return (p.outsideCell()) ? std::nullopt : std::make_optional<double>(drifted_time);
}

void ModelSimulator::runPhononByPhonon(double t_eq) {
    std::vector<std::unique_ptr<Phonon>> phonons;
    phonons.reserve(total_phonons_);
    for (BuilderObj& builderObj : phonon_builders_) {
        std::visit([&] (auto& builder) {
            while (builder.hasPhonons()) {
                phonons.emplace_back(std::make_unique<Phonon>(builder(t_eq)));
            }
        }, builderObj);
    }
    // Shuffling the phonons may help with the performance bottleneck when recording a phonon's contribution
    // to a sensor's energy/flux at each measurement step (Sensor::updateHeatParams)
    std::shuffle(std::begin(phonons), std::end(phonons), std::random_device());
    std::for_each(std::execution::par, std::begin(phonons), std::end(phonons), [&](auto& phonon) {
        simulatePhonon(std::move(*phonon), step_times_.size());
    });
}

void ModelSimulator::runUsingBuilders(double t_eq) {
    std::for_each(std::execution::par, std::begin(phonon_builders_), std::end(phonon_builders_), [&](BuilderObj& builderObj) {
        std::visit([&] (auto& builder) {
            while (builder.hasPhonons()) {
                simulatePhonon(builder(t_eq), step_times_.size());
            }
        }, builderObj);
    });
}


