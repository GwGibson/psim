#include <execution> // execution policy

#include "ModelSimulator.h"
#include "Geometry.h"
#include "Cell.h"
#include "Utils.h"

using Point = Geometry::Point;
using Line = Geometry::Line;
using Polar = Material::Polar;

namespace {
    static constexpr double SCALING_FACTOR{1e9}; // Factor to scale the scattering time -> related to GEOEPS in geometry
    static constexpr std::size_t PHONON_CUTOFF{5'000'000}; // Switch from individual phonons to builders at this point
    static constexpr std::size_t BUILDER_MAX_PHONONS{100'000};
    static constexpr std::size_t MAX_COLLISIONS{100}; // Prevent phonons from endlessly bouncing in tight corners
}

ModelSimulator::ModelSimulator(std::size_t measurement_steps, double simulation_time)
    : step_time_{simulation_time / measurement_steps}
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
            // Using max_phonons/2 since init_phonons generally take longer to simulate compared to emitted phonons
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
                    phonon_builders_.emplace_back(SurfaceOriginBuilder{cell, es, BUILDER_MAX_PHONONS});
                    emit_phonons -= BUILDER_MAX_PHONONS;
                }
                phonon_builders_.emplace_back(SurfaceOriginBuilder{cell, es, emit_phonons});
            }
        }
    } if (cb.hasPhonons()) {
        phonon_builders_.emplace_back(std::move(cb));
    }
}

std::optional<double> ModelSimulator::nextImpact(Phonon& p, double time) const noexcept {
    // Get some necessary information
    const auto cell = p.getCell();
    const auto& [px, py] = p.getPosition();
    const auto& [vx, vy] = p.getVelVector();
    const Point start_point{px, py};
    const Point end_point{px+time*vx, py+time*vy};
    const auto boundaryLines = cell->getBoundaryLines();
    if (start_point == end_point) {
        return std::nullopt;
    }
    const Line phonon_path{start_point, end_point};

    // Find the nearest impact time and corresponding impact point
    double nearest_impact_time = time;
    std::optional<Point> impact_point = std::nullopt;
    for (const auto& line : boundaryLines) {
        // If there is a point of intersection that is not the start point
        if (const auto poi = line.getIntersection(phonon_path); poi && !(poi == start_point)) {
            // If the time taken to hit this POI is <= previous shortest time -> store POI and time
            const auto impact_time = (vx == 0.) ? ((*poi).y - start_point.y) / vy : ((*poi).x - start_point.x) / vx;
            if (impact_time <= nearest_impact_time) {
                nearest_impact_time = impact_time;
                impact_point = poi;
            }
        }
    }
    if (impact_point) {
        p.setPosition((*impact_point).x, (*impact_point).y);
        cell->handleSurfaceCollision(p, *impact_point, step_time_);
        return std::make_optional(nearest_impact_time);
    }
    return std::nullopt;
}

void ModelSimulator::scatter(Phonon& p, const std::array<double, 3>& relax_rates) noexcept {
    const auto [tau_N_inv, tau_U_inv, tau_I_inv] = relax_rates;
    const double tau_inv = std::accumulate(std::cbegin(relax_rates), std::cend(relax_rates), 0.);
    const double rand = Utils::urand();
    if (rand <= (tau_N_inv + tau_U_inv) / tau_inv) { // Not an impurity scatter
        // Resample the new phonon (freq, vel & polarization)
        p.getCell()->scatterUpdate(p);
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
    std::array<double, 3> relax_rates{};
    double time_to_scatter = 0.;
    double time_to_measurement = 0.;

    auto get_scatter_info = [&step](const Phonon& p) {
        const auto& cell = p.getCell();
        auto scatter_rates = cell->getMaterial().scatteringRates(cell->getSteadyTemp(step), p.getFreq(), p.getPolar());
        return std::make_pair(scatter_rates,
                              SCALING_FACTOR * -log(Utils::urand()) / std::accumulate(std::cbegin(scatter_rates),
                                                                                      std::cend(scatter_rates), 0.));
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
        auto drift_time = std::min(time_to_scatter, time_to_measurement);
        const auto sensor_id = p.getCell()->getSensorID();
        const auto drifted_time = handleImpacts(p, drift_time, sensor_id);

        if (drifted_time) { // If the phonon had surface collisions
            if (p.getCell()->getSensorID() != sensor_id) { // If the phonon has transitioned to a new sensor area
                // Adjust drift_time to reflect there may be additional impacts but, first we need to find
                // a new scattering time before continuing
                drift_time = *drifted_time;
            }
            p.drift(drift_time-*drifted_time);
            phonon_age += drift_time;
            time_to_measurement -= drift_time;
            time_to_scatter -= drift_time;
            if (time_to_measurement == 0.) {
                if (++step < measurement_steps) { // Simulation time has not been exceeded
                    p.setLifeStep(step);
                    p.getCell()->updateHeatParams(p, step);
                } else { // Exceeds simulation time
                    phonon_alive = false;
                }
            } else if (time_to_scatter == 0.) {
                scatter(p, relax_rates);
            } else { // This is the condition when the phonon transitions to a new sensor area
                // TODO: If material is also different - reset phonon based on new material properties
                time_to_scatter = 0.; // Reset scattering time based on new sensor area properties
            }
        } else { // Phonon made impact with an emitting surface (left system)
            phonon_alive = false;
        }
    }
}

std::optional<double> ModelSimulator::handleImpacts(Phonon& p, double drift_time, std::size_t sensor_id) const noexcept {
    auto impact_time = nextImpact(p, drift_time);
    double drifted_time = 0.;
    std::size_t collision_counter = 0;
    while (impact_time) {
        // If phonon is stuck, move it to a random location in the cell
        if (++collision_counter > MAX_COLLISIONS) {
            const auto [x, y] = p.getCell()->getRandPoint(Utils::urand(), Utils::urand());
            p.setPosition(x, y);
        }
        // Impact with an emitting surface will set the phonon cell to nullptr
        if (const auto& cell = p.getCell(); cell != nullptr) {
            drifted_time += *impact_time;
            // If the phonon has changed sensor areas - return immediately as scatter time must be reset
            const auto cur_sensor_id = cell->getSensorID();
            if (sensor_id != cur_sensor_id ) {
                return std::make_optional<double>(drifted_time);
            }
            sensor_id = cur_sensor_id;
            impact_time = nextImpact(p, drift_time - drifted_time);
        } else {
            return std::nullopt;
        }
    }
    return std::make_optional<double>(drifted_time);
}

void ModelSimulator::runPhononByPhonon(double t_eq) {
    std::vector<std::unique_ptr<Phonon>> phonons;
    phonons.reserve(total_phonons_);
    for (BuilderObj& builderObj : phonon_builders_) {
        std::visit([&] (auto& builder) {
            while (auto phonon = builder(t_eq)) {
                phonons.emplace_back(std::make_unique<Phonon>(*phonon));
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
            while (auto phonon = builder(t_eq)) {
                simulatePhonon(std::move(*phonon), step_times_.size());
            }
        }, builderObj);
    });
}