#include <algorithm>
#include <cmath>
#include <execution>

#include "Model.h"

#include <iostream> // Testing

namespace {
    // Maximum number of simulations resets
    static constexpr std::size_t MAX_ITERS{10};
    // Threshold (percentage*100) of sensors that must be stable for the system to be considered stable
    static constexpr std::size_t RESET_THRESHOLD{90};
    // Threshold (percentage*1000) that t_eq must be within between runs for the system to be considered stable
    // 5 means t_eq must be within 0.5% of previous t_eq
    static constexpr std::size_t TEQ_THRESHOLD{5};
}

using Point = Geometry::Point;
using Line = Geometry::Line;

Model::Model(std::size_t num_cells, std::size_t num_sensors, std::size_t measurement_steps, std::size_t num_phonons,
             double simulation_time, double t_eq)
        : num_cells_{num_cells},
          measurement_steps_{measurement_steps},
          simulation_time_{simulation_time},
          num_phonons_{num_phonons},
          t_eq_{t_eq},
          simulator_{measurement_steps, simulation_time},
          outputManager_{},
          interpreter_{},
          addMeasurementMutex_{std::make_unique<std::mutex>()}
{
    cells_.reserve(num_cells);
    sensors_.reserve(num_sensors);
}

void Model::setSimulationType(SimulationType type, std::size_t step_interval) {
    sim_type_ = type;
    // Throwing here so the user doesn't run a full simulation only to realize after that they
    // are using the incorrect settings.
    if ( (type == SimulationType::Transient || type == SimulationType::Periodic) && step_interval == 0) {
        throw std::runtime_error(std::string("Step interval of 0 is invalid for transient and periodic simulations.\n"));
    }
    if (type == SimulationType::SteadyState && step_interval > 0) {
        throw std::runtime_error(std::string("Step interval > 0 is invalid for steady-state simulations\n"));
    }
    outputManager_.setStepInterval(step_interval);
}

void Model::addMaterial(const std::string& material_name, Material material) {
    const auto exists = materials_.find(material_name);
    if (exists != std::end(materials_)) {
        throw std::runtime_error(std::string("A duplicate material name was detected.\n"));
    }
    materials_.emplace(material_name, material);
}

// Assumes the sensor material exists in the materials_ map
void Model::addSensor(std::size_t ID, const std::string& material_name, double t_init, Sensor::SimulationType type) {
    auto sensor = std::find_if(std::begin(sensors_), std::end(sensors_), [&ID](auto& s){
        return s.getID() == ID;
    });
    if (sensor == std::end(sensors_)) {
        sensors_.emplace_back(Sensor{ID, materials_.at(material_name), measurement_steps_, t_init, type});
    } else {
        throw std::runtime_error(std::string("Sensor with this ID already exists\n"));
    }
}

void Model::addCell(Geometry::Triangle&& triangle, std::size_t sensor_ID, double spec) {
    if (cells_.size() >= num_cells_) {
        throw std::runtime_error(std::string("Too many cells\n"));
    }
    cells_.emplace_back(Cell{std::move(triangle), getSensor(sensor_ID), spec});
    Cell& inc_cell = cells_.back();
    std::size_t identical_cells = 0; // Should only be 1 identical cell (check against itself)
    for (auto& cell : cells_) {
        if (inc_cell != cell) {
            // Throws if the incoming cell is incompatible with any existing cells
            // TODO: Misses cases when all 3 points of a cell are on the edges of another cell
            inc_cell.validate(cell);
            inc_cell.findTransitionSurface(cell);
        } else if (++identical_cells == 2) {
            throw std::runtime_error(std::string("Duplicate cell detected.\n"));
        }
    }
}

// Not used if the python interface is used
void Model::addCell(Point&& p1, Point&& p2, std::size_t sensor_ID, double spec) {
    if (p1.x == p2.x || p1.y == p2.y) {
        throw std::runtime_error(std::string("These points do not specify a rectangle\n"));
    }
    addCell(Geometry::Triangle{p1, Point(p1.x, p2.y), Point(p2.x, p1.y)}, sensor_ID, spec);
    addCell(Geometry::Triangle{p2, Point(p2.x, p1.y), Point(p1.x, p2.y)}, sensor_ID, spec);
}

// TODO: Fix emit surface placement failing when incoming surface is not exact in rare circumstances
bool Model::setEmitSurface(const Point& p1, const Point& p2, double temp, double duration, double start_time) {

    if ( (start_time < 0. || start_time >= simulation_time_) || (duration < 0. || duration > simulation_time_-start_time) ) {
        throw std::runtime_error(std::string("Transient Surface start_time or duration specifications are invalid.\n"));
    }
    if ( (start_time > 0. || duration < simulation_time_) && sim_type_ != SimulationType::Transient) {
        throw std::runtime_error(std::string("Cannot add a transient surface to a non transient simulation.\n"));
    }
    for (auto& cell : cells_) {
        if (cell.setEmitSurface(Line{p1, p2}, temp, duration, start_time)) {
            return true;
        }
    }
    return false;
}

void Model::runSimulation() {
    // TODO: could run checks here to verify there is at least 1 sensor/cell etc.
    setTemperatureBounds();

    auto refresh = [this]() {
        auto total_energy = getTotalInitialEnergy();
        double energy_per_phonon = total_energy / num_phonons_;
        interpreter_.setParams(t_eq_, energy_per_phonon);
        return std::pair{total_energy, energy_per_phonon};
    };

    auto [total_energy, energy_per_phonon] = refresh();

    std::size_t iter = 0;
    bool reset_required = true;
    while (reset_required && ++iter <= MAX_ITERS) {
        simulator_.initPhononBuilders(cells_, t_eq_, energy_per_phonon);
        simulator_.runSimulation(t_eq_);

        // Check if sensor temperatures are stable
        if (const auto new_t_eq = resetRequired(); new_t_eq && iter < MAX_ITERS) {
            reset();
            t_eq_ = *new_t_eq;
            std::cout << "system not stable\n";
            std::cout << "updated t_eq: " << t_eq_ << '\n';
        } else {
            reset_required = false;
        }

        std::tie(total_energy, energy_per_phonon) = refresh();

    }
    if (iter >= MAX_ITERS) { std::cout << "System did not stabilize!!\n"; } // Should log this or include in output title
    storeResults();
}

void Model::exportResults(const fs::path& filepath) const {
    outputManager_.steadyStateExport(filepath);
    // If transient or periodic simulation -> do periodic export (This is useless for ss simulation that require resets)
    if (sim_type_ != SimulationType::SteadyState) {
        outputManager_.periodicExport(filepath);
    }
}

Sensor& Model::getSensor(std::size_t ID) {
    auto sensor = std::find_if(std::begin(sensors_), std::end(sensors_), [&ID](auto& s){
        return s.getID() == ID;
    });
    if (sensor == std::end(sensors_)) {
        throw std::runtime_error(std::string("Sensor does not exist\n"));
    }
    return *sensor;
}

double Model::getTotalInitialEnergy() const noexcept {
    return std::transform_reduce(std::execution::seq, std::cbegin(cells_), std::cend(cells_), 0., std::plus{},
                               [&](const auto& cell) { return cell.getInitEnergy(t_eq_) + cell.getEmitEnergy(t_eq_); });
}

void Model::setTemperatureBounds() noexcept {
    constexpr double temp_eps{10.};
    std::vector<double> temperatures;
    for (const auto& cell : cells_) {
        temperatures.push_back(cell.getInitTemp());
        for (const auto& surface : cell.getBoundaries()) {
            for (const auto& es : surface.getEmitSurfaces()) {
                temperatures.push_back(es.getTemp());
            }
        }
    }
    const auto [min, max] = std::minmax_element(std::cbegin(temperatures), std::cend(temperatures));
    interpreter_.setBounds(*max+temp_eps, std::max(*min-temp_eps, 0.));
}

double Model::avgTemp() const {
    const double total_area = std::transform_reduce(std::execution::seq, std::cbegin(sensors_), std::cend(sensors_), 0.,
                                                    std::plus{}, [](const auto& sensor) { return sensor.getArea(); });
    return std::transform_reduce(std::execution::seq, std::cbegin(sensors_), std::cend(sensors_), 0., std::plus{},
                            [&](const auto& sensor) { return sensor.getSteadyTemp() * sensor.getArea() / total_area; });
}

void Model::storeResults() noexcept {
    std::for_each(std::execution::par, std::cbegin(sensors_), std::cend(sensors_), [&](const auto& sensor) {
        auto measurement = interpreter_.scaleHeatParams(sensor);
        std::scoped_lock lg(*addMeasurementMutex_);
        outputManager_.addMeasurement(std::move(measurement));
    });
    outputManager_.sortMeasurements();
}

std::optional<double> Model::resetRequired() noexcept {
    auto t_diff = [](const auto& t_final, const auto& t_init) {
        return std::fabs(t_final-t_init)/t_init * 1000 > TEQ_THRESHOLD;
    };

    const std::size_t total_sensors = sensors_.size();
    std::size_t stable_sensors = 0;
    for (auto& sensor : sensors_) {
        if (sim_type_ != SimulationType::Transient) {
            // The average temperature of the last 10% of measurement steps
            if (sensor.resetRequired(interpreter_.getFinalTemp(sensor))) { ++stable_sensors; }
        } else {
            // Get the temperature at each measurement step
            if (sensor.resetRequired(0., interpreter_.getFinalTemps(sensor))) { ++stable_sensors; }
        }
    }
    std::cout << "Stable sensors: " << stable_sensors << '\n';
    const auto new_t_eq = (t_eq_ == 0. || sim_type_ == SimulationType::Transient) ? t_eq_ : avgTemp();
    // if either the sensors or t_eq are not stable, return true indicated a reset is required
    return (stable_sensors * 100 / total_sensors < RESET_THRESHOLD) || t_diff(new_t_eq, t_eq_) ? std::make_optional(new_t_eq)
                                                                                               : std::nullopt;
}

void Model::reset() noexcept {
    simulator_.reset();
    for (auto& sensor : sensors_) { sensor.reset(); }
}