#include <numeric>
#include <algorithm>
#include <cmath>
#include <execution>

#include "SensorInterpreter.h"

namespace {
    // For the numerical inversions
    constexpr double EPS{.0001};
    constexpr std::size_t MAX_ITERS{40};
    // If sensor temperature is greater than .5% -> needs to be reset - check flux as well?
    constexpr double RESET_THRESHOLD{0.005};
}

void SensorInterpreter::setParams(double t_eq, double eff_energy) noexcept {
    t_eq_ = t_eq;
    eff_energy_ = eff_energy;
}

SensorMeasurements SensorInterpreter::scaleHeatParams(const Sensor& sensor) const noexcept {
    SensorMeasurements sm;
    const auto& flux = sensor.getFluxes();
    const auto num_measurements = flux.size();
    sm.id = sensor.getID();
    sm.final_temps = findTemperature(sensor);
    sm.final_temps.front() = sensor.getInitTemp();
    sm.final_fluxes.resize(num_measurements);
    const auto f_factor = eff_energy_ / sensor.getArea();

    std::transform(std::cbegin(flux), std::cend(flux), std::begin(sm.final_fluxes), [&f_factor](auto elem) {
        const auto [vx, vy] = elem;
        return std::array{vx * f_factor, vy * f_factor};
    });

    // Store steady state results in Sensor object - may want to consider moving this back to the output manager class
    auto avgAndStdError = [](const std::vector<double>& data) {
        const std::size_t size = data.size();
        const auto avg = std::reduce(std::execution::seq, std::cbegin(data), std::cend(data), 0.) / size;
        const auto std_dev = std::sqrt(std::transform_reduce(std::execution::seq, std::cbegin(data), std::cend(data), 0., std::plus{},
                                                         [&avg](const auto& val) { return (avg - val) * (avg - val); }) / size);
        return std::pair<double, double>(avg, std_dev/std::sqrt(size));
    };

    std::vector<double> t_data; t_data.reserve(num_measurements);
    std::vector<double> fx_data; fx_data.reserve(num_measurements);
    std::vector<double> fy_data; fy_data.reserve(num_measurements);

    std::transform(std::cbegin(sm.final_temps), std::cend(sm.final_temps), std::back_inserter(t_data),
                   [](const auto temp) { return temp; });
    std::transform(std::cbegin(sm.final_fluxes), std::cend(sm.final_fluxes), std::back_inserter(fx_data),
                   [](const auto& flux) { return flux[0]; });
    std::transform(std::cbegin(sm.final_fluxes), std::cend(sm.final_fluxes), std::back_inserter(fy_data),
                   [](const auto& flux) { return flux[1]; });

    std::tie(sm.t_steady, sm.std_t_steady) = avgAndStdError(t_data);
    std::tie(sm.x_flux, sm.std_x_flux) = avgAndStdError(fx_data);
    std::tie(sm.y_flux, sm.std_y_flux) = avgAndStdError(fy_data);
    return sm;
}

double SensorInterpreter::getFinalTemp(const Sensor& sensor, std::size_t start_step) const noexcept {
    // TODO: it is probably a user input error if no cells are linked to a given sensor
    if (sensor.getArea() == 0.) { return 0.; }
    const auto steps = sensor.getEnergies().size();
    const auto temps = findTemperature(sensor, start_step);
    return std::accumulate(std::cbegin(temps), std::cend(temps), 0.) / (steps - start_step);
}

std::vector<double> SensorInterpreter::getFinalTemps(const Sensor& sensor) const noexcept {
    return findTemperature(sensor); // Transient simulations need the temperature at each measurement step
}

std::vector<double> SensorInterpreter::findTemperature(const Sensor& sensor, std::size_t start_step) const noexcept {
    const auto& energies = sensor.getEnergies();
    auto start = std::cbegin(energies)+start_step; auto end = std::cend(energies);
    std::vector<double> temps(std::distance(start, end));
    const auto& material = sensor.getMaterial();
    const auto area = sensor.getArea();

    auto inversion = [&](double current_energy, bool pseudo=false) {
        double temp = 0., de, ub = ub_, lb = lb_;
        std::size_t iter = 0;
        while ( (ub - lb >= EPS) && (++iter != MAX_ITERS) ) {
            temp = (ub + lb) / 2.;
            de = (material.theoreticalEnergy(temp, pseudo) * area) - current_energy;
            (de < 0.) ? lb = temp : ub = temp;
        }
        return temp;
    };

    std::transform(std::execution::seq, start, end, std::begin(temps), [&, index=0](const auto& energy_units) mutable {
        // Multiply the number of energy units by the phonon effective energy to get the total energy at each measurement step
        const double energy = eff_energy_ * energy_units;
        if (t_eq_ != 0.) { // Do approximation to find the temperature
            // If it is a steady-state simulation, the index will be disregarded when finding the heat capacity
            return energy / (area * sensor.getHeatCapacity(index++)) + t_eq_;
        }
        else { // Do numerical inversion
            return inversion(energy);
        }
    });

    return temps;
}


