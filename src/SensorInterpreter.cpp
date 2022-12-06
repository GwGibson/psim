#include <numeric>
#include <algorithm>
#include <cmath>
#include <execution>

#include "SensorInterpreter.h"

namespace {
    // For the numerical inversions
    static constexpr double EPS{.0001};
    static constexpr std::size_t MAX_ITERS{40};
    // If sensor temperature is greater than .5% -> needs to be reset - check flux as well?
    static constexpr double RESET_THRESHOLD{0.005};
    // Percentage of measurement steps that will be used for steady state calculations
    static constexpr double SS_STEPS_PERCENT{0.1};
}

SensorMeasurements SensorInterpreter::scaleHeatParams(const Sensor& sensor) const noexcept {
    SensorMeasurements sm;
    // Scale the results
    const auto& energies = sensor.getEnergies();
    const auto& flux = sensor.getFluxes();
    const auto num_measurements = energies.size();
    const auto ss_steps = static_cast<int>(num_measurements * SS_STEPS_PERCENT);
    const auto start = static_cast<int>(num_measurements - ss_steps);
    sm.id = sensor.getID();
    sm.final_temps = findTemperature(sensor, std::cbegin(energies), std::cend(energies));
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

    std::vector<double> t_data; t_data.reserve(ss_steps);
    std::vector<double> fx_data; fx_data.reserve(ss_steps);
    std::vector<double> fy_data; fy_data.reserve(ss_steps);

    std::transform(std::cbegin(sm.final_temps) + start, std::cend(sm.final_temps), std::back_inserter(t_data),
                   [](const auto temp) { return temp; });
    std::transform(std::cbegin(sm.final_fluxes) + start, std::cend(sm.final_fluxes), std::back_inserter(fx_data),
                   [](const auto& flux) { return flux[0]; });
    std::transform(std::cbegin(sm.final_fluxes) + start, std::cend(sm.final_fluxes), std::back_inserter(fy_data),
                   [](const auto& flux) { return flux[1]; });

    std::tie(sm.t_steady, sm.std_t_steady) = avgAndStdError(t_data);
    std::tie(sm.x_flux, sm.std_x_flux) = avgAndStdError(fx_data);
    std::tie(sm.y_flux, sm.std_y_flux) = avgAndStdError(fy_data);
    return sm;
}

double SensorInterpreter::getFinalTemp(const Sensor& sensor) const noexcept {
    // TODO: it is probably a user input error if no cells are linked to a given sensor
    if (sensor.getArea() == 0.) { return 0.; }
    const auto& inc_energy = sensor.getEnergies();
    const auto total_steps = inc_energy.size();
    const auto steps = static_cast<std::size_t>(total_steps*SS_STEPS_PERCENT + 1); // Use 10% of steps to determine steady state temp
    const auto start = static_cast<int>(total_steps - steps);
    const auto temps = findTemperature(sensor, std::cbegin(inc_energy)+start, std::cend(inc_energy));

    return std::accumulate(std::cbegin(temps), std::cend(temps), 0.) / steps;
}

std::vector<double> SensorInterpreter::getFinalTemps(const Sensor& sensor) const noexcept {
    const auto& inc_energy = sensor.getEnergies();
    return findTemperature(sensor, std::cbegin(inc_energy), std::cend(inc_energy));
}

std::vector<double> SensorInterpreter::findTemperature(const Sensor& sensor, std::vector<int>::const_iterator start,
                                                      std::vector<int>::const_iterator end) const noexcept {
    // Temperatures stored here - this vector is the return vector - gives the local temperature at each requested measurement step
    std::vector<double> temps(std::distance(start, end));
    if (t_eq_ != 0.) { // Do approximation
        std::transform(std::execution::par, start, end, std::begin(temps), [&, index=0](const auto& energy_units) mutable {
            // Multiply the number of energy units by the phonon effective energy to get the total energy at each measurement step
            const double energy = eff_energy_ * energy_units;
            // If it is a steady-state simulation, the index will be disregarded when finding the heat capacity
            return energy / (sensor.getArea() * sensor.getHeatCapacity(index++)) + t_eq_;
        });
    } else { // Do numerical inversion
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

        std::transform(std::execution::par, start, end, std::begin(temps), [&](const auto& energy_units) {
            // Get the local temperature by inverting the sum of the energy_units in each energy array
            return inversion(eff_energy_ * energy_units);
        });
    }
    return temps;
}

