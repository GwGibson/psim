#ifndef GEOMETRY_SENSORINTERPRETER_H
#define GEOMETRY_SENSORINTERPRETER_H

#include "Sensor.h"

class SensorInterpreter {
public:
    SensorInterpreter() = default;

    // The temperature bounds and t_eq / eff_energy parameters cannot be set in the constructor
    // as they are contingent on other aspects of the model and cannot be set until the model is ready
    // to start the simulation
    void setBounds(double lb, double ub) noexcept { lb_ = lb; ub_ = ub; }
    void setParams(double t_eq, double eff_energy) noexcept;
    [[nodiscard]] SensorMeasurements scaleHeatParams(const Sensor& sensor) const noexcept;
    [[nodiscard]] double getFinalTemp(const Sensor& sensor, std::size_t start_step) const noexcept;
    [[nodiscard]] std::vector<double> getFinalTemps(const Sensor& sensor) const noexcept; // For transient simulations

private:
    double ub_;         // Adjusted in the Model class based on the emitting max emitting surface temperature
    double lb_;         // Same here but min emitting surface temperature
    double t_eq_;       // Adjusted in the Model class
    double eff_energy_; // Adjusted in the Model class once the effective energy is calculated

    [[nodiscard]] std::vector<double> findTemperature(const Sensor& sensor, std::size_t start_step=0) const noexcept;
};


#endif //GEOMETRY_SENSORINTERPRETER_H
