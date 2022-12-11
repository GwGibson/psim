#ifndef GEOMETRY_SENSORINTERPRETER_H
#define GEOMETRY_SENSORINTERPRETER_H

#include "Sensor.h"

class SensorInterpreter {
public:
    SensorInterpreter() = default;

    void setBounds(double ub, double lb) noexcept { ub_ = ub; lb_ = lb; }
    void setParams(double t_eq, double eff_energy) { t_eq_ = t_eq; eff_energy_ = eff_energy; }
    [[nodiscard]] SensorMeasurements scaleHeatParams(const Sensor& sensor) const noexcept;
    [[nodiscard]] double getFinalTemp(const Sensor& sensor) const noexcept;
    [[nodiscard]] std::vector<double> getFinalTemps(const Sensor& sensor) const noexcept; // For transient simulations

private:
    double ub_{1000.};
    double lb_{0.};
    double t_eq_{0.};
    double eff_energy_{0.};

    [[nodiscard]] std::vector<double> findTemperature(const Sensor& sensor, std::vector<int>::const_iterator start,
                                                      std::vector<int>::const_iterator end) const noexcept;
};


#endif //GEOMETRY_SENSORINTERPRETER_H
