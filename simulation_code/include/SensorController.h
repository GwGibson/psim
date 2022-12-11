#ifndef GEOMETRY_SENSORCONTROLLER_H
#define GEOMETRY_SENSORCONTROLLER_H

#include "Material.h"

class SensorController {
public:
    SensorController(const Material& material, double t_init, std::size_t num_measurements=0);
    virtual ~SensorController() = default;

    [[nodiscard]] const Material& getMaterial() const noexcept { return material_; }
    [[nodiscard]] double getHeatCapacityAtFreq(std::size_t freq_index) const noexcept;
    [[nodiscard]] virtual double getHeatCapacity(std::size_t step) const noexcept = 0;
    [[nodiscard]] virtual double getInitTemp() const noexcept = 0;
    [[nodiscard]] virtual double getSteadyTemp(std::size_t step) const noexcept = 0;

    void initialUpdate(Phonon& p, const Material::Table* table) const noexcept;
    virtual void scatterUpdate(Phonon& p) const noexcept;
    /**
     * Returns true if the input sensor's temperature has not significantly changed over the course of the simulation
     * @param t_final - The temperature at the end of the run
     * @param final_temps - A vector of temperature from the previous measurement steps -> for transient simulations
     * @return - true if the temp of this sensor is unstable (final temp not within some percentage of initial temp)
     */
    [[nodiscard]] virtual bool resetRequired(double t_final, std::vector<double>&&) noexcept;
    virtual void reset() noexcept = 0;
protected:
    const Material& material_;
    const double t_init_;
    const std::size_t num_measurements_; // For transient surfaces only

    double t_steady_{0.}; // Steady state temperature of the cell. Used to set the energy tables & heat_capacity_
    double heat_capacity_{0.};
    Material::Table base_table_;
    Material::Table scatter_table_;

    // Transient sensor containers -> not needed for steady state or periodic simulations
    std::vector<Material::Table> scatter_tables_; // Transient controllers need a scatter table for each measurement step
    std::vector<double> heat_capacities_;
    std::vector<double> steady_temps_;
};

class SteadyStateController : public SensorController {
public:
    using SensorController::SensorController;

    [[nodiscard]] double getHeatCapacity(std::size_t) const noexcept override { return heat_capacity_; }
    [[nodiscard]] double getInitTemp() const noexcept override { return t_steady_; }
    [[nodiscard]] double getSteadyTemp(std::size_t) const noexcept override { return t_steady_; }

    void reset() noexcept override;
};

class PeriodicController : public SensorController {
public:
    using SensorController::SensorController;

    [[nodiscard]] double getHeatCapacity(std::size_t) const noexcept override { return heat_capacity_; }
    [[nodiscard]] double getInitTemp() const noexcept override { return t_init_; }
    [[nodiscard]] double getSteadyTemp(std::size_t) const noexcept override { return t_steady_; }

    void reset() noexcept override;
};

// May want to average results of final 10% of runs or something like this
class TransientController : public SensorController {
public:
    using SensorController::SensorController;

    [[nodiscard]] double getHeatCapacity(std::size_t step) const noexcept override { return heat_capacities_[step]; }
    [[nodiscard]] double getInitTemp() const noexcept override { return t_init_; }
    // 0 if no step and steady_temps_[step] if step specified (t_eq update pointless)
    [[nodiscard]] double getSteadyTemp(std::size_t step) const noexcept override;

    void scatterUpdate(Phonon& p) const noexcept override;
    [[nodiscard]] bool resetRequired(double, std::vector<double>&& final_temps) noexcept override;
    void reset() noexcept override;
private:
};

#endif //GEOMETRY_SENSORCONTROLLER_H
