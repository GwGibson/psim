#ifndef GEOMETRY_SENSOR_H
#define GEOMETRY_SENSOR_H

#include <mutex>
#include <memory>

#include "Material.h"
#include "SensorController.h"

struct SensorMeasurements;

class Sensor {
public:
    enum class SimulationType {SteadyState, Periodic, Transient};
    Sensor(std::size_t ID, const Material& material, SimulationType type, std::size_t num_measurements, double t_init);

    void initialUpdate(Phonon& p, const Material::Table& table) const noexcept;
    void initialUpdate(Phonon& p) const noexcept;
    void scatterUpdate(Phonon& p) const noexcept;
    void addToArea(double area) noexcept { area_covered_ += area; }
    // final temps vector is only needed for transient simulations -> included for all calls to keep a common interface
    [[nodiscard]] bool resetRequired(double t_final, std::vector<double>&& final_temps={}) noexcept;

    [[nodiscard]] std::size_t getID() const noexcept { return ID_; }
    [[nodiscard]] const Material& getMaterial() const noexcept { return controller_->getMaterial(); }
    [[nodiscard]] double getHeatCapacity(std::size_t step=0) const noexcept { return controller_->getHeatCapacity(step); }
    [[nodiscard]] double getHeatCapacityAtFreq(std::size_t freq_index) const noexcept { return controller_->getHeatCapacityAtFreq(freq_index); }
    [[nodiscard]] double getInitTemp() const noexcept { return controller_->getInitTemp(); }
    [[nodiscard]] double getSteadyTemp(std::size_t step=0) const noexcept { return controller_->getSteadyTemp(step); }
    [[nodiscard]] double getArea() const noexcept { return area_covered_; }
    [[nodiscard]] const std::vector<std::array<double, 2>>& getFluxes() const noexcept { return inc_flux_; }
    [[nodiscard]] const std::vector<int>& getEnergies() const noexcept { return inc_energy_; }

    /**
     * Updates the heat parameters (inc_energy_ & inc_flux_) at the given measurement step
     * @param Phonon - The phonon that is transferring heat/flux to the system
     * @param step - The measurement step to update
     */
    void updateHeatParams(const Phonon& p, std::size_t step) noexcept;
    void reset() noexcept;
    void updateTables() { controller_->updateTables(); }
private:
    const std::size_t ID_;
    std::unique_ptr<SensorController> controller_;
    double area_covered_{0.};

    std::vector<int> inc_energy_;
    std::vector<std::array<double, 2>> inc_flux_;
    std::unique_ptr<std::mutex> updateMutex_;
};

struct SensorMeasurements {
    std::size_t id{0};
    double t_steady{0.};
    double std_t_steady{0.};
    double x_flux{0.};
    double std_x_flux{0.};
    double y_flux{0.};
    double std_y_flux{0.};
    // For animations and full analyses
    std::vector<double> final_temps;
    std::vector<std::array<double, 2>> final_fluxes;
};

#endif //GEOMETRY_SENSOR_H