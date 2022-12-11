#ifndef GEOMETRY_MODEL_H
#define GEOMETRY_MODEL_H

#include <filesystem>
#include <unordered_map>

#include "Cell.h"
#include "ModelSimulator.h"
#include "SensorInterpreter.h"
#include "OutputManager.h"

/**
 * The Model class primarily controls the geometrical aspects of the simulation.
 * Materials must be added first with the addMaterials function. Then the sensors must be added as each
 * sensor is linked to a material. Then the cells can be added as they are linked to sensors.
 * If a cell contains emitting surfaces, they must be added immediately after the cell is added.
 * The class will throw if components are not added in this sequence.
 */
class Model {
public:
    using Point = Geometry::Point;
    using SimulationType = Sensor::SimulationType;

    /**
     * @param num_cells - The number of cells that the model will contain.
     * @param num_sensors - The number of sensors should be <= num_cells
     * @param measurement_steps - The total number of measurement steps (temp/flux recordings)
     * @param num_phonons - The total number of phonons to simulate. More phonons gives more precision but
     * increases simulation time.
     * @param simulation_time - The duration of the simulation
     * @param t_eq - The linearization (equilibrium) temperature of the system - 0 indicates a full simulation
     */
    Model(std::size_t num_cells, std::size_t num_sensors, std::size_t measurement_steps, std::size_t num_phonons,
          double simulation_time, double t_eq);
    /**
     * @param type - The type of simulation. Default is steady state.
     * @param step_interval - For periodic and transient simulations only. The distance between steps for which
     * measurements are recorded.
     * Throws if the step interval is not valid for the given simulation type. This throw can be handled but an error
     * here likely indicates the user is not performing the simulation they think they are.
     */
    void setSimulationType(SimulationType type, std::size_t step_interval=0);
    /**
     * Adds a Material object to the model. A Material is linked to sensor objects which it turn specifies how the
     * cells attached to those sensors should act.
     * @param material_name - A string representing the name of the material
     * @param material - An object that specifies the property of the material
     */
    void addMaterial(const std::string& material_name, Material material);
    /**
     * Adds a Sensor object to the model. Each sensor is linked to a number of cells so that not every cell
     * must independently track its own temperature and flux values in addition to scattering considerations.
     * @param ID - The ID of the sensor, each sensor needs a unique ID. Throws if the ID already exists in the model
     * @param material_name - The name of the material this sensor is linked to. Throws if the material does not exist
     * @param t_init - The initial temperature of the cells linked to the sensor
     * @param type - The type of simulation - steady-state periodic or transient
     */
    void addSensor(std::size_t ID, const std::string& material_name, double t_init, Sensor::SimulationType type);
    /**
     * Adds a triangular cell to the system.
     * Ensures new cells are not contained in existing cells and existing cells do not contain the new cell.
     * Also verifies that the new cell does not intersect any existing cells. Will throw if any of these occur.
     * @param triangle - The underlying shape of the cell
     * @param sensor_ID - The sensor this cell is linked to -> will throw if the sensor does not exist
     * @param spec - The specularity of the cell's boundary surfaces [0-1]
     */
    void addCell(Geometry::Triangle&& triangle, std::size_t sensor_ID, double spec=1.);
    /**
     * This not not used as the python interface takes care of this
     *
     * Adds a rectangular cell to the system. This will count as 2 cells as the surface is broken
     * down into two right angled triangles.
     * @param p1 - The lower left point of the rectangle
     * @param p2 - The upper right point of the rectangle
     * @param sensor_ID - The sensor this cell is linked to -> will throw if the sensor does not exist
     * @param spec - The specularity of the cell's boundary surfaces [0-1]
     */
    void addCell(Point&& p1, Point&& p2, std::size_t sensor_ID, double spec=1.);
    /**
     * Tries to set the line given by p1 and p2 to an emitting surface. A transient surface is set if a start time
     * and duration are specified. Otherwise, the surface emits from phonons from time 0 until the end of the
     * simulation - non transient cases. If a transient surface is specified, the simulation type must be transient.
     *
     * Throws if the input line overlaps with an existing transition or existing emitting surface or if start_times
     * < 0 or >= simulation_time_ and if duration is < 0 or >= simulation_time_-start_time. Also throws if a transient surface
     * is specified but the simulation type is not transient.
     * @param p1 - p1 of a line segment
     * @param p2 - p2 of a line segment
     * @param temp - The temperature of the emitting surface
     * @param duration - The duration of time the transient surface will emit phonons
     * @param start_time - The amount of time into the simulation at which this surface starts emitting phonons
     * @return - True if the surface is successfully added and false if the incoming line segment cannot be turned
     * cannot be found in the existing system.
     */
    [[nodiscard]] bool setEmitSurface(const Point& p1, const Point& p2, double temp, double duration, double start_time);

    void runSimulation();
    void exportResults(const fs::path& filepath) const;
private:
    SimulationType sim_type_{SimulationType::SteadyState};
    const std::size_t num_cells_;
    const std::size_t measurement_steps_;
    const double simulation_time_;
    const std::size_t num_phonons_;
    double t_eq_{0.}; // Changes as the system evolves between runs

    ModelSimulator simulator_;
    OutputManager outputManager_;
    SensorInterpreter interpreter_;
    std::unique_ptr<std::mutex> addMeasurementMutex_;

    std::vector<Cell> cells_;
    std::vector<Sensor> sensors_;
    std::unordered_map<std::string, Material> materials_;

    /**
     * Return a reference to the sensor with the input ID. If the sensor does not exist, an exception is thrown
     * @param ID - The sensor ID
     * @return A reference to a sensor object
     */
    [[nodiscard]] Sensor& getSensor(std::size_t ID);
    [[nodiscard]] double getTotalInitialEnergy() const noexcept;
    /**
     * Find the maximum and minimum possible temperatures of the system. These are used as bounds for the numerical
     * inversions.
     */
    void setTemperatureBounds() noexcept;
    [[nodiscard]] double avgTemp() const;
    void storeResults() noexcept;
    /**
     * @return - The new t_eq if the system needs to be reset and re-run - less than 90% of the model's sensor
     * temperatures are stable (90%) and std::nullopt otherwise
     */
    [[nodiscard]] std::optional<double> resetRequired() noexcept;
    void reset() noexcept;
};


#endif //GEOMETRY_MODEL_H
