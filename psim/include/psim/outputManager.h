#ifndef PSIM_OUTPUTMANAGER_H
#define PSIM_OUTPUTMANAGER_H

#include "utils.h"
#include <filesystem>
#include <unordered_map>
#include <vector>

struct SensorMeasurements;
namespace fs = std::filesystem;

class OutputManager {
public:
    OutputManager() = default;
    /**
     * Exports results from each measurement step so the evolution of the system can be visualized.
     * Step intervals of 1 will write every measurement to file. If there are 100 measurements,
     * then there will be 100 * num_sensors entries written. Step intervals greater than 1 will write the average of
     * the measurements to file. Step intervals of 10 with 100 measurements -> 10 * num_sensors entries. First entry is
     * the avg of steps 1-10, second is avg of steps 11-20, etc.
     * @param filepath - path and filename where the results will be written - existing file will be overwritten
     * @param time - The time taken to run the simulation.
     * @param num_runs - The number of runs the simulation was executed for.
     * @param type - The type of simulation that was run.
     */
    void exportResults(const fs::path& filepath, double time, std::size_t num_runs, SimulationType type) const;
    void addMeasurement(std::size_t runId, SensorMeasurements&& measurement) noexcept;
    void sortMeasurements(std::size_t runId) noexcept;
    void setStepInterval(std::size_t interval) noexcept {
        step_interval_ = interval;
    }

private:
    std::size_t step_interval_{ 1 };
    std::unordered_map<std::size_t, std::vector<SensorMeasurements>> measurements_;

    std::vector<SensorMeasurements> calculateAverages() const;
    void steadyStateExport(std::ofstream& output) const;
    void periodicExport(std::ofstream& output) const;
    [[nodiscard]] static std::filesystem::path adjustPath(const fs::path& filepath, const std::string& prepend);
    [[nodiscard]] static std::string getCurrentDateTime();
};


#endif// PSIM_OUTPUTMANAGER_H
