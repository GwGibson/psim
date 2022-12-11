#ifndef GEOMETRY_OUTPUTMANAGER_H
#define GEOMETRY_OUTPUTMANAGER_H

#include <vector>
#include <string>
#include <ostream>
#include <filesystem>

namespace fs = std::filesystem;

struct SensorMeasurements;

class OutputManager {
public:
    OutputManager() = default;
    void steadyStateExport(const fs::path& filepath) const;
    /**
     * Exports results from each measurement step so the evolution of the system can be visualized.
     * Step intervals of 1 will write every measurement to file. If there are 100 measurements,
     * then there will bee 100 * num_sensors entries written. Step intervals greater than 1 will write the average of
     * the measurements to file. Step intervals of 10 with 100 measurements -> 10 * num_sensors entries. First entry is
     * the avg of steps 1-10, second is avg of steps 11-20, etc.
     * @param filepath - path and filename where the results will be written - existing file will be overwritten
     * @param step_intervals - The distance between steps. If > 1, the average of the measurements is used.
     */
    void periodicExport(const fs::path& filepath) const;
    void addMeasurement(SensorMeasurements&& measurement) noexcept;
    void sortMeasurements() noexcept;
    void setStepInterval(std::size_t interval) noexcept { step_interval_ = interval; }
private:
    std::size_t step_interval_{1};
    std::vector<SensorMeasurements> measurements_;

    [[nodiscard]] static std::filesystem::path adjustPath(const fs::path& filepath, const std::string& prepend);
    [[nodiscard]] static std::string getCurrentDateTime();
};


#endif //GEOMETRY_OUTPUTMANAGER_H




