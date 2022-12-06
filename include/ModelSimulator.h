#ifndef GEOMETRY_MODELSIMULATOR_H
#define GEOMETRY_MODELSIMULATOR_H

#include <vector>
#include <variant>

#include "PhononBuilder.h"

class Cell;

class ModelSimulator {
public:
    using BuilderObj = std::variant<CellOriginBuilder, SurfaceOriginBuilder>;
    ModelSimulator(std::size_t measurement_steps, double simulation_time);

    void runSimulation(double t_eq);
    void initPhononBuilders(std::vector<Cell>& cells, double t_eq, double eff_energy) noexcept;
    [[nodiscard]] std::optional<double> nextImpact(Phonon& p, double time) const noexcept;
    void reset() noexcept { total_phonons_ = 0; phonon_builders_.clear(); }
private:
    std::vector<BuilderObj> phonon_builders_;
    std::vector<double> step_times_;
    const double step_time_;

    std::size_t total_phonons_{0};

    void runPhononByPhonon(double t_eq);
    void runUsingBuilders(double t_eq);
    static void scatter(Phonon& p, const std::array<double, 3>& relax_rates) noexcept;
    void simulatePhonon(Phonon&& p, std::size_t measurement_steps) const;
    std::optional<double> handleImpacts(Phonon& p, double drift_time, std::size_t sensor_id) const noexcept;
};

#endif //GEOMETRY_MODELSIMULATOR_H
