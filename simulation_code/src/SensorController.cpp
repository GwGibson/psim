#include <tuple> // tie
#include <algorithm> // generate_n
#include <cmath> // fabs

#include "SensorController.h"

namespace {
    // Sensor temperature must be within this percentage to be considered stable -> 0.005 = 0.5%
    constexpr double RESET_THRESHOLD{0.005};
    // Temperature at each measurement step must be within this percentage to be considered stable
    // compared to the temperature at that step on the previous simulation iteration
    constexpr double TRANSIENT_RESET_THRESHOLD{0.02};
}

SensorController::SensorController(const Material& material, double t_init, std::size_t num_measurements)
    : material_{material},
      t_init_{t_init},
      num_measurements_{num_measurements},
      t_steady_{t_init}
{}

double SensorController::getHeatCapacityAtFreq(std::size_t freq_index) const noexcept {
    return (freq_index == 0) ? (*base_table_)[0].first : (*base_table_)[freq_index].first - (*base_table_)[freq_index-1].first;
}

void SensorController::initialUpdate(Phonon& p, const Material::Table& table) const noexcept {
    const auto& [index, polar] = Material::freqIndex(table);
    p.scatterUpdate(index, material_.getFreq(index), material_.getVel(index, polar), polar);
}

void SensorController::initialUpdate(Phonon& p) const noexcept {
    const auto& [index, polar] = Material::freqIndex(*base_table_);
    p.scatterUpdate(index, material_.getFreq(index), material_.getVel(index, polar), polar);
}

void SensorController::updateTables() {
    base_table_ = material_.baseTable(t_init_);
    heat_capacity_ = material_.baseEnergy(t_init_);
    scatter_table_ = material_.scatterTable(t_init_);
    if (num_measurements_ > 0) { // Set up vectors for a transient simulation
        scatter_tables_.resize(num_measurements_);
        heat_capacities_.resize(num_measurements_);
        steady_temps_.resize(num_measurements_);
        std::generate_n(std::begin(scatter_tables_), num_measurements_, [this](){ return scatter_table_; });
        std::generate_n(std::begin(heat_capacities_), num_measurements_, [this](){ return heat_capacity_; });
        std::generate_n(std::begin(steady_temps_), num_measurements_, [this](){ return t_init_; });
    }
}

void SensorController::scatterUpdate(Phonon& p) const noexcept {
    const auto& [index, polar] = Material::freqIndex(*scatter_table_);
    p.scatterUpdate(index, material_.getFreq(index), material_.getVel(index, polar), polar);
}

bool SensorController::resetRequired(double t_final, std::vector<double>&&) noexcept {
    const auto t_diff = std::fabs(t_final - t_steady_);
    const bool temp_stable = t_diff/t_steady_ <= RESET_THRESHOLD || t_diff < 1.;
    t_steady_ = t_final;
    return temp_stable;
}

void SteadyStateController::reset() noexcept {
    // Update tables and heat capacity using the steady state temperature at the end of the current run
    base_table_ = material_.baseTable(t_steady_);
    heat_capacity_ = material_.baseEnergy(t_steady_);
    scatter_table_ = material_.scatterTable(t_steady_);
}

void PeriodicController::reset() noexcept {
    // Do not update the heat capacity -> this will have no effect for full simulations but will
    // further limit the temperature ranges of approximation simulations.
    base_table_ = material_.baseTable(t_steady_);
    scatter_table_ = material_.scatterTable(t_steady_);
}

double TransientController::getSteadyTemp(std::size_t step) const noexcept {
    return (step == 0) ? t_init_ : steady_temps_[step];
}

void TransientController::scatterUpdate(Phonon &p) const noexcept {
    // Get the scatter at the correct measurement step
    const auto& [index, polar] = Material::freqIndex(*scatter_tables_[p.getLifeStep()]);
    p.scatterUpdate(index, material_.getFreq(index), material_.getVel(index, polar), polar);
}

bool TransientController::resetRequired(double, std::vector<double>&& final_temps) noexcept {
    const bool temp_stable = std::equal(std::cbegin(steady_temps_), std::cend(steady_temps_), std::cbegin(final_temps),
                   [](double t1, double t2) { return std::fabs(t2 - t1)/t1 <= TRANSIENT_RESET_THRESHOLD; });
    steady_temps_ = final_temps;
    return temp_stable;
}

void TransientController::reset() noexcept {
    std::transform(std::cbegin(steady_temps_), std::cend(steady_temps_), std::begin(heat_capacities_),
                   [this](double temp){ return material_.baseEnergy(temp); });
    std::transform(std::cbegin(steady_temps_), std::cend(steady_temps_), std::begin(scatter_tables_),
                   [this](double temp){ return material_.scatterTable(temp); });
}




