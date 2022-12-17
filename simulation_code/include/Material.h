#ifndef GEOMETRY_MATERIAL_H
#define GEOMETRY_MATERIAL_H

#include <array>
#include <vector>

#include "Phonon.h"

struct RelaxationData;
struct DispersionData;
struct TableData;

class Material {
public:
    static constexpr std::size_t NUM_FREQ_BINS = 1000; // Used on one occasion outside this class

    using Array = std::array<double, NUM_FREQ_BINS>;
    using Table = std::array<std::pair<double, double>, NUM_FREQ_BINS>;
    using Polar = Phonon::Polarization;

    Material(std::size_t id, const DispersionData& disp_data, const RelaxationData& relax_data);

    void setFullSimulation() noexcept { full_simulation_ = true; }

    [[nodiscard]] std::size_t id() const noexcept { return id_; }
    [[nodiscard]] double max_freq_la() const noexcept { return w_max_la_; }
    [[nodiscard]] double max_freq_ta() const noexcept { return w_max_ta_; }
    [[nodiscard]] Phonon::RelaxRates relaxRates(double temp, double freq, Polar polarization) const noexcept;
    [[nodiscard]] Phonon::RelaxRates relaxRates(std::size_t freq_index, Polar polarization, double temp) const;
    [[nodiscard]] static std::pair<std::size_t, Polar> freqIndex(const Table& dist) noexcept;

    [[nodiscard]] const Array& getFrequencies() const noexcept { return frequencies_; }
    [[nodiscard]] double getFreq(std::size_t index) const noexcept;
    [[nodiscard]] double getVel(std::size_t index, Polar polar) const noexcept;
    [[nodiscard]] const Table* baseTable(double temp) const { return baseData(temp).first; }
    [[nodiscard]] double baseEnergy(double temp) const { return baseData(temp).second; }
    [[nodiscard]] const Table* emitTable(double temp) const { return emitData(temp).first; }
    // Needed to determine the amount of energy emitted by emitting surfaces
    [[nodiscard]] double emitEnergy(double temp) const { return emitData(temp).second; }
    [[nodiscard]] const Table* scatterTable(double temp) const { return scatterData(temp).first; }
    // Needed for second numerical inversion when doing a full simulation
    [[nodiscard]] double scatterEnergy(double temp) const { return scatterData(temp).second; }
    // If pseudo=true -> scales results by the relaxation rates (returns scatterEnergy(temp))
    [[nodiscard]] double theoreticalEnergy(double temp, bool pseudo=false) const noexcept;
    void initializeTables(double low_temp, double high_temp, float temp_interval);
private:
    const std::size_t id_;
    const double b_l_;
    const double b_tn_;
    const double b_tu_;
    const double b_i_;
    const double w_; // Cutoff for TA Umklapp scattering

    const double w_max_la_; // Maximum frequency of LA phonons that can exist in this material
    const double w_max_ta_;

    const double freq_width_; // Width of frequency bin.
    bool full_simulation_{false};

    Array frequencies_{0.};
    Array densities_la_{0.};
    Array densities_ta_{0.};
    Array velocities_la_{0.};
    Array velocities_ta_{0.};

    std::vector<double> temps_;
    std::vector<TableData> base_tables_;
    std::vector<TableData> emit_tables_;
    std::vector<TableData> scatter_tables_;

    [[nodiscard]] static double getK(double freq, std::array<double, 3> coeffs);
    [[nodiscard]] static double getGv(double freq, std::array<double, 3> coeffs);

    [[nodiscard]] std::pair<const Table*, double> baseData(double temp) const;
    [[nodiscard]] std::pair<const Table*, double> emitData(double temp) const;
    [[nodiscard]] std::pair<const Table*, double> scatterData(double temp) const;

    [[nodiscard]] std::pair<Table, double> cumulDistEmit(Array&& la_dist, Array&& ta_dist) const;
    [[nodiscard]] std::pair<Table, double> cumulDistScatter(Array&& la_dist, Array&& ta_dist, double temp) const;
    [[nodiscard]] static Table buildCumulDist(const Array& t1, const Array& t2);
    [[nodiscard]] Array phononDist(double temp, Polar polarization) const;

    [[nodiscard]] double tauNInv(double temp, double freq, Polar polarization) const noexcept;
    [[nodiscard]] double tauUInv(double temp, double freq, Polar polarization) const noexcept;
    [[nodiscard]] double tauIInv(double freq) const noexcept; // Impurity scattering

    static void verifyInput(double low_temp, double high_temp, float interval_size);
    [[nodiscard]] std::size_t getTempIndex(double temp) const noexcept;
};

// Should extend this to a class
struct RelaxationData {
    double b_l {0.};
    double b_tn {0.};
    double b_tu {0.};
    double b_i {0.};
    double w {0.};
};

// Only handles quadratic data - consider extending to a class as well
struct DispersionData {
    std::array<double, 3> LA_data {0., 0., 0.,};
    std::array<double, 3> TA_data {0., 0., 0.,};
    //std::Array<double, 3> LO_data {0., 0., 0.,};
    //std::Array<double, 3> TO_data {0., 0., 0.,};
    double w_max_la {0.};
    double w_max_ta {0.};
};

struct TableData {
    TableData(Material::Table table, double cumul_sum)
    : table{std::move(table)},
      cumul_sum{cumul_sum}
    {}

    std::array<std::pair<double, double>, Material::NUM_FREQ_BINS> table;
    double cumul_sum {0.};
};

#endif //GEOMETRY_MATERIAL_H
