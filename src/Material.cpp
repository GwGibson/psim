#define _USE_MATH_DEFINES   // Allow for M_PI
#include <cmath>            // M_PI & trig functions
#include <numeric>          // accumulate
#include <utility>          // std::as_const
#include <algorithm>        // generate

#include "Material.h"
#include "Utils.h"

namespace {
    static constexpr double HBAR = 1.054517e-34;
    static constexpr double BOLTZ = 1.38065e-23;
}

using Array = Material::Array;
using Table = Material::Table;
using Polar = Material::Polar;

Material::Material(std::size_t id, const DispersionData& disp_data, const RelaxationData& relax_data) :
        id_{id},
        b_l_{relax_data.b_l},
        b_tn_{relax_data.b_tn},
        b_tu_{relax_data.b_tu},
        b_i_{relax_data.b_i},
        w_{relax_data.w},
        w_max_la_{disp_data.w_max_la},
        w_max_ta_{disp_data.w_max_ta},
        freq_width_{std::max(w_max_la_, w_max_ta_) / NUM_FREQ_BINS}
{
    const auto& [LA_coeffs, TA_coeffs, _unused, _unused1] = disp_data;
    // Build frequency Table
    std::generate(std::begin(frequencies_), std::end(frequencies_),
                  [n = -1., this]() mutable { return (n += 2) * freq_width_ / 2.; });

    // Build velocity and state density tables
    std::size_t index = 0;
    for (const auto& freq : frequencies_) {
        const double la_gv = getGv(freq, LA_coeffs);
        velocities_la_[index] = la_gv;
        densities_la_[index] = pow(getK(freq, LA_coeffs), 2) / 2. / pow(M_PI, 2) / la_gv;
        const double ta_gv = getGv(freq, TA_coeffs);
        if (!std::isnan(ta_gv)) {
            velocities_ta_[index] = ta_gv;
            densities_ta_[index] = pow(getK(freq, TA_coeffs), 2) / 2. / pow(M_PI, 2) / ta_gv;
        }
        ++index;
    }
}

// Returns the scattering rates in the following order [N_rate, U_rate, I_rate]
std::array<double, 3> Material::scatteringRates(double temp, double freq, Polar polarization) const noexcept {
    return std::array {tauNInv(temp, freq, polarization), tauUInv(temp, freq, polarization), tauIInv(freq)};
}

std::pair<std::size_t, Polar> Material::freqIndex(const Table& dist) noexcept {
    std::size_t low = 0;
    std::size_t high = dist.size()-1;
    auto mid = static_cast<std::size_t>(low + (high - low)/2.);
    const double rand = Utils::urand();
    // Perform bisection search
    while (high - low > 1) {
        rand < dist[mid].first ? high = mid : low = mid;
        mid = floor((low + high)/2);
    }
    return {high, Utils::urand() <= dist[high].second ? Polar::LA : Polar::TA};
}

std::size_t Material::getFreqIndex(double frequency) const noexcept {
    return std::distance(std::cbegin(frequencies_), std::max_element(std::cbegin(frequencies_),
                                                    std::find_if(std::cbegin(frequencies_), std::cend(frequencies_),
                                                    [&frequency](double freq) { return freq >= frequency;})));
}

double Material::getFreq(std::size_t index) const noexcept {
    // TODO: This causes full simulations (teq == 0) to give bizarre results. Unsure why this happens.
    //return frequencies_[index] + (2. * Utils::urand() - 1.) * freq_width_ / 2.;
    return frequencies_[index];
}

double Material::getVel(std::size_t index, Polar polar) const noexcept {
    return (polar == Polar::LA) ? velocities_la_[index] : velocities_ta_[index];
}

double Material::getK(double freq, std::array<double, 3> coeffs) {
    const double d = pow(coeffs[1], 2) - 4.*coeffs[0]*(coeffs[2]-freq);
    const double a = (-coeffs[1]-sqrt(d)) / (2.*coeffs[0]);
    const double b = (-coeffs[1]+sqrt(d)) / (2.*coeffs[0]);
    return (a < b) ? a : b;
}

double Material::getGv(double freq, std::array<double, 3> coeffs) {
    return 2. * coeffs[0] * getK(freq, coeffs) + coeffs[1];
}

std::pair<Table, double> Material::baseData(double temp) const {
    Array la_base = phononDist(temp, Polar::LA);
    Array ta_base = phononDist(temp, Polar::TA);
    const auto heat_capacity = std::accumulate(std::cbegin(la_base), std::cend(la_base), 0.)
                             + std::accumulate(std::cbegin(ta_base), std::cend(ta_base), 0.);
    return {buildCumulDist(la_base, ta_base), heat_capacity};
}

Table Material::emitTable(double temp) const {
    return emitData(temp).first;
}

double Material::emitEnergy(double temp) const {
    return emitData(temp).second;
}

Table Material::scatterTable(double temp) const {
    return scatterData(temp).first;
}

double Material::scatterEnergy(double temp) const {
    return scatterData(temp).second;
}

double Material::theoreticalEnergy(double temp, bool pseudo) const noexcept {
    return (pseudo) ? scatterEnergy(temp) : baseData(temp).second;
}

std::pair<Table, double> Material::emitData(double temp) const {
    return cumulDistEmit(phononDist(temp, Polar::LA), phononDist(temp, Polar::TA));
}

std::pair<Table, double> Material::scatterData(double temp) const {
    return cumulDistScatter(phononDist(temp, Polar::LA), phononDist(temp, Polar::TA), temp);
}

std::pair<Table, double> Material::cumulDistEmit(Array&& la_dist, Array&& ta_dist) const {
    auto transform = [](auto& dist, const auto& velocities){
        std::transform(std::cbegin(dist), std::cend(dist), std::cbegin(velocities),
                       std::begin(dist), std::multiplies<>());
    };
    transform(la_dist, velocities_la_);
    transform(ta_dist, velocities_ta_);
    return {buildCumulDist(la_dist, ta_dist), std::accumulate(std::cbegin(la_dist), std::cend(la_dist), 0.)
                                            + std::accumulate(std::cbegin(ta_dist), std::cend(ta_dist), 0.) };
}

std::pair<Table, double> Material::cumulDistScatter(Array&& la_dist, Array&& ta_dist, double temp) const {
    auto transform = [&, this](auto& dist, const auto& polar){
        std::transform(std::cbegin(dist), std::cend(dist), std::cbegin(frequencies_), std::begin(dist),
                       [&, this](auto& elem, const auto& freq) {
            const auto scatter_rates = scatteringRates(temp, freq, polar);
            return elem * std::accumulate(std::cbegin(scatter_rates), std::cend(scatter_rates), 0.);
        });
    };
    transform(la_dist, Polar::LA);
    transform(ta_dist, Polar::TA);
    return {buildCumulDist(la_dist, ta_dist), std::accumulate(std::cbegin(la_dist), std::cend(la_dist), 0.)
                                            + std::accumulate(std::cbegin(ta_dist), std::cend(ta_dist), 0.) };
}

Table Material::buildCumulDist(const Array& t1, const Array& t2) {
    Table cumul_dist;
    double const cumul_sum = std::accumulate(std::cbegin(t1), std::cend(t1), 0.)
                                 + std::accumulate(std::cbegin(t2), std::cend(t2), 0.);
    // Algorithmic way to do this?
    cumul_dist[0] = { (t1[0] + t2[0]) / cumul_sum, t1[0] / (t1[0] + t2[0]) };
    for (std::size_t i = 1; i != NUM_FREQ_BINS; ++i) {
        cumul_dist[i] = { cumul_dist[i-1].first + (t1[i] + t2[i]) / cumul_sum,
                          t1[i] / (t1[i] + t2[i]) };
    }
    return cumul_dist;
}

// Returns an array where each entry represents the number of phonons that a present in the corresponding
// frequency bin. The sum of all the entries is the total number of phonons per volume in this material.
Array Material::phononDist(double temp, Polar polarization) const {
    // Build bose-einstein Array (including (HBAR * freq) in original distribution - hw/(exp(hw/(k_b*T)-1)
    Array bose_einstein;
    const auto& densities = (polarization == Polar::LA) ? densities_la_ : densities_ta_;
    // This will need revision if the optical branches are added
    if (fullSimulation_) { // Simulate all phonons and use BE for phonon distributions
        std::transform(std::cbegin(frequencies_), std::cend(frequencies_), std::cbegin(densities),
                       std::begin(bose_einstein),
                       [&, const_calc = HBAR / (BOLTZ * temp)](const auto freq, const auto density) {
                           const auto dist = freq * HBAR / (exp(const_calc*freq)-1.) * freq_width_ * density;
                           return (polarization == Polar::LA) ? dist : 2. * dist;
                       });
    } else { // Only simulate phonons that deviate from t_eq & use BE derivative for phonon distributions
        std::transform(std::cbegin(frequencies_), std::cend(frequencies_), std::cbegin(densities),
                       std::begin(bose_einstein),
                       [&, const_calc = HBAR / (BOLTZ * temp)](const auto freq, const auto density) {
                           const auto dist = pow(const_calc * freq, 2) * BOLTZ * exp(const_calc * freq) /
                                             pow((exp(const_calc * freq) - 1.), 2) *
                                             freq_width_ * density;
                           return (polarization == Polar::LA) ? dist : 2. * dist;
                       });
    }
    return bose_einstein;
}

double Material::tauNInv(double temp, double freq, Polar polarization) const noexcept {
    const auto tau_inv = [&](){
        switch(polarization) {
            case Polar::LA:
                return b_l_ * freq * freq * pow(temp, 3);
            case Polar::TA:
                if (freq < w_) {
                    return b_tn_ * freq * pow(temp, 4);
                }
                return 0.;
            default:
                return 0.;
        }
    }();
    return tau_inv;
}

double Material::tauUInv(double temp, double freq, Polar polarization) const noexcept {
    const auto tau_inv = [&](){
        switch(polarization) {
            case Polar::LA:
                return b_l_ * freq * freq * pow(temp, 3);
            case Polar::TA:
                if (freq >= w_) {
                    return b_tu_ * freq * freq / sinh(HBAR * freq / (temp * BOLTZ));
                }
                return 0.;
            default:
                return 0.;
        }
    }();
    return tau_inv;
}

// Impurity Scattering
double Material::tauIInv(double freq) const noexcept {
    return b_i_ * pow(freq, 4);
}





