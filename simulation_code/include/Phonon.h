#ifndef GEOMETRY_PHONON_H
#define GEOMETRY_PHONON_H

#include <optional>
#include <ostream>

#include "Geometry.h"

class Cell;

class Phonon {
public:
    constexpr std::size_t NUM_RELAX_RATES = 3;
    enum class Polarization { LA, TA };

    using RelaxRates = std::array<double, NUM_RELAX_RATES>;

    Phonon(signed char sign, double lifetime, Cell* cell);

    [[nodiscard]] signed char getSign() const noexcept { return sign_; }
    [[nodiscard]] std::pair<double, double> getPosition() const noexcept { return {px_, py_}; }
    [[nodiscard]] std::pair<double, double> getVelVector() const noexcept { return {dx_ * velocity_, dy_ * velocity_}; }
    [[nodiscard]] std::pair<double, double> getDirection() const noexcept { return {dx_, dy_}; }
    [[nodiscard]] std::size_t getFreqIndex() const noexcept { return freq_index_; }
    [[nodiscard]] double getFreq() const noexcept { return freq_; }
    [[nodiscard]] Polarization getPolar() const noexcept { return polar_; }
    [[nodiscard]] double getLifetime() const noexcept { return lifetime_; }
    [[nodiscard]] std::size_t getLifeStep() const noexcept { return lifestep_; }
    [[nodiscard]] bool outsideCell() const noexcept { return cell_ == nullptr; }

    void scatterUpdate(std::size_t freq_index, double freq, double velocity, Polarization polar) noexcept;
    void setPosition(double px, double py) noexcept { px_ = px; py_ = py; }
    void setDirection(double dx, double dy) noexcept { dx_ = dx;  dy_ = dy; }
    void setCell(Cell* cell) noexcept { cell_ = cell; }
    void setLifeStep(std::size_t step) { lifestep_ = step; }
    void scatterUpdate();
    void drift(double time) noexcept;
    void setRandDirection(double r1, double r2) noexcept;

    // All these methods will throw if the phonon is not in a cell (cell_ == nullptr)
    [[nodiscard]] std::size_t getCellSensorID() const;
    [[nodiscard]] std::size_t getCellMaterialID() const;
    [[nodiscard]] double getCellHeatCapacityAtFreq(std::size_t index) const;
    [[nodiscard]] RelaxRates getRelaxRates(std::size_t step) const;
    [[nodiscard]] std::array<Geometry::Line, 3> getCellBoundaryLines() const;
    void handleSurfaceCollision(const Geometry::Point& poi, double step_time);
    void updateCellHeatParams(std::size_t step) const;
    void setRandPoint(double r1, double r2);

    friend std::ostream &operator<<(std::ostream& os, const Phonon& phonon) {
        os << "px: " << phonon.px_ << " py: " << phonon.py_
           << "\nvx_: " << phonon.dx_ * phonon.velocity_ << " vy_: " << phonon.dy_ * phonon.velocity_;
        return os;
    }
private:
    const signed char sign_; // Determines how the phonon energy/flux is handled (-1, 1)
    const double lifetime_;
    std::size_t lifestep_{0};

    double px_{0.};
    double py_{0.};
    double dx_{0.};
    double dy_{0.};

    std::size_t freq_index_{0};
    double freq_{0.};
    double velocity_{0.};
    Polarization polar_{Polarization::LA};

    Cell* cell_{nullptr};
};

#endif //GEOMETRY_PHONON_H
