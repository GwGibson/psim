#ifndef GEOMETRY_PHONONBUILDER_H
#define GEOMETRY_PHONONBUILDER_H

#include <stack>

#include "Phonon.h"

class EmitSurface;
class Cell;

// TODO: is abstract class necessary here - using variants in ModelSimulator
class PhononBuilder {
public:
    PhononBuilder() = default;
    virtual ~PhononBuilder() = default;
    /**
     * Generates a phonon in accordance to the builder specifications
     * @param t_eq - Equilibrium temperature of the system
     * @return A phonon that is created according to the builder specifications
     * or nullopt if no phonons remain in the builder
     */
    [[nodiscard]] virtual std::optional<Phonon> operator()(double t_eq) noexcept = 0;
    [[nodiscard]] virtual bool hasPhonons() const noexcept = 0;

    [[nodiscard]] std::size_t totalPhonons() const noexcept { return total_phonons_; }
protected:
    std::size_t total_phonons_{0};
};

class CellOriginBuilder : public PhononBuilder {
public:
    CellOriginBuilder() = default;

    [[nodiscard]] std::optional<Phonon> operator()(double t_eq) noexcept override;
    [[nodiscard]] bool hasPhonons() const noexcept override { return !cells_.empty(); }

    void addCellPhonons(Cell* cell, std::size_t num_phonons) noexcept;
private:
    std::stack<std::pair<Cell*, std::size_t>> cells_;
};

class SurfaceOriginBuilder : public PhononBuilder {
public:
    SurfaceOriginBuilder(Cell& cell, const EmitSurface& surface, std::size_t num_phonons);

    [[nodiscard]] std::optional<Phonon> operator()(double t_eq) noexcept override;
    [[nodiscard]] bool hasPhonons() const noexcept override { return total_phonons_ > 0; }
private:
    Cell& cell_;
    const EmitSurface& surface_;
};

#endif //GEOMETRY_PHONONBUILDER_H
