#include "PhononBuilder.h"
#include "Utils.h"
#include "Cell.h"

std::optional<Phonon> CellOriginBuilder::operator()(double t_eq) noexcept {
    if (!hasPhonons()) {
        return std::nullopt;
    }
    auto& [cell, phonons] = cells_.top();
    const int sign = (cell->getInitTemp() > t_eq) ? 1 : -1;
    Phonon p{sign, 0., cell};
    cell->initialUpdate(p, nullptr); // Use the base_table_ in the cell's sensor
    const auto& [px, py] = cell->getRandPoint(Utils::urand(), Utils::urand());
    p.setPosition(px, py);
    p.setRandDirection(Utils::urand(), Utils::urand());
    if (--phonons == 0) {
        cells_.pop();
    }
    return p;
}

void CellOriginBuilder::addCellPhonons(Cell* cell, std::size_t num_phonons) noexcept {
    if (num_phonons > 0) {
        total_phonons_ += num_phonons;
        cells_.emplace(std::pair{cell, num_phonons});
    }
}

SurfaceOriginBuilder::SurfaceOriginBuilder(Cell& cell, const EmitSurface& surface, std::size_t num_phonons)
    : cell_{cell}, surface_{surface}
{
    total_phonons_ += num_phonons;
}

std::optional<Phonon> SurfaceOriginBuilder::operator()(double t_eq) noexcept {
    if (!hasPhonons()) {
        return std::nullopt;
    }
    --total_phonons_;
    const int sign = (surface_.getTemp() > t_eq) ? 1 : -1;
    Phonon p{sign, surface_.getPhononTime(), &cell_};
    cell_.initialUpdate(p, surface_.getTable());
    const auto& [px, py] = surface_.getRandPoint(Utils::urand()); // Use the surface's emitting table
    p.setPosition(px, py);
    surface_.redirectPhonon(p);
    return p;
}