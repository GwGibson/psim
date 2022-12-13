#include "PhononBuilder.h"
#include "Utils.h"
#include "Cell.h"

Phonon CellOriginBuilder::operator()(double t_eq) noexcept {
    auto& [cell, phonons] = cells_.top();
    const signed char sign = (cell->getInitTemp() > t_eq) ? 1 : -1;
    Phonon p{sign, 0., cell};
    cell->initialUpdate(p); // Use the base_table_ in the cell's sensor
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

Phonon SurfaceOriginBuilder::operator()(double t_eq) noexcept {
    --total_phonons_;
    const signed char sign = (surface_.getTemp() > t_eq) ? 1 : -1;
    Phonon p{sign, surface_.getPhononTime(), &cell_};
    cell_.initialUpdate(p, surface_.getTable()); // Phonon Freq, Velocity & Polarization set here
    const auto& [px, py] = surface_.getRandPoint(Utils::urand()); // Use the surface's emitting table
    p.setPosition(px, py);
    surface_.redirectPhonon(p); // Phonon direction is set in a biased manner based on the surface orientation in space
    return p;
}

Phonon PhasorBuilder::operator()(double t_eq) noexcept {
    const auto& [nx, ny] = surface_.getNormal();
    Phonon p = SurfaceOriginBuilder::operator()(t_eq);
    // Not ideal, but it works
    p.scatterUpdate(1,1,1000,Phonon::Polarization::LA); // All irrelevant except velocity
    p.setDirection(nx,ny);
    return p;
}
