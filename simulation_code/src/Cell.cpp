#include <numeric>
#include <cmath>
#include <execution>
#include <sstream> // for exception message

#include "Cell.h"

using Line = Geometry::Line;
using Point = Geometry::Point;

// Assumes triangle has no intersecting surfaces - this is handled by the triangle class
Cell::Cell(Cell::Triangle cell, Sensor& sensor, double spec)
    : cell_{std::move(cell)}, sensor_{sensor}, boundaries_{buildCompositeSurfaces(spec)}
{
    sensor_.addToArea(getArea());
}

// Checks if an incoming cell is valid if this cell already exists in the model
void Cell::validate(const Cell& other) const {
    if (cell_.intersects(other.cell_)) {
        throw IntersectError(this->cell_, other.cell_);
    }
    if (cell_.contains(other.cell_)) {
        throw OverlapError(this->cell_, other.cell_);
    }
    if (other.cell_.contains(cell_)) {
        throw OverlapError(other.cell_, this->cell_);
    }
}

// Sets an emitting surface on the edge of this cell based on the coordinates of line
// returns false if the incoming line is not a valid location in this cell (the line is not part of the cell surfaces)
bool Cell::setEmitSurface(const Line& line, double temp, double duration, double start_time) {
    const int norm_sign = (cell_.isClockwise()) ? 1 : -1;
    for (auto& boundary : boundaries_) {
        if (boundary.addEmitSurface(line, *this, getMaterial(), temp, norm_sign, duration, start_time)) {
            return true;
        }
    }
    return false;
}

// Gets this initial amount of energy in this cell based on the cell starting temperature and size
double Cell::getInitEnergy(double t_eq) const noexcept {
    const double init_energy = getArea() * sensor_.getHeatCapacity();
    return (t_eq == 0.) ? init_energy : init_energy * std::fabs(getInitTemp() - t_eq);
}

// Gets the amount of energy that is emitted from surfaces contained in this cell over the duration of the simulation
// Will be 0 if the cell has no emitting/isothermal surfaces
double Cell::getEmitEnergy(double t_eq) const noexcept {
    double emit_energy = 0.;
    for (const auto& surface : getBoundaries()) {
        const auto& emit_surfaces = surface.getEmitSurfaces();
        // sum the energy of all the emitting surfaces contained within this cell
        emit_energy += std::transform_reduce(std::execution::seq, std::cbegin(emit_surfaces),
                                             std::cend(emit_surfaces), 0., std::plus{}, [&](const auto& sub_surface) {
            const auto temp = sub_surface.getTemp();
            const double emit_energy = sub_surface.getSurfaceLine().length * sub_surface.getEmitDuration()
                                     * getMaterial().emitEnergy(temp) / 4.;
            return (t_eq == 0.) ? emit_energy : emit_energy * std::fabs(temp - t_eq);
        });
    }
    return emit_energy;
}

// Updates the tables the emitting/isothermal surfaces use to generate phonons. Since the emitting surfaces
// are at a constant temperature for the duration of the simulation, this is only called once during the initialization
// phase and the tables are generated using the emitting temperature of whichever emitting surfaces are contained
// in the cell
void Cell::updateEmitTables() noexcept {
    for (auto& surface : boundaries_) {
        surface.updateEmitSurfaceTables();
    }
}

// The incoming phonon contributes its energy & flux to the sensor that is linked to this cell. The step refers to the
// time frame where this contribution occurs
void Cell::updateHeatParams(const Phonon& p, std::size_t step) noexcept {
    sensor_.updateHeatParams(p, step);
}

// Assumes cells can only have a single transition surface between them - this will likely remain a constraint that
// can be worked around with an individual boundary surface placement feature
void Cell::findTransitionSurface(Cell& other) {
    const auto& lines = getBoundaryLines();
    const auto& other_lines = other.getBoundaryLines();

    for (const auto& l1 : lines) {
        for (const auto& l2 : other_lines) {
            if (l1.contains(l2)) { // The only potential transition area is found
                if (!setTransitionSurface(l2, other) || !other.setTransitionSurface(l2, *this)) {
                    throw CompositeSurfaceError(l1, l2);
                }
            } else if (l2.contains(l1)) {
                if (!setTransitionSurface(l1, other) || !other.setTransitionSurface(l1, *this)) {
                    throw CompositeSurfaceError(l2, l1);
                }
            }
        }
    }
}

// Handles the interaction of a phonon colliding with the surface in which poi is located.
// The step_time is only needed for transient simulations. It is used to check whether a surface
// is currently acting as an emitting surface or whether it is currently acting as a boundary surface.
void Cell::handleSurfaceCollision(Phonon& p, const Point& poi, double step_time) const noexcept {
    for (const auto& boundary : boundaries_) {
        if (boundary.contains(poi)) {
            boundary.handlePhonon(p, poi, step_time);
            return;
        }
    }
}

std::array<Line, 3> Cell::getBoundaryLines() const noexcept {
    return { boundaries_[0].getSurfaceLine(), boundaries_[1].getSurfaceLine(), boundaries_[2].getSurfaceLine() };
}

std::array<CompositeSurface, 3> Cell::buildCompositeSurfaces(double spec) noexcept {
    if (spec < 0.) { spec = 0.; }
    else if (spec > 1.) { spec = 1.; }
    const int norm_sign = (cell_.isClockwise()) ? 1 : -1;
    return { CompositeSurface{Surface(Line{cell_.p1, cell_.p2}, *this, spec, norm_sign)},
             CompositeSurface{Surface(Line{cell_.p2, cell_.p3}, *this, spec, norm_sign)},
             CompositeSurface{Surface(Line{cell_.p3, cell_.p1}, *this, spec, norm_sign)}
           };
}

bool Cell::setTransitionSurface(const Line& line, Cell& cell) {
    const int norm_sign = (cell_.isClockwise()) ? 1 : -1;
    for (auto& boundary : boundaries_) {
        if (boundary.addTransitionSurface(line, cell, norm_sign)) {
            return true;
        }
    }
    return false;
}

std::ostream& operator<<(std::ostream& os, const Cell& cell) {
    os << "(" << cell.sensor_.getID() << ") " << cell.cell_ << '\n' <<
                                cell.boundaries_[0] << '\n' <<
                                cell.boundaries_[1] << '\n' <<
                                cell.boundaries_[2] << '\n';
    return os;
}

bool Cell::operator==(const Cell& rhs) const {
    return cell_ == rhs.cell_;
}

bool Cell::operator!=(const Cell& rhs) const {
    return !(rhs == *this);
}

IntersectError::IntersectError(Geometry::Triangle existing, Geometry::Triangle incoming)
    : existing_{std::move(existing)}, incoming_{std::move(incoming)}
{
    std::ostringstream os;
    os << "Incoming " << incoming_ << "\nintersects\nExisting" << existing_ << '\n';
    setMessage(os.str());
}

OverlapError::OverlapError(Geometry::Triangle bigger, Geometry::Triangle smaller)
    : bigger_{std::move(bigger)}, smaller_{std::move(smaller)}
{
    std::ostringstream os;
    os << smaller_ << "\nis contained within\n" << bigger_ << '\n';
    setMessage(os.str());
}
