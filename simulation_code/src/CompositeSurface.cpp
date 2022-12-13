#include <algorithm> // std::sort
#include <sstream> // for exception message

#include "CompositeSurface.h"

using Line = Geometry::Line;

CompositeSurface::CompositeSurface(Surface&& main_surface) : main_surface_{main_surface} {}


void CompositeSurface::updateEmitSurfaceTables() noexcept {
    for (auto& es : emit_sub_surfaces_) {
        es.updateTable();
    }
}

bool CompositeSurface::addEmitSurface(const Line& surface_line, Cell& cell, const Material& mat, double temp,
                                      int norm_sign, double duration, double start_time) {
    if (verifySurfaceLine(surface_line)) {
        emit_sub_surfaces_.emplace_back(EmitSurface{surface_line, cell, main_surface_.getSpecularity(),
                                                     norm_sign, mat, temp, duration, start_time});
        // Ensure normal is correct
        emit_sub_surfaces_.back().setNormal(main_surface_.getNormal());
        return true;
    }
    return false;
}

bool CompositeSurface::addTransitionSurface(const Line& surface_line, Cell& cell, int norm_sign) {
    if (verifySurfaceLine(surface_line)) {
        // Transition surfaces cause the phonon to diffusely (spec=0) backscatter when they lack a corresponding
        // state in the new material (frequency too high)
        transition_sub_surfaces_.emplace_back(TransitionSurface{surface_line, cell, 0., norm_sign});
        return true;
    }
    return false;
}

void CompositeSurface::handlePhonon(Phonon& p, const Point& poi, double step_time) const noexcept {
    // Check to see if the phonon interacts with sub-surfaces
    for (const auto& ts : transition_sub_surfaces_) {
        if (ts.getSurfaceLine().contains(poi)) {
            ts.handlePhonon(p);
            return;
        }
    }
    for (const auto& es : emit_sub_surfaces_) {
        if (es.getSurfaceLine().contains(poi)) {
            es.handlePhonon(p, step_time);
            return;
        }
    }
    // If not, it must have hit the main surface
    main_surface_.boundaryHandlePhonon(p);
}

bool CompositeSurface::verifySurfaceLine(const Line& surface_line) const {
    const auto &main = main_surface_.getSurfaceLine();
    // If the incoming surface is contained within the primary surface
    if (main.contains(surface_line)) {
        for (const auto& ts : transition_sub_surfaces_) {
            // If the incoming subsurface overlaps with an existing subsurface -> invalid user specifications
            if (const auto& existing = ts.getSurfaceLine(); surface_line.overlaps(existing)) {
                throw CompositeSurfaceError(existing, surface_line);
            }
        }
        for (const auto& es : emit_sub_surfaces_) {
            // If the incoming subsurface overlaps with an existing subsurface -> invalid user specifications
            if (const auto& existing = es.getSurfaceLine(); surface_line.overlaps(existing)) {
                throw CompositeSurfaceError(existing, surface_line);
            }
        }
        return true;
    }
    return false;
}

std::ostream& operator<<(std::ostream& os, const CompositeSurface& surface) {
    os << "Main Surface: " << surface.main_surface_.getSurfaceLine() << '\n';
    for (const auto& ts : surface.transition_sub_surfaces_) {
        os << '\t' << "Transition Surface: " << ts.getSurfaceLine() << '\n';
    }
    for (const auto& es : surface.emit_sub_surfaces_) {
        os << '\t' << "Emit Surface: " << es.getSurfaceLine() << '\n';
    }
    return os;
}

CompositeSurfaceError::CompositeSurfaceError(Line main, Line inc)
: main_{std::move(main)}, inc_{std::move(inc)}
{
    std::ostringstream os;
    os << "An existing surface conflicts with the location of the incoming surface.\n";
    os << "Existing surface points are: " << main_.p1 << ' ' << main_.p2 << '\n';
    os << "Incoming surface points are: " << inc_.p1 << ' ' << inc_.p2 << '\n';
    setMessage(os.str());
}