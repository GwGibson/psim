#define _USE_MATH_DEFINES   // Allow for M_PI
#include <cmath>            // M_PI & trig functions

#include "Surface.h"
#include "Utils.h"
#include "Cell.h"

using Point = Geometry::Point;
using Polar = Phonon::Polarization;
using Utils::urand;

Surface::Surface(Line surface_line, Cell& cell, double specularity, int norm_sign)
    : surface_line_{std::move(surface_line)},
      cell_{cell},
      normal_{surface_line.normal(norm_sign)},
      specularity_{specularity}
{}

void Surface::redirectPhonon(Phonon& p) const noexcept {
    const auto& [nx, ny] = normal_;
    const auto rand = urand();
    const auto new_dx = sqrt(rand);
    const auto new_dy = sqrt(1.-rand) * cos(2. * M_PI * urand());
    p.setDirection(nx*new_dx-ny*new_dy, ny*new_dx+nx*new_dy);
}

void Surface::boundaryHandlePhonon(Phonon& p) const noexcept {
    if (specularity_ == 1. || urand() < specularity_) { // Reflective Scatter
        const auto& [nx, ny] = normal_;
        const auto& [dx, dy] = p.getDirection();
        const auto new_dx = -dx*nx-dy*ny;
        const auto new_dy = -dx*ny+dy*nx;
        p.setDirection(nx*new_dx-ny*new_dy, ny*new_dx+nx*new_dy);
    } else { // Diffuse Scatters
        redirectPhonon(p);
    }
}

EmitSurface::EmitSurface(Line surface_line, Cell& cell, double specularity, int norm_sign,
                         const Material& mat, double temp, double duration, double start_time)
    : Surface{std::move(surface_line), cell, specularity, norm_sign},
      temp_{temp}, emit_table_{mat.emitTable(temp)}, duration_{duration}, start_time_{start_time} {}

void EmitSurface::handlePhonon(Phonon& p, double step_time) const noexcept {
    const auto phonon_time = p.getLifeStep() * step_time;
    (phonon_time < start_time_ || phonon_time + step_time > start_time_ + duration_) ? Surface::boundaryHandlePhonon(p)
                                                                                     : p.setCell(nullptr);
}

double EmitSurface::getPhononTime() const noexcept {
    return start_time_ + duration_ * Utils::urand();
}

void TransitionSurface::handlePhonon(Phonon& p) const noexcept {
    // Material is the same between sensor areas
    if (p.getCell()->getMaterialID() == cell_.getMaterialID()) {
        p.setCell(&cell_);
    } else { // Phonon is passing from one material to another
        const auto& material = cell_.getMaterial();
        // Find maximum frequency allowable in the new material
        const auto max_freq = [&p, &material](){
            switch(p.getPolar()) {
                case Polar::LA:
                    return material.max_freq_la();
                case Polar::TA:
                    return material.max_freq_ta();
                default:
                    return std::max(material.max_freq_la(), material.max_freq_ta());
            }
        }();
        // TODO: Finish
        // Lazy evaluate whether the phonon is transmitted
        auto isTransmitted = [this](const Phonon& p){
            const auto freq_index = p.getFreqIndex();
            const auto c1 = cell_.getHeatCapacityAtFreq(freq_index);
            const auto c2 = p.getCell()->getHeatCapacityAtFreq(freq_index);
            return false;
        };

        // Transition surfaces cause the phonon to diffusely (spec=0) backscatter when they lack a corresponding
        // state in the new material (frequency too high)
        if (p.getFreq() > max_freq) { // Or transmission fails
            redirectPhonon(p); // Backscatters into the same cell
        } else {
            // TODO: if phonon cell can pass into new material -> work to do
            // Material interface - hard part

            p.setCell(&cell_);
        }
    }
}

