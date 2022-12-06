#define _USE_MATH_DEFINES   // Allow for M_PI
#include <cmath>            // M_PI & trig functions

#include "Phonon.h"

Phonon::Phonon(int sign, double lifetime, Cell *cell)
    : sign_{sign}, lifetime_{lifetime}, cell_{cell} {}

void Phonon::drift(double time) noexcept {
    const auto factor = velocity_ * time;
    px_ += dx_ * factor;
    py_ += dy_ * factor;
}

void Phonon::scatterUpdate(std::size_t freq_index, double freq, double velocity, Phonon::Polarization polar) noexcept {
    freq_index_ = freq_index;
    freq_ = freq;
    velocity_ = velocity;
    polar_ = polar;
}

void Phonon::setRandDirection(double r1, double r2) noexcept {
    dx_ = 2.*r1-1.;
    dy_ = sqrt(1.-dx_*dx_) * cos(2.*M_PI*r2);
}
