#define _USE_MATH_DEFINES   // Allow for M_PI
#include <cmath>            // M_PI & trig functions

#include "Phonon.h"
#include "Cell.h"

Phonon::Phonon(signed char sign, double lifetime, Cell *cell)
    : sign_{sign}, lifetime_{lifetime}, cell_{cell} {}



void Phonon::scatterUpdate(std::size_t freq_index, double freq, double velocity, Phonon::Polarization polar) noexcept {
    freq_index_ = freq_index;
    freq_ = freq;
    velocity_ = velocity;
    polar_ = polar;
}

void Phonon::drift(double time) noexcept {
    const auto factor = velocity_ * time;
    px_ += dx_ * factor;
    py_ += dy_ * factor;
}

void Phonon::setRandDirection(double r1, double r2) noexcept {
    dx_ = 2.*r1-1.;
    dy_ = sqrt(1.-dx_*dx_) * cos(2.*M_PI*r2);
}

// TODO: Add better/reusable exceptions for the methods below
std::size_t Phonon::getCellSensorID() const {
    if (cell_ == nullptr) {
        throw std::runtime_error(std::string("Cannot get the sensorID if the phonon it is not in a cell.\n"));
    }
    return cell_->getSensorID();
}

std::size_t Phonon::getCellMaterialID() const {
    if (cell_ == nullptr) {
        throw std::runtime_error(std::string("Cannot get the materialID if the phonon it is not in a cell.\n"));
    }
    return cell_->getMaterialID();
}

double Phonon::getCellHeatCapacityAtFreq(std::size_t index) const {
    if (cell_ == nullptr) {
        throw std::runtime_error(std::string("Cannot get the heat capacity if the phonon it is not in a cell.\n"));
    }
    return cell_->getHeatCapacityAtFreq(index);
}

Phonon::RelaxRates Phonon::getRelaxRates(std::size_t step) const {
    if (cell_ == nullptr) {
        throw std::runtime_error(std::string("Cannot get the scattering rates if the phonon it is not in a cell.\n"));
    }
    return cell_->getMaterial().relaxRates(cell_->getSteadyTemp(step), freq_, polar_);
}

std::array<Geometry::Line, 3> Phonon::getCellBoundaryLines() const {
    if (cell_ == nullptr) {
        throw std::runtime_error(std::string("Cannot get the cell boundary lines if the phonon it is not in a cell.\n"));
    }
    return cell_->getBoundaryLines();
}

void Phonon::handleSurfaceCollision(const Geometry::Point& poi, double step_time) {
    if (cell_ == nullptr) {
        throw std::runtime_error(std::string("Cannot handle surface collisions if the phonon it is not in a cell.\n"));
    }
    cell_->handleSurfaceCollision(*this, poi, step_time);
}

void Phonon::scatterUpdate() {
    if (cell_ == nullptr) {
        throw std::runtime_error(std::string("Cannot scatter a phonon if is not in a cell.\n"));
    }
    cell_->scatterUpdate(*this);
}

void Phonon::updateCellHeatParams(std::size_t step) const {
    if (cell_ == nullptr) {
        throw std::runtime_error(std::string("Cannot update cell heat params if the phonon is not in a cell.\n"));
    }
    cell_->updateHeatParams(*this, step);
}

void Phonon::setRandPoint(double r1, double r2) {
    if (cell_ == nullptr) {
        throw std::runtime_error(std::string("Cannot set a phonon to a random location when it is not in a cell.\n"));
    }
    const auto [x, y] = cell_->getRandPoint(r1, r2);
    setPosition(x, y);
}



