#include <fstream>
#include <iostream>

#include "InputManager.h"
#include "Json.h"

using json = nlohmann::json;
using Point = Geometry::Point;
using Triangle = Geometry::Triangle;
using SimulationType = Sensor::SimulationType;

std::optional<Model> InputManager::deserialize(const std::filesystem::path& filepath) {

    // ************************ Helper methods **********************************
    auto buildModel = [](std::size_t num_cells, std::size_t num_sensors, const auto& s_data) {
        const bool phasor_sim = s_data.at("phasor_sim").dump() == "true";
        return Model(num_cells, num_sensors, s_data.at("num_measurements"),
                     s_data.at("num_phonons"), s_data.at("sim_time"), s_data.at("t_eq"), phasor_sim);
    };

    auto addMaterial = [](auto& model, double t_eq, const auto& m_data, std::size_t id) {
        const auto& jd_data = m_data.at("d_data");
        const auto& jr_data = m_data.at("r_data");

        auto la_data = static_cast<std::array<double, 3>>(jd_data.at("la_data"));
        auto ta_data = static_cast<std::array<double, 3>>(jd_data.at("ta_data"));
        DispersionData d_data { la_data, ta_data, jd_data.at("max_freq_la"), jd_data.at("max_freq_ta") };
        RelaxationData r_data { jr_data.at("b_l"), jr_data.at("b_tn"), jr_data.at("b_tu"), jr_data.at("b_i"),
                                jr_data.at("w") };
        auto material = Material(id, d_data, r_data);
        if (t_eq == 0.) { material.setFullSimulation(); }
        model.addMaterial(m_data.at("name"), material);
    };

    auto getPoint = [](const auto& p_data) {
        return Point{p_data.at("x"), p_data.at("y")};
    };

    auto addCell = [&](auto& model, const auto& c_data) {
        const auto& t_data = c_data.at("triangle");
        Triangle t{ getPoint(t_data.at("p1")), getPoint(t_data.at("p2")), getPoint(t_data.at("p3")) };
        model.addCell(std::move(t), c_data.at("sensorID"), c_data.at("specularity"));
    };
    // ************************ End Helper methods **********************************

    if (is_regular_file(filepath)) {
        if (std::ifstream file(filepath.string().data()); file.is_open()) {
            json jdata;
            try {
                file >> jdata;

                const auto& settings = jdata.at("settings");
                const auto type = [&settings](){
                    switch (static_cast<int>(settings.at("sim_type"))) {
                        case 1:
                            return SimulationType::Periodic;
                        case 2:
                            return SimulationType::Transient;
                        default:
                            return SimulationType::SteadyState;
                    }
                }();

                const auto& sensors = jdata.at("sensors");
                const auto& cells = jdata.at("cells");
                const auto t_eq = settings.at("t_eq");
                // Generate the base model
                Model model = buildModel(cells.size(), sensors.size(), settings);
                // Set the model simulation type
                (type == SimulationType::SteadyState) ? model.setSimulationType(type)
                                                      : model.setSimulationType(type, settings.at("step_interval"));

                const auto& materials = jdata.at("materials");
                std::size_t material_id = 0;
                // Incorporate material data into the model
                for (const auto& m_data : materials) {
                    addMaterial(model, t_eq, m_data, material_id++);
                }
                // Add sensors to the model
                for (const auto& s_data : sensors) {
                    model.addSensor(s_data.at("id"), s_data.at("material"), s_data.at("t_init"), type);
                }
                // Add cells to the model
                for (const auto& c_data : cells) {
                    addCell(model, c_data);
                }
                const auto& surfaces = jdata.at("emit_surfaces");
                // Add the emitting surfaces to the model
                for (const auto& s_data : surfaces) {
                    // TODO: better exception here
                    if (!model.setEmitSurface(getPoint(s_data.at("p1")), getPoint(s_data.at("p2")),
                                              s_data.at("temp"), s_data.at("duration"), s_data.at("start_time"))) {
                        throw std::runtime_error(std::string("Unable to add emitting surface.\n"));
                    }
                }
                return model;
            } catch (const std::exception& e) {
                std::cerr << e.what() << '\n';
            }
        } else {
            std::cerr << "Error opening file at " << filepath << '\n';
        }
    } else {
        std::cerr << "File at " << filepath << " is not a regular file.\n";
    }
    return std::nullopt;
}