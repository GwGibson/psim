#ifndef GEOMETRY_INPUTMANAGER_H
#define GEOMETRY_INPUTMANAGER_H

#include "Model.h"

class InputManager {
public:
    InputManager() = delete;
    /***
    * @param filepath - Filepath to the JSON file containing the simulation settings
    * @return A model object that is prepared to be simulated or nullopt if there is an error generating the model
    */
    [[nodiscard]] static std::optional<Model> deserialize(const std::filesystem::path& filepath);
};


#endif //GEOMETRY_INPUTMANAGER_H
