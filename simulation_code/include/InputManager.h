#ifndef GEOMETRY_INPUTMANAGER_H
#define GEOMETRY_INPUTMANAGER_H

#include "Model.h"

class InputManager {
public:
    InputManager() = delete;
    [[nodiscard]] static std::optional<Model> deserialize(const std::filesystem::path& filepath);
};


#endif //GEOMETRY_INPUTMANAGER_H
