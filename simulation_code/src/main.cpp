#include "Timer.h"
#include "InputManager.h"

int main(int argc, char* argv[]) {
    if (argc > 1) {
        try {
            std::vector<std::string> filenames(argv+1, argv + argc);
            for (const auto& filename : filenames) {
                const std::filesystem::path filepath = filename;
                if (auto model = InputManager::deserialize(filepath); model) {
                    auto &m = *model;
                    Timer timer;
                    m.runSimulation();
                    timer.time();
                    m.exportResults(filepath);
                } else {
                    std::cerr << "There was an error reading the data from the file at " << filepath << '\n';
                }
            }
        } catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
        }
    } else {
        std::cout << "Need filenames\n";
    }

    std::cout << "done\n";
}

