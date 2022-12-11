#ifndef GEOMETRY_TIMER_H
#define GEOMETRY_TIMER_H

#include <iostream>
#include <string>
#include <chrono>
#include <functional>

class Timer {
public:
    Timer() : last(std::chrono::steady_clock::now()) {}

    double get_time_diff() {
        auto now = std::chrono::steady_clock::now();
        std::chrono::duration<double> diff = now - last;
        last = std::chrono::steady_clock::now();
        return diff.count();
    }

    void time (const std::string& msg = "Time Taken: ") {
        auto now = std::chrono::steady_clock::now();
        std::chrono::duration<double> diff = now - last;
        std::cout << msg << diff.count() << "[s]\n";
        last = std::chrono::steady_clock::now();
    }
private:
    std::chrono::steady_clock::time_point last;
};

template <typename Time = std::chrono::microseconds, typename Clock = std::chrono::high_resolution_clock>
struct FuncTimer {
    template<typename F, typename... Args>
    static Time duration(F&& f, Args... args) {
        auto start = Clock::now();
        std::invoke(std::forward<F>(f), std::forward<Args>(args)...);
        auto end = Clock::now();
        return std::chrono::duration_cast<Time>(end - start);
    }
};
#endif //GEOMETRY_TIMER_H
