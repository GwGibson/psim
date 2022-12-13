#ifndef GEOMETRY_UTILS_H
#define GEOMETRY_UTILS_H

#include <array>
#include <random>

namespace Utils {
    // Generates a random number from a uniform distribution over [0,1].
    inline double urand() noexcept {
        static std::random_device rd;
        static thread_local std::mt19937 generator(10);
        std::uniform_real_distribution<double> dist(0., std::nextafter(1.0, 2.0));
        return dist(generator);
    }

    // Allows usage of enum classes as integral values similar to unscoped enums. Enum classes seems to
    // require a static_cast<std::size_t> so this save having to type that everytime.
    // Can use eType(enum) instead of std::static_cast<std::size_t>(enum).
    // Must be constexpr to work with std::get when using variants
    template<typename E>
    constexpr auto toInteger(E enumerator) noexcept { return static_cast<std::underlying_type_t<E>>(enumerator); }

    // Custom to function to emulate python zip functionality
    // General implementation
    template<typename In1, typename In2, typename Out>
    void zip(In1 begin1, In2 end1, In2 begin2, In2 end2, Out result) {
        auto it1 {begin1}; auto it2 {begin2};
        while (it1 != end1 && it2 != end2) {
            result++ = {*it1++, *it2++};
        }
    }

    // Specifically for vectors
    template <typename T, typename U>
    std::vector<std::pair<T, U>> zip (const std::vector<T>& r1, const std::vector<U> r2) {
        std::vector<std::pair<T,U>> result;
        zip(std::cbegin(r1), std::cend(r1), std::cbegin(r2), std::cend(r2), std::back_inserter(result));
        return result;
    }

    // Specifically for arrays
    template <typename T, typename U, std::size_t V>
    std::array<std::pair<T, U>, V> zip (const std::array<T, V>& r1, const std::array<U, V>& r2) {
        std::array<std::pair<T,U>, V> result;
        std::size_t index = 0;
        auto it1 {std::cbegin(r1)}; auto it2 {std::cbegin(r2)};
        while (it1 != std::cend(r1) && it2 != std::cend(r2)) {
            result[index++] = {*it1++, *it2++};
        }
        return result;
    }
}

#endif //GEOMETRY_UTILS_H
