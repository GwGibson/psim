#ifndef GEOMETRY_GEOMETRY_H
#define GEOMETRY_GEOMETRY_H

#include <optional>
#include <ostream>
#include <array>

namespace Geometry {

    struct Vector2D {
        constexpr Vector2D(double x, double y) noexcept: x{x}, y{y} {}
        double x;
        double y;
    };

    struct Point {
        constexpr Point(double x, double y) noexcept : x{x}, y{y} {}
        double x;
        double y;

        bool operator==(const Point& rhs) const;
    };
    std::ostream& operator<<(std::ostream& os, const Point& point);
    Point operator-(const Point& lhs, const Point& rhs);

    using PointPair = std::pair<Point, Point>;

    struct Line {
        Line(Point p1, Point p2);
        const Point p1;
        const Point p2;

        // maybe only store the length and bounding box
        const double slope;
        const double intercept;
        const PointPair boundingBox;
        const double length;

        [[nodiscard]] PointPair getPoints() const noexcept { return {p1, p2}; }

        [[nodiscard]] bool overlaps(const Line& other) const noexcept;
        [[nodiscard]] bool contains(const Line& other) const noexcept;
        [[nodiscard]] bool contains(const Point& point) const noexcept;
        [[nodiscard]] bool intersects(const Line& other) const noexcept;
        // norm_sign controls the direction of the normal vector. It can be used
        // to ensure the normal vector always points into a surface.
        // Should be 1 if the polygon is constructed in a clockwise fashion and -1 otherwise
        [[nodiscard]] Vector2D normal(int norm_sign=1) const noexcept;
        // Returns the point of intersection if this line intersects with the input line.
        // Return std::nullopt if the lines do not intersect
        [[nodiscard]] std::optional<Point> getIntersection(const Line& other) const noexcept;
        [[nodiscard]] Point getRandPoint(double r1) const noexcept;

        bool operator==(const Line& rhs) const { return p1 == rhs.p1 && p2 == rhs.p2; }
        bool operator>(const Line& rhs) const { return length > rhs.length; }
    };
    std::ostream& operator<<(std::ostream& os, const Line& line);

    struct Triangle {
        Triangle(Point p1, Point p2, Point p3);
        const Point p1;
        const Point p2;
        const Point p3;

        [[nodiscard]] std::array<Line, 3> lines() const noexcept;

        [[nodiscard]] double area() const noexcept;
        [[nodiscard]] bool intersects(const Triangle& other) const noexcept;
        // Returns true if any point of the incoming triangle is contained in this triangle
        [[nodiscard]] bool contains(const Triangle& other) const noexcept;
        // Return false if the point is on an edge of the triangle
        [[nodiscard]] bool contains(const Point& p) const noexcept;
        [[nodiscard]] bool isClockwise() const noexcept;
        [[nodiscard]] Point getRandPoint(double r1, double r2) const noexcept;

        bool operator==(const Triangle& rhs) const;
    };
    std::ostream& operator<<(std::ostream& os, const Triangle& triangle);
}

#include <string>
class ShapeError : public std::exception {
public:
    [[nodiscard]] const char* what() const noexcept override { return message_.c_str(); }
protected:
    void setMessage(std::string_view message) { message_ = message; }
private:
    std::string message_;
};

class LineError : public ShapeError {
public:
    explicit LineError(Geometry::Line line);
private:
    const Geometry::Line line_;
};

class TriangleError : public ShapeError {
public:
    explicit TriangleError(Geometry::Triangle triangle);
private:
    const Geometry::Triangle triangle_;
};

#endif //GEOMETRY_GEOMETRY_H
