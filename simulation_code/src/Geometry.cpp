#include <cmath>
#include <algorithm> // minmax
#include <array>

#include "Geometry.h"

using Vector2D = Geometry::Vector2D;
using Point = Geometry::Point;
using Line = Geometry::Line;
using Triangle = Geometry::Triangle;
using PointPair = Geometry::PointPair;

namespace { constexpr double GEOEPS {10e-9}; }

[[nodiscard]] static PointPair findBoundingBox(const Line& line) noexcept;
[[nodiscard]] static bool boxesIntersect(const PointPair& lhs, const PointPair& rhs) noexcept;
[[nodiscard]] static bool isPointOnLine(const Line& line, const Point& p) noexcept;
[[nodiscard]] static bool isPointRightOfLine(const Line& line, const Point& p) noexcept;
[[nodiscard]] static bool doesSegmentCrossLine(const Line& l1, const Line& l2) noexcept;
[[nodiscard]] static double crossProduct(const Point& p1, const Point& p2) noexcept;
[[nodiscard]] static double dotProduct(const Point& p1, const Point& p2) noexcept;
[[nodiscard]] static double getLineLength(const Line& line) noexcept;
[[nodiscard]] static double getSlope(const Point& p1, const Point& p2) noexcept;
[[nodiscard]] static double getIntercept(double x, double y, double slope) noexcept;

std::ostream& Geometry::operator<<(std::ostream& os, const Point& point) {
    os << "Point (" << point.x << ", " << point.y << ')';
    return os;
}

Point Geometry::operator-(const Point& lhs, const Point& rhs) {
    return {lhs.x - rhs.x, lhs.y - rhs.y};
}

bool Geometry::Point::operator==(const Point& rhs) const {
    return (x <= rhs.x + GEOEPS && x >= rhs.x - GEOEPS) &&
           (y <= rhs.y + GEOEPS && y >= rhs.y - GEOEPS);
}

bool Point::operator!=(const Point &rhs) const {
    return !(rhs == *this);
}

Line::Line(Point p1, Point p2)
    : p1{p1}, p2{p2},
      slope{getSlope(p1, p2)}, intercept{getIntercept(p1.x, p1.y, slope)},
      boundingBox{findBoundingBox(*this)},
      length{getLineLength(*this)}
{
    if (p1 == p2) {
        throw LineError(*this);
    }
}

// Assumes the lines are parallel. Returns true if the overlap at more than a single point.
bool Line::overlaps(const Line& other) const noexcept {
    if (*this == other) {
        return true;
    }

    auto minmax_point = [](const auto& pt1, const auto& pt2){
        return (pt1.x + pt1.y - (pt2.x + pt2.y) >= 0.) ? std::pair{pt2, pt1} : std::pair{pt1, pt2};
    };

    bool overlaps = false;
    const auto& [pa, pb] = minmax_point(p1, p2);
    const auto& [pc, pd] = minmax_point(other.p1, other.p2);
    if (!(pd == pa || pb == pc)) { // If the lines do not overlap at a single point
        if ((pd.x + pd.y - (pa.x + pa.y) >= 0.) && (pb.x + pb.y - (pc.x + pc.y) >= 0.)) { // Use GEOEPS here?
            overlaps = true;
        }
    }
    return overlaps;
}

// True if this line segment contains the other line segment (<=) - other line segment fully contained in this one
bool Line::contains(const Line& other) const noexcept {
    const auto& [p3, p4] = other.getPoints();
    if (isPointOnLine(*this, p3) && isPointOnLine(*this, p4) && length >= other.length) {
        const auto& [l1_x_min, l1_x_max] = std::minmax(p1.x, p2.x);
        const auto& [l2_x_min, l2_x_max] = std::minmax(p3.x, p4.x);
        const auto& [l1_y_min, l1_y_max] = std::minmax(p1.y, p2.y);
        const auto& [l2_y_min, l2_y_max] = std::minmax(p3.y, p4.y);
        return (l1_x_max >= l2_x_max - GEOEPS &&
                l1_x_min <= l2_x_min + GEOEPS &&
                l1_y_max >= l2_y_max - GEOEPS &&
                l1_y_min <= l2_y_min + GEOEPS);
    }
    return false;
}

bool Line::contains(const Point& point) const noexcept {
    return isPointOnLine(*this, point);
}

bool Line::intersects(const Line& other) const noexcept {
    return boxesIntersect(boundingBox, other.boundingBox) &&
           doesSegmentCrossLine(*this, other) &&
           doesSegmentCrossLine(other, *this);
}

Vector2D Line::normal(int norm_sign) const noexcept {
    (norm_sign >= 0) ? norm_sign = 1 : norm_sign = -1;
    return {norm_sign*(p2.y - p1.y)/length, -norm_sign*(p2.x - p1.x)/length};
}

// If a line is parallel and overlaps another line, this acts as though there is no
// intersection point
std::optional<Point> Line::getIntersection(const Line& other) const noexcept {
    std::optional<Point> poi = std::nullopt;

    auto equals = [](double a, double b) {
        return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * GEOEPS);
    };

    auto handle_vert = [&](const auto& pt1, const auto& pt3, const auto& pt4) {
        if (equals(pt3.x,pt4.x)) { // both lines vertical
            poi = std::nullopt;
        } else { // line from pt3 to pt4 is not vertical
            const auto m = getSlope(pt3, pt4);
            const auto b = getIntercept(pt3.x, pt3.y, m);
            poi = {pt1.x, m*pt1.x + b};
        }
    };
    if (intersects(other)) { // If the line segments intersect, find the poi
        const auto& [t_p1, t_p2] = getPoints();
        const auto& [o_p3, o_p4] = other.getPoints();
        if (equals(t_p1.x, t_p2.x)) { // l1 is a vertical line - needs special consideration
            handle_vert(t_p1, o_p3, o_p4);
        } else if (equals(o_p3.x, o_p4.x)) { // l2 is a vertical line
            handle_vert(o_p3, t_p1, t_p2);
        } else {
            if (equals(slope, other.slope)) { // Lines are parallel (by convention, no intersection - will this be a problem?)
                poi = std::nullopt;
            } else { // Only 1 point of intersection -> should be the most common case
                const double x = (other.intercept - intercept)/(slope-other.slope);
                poi = {x, slope*x+intercept};
            }
        }
    }
    return poi;
}

Point Line::getRandPoint(double r1) const noexcept {
    const double r2 = 1. - r1;
    return {p1.x*r1+p2.x*r2, p1.y*r1+p2.y*r2};
}

std::ostream& Geometry::operator<<(std::ostream& os, const Line& line) {
    os << "Line Segment: [" << line.p1 << ", " << line.p2 << ']';
    return os;
}


Triangle::Triangle(Point p1, Point p2, Point p3)
    : p1{p1}, p2{p2}, p3{p3}
{
    if  ( (fabs(getSlope(p1, p2)) == fabs(getSlope(p2, p3))  && fabs(getSlope(p2, p3)) == fabs(getSlope(p3, p1)))
        || (p1 == p2) || (p2 == p3) || (p3 == p1) ) {
        throw TriangleError(*this);
    }
}

[[nodiscard]] std::array<Line, 3> Triangle::lines() const noexcept {
    return std::array{Line(p1, p2), Line(p2, p3), Line(p3,p1)};
}

// If this triangle intersects the other triangle. Intersections at end points
// or along parallel lines are not counted. This is used to build the model and ensure
// cell placement is valid.
bool Geometry::Triangle::intersects(const Triangle &other) const noexcept {
    const auto intersect = [](const Line& line1, const Line& line2){
        auto intersects = false;
        const auto poi = line1.getIntersection(line2); // This function returns std::nullopt for overlapping || lines
        if (poi) {
            // Don't count intersections at the end points
            if (!(*poi == line1.p1 || *poi == line1.p2 || *poi == line2.p1 || *poi == line2.p2)) {
                intersects = true;
            }
        }
        return intersects;
    };
    const auto& this_lines = lines();
    const auto& other_lines = other.lines();
    // Check if any of the sides intersect
    for (const auto& line1 : this_lines) {
        if (std::any_of(std::cbegin(other_lines), std::cend(other_lines),
                        [&](const auto& line2){ return intersect(line1, line2); })) {
            return true;
        }
    }
    return false;
}

bool Triangle::contains(const Point& p) const noexcept {
    auto getArea = [](auto pt1, auto pt2, auto pt3) {
        try {
            auto t = Triangle(pt1,pt2,pt3);
            return t.area();
        } catch (const TriangleError& e) {
            return 0.;
        }
    };

    // bp1, bp2, d00, d01, d11 & denom can all be cached if necessary
    const auto bp1 = p2-p1;
    const auto bp2 = p3-p1;
    const auto d00 = dotProduct(bp1, bp1);
    const auto d01 = dotProduct(bp1, bp2);
    const auto d11 = dotProduct(bp2, bp2);
    const auto denom = d00 * d11 - d01 * d01;
    // These values cannot be cached
    const auto bp3 = p - p1;
    const auto d20 = dotProduct(bp3, bp1);
    const auto d21 = dotProduct(bp3, bp2);
    const auto u = (d11 * d20 - d01 * d21) / denom;
    const auto v = (d00 * d21 - d01 * d20) / denom;
    if (u >= GEOEPS && u <= 1. + GEOEPS &&
    v >= GEOEPS && v <= 1. + GEOEPS &&
    u + v <= 1. + GEOEPS) {
        if ((area() - getArea(p, p2, p3) + getArea(p1, p, p3) + getArea(p1, p2, p)) < GEOEPS) {
            return true;
        }

    }
    return false;
}

bool Geometry::Triangle::isClockwise() const noexcept {
    const auto& [p1_x, p1_y] = std::pair{p2.x - p1.x, p2.y + p1.y};
    const auto& [p2_x, p2_y] = std::pair{p3.x - p2.x, p3.y + p2.y};
    const auto& [p3_x, p3_y] = std::pair{p1.x - p3.x, p1.y + p3.y};
    return p1_x * p1_y + p2_x * p2_y + p3_x * p3_y >= 0.;
}
#if 0
Point Triangle::getRandPoint(double r1, double r2) const noexcept {
    const double root_r1 = sqrt(r1);
    const auto f1 = 1.-root_r1;
    const auto f2 = root_r1 * (1.-r2);
    const auto f3 = root_r1 * r2;
    const auto x = f1*p1.x + f2*p2.x + f3*p3.x;
    const auto y = f1*p1.y + f2*p2.y + f3*p3.y;
    return {x,y};
}
#else
Point Triangle::getRandPoint(double r1, double r2) const noexcept {
    if (r1 + r2 > 1.) {
        r1 = 1 - r1;
        r2 = 1 - r2;
    }
    const auto x = p1.x + (p2.x - p1.x) * r1 + (p3.x - p1.x) * r2;
    const auto y = p1.y + (p2.y - p1.y) * r1 + (p3.y - p1.y) * r2;
    return {x,y};
}
#endif
bool Triangle::operator==(const Triangle& rhs) const {
    auto contains = [](const auto& t, const auto& p) {
        return (p == t.p1 || p == t.p2 || p == t.p3);
    };

    return contains(rhs, p1) &&
           contains(rhs, p2) &&
           contains(rhs, p3);
}

double Triangle::area() const noexcept {
    const auto& [l1, l2, l3] = lines();
    const auto a = l1.length;
    const auto b = l2.length;
    const auto c = l3.length;
    const auto p = (a+b+c)/2.;
    return sqrt(p*(p-a)*(p-b)*(p-c));
}

// TODO: This misses some cases -> Example: [T1] [(0,0), (10,0), (0, 10)] & [T2] [(0,1), (10,0), (0,10)]
// TODO: Catch cases where all 3 points are on the edges of existing triangle (2 is ok - transition surface)
// Returns true if the other triangle is at all (partially or fully) contained within this triangle
bool Triangle::contains(const Triangle& other) const noexcept {
    return contains(other.p1) || contains(other.p2) || contains(other.p3);
}

std::ostream& Geometry::operator<<(std::ostream& os, const Triangle& triangle) {
    os << "Triangle: [" << triangle.p1 << ", " << triangle.p2 << ", " << triangle.p3 << "]";
    return os;
}

PointPair findBoundingBox(const Line& line) noexcept {
    const auto& [p1, p2] = line.getPoints();
    // Bottom left point of bounding box - assume it is p1
    auto bl_x = p1.x;
    auto bl_y = p1.y;
    // Top right point of bounding box - assume it is p2
    auto tr_x = p2.x;
    auto tr_y = p2.y;
    if (p2.x <= p1.x) {
        std::swap(bl_x, tr_x);
    }
    if (p2.y <= p1.y) {
        std::swap(bl_y, tr_y);
    }
    return { {bl_x, bl_y}, {tr_x, tr_y} };
}

bool boxesIntersect(const PointPair& lhs, const PointPair& rhs) noexcept {
    const auto& [bl_1, tr_1] = lhs;
    const auto& [bl_2, tr_2] = rhs;
    // If all these conditions are fulfilled, the boxes intersect, otherwise they do not
    return bl_1.x <= tr_2.x - GEOEPS && // first box bottom left x <= second box top right x
           tr_1.x >= bl_2.x + GEOEPS && // first box top right x >= second box bottom left x
           bl_1.y <= tr_2.y - GEOEPS && // first box bottom left y <= second box top right y
           tr_1.y >= bl_2.y + GEOEPS;   // first box top right y >= second box bottom left y
}

// Returns if the point is on the line (not just the line segment)
bool isPointOnLine(const Line& line, const Point& p) noexcept {
    const auto& [p1, p2] = line.getPoints();
    return fabs(crossProduct( {p2.x - p1.x, p2.y - p1.y}, {p.x - p1.x, p.y - p1.y} )) < GEOEPS;
}

bool isPointRightOfLine(const Line &line, const Point& p) noexcept {
    const auto& [p1, p2] = line.getPoints();
    return crossProduct( {p2.x - p1.x, p2.y - p1.y}, {p.x - p1.x, p.y - p1.y} ) < 0.;
}

bool doesSegmentCrossLine(const Line& l1, const Line& l2) noexcept {
    return isPointOnLine(l1, l2.p1) ||
           isPointOnLine(l1, l2.p2) ||
           (isPointRightOfLine(l1, l2.p1) ^ isPointRightOfLine(l1, l2.p2));
}

double crossProduct(const Point& p1, const Point& p2) noexcept {
    return p1.x * p2.y - p2.x * p1.y;
}

double dotProduct(const Point& p1, const Point& p2) noexcept {
    return p1.x * p2.x + p1.y * p2.y;
}

double getLineLength(const Line& line) noexcept {
    const auto p1 = line.p1;
    const auto p2 = line.p2;
    return sqrt((p2.x-p1.x)*(p2.x-p1.x)+(p2.y-p1.y)*(p2.y-p1.y));
}

// Returns 0 for a vertical line!
double getSlope(const Point& p1, const Point& p2) noexcept {
    return (p1.x == p2.x) ? 0. : (p1.y - p2.y) / (p1.x - p2.x);
}

double getIntercept(double x, double y, double slope) noexcept {
    return y - slope * x;
}

#include <sstream>

LineError::LineError(Geometry::Line line) : line_{std::move(line)} {
    std::ostringstream os;
    os << line_.p1 << ' ' << line_.p2;
    setMessage("Cannot create a line using 2 identical points -> " + os.str());
}

TriangleError::TriangleError(Geometry::Triangle triangle) : triangle_{std::move(triangle)} {
    std::ostringstream os;
    os << triangle_;
    setMessage("These 3 points do not allow for a valid triangle -> " + os.str());
}
