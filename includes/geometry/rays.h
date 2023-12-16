#ifndef RAYS
#define RAYS

#include <iostream>
#include "vectors.h"
#include "points.h"

using namespace std;

/**
 * @brief A class for 3D rays.
*/
class Ray {
public:
    Point3f origin;
    Vector3f direction;
    mutable float t_max;
    float time;

    /**
     * @brief Default constructor. Initializes t_max to infinity and time to 0.
     */
    Ray() : t_max(INFINITY), time(0) {}

    /**
     * @brief Constructor with origin, direction, t_max, and time.
     * @param origin The origin of the ray.
     * @param direction The direction of the ray.
     * @param t_max The maximum value of t for which the ray is defined.
     * @param time The time at which the ray exists.
     */
    Ray(const Point3f& origin, const Vector3f& direction, float t_max = INFINITY, float time = 0) :
        origin(origin), direction(direction), t_max(t_max), time(time) {}

    /**
     * @brief Computes the point at a given t value along the ray.
     * @param t The t value.
     * @return The point at t along the ray.
     */
    Point3f operator()(float t);

    /**
     * @brief Outputs the ray to an ostream.
     * @param os The ostream.
     * @param r The ray.
     * @return The ostream.
     */
    friend ostream& operator<<(ostream& os, const Ray& r);
};

/**
 * @brief A class for rays with differentials.
 * @details This class is used for antialiasing.
 * @see Ray
*/
class RayDifferential : public Ray {
public:
    /**
     * @brief Default constructor. Initializes has_differentials to false.
     */
    RayDifferential() : has_differentials(false) {}

    /**
     * @brief Constructor with origin, direction, t_max, and time. Initializes has_differentials to false.
     * @param origin The origin of the ray.
     * @param direction The direction of the ray.
     * @param t_max The maximum value of t for which the ray is defined.
     * @param time The time at which the ray exists.
     */
    RayDifferential(const Point3f& origin, const Vector3f& direction, float t_max = INFINITY, float time = 0) :
        Ray(origin, direction, t_max, time), has_differentials(false) {}

    /**
     * @brief Constructor that initializes a RayDifferential from a Ray. Initializes has_differentials to false.
     * @param ray The Ray to convert.
     */
    RayDifferential(const Ray& ray) : Ray(ray), has_differentials(false) {}

    /**
     * @brief Scales the differentials of the ray.
     * @param s The scale factor.
     */
    void scale_differentials(float s);

    bool has_differentials; ///< Indicates whether the ray has differentials.
    Point3f rx_origin, ry_origin; ///< The x and y origins of the differentials.
    Vector3f rx_direction, ry_direction; ///< The x and y directions of the differentials.
};

#endif // !RAYS