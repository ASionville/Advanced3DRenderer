#include <iostream>
#include "vectors.h"
#include "points.h"
#include "rays.h"
#include "transforms.h"
#include "bboxes.h"
#include "animated_transforms.h"
#include "../utils/maths.h"

using namespace std;

/**
 * Calculates the dot product between two quaternions.
 *
 * @param q1 The first quaternion.
 * @param q2 The second quaternion.
 * @return The dot product of the two quaternions.
 */
inline float dot(const Quaternion& q1, const Quaternion& q2) {
    return dot(q1.v, q2.v) + q1.w * q2.w;
}

/**
 * @brief Normalizes a quaternion.
 * 
 * @param q The quaternion to be normalized.
 * @return The normalized quaternion.
 */
inline Quaternion normalize(const Quaternion& q) {
    return q / sqrt(dot(q, q));
}

/**
 * Converts a Quaternion to a Transform.
 *
 * @param q The Quaternion to convert.
 * @return The resulting Transform.
 */
Transform quaternion_to_transform(const Quaternion& q) {
    float x = q.v.x;
    float y = q.v.y;
    float z = q.v.z;
    float w = q.w;

    float x2 = x * x;
    float y2 = y * y;
    float z2 = z * z;
    float w2 = w * w;

    float xy = x * y;
    float xz = x * z;
    float xw = x * w;
    float yz = y * z;
    float yw = y * w;
    float zw = z * w;

    Matrix4x4 m = Matrix4x4();
    m.m[0][0] = 1 - 2 * (y2 + z2);
    m.m[0][1] = 2 * (xy + zw);
    m.m[0][2] = 2 * (xz - yw);

    m.m[1][0] = 2 * (xy - zw);
    m.m[1][1] = 1 - 2 * (x2 + z2);
    m.m[1][2] = 2 * (yz + xw);

    m.m[2][0] = 2 * (xz + yw);
    m.m[2][1] = 2 * (yz - xw);
    m.m[2][2] = 1 - 2 * (x2 + y2);

    return Transform(m, transpose(m));
}

/**
 * Performs spherical linear interpolation between two quaternions.
 *
 * @param q1 The starting quaternion.
 * @param q2 The ending quaternion.
 * @param t The interpolation factor (between 0 and 1).
 * @return The interpolated quaternion.
 */
Quaternion slerp(const Quaternion& q1, const Quaternion& q2, float t) {
    float cos_theta = dot(q1, q2);
    if (cos_theta > 0.9995f) {
        return normalize(q1 * (1 - t)  + q2 * t);
    } else {
        float theta = acos(clamp(cos_theta, -1, 1));
        float theta_prime = theta * t;
        Quaternion q_perp = normalize(q2 - q1 * cos_theta);
        return q1 * cos(theta_prime) + q_perp * sin(theta_prime);
    }
}

/**
 * Transforms a 3D point using the specified transformation matrix.
 *
 * @param point The 3D point to be transformed.
 * @param t The transformation matrix.
 * @return The transformed 3D point.
 */
Point3f transform(const Point3f& point, const Transform& t) {
    float x = point.x;
    float y = point.y;
    float z = point.z;

    float xp = t.matrix.m[0][0] * x + t.matrix.m[0][1] * y + t.matrix.m[0][2] * z + t.matrix.m[0][3];
    float yp = t.matrix.m[1][0] * x + t.matrix.m[1][1] * y + t.matrix.m[1][2] * z + t.matrix.m[1][3];
    float zp = t.matrix.m[2][0] * x + t.matrix.m[2][1] * y + t.matrix.m[2][2] * z + t.matrix.m[2][3];
    float wp = t.matrix.m[3][0] * x + t.matrix.m[3][1] * y + t.matrix.m[3][2] * z + t.matrix.m[3][3];

    if (wp == 1) {
        return Point3f(xp, yp, zp);
    } else {
        return Point3f(xp, yp, zp) / wp;
    }
}

/**
 * Transforms a 3D vector using the specified transformation.
 *
 * @param vector The vector to be transformed.
 * @param t The transformation to apply.
 * @return The transformed vector.
 */
Vector3f transform(const Vector3f& vector, const Transform& t) {
    float x = vector.x;
    float y = vector.y;
    float z = vector.z;

    return Vector3f(t.matrix.m[0][0] * x + t.matrix.m[0][1] * y + t.matrix.m[0][2] * z,
                    t.matrix.m[1][0] * x + t.matrix.m[1][1] * y + t.matrix.m[1][2] * z,
                    t.matrix.m[2][0] * x + t.matrix.m[2][1] * y + t.matrix.m[2][2] * z);
}

/**
 * Transforms a ray using the given transformation matrix.
 *
 * @param ray The ray to be transformed.
 * @param t The transformation matrix.
 * @return The transformed ray.
 */
Ray transform(const Ray& ray, const Transform& t) {
    Point3f new_origin = transform(ray.origin, t);
    Vector3f new_direction = transform(ray.direction, t);
    return Ray(new_origin, new_direction, ray.t_max, ray.time);
}

/**
 * Transforms a RayDifferential object using the specified Transform.
 *
 * @param ray The RayDifferential object to be transformed.
 * @param t The Transform to apply to the ray.
 * @return The transformed RayDifferential object.
 */
RayDifferential transform(const RayDifferential& ray, const Transform& t) {
    Ray tr = transform(Ray(ray), t);
    RayDifferential ret(tr.origin, tr.direction, tr.t_max, tr.time);
    ret.has_differentials = ray.has_differentials;
    ret.rx_origin = transform(ray.rx_origin, t);
    ret.ry_origin = transform(ray.ry_origin, t);
    ret.rx_direction = transform(ray.rx_direction, t);
    ret.ry_direction = transform(ray.ry_direction, t);
    return ret;
}

/**
 * Transforms a bounding box using a given transformation matrix.
 *
 * @param b The original bounding box.
 * @param t The transformation matrix.
 * @return The transformed bounding box.
 */
BBox3f transform(const BBox3f& b, const Transform& t) {
    BBox3f ret(transform(Point3f(b.p_min.x, b.p_min.y, b.p_min.z), t));
    ret = union_bbox(ret, transform(Point3f(b.p_max.x, b.p_min.y, b.p_min.z), t));
    ret = union_bbox(ret, transform(Point3f(b.p_min.x, b.p_max.y, b.p_min.z), t));
    ret = union_bbox(ret, transform(Point3f(b.p_min.x, b.p_min.y, b.p_max.z), t));
    ret = union_bbox(ret, transform(Point3f(b.p_min.x, b.p_max.y, b.p_max.z), t));
    ret = union_bbox(ret, transform(Point3f(b.p_max.x, b.p_max.y, b.p_min.z), t));
    ret = union_bbox(ret, transform(Point3f(b.p_max.x, b.p_min.y, b.p_max.z), t));
    ret = union_bbox(ret, transform(Point3f(b.p_max.x, b.p_max.y, b.p_max.z), t));
    return ret;
}

BBox3f motion_bounds(AnimatedTransform& t, const BBox3f& b) {
    // Copy the start and end transforms
    Transform t_start = *t.start_transform;
    Transform t_end = *t.end_transform;

    if (!t.actually_animated) {
        return transform(b, t_start);
    }
    if (t.has_rotation == false) {
        return union_bbox(transform(b, t_start), transform(b, t_end));
    }

    // Compute the bounding box of the object's motion
    BBox3f bounds;
    for (int corner = 0; corner < 8; corner++) {
        bounds = union_bbox(bounds, transform(b.corner(corner), t_start));
        bounds = union_bbox(bounds, transform(b.corner(corner), t_end));
    }
