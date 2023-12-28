#ifndef ANIMATED_TRANSFORMS
#define ANIMATED_TRANSFORMS

#include <iostream>
#include "vectors.h"
#include "points.h"
#include "rays.h"
#include "transforms.h"
#include "bboxes.h"
#include "../utils/maths.h"

using namespace std;

/**
 * @brief Represents a quaternion in 3D space.
 */
class Quaternion {
    public :
        Vector3f v;
        float w;

        /**
         * @brief Default constructor.
         */
        Quaternion() {
            v = Vector3f();
            w = 1;
        }

        /**
         * @brief Constructor.
         * 
         * @param matrix The matrix to be converted to a quaternion.
         */
        Quaternion(const Matrix4x4& matrix) {
            float trace = matrix.m[0][0] + matrix.m[1][1] + matrix.m[2][2];
            if (trace > 0) {
                float s = 0.5f / sqrt(trace + 1.0f);
                w = 0.25f / s;
                v.x = (matrix.m[2][1] - matrix.m[1][2]) * s;
                v.y = (matrix.m[0][2] - matrix.m[2][0]) * s;
                v.z = (matrix.m[1][0] - matrix.m[0][1]) * s;
            } else {
                if (matrix.m[0][0] > matrix.m[1][1] && matrix.m[0][0] > matrix.m[2][2]) {
                    float s = 2.0f * sqrt(1.0f + matrix.m[0][0] - matrix.m[1][1] - matrix.m[2][2]);
                    w = (matrix.m[2][1] - matrix.m[1][2]) / s;
                    v.x = 0.25f * s;
                    v.y = (matrix.m[0][1] + matrix.m[1][0]) / s;
                    v.z = (matrix.m[0][2] + matrix.m[2][0]) / s;
                } else if (matrix.m[1][1] > matrix.m[2][2]) {
                    float s = 2.0f * sqrt(1.0f + matrix.m[1][1] - matrix.m[0][0] - matrix.m[2][2]);
                    w = (matrix.m[0][2] - matrix.m[2][0]) / s;
                    v.x = (matrix.m[0][1] + matrix.m[1][0]) / s;
                    v.y = 0.25f * s;
                    v.z = (matrix.m[1][2] + matrix.m[2][1]) / s;
                } else {
                    float s = 2.0f * sqrt(1.0f + matrix.m[2][2] - matrix.m[0][0] - matrix.m[1][1]);
                    w = (matrix.m[1][0] - matrix.m[0][1]) / s;
                    v.x = (matrix.m[0][2] + matrix.m[2][0]) / s;
                    v.y = (matrix.m[1][2] + matrix.m[2][1]) / s;
                    v.z = 0.25f * s;

                }
            }

            float length = sqrt(dot(v, v) + w * w);
            v /= length;
            w /= length;

            this->v = v;
            this->w = w;
        }

        /**
         * @brief Constructor.
         *
         * @param v The vector part of the quaternion.
         * @param w The scalar part of the quaternion.
         */
        Quaternion(const Vector3f& v, float w) {
            this->v = v;
            this->w = w;
        }

        /**
         * @brief Constructor.
         *
         * @param v The vector part of the quaternion.
         */
        Quaternion(const Vector3f& v) {
            this->v = v;
            this->w = 1;
        }

        /**
         * @brief Constructor.
         * 
         * @param transform The transform to be converted to a quaternion.
         */
        Quaternion(const Transform& transform) {
            float trace = transform.matrix.m[0][0] + transform.matrix.m[1][1] + transform.matrix.m[2][2];
            if (trace > 0) {
                float s = 0.5f / sqrt(trace + 1.0f);
                w = 0.25f / s;
                v.x = (transform.matrix.m[2][1] - transform.matrix.m[1][2]) * s;
                v.y = (transform.matrix.m[0][2] - transform.matrix.m[2][0]) * s;
                v.z = (transform.matrix.m[1][0] - transform.matrix.m[0][1]) * s;
            } else {
                if (transform.matrix.m[0][0] > transform.matrix.m[1][1] && transform.matrix.m[0][0] > transform.matrix.m[2][2]) {
                    float s = 2.0f * sqrt(1.0f + transform.matrix.m[0][0] - transform.matrix.m[1][1] - transform.matrix.m[2][2]);
                    w = (transform.matrix.m[2][1] - transform.matrix.m[1][2]) / s;
                    v.x = 0.25f * s;
                    v.y = (transform.matrix.m[0][1] + transform.matrix.m[1][0]) / s;
                    v.z = (transform.matrix.m[0][2] + transform.matrix.m[2][0]) / s;
                } else if (transform.matrix.m[1][1] > transform.matrix.m[2][2]) {
                    float s = 2.0f * sqrt(1.0f + transform.matrix.m[1][1] - transform.matrix.m[0][0] - transform.matrix.m[2][2]);
                    w = (transform.matrix.m[0][2] - transform.matrix.m[2][0]) / s;
                    v.x = (transform.matrix.m[0][1] + transform.matrix.m[1][0]) / s;
                    v.y = 0.25f * s;
                    v.z = (transform.matrix.m[1][2] + transform.matrix.m[2][1]) / s;
                } else {
                    float s = 2.0f * sqrt(1.0f + transform.matrix.m[2][2] - transform.matrix.m[0][0] - transform.matrix.m[1][1]);
                    w = (transform.matrix.m[1][0] - transform.matrix.m[0][1]) / s;
                    v.x = (transform.matrix.m[0][2] + transform.matrix.m[2][0]) / s;
                    v.y = (transform.matrix.m[1][2] + transform.matrix.m[2][1]) / s;
                    v.z = 0.25f * s;

                }
            }

            float length = sqrt(dot(v, v) + w * w);
            v /= length;
            w /= length;

            this->v = v;
            this->w = w;
        }

        /**
         * @brief Addition operator for Quaternions.
         * 
         * @param q The Quaternion to be added.
         * @return The result of the addition.
         */
        Quaternion operator+(const Quaternion& q) const {
            return Quaternion(v + q.v, w + q.w);
        }

        /**
         * @brief In-place addition operator for Quaternions.
         * 
         * @param q The Quaternion to be added.
         * @return The result of the addition.
         */
        Quaternion& operator+=(const Quaternion& q) {
            v += q.v;
            w += q.w;
            return *this;
        }

        /**
         * @brief Subtraction operator for Quaternions.
         * 
         * @param q The Quaternion to be subtracted.
         * @return The result of the subtraction.
         */
        Quaternion operator-(const Quaternion& q) const {
            return Quaternion(v - q.v, w - q.w);
        }

        /**
         * @brief In-place subtraction operator for Quaternions.
         * 
         * @param q The Quaternion to be subtracted.
         * @return The result of the subtraction.
         */
        Quaternion& operator-=(const Quaternion& q) {
            v -= q.v;
            w -= q.w;
            return *this;
        }

        /**
         * @brief Negation operator for Quaternions.
         * 
         * @return The negated Quaternion.
         */
        Quaternion operator-() const {
            return Quaternion(-v, -w);
        }

        /**
         * @brief Multiplication operator for Quaternions.
         * 
         * @param s The scalar to be multiplied.
         * @return The result of the multiplication.
         */
        Quaternion operator*(float s) const {
            return Quaternion(v * s, w * s);
        }

        /**
         * @brief In-place multiplication operator for Quaternions.
         * 
         * @param s The scalar to be multiplied.
         * @return The result of the multiplication.
         */
        Quaternion& operator*=(float s) {
            v *= s;
            w *= s;
            return *this;
        }

        /**
         * @brief Division operator for Quaternions.
         * 
         * @param s The scalar to be divided.
         * @return The result of the division.
         */
        Quaternion operator/(float s) const {
            return Quaternion(v / s, w / s);
        }

        /**
         * @brief In-place division operator for Quaternions.
         * 
         * @param s The scalar to be divided.
         * @return The result of the division.
         */
        Quaternion& operator/=(float s) {
            v /= s;
            w /= s;
            return *this;
        }

        /**
         * @brief Outputs the quaternion to an ostream.
         * @param os The ostream.
         * @param r The quaternion.
         * @return The ostream.
         */
        friend ostream& operator<<(ostream& os, const Quaternion& q) {
            os << "Quaternion(" << q.v << ", " << q.w << ")";
            return os;
        }
};

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
inline Quaternion normalize(const Quaternion& q);

/**
 * Converts a Quaternion to a Transform.
 *
 * @param q The Quaternion to convert.
 * @return The resulting Transform.
 */
Transform quaternion_to_transform(const Quaternion& q);

/**
 * Performs spherical linear interpolation between two quaternions.
 *
 * @param q1 The starting quaternion.
 * @param q2 The ending quaternion.
 * @param t The interpolation factor (between 0 and 1).
 * @return The interpolated quaternion.
 */
Quaternion slerp(const Quaternion& q1, const Quaternion& q2, float t);




/**
 * @brief Represents an animated transformation between two Transform objects.
 * 
 * The AnimatedTransform class stores the start and end Transform objects, along with the corresponding start and end times.
 * It provides methods to decompose the transformation matrix into translation, rotation, and scaling components.
 */
class AnimatedTransform {
    public:
        const Transform* start_transform;
        const Transform* end_transform;
        const float start_time;
        const float end_time;
        bool actually_animated;
        Matrix4x4 scaling[2];
        Vector3f translation[2];
        Quaternion rotation[2];
        bool has_rotation;

        AnimatedTransform(const Transform* start_transform, float start_time, const Transform* end_transform, float end_time) :
            start_transform(start_transform),
            end_transform(end_transform),
            start_time(start_time),
            end_time(end_time),
            actually_animated(*start_transform != *end_transform),
            has_rotation(dot(rotation[0], rotation[1]) < 0.9995f) {
                if (has_rotation) {
                    rotation[1] = -rotation[1];
                }
                decompose(start_transform->matrix, &translation[0], &rotation[0], &scaling[0]);
                decompose(end_transform->matrix, &translation[1], &rotation[1], &scaling[1]);
            }

        void decompose(Matrix4x4 transform_matrix, Vector3f* translation, Quaternion* rotation, Matrix4x4* scaling) const {
            // Extract translation
            translation->x = transform_matrix.m[0][3];
            translation->y = transform_matrix.m[1][3];
            translation->z = transform_matrix.m[2][3];

            // Remove translation from matrix
            Matrix4x4 mat = transform_matrix;
            for (int i = 0; i < 3; i++) {
                mat.m[i][3] = 0;
            }
            mat.m[3][3] = 1;

            // Extract rotation
            float norm;
            int count = 0;
            Matrix4x4 rot = mat;
            do {
                Matrix4x4 rot_next;
                Matrix4x4 rot_inv_transpose = invert(transpose(rot));
                for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < 4; i++) {
                        rot_next.m[i][j] = 0.5f * (rot.m[i][j] + rot_inv_transpose.m[i][j]);
                    }
                }

                norm = 0;
                for (int i = 0; i < 3; i++) {
                    float n =   abs(rot.m[i][0] - rot_next.m[i][0]) +
                                abs(rot.m[i][1] - rot_next.m[i][1]) +
                                abs(rot.m[i][2] - rot_next.m[i][2]);
                    norm = max(norm, n);
                }

                rot = rot_next;
            } while (++count < 100 && norm > 0.0001f);

            *rotation = Quaternion(rot);

            // Extract scaling
            *scaling = invert(rot) * mat;
        }

        void interpolate(float time, Transform* t) {
            // Verify that the transformation is actually animated
            if (!actually_animated || time <= start_time) {
                *t = *start_transform;
                return;
            } else if (time >= end_time) {
                *t = *end_transform;
                return;
            }

            float dt = (time - start_time) / (end_time - start_time);

            // Interpolate translation
            Vector3f trans = (1 - dt) * translation[0] + dt * translation[1];

            // Interpolate rotation
            Quaternion rot = slerp(rotation[0], rotation[1], dt);

            // Interpolate scaling
            Matrix4x4 scale;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 4; j++) {
                    scale.m[i][j] = lerp(dt, scaling[0].m[i][j], scaling[1].m[i][j]);
                }
            }

            // Interpolated = translation * rotation * scaling
            *t = translate(trans) * quaternion_to_transform(rot) * Transform(scale, invert(scale));
        }
};

/**
 * Transforms a 3D point using the specified transformation matrix.
 *
 * @param point The 3D point to be transformed.
 * @param t The transformation matrix.
 * @return The transformed 3D point.
 */
Point3f transform(const Point3f& point, const Transform& t);

/**
 * Transforms a 3D vector using the specified transformation.
 *
 * @param vector The vector to be transformed.
 * @param t The transformation to apply.
 * @return The transformed vector.
 */
Vector3f transform(const Vector3f& vector, const Transform& t);

/**
 * Transforms a ray using the given transformation matrix.
 *
 * @param ray The ray to be transformed.
 * @param t The transformation matrix.
 * @return The transformed ray.
 */
Ray transform(const Ray& ray, const Transform& t);

/**
 * Transforms a RayDifferential object using the specified Transform.
 *
 * @param ray The RayDifferential object to be transformed.
 * @param t The Transform to apply to the ray.
 * @return The transformed RayDifferential object.
 */
RayDifferential transform(const RayDifferential& ray, const Transform& t);

#endif // !ANIMATED_TRANSFORMS