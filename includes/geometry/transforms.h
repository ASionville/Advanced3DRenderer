#ifndef TRANSFORMS
#define TRANSFORMS

#include <iostream>
#include <cassert>
#include "vectors.h"
#include "points.h"
#include "normals.h"
#include "rays.h"
#include "bboxes.h"
#include "../utils/maths.h"

using namespace std;

/**
 * @brief Represents a 4x4 matrix used for transformations in 3D space.
 */
class Matrix4x4 {
    public:
        float m[4][4];

        /**
         * @brief Default constructor. Initializes the matrix to the identity matrix. 
         */
        Matrix4x4() {
            m[0][0] = 1; m[0][1] = 0; m[0][2] = 0; m[0][3] = 0;
            m[1][0] = 0; m[1][1] = 1; m[1][2] = 0; m[1][3] = 0;
            m[2][0] = 0; m[2][1] = 0; m[2][2] = 1; m[2][3] = 0;
            m[3][0] = 0; m[3][1] = 0; m[3][2] = 0; m[3][3] = 1;
        }

        /**
         * @brief Constructor with a 2D array of floats.
         * @param mat The 2D array of floats.
         */
        Matrix4x4(float mat[4][4]) {
            memcpy(m, mat, 16 * sizeof(float));
        }

        /**
         * @brief Constructor with 16 floats.
         * @param tij The 16 floats of the matrix written as M_ij
         */
        Matrix4x4(float t00, float t01, float t02, float t03,
                  float t10, float t11, float t12, float t13,
                  float t20, float t21, float t22, float t23,
                  float t30, float t31, float t32, float t33) {
            m[0][0] = t00; m[0][1] = t01; m[0][2] = t02; m[0][3] = t03;
            m[1][0] = t10; m[1][1] = t11; m[1][2] = t12; m[1][3] = t13;
            m[2][0] = t20; m[2][1] = t21; m[2][2] = t22; m[2][3] = t23;
            m[3][0] = t30; m[3][1] = t31; m[3][2] = t32; m[3][3] = t33;
        }

        /**
         * @brief Outputs the matrix to an ostream.
         * @param os The ostream.
         * @param mat The matrix.
         * @return The ostream.
         */
        friend ostream& operator<<(ostream& os, const Matrix4x4& mat) {
            os << "Matrix4x4(" << mat.m[0][0] << ", " << mat.m[0][1] << ", " << mat.m[0][2] << ", " << mat.m[0][3] << ", "
                               << mat.m[1][0] << ", " << mat.m[1][1] << ", " << mat.m[1][2] << ", " << mat.m[1][3] << ", "
                               << mat.m[2][0] << ", " << mat.m[2][1] << ", " << mat.m[2][2] << ", " << mat.m[2][3] << ", "
                               << mat.m[3][0] << ", " << mat.m[3][1] << ", " << mat.m[3][2] << ", " << mat.m[3][3] << ")";
            return os;
        }

        /**
         * @brief Access operator.
         * @param i The row index.
         * @param j The column index.
         * @return The element at the given row and column.
         */
        float operator()(int i, int j) const {
            assert(i >= 0 && i <= 3);
            assert(j >= 0 && j <= 3);
            return m[i][j];
        }

        /**
         * @brief Non-const access operator.
         * @param i The row index.
         * @param j The column index.
         * @return The element at the given row and column.
         */
        float& operator()(int i, int j) {
            assert(i >= 0 && i <= 3);
            assert(j >= 0 && j <= 3);
            return m[i][j];
        }

        /**
         * @brief Negation operator.
         * @return The negation of this matrix.
         */
        Matrix4x4 operator-() const {
            return Matrix4x4(-m[0][0], -m[0][1], -m[0][2], -m[0][3],
                             -m[1][0], -m[1][1], -m[1][2], -m[1][3],
                             -m[2][0], -m[2][1], -m[2][2], -m[2][3],
                             -m[3][0], -m[3][1], -m[3][2], -m[3][3]);
        }

        /**
         * @brief Addition operator with a matrix.
         * @param mat The matrix to add.
         * @return The result of the addition.
         */
        Matrix4x4 operator+(const Matrix4x4& mat) const {
            return Matrix4x4(m[0][0] + mat.m[0][0], m[0][1] + mat.m[0][1], m[0][2] + mat.m[0][2], m[0][3] + mat.m[0][3],
                             m[1][0] + mat.m[1][0], m[1][1] + mat.m[1][1], m[1][2] + mat.m[1][2], m[1][3] + mat.m[1][3],
                             m[2][0] + mat.m[2][0], m[2][1] + mat.m[2][1], m[2][2] + mat.m[2][2], m[2][3] + mat.m[2][3],
                             m[3][0] + mat.m[3][0], m[3][1] + mat.m[3][1], m[3][2] + mat.m[3][2], m[3][3] + mat.m[3][3]);
        }

        /**
         * @brief In-place addition operator with a matrix.
         * @param mat The matrix to add.
         * @return This matrix after the addition.
         */
        Matrix4x4& operator+=(const Matrix4x4& mat) {
            m[0][0] += mat.m[0][0]; m[0][1] += mat.m[0][1]; m[0][2] += mat.m[0][2]; m[0][3] += mat.m[0][3];
            m[1][0] += mat.m[1][0]; m[1][1] += mat.m[1][1]; m[1][2] += mat.m[1][2]; m[1][3] += mat.m[1][3];
            m[2][0] += mat.m[2][0]; m[2][1] += mat.m[2][1]; m[2][2] += mat.m[2][2]; m[2][3] += mat.m[2][3];
            m[3][0] += mat.m[3][0]; m[3][1] += mat.m[3][1]; m[3][2] += mat.m[3][2]; m[3][3] += mat.m[3][3];
            return *this;
        }

        /**
         * @brief Subtraction operator with a matrix.
         * @param mat The matrix to subtract.
         * @return The result of the subtraction.
         */
        Matrix4x4 operator-(const Matrix4x4& mat) const {
            return Matrix4x4(m[0][0] - mat.m[0][0], m[0][1] - mat.m[0][1], m[0][2] - mat.m[0][2], m[0][3] - mat.m[0][3],
                             m[1][0] - mat.m[1][0], m[1][1] - mat.m[1][1], m[1][2] - mat.m[1][2], m[1][3] - mat.m[1][3],
                             m[2][0] - mat.m[2][0], m[2][1] - mat.m[2][1], m[2][2] - mat.m[2][2], m[2][3] - mat.m[2][3],
                             m[3][0] - mat.m[3][0], m[3][1] - mat.m[3][1], m[3][2] - mat.m[3][2], m[3][3] - mat.m[3][3]);
        }

        /**
         * @brief In-place subtraction operator with a matrix.
         * @param mat The matrix to subtract.
         * @return This matrix after the subtraction.
         */
        Matrix4x4& operator-=(const Matrix4x4& mat) {
            m[0][0] -= mat.m[0][0]; m[0][1] -= mat.m[0][1]; m[0][2] -= mat.m[0][2]; m[0][3] -= mat.m[0][3];
            m[1][0] -= mat.m[1][0]; m[1][1] -= mat.m[1][1]; m[1][2] -= mat.m[1][2]; m[1][3] -= mat.m[1][3];
            m[2][0] -= mat.m[2][0]; m[2][1] -= mat.m[2][1]; m[2][2] -= mat.m[2][2]; m[2][3] -= mat.m[2][3];
            m[3][0] -= mat.m[3][0]; m[3][1] -= mat.m[3][1]; m[3][2] -= mat.m[3][2]; m[3][3] -= mat.m[3][3];
            return *this;
        }

        /**
         * @brief Multiplication operator with a scalar.
         * @param s The scalar to multiply with.
         * @return The result of the multiplication.
         */
        Matrix4x4 operator*(float s) const {
            return Matrix4x4(m[0][0] * s, m[0][1] * s, m[0][2] * s, m[0][3] * s,
                             m[1][0] * s, m[1][1] * s, m[1][2] * s, m[1][3] * s,
                             m[2][0] * s, m[2][1] * s, m[2][2] * s, m[2][3] * s,
                             m[3][0] * s, m[3][1] * s, m[3][2] * s, m[3][3] * s);
        }

        /**
         * @brief In-place multiplication operator with a scalar.
         * @param s The scalar to multiply with.
         * @return This matrix after the multiplication.
         */
        Matrix4x4& operator*=(float s) {
            m[0][0] *= s; m[0][1] *= s; m[0][2] *= s; m[0][3] *= s;
            m[1][0] *= s; m[1][1] *= s; m[1][2] *= s; m[1][3] *= s;
            m[2][0] *= s; m[2][1] *= s; m[2][2] *= s; m[2][3] *= s;
            m[3][0] *= s; m[3][1] *= s; m[3][2] *= s; m[3][3] *= s;
            return *this;
        }

        /**
         * @brief Division operator with a scalar.
         * @param s The scalar to divide with.
         * @return The result of the division.
         */
        Matrix4x4 operator/(float s) const {
            float inv = 1 / s;
            return Matrix4x4(m[0][0] * inv, m[0][1] * inv, m[0][2] * inv, m[0][3] * inv,
                             m[1][0] * inv, m[1][1] * inv, m[1][2] * inv, m[1][3] * inv,
                             m[2][0] * inv, m[2][1] * inv, m[2][2] * inv, m[2][3] * inv,
                             m[3][0] * inv, m[3][1] * inv, m[3][2] * inv, m[3][3] * inv);
        }

        /**
         * @brief In-place division operator with a scalar.
         * @param s The scalar to divide with.
         * @return This matrix after the division.
         */
        Matrix4x4& operator/=(float s) {
            float inv = 1 / s;
            m[0][0] *= inv; m[0][1] *= inv; m[0][2] *= inv; m[0][3] *= inv;
            m[1][0] *= inv; m[1][1] *= inv; m[1][2] *= inv; m[1][3] *= inv;
            m[2][0] *= inv; m[2][1] *= inv; m[2][2] *= inv; m[2][3] *= inv;
            m[3][0] *= inv; m[3][1] *= inv; m[3][2] *= inv; m[3][3] *= inv;
            return *this;
        }

        /**
         * @brief Multiplication operator with a matrix.
         * @param mat The matrix to multiply with.
         * @return The result of the multiplication.
         */
        Matrix4x4 operator*(const Matrix4x4& mat) const {
            Matrix4x4 result;
            for (int i = 0; i < 4; i++) {
                for (int j = 0; i < 4; i++) {
                    result.m[i][j] = m[i][0] * mat.m[0][j] +
                                     m[i][1] * mat.m[1][j] +
                                     m[i][2] * mat.m[2][j] +
                                     m[i][3] * mat.m[3][j];
                }
            }
            return result;
        }

        /**
         * @brief In-place multiplication operator with a matrix.
         * @param mat The matrix to multiply with.
         * @return This matrix after the multiplication.
         */
        Matrix4x4& operator*=(const Matrix4x4& mat) {
            Matrix4x4 result;
            for (int i = 0; i < 4; i++) {
                for (int j = 0; i < 4; i++) {
                    result.m[i][j] = m[i][0] * mat.m[0][j] +
                                     m[i][1] * mat.m[1][j] +
                                     m[i][2] * mat.m[2][j] +
                                     m[i][3] * mat.m[3][j];
                }
            }
            memcpy(m, result.m, 16 * sizeof(float));
            return *this;
        }
};

/**
 * Transposes a 4x4 matrix.
 *
 * @param m The input matrix to be transposed.
 * @return The transposed matrix.
 */
Matrix4x4 transpose(const Matrix4x4& m);

/**
 * @brief Calculates the inverse of a 4x4 matrix.
 * @details Uses Gauss-Jordan with pivoting
 *
 * @param m The input matrix to invert.
 * @return The inverted matrix.
 */
Matrix4x4 invert(const Matrix4x4& m);




/**
 * @brief A class representing 3D transformations
 * @details Uses homogenous coordinates to represent 3D transformations with a 4x4 matrix.
*/
class Transform {
    public:
        Matrix4x4 matrix;
        Matrix4x4 inv_matrix;

        /**
         * @brief Default constructor. Initializes the matrix to the identity matrix.
         * @detais We use the default constructor of Matrix4x4.
         */
        Transform() { }

        /**
         * @brief Constructor with a 4x4 matrix of floats.
         * @param mat The 4x4 matrix of floats.
         */
        Transform(const float mat[4][4]) {
            matrix = Matrix4x4(mat[0][0], mat[0][1], mat[0][2], mat[0][3],
                               mat[1][0], mat[1][1], mat[1][2], mat[1][3],
                               mat[2][0], mat[2][1], mat[2][2], mat[2][3],
                               mat[3][0], mat[3][1], mat[3][2], mat[3][3]);
            inv_matrix = invert(matrix);
        }

        /**
         * @brief Constructor with a 4x4 matrix and its inverse.
         * @param mat The 4x4 matrix
         * @param inv_mat The inverse of the 4x4 matrix
         */
        Transform(const Matrix4x4& mat, const Matrix4x4& inv_mat) :
            matrix(mat), inv_matrix(inv_mat) { }

        /**
         * @brief Inverts the given Transform.
         * 
         * @param t The Transform to invert.
         * @return The inverted Transform.
         */
        friend Transform invert(const Transform& t) {
            return Transform(t.inv_matrix, t.matrix);
        }

        /**
         * @brief Transposes the given Transform.
         * 
         * @param t The Transform to be transposed.
         * @return The transposed Transform.
         */
        friend Transform transpose(const Transform& t) {
            return Transform(t.matrix, t.inv_matrix);
        }

        /**
         * @brief Outputs the Transform to an ostream.
         * @param os The ostream.
         * @param t The Transform.
         * @return The ostream.
         */
        friend ostream& operator<<(ostream& os, const Transform& t) {
            os << "Transform(" << t.matrix << ", " << t.inv_matrix << ")";
            return os;
        }

        Transform& operator*(const Transform& t2) const {
            return Transform(matrix * t2.matrix, t2.inv_matrix * inv_matrix);
        }

        /**
         * Checks if the Transform has scaling applied.
         *
         * @return True if scaling is applied, false otherwise.
         */
        bool has_scale() const {
            // Get the transformed basis vectors
            Vector3f la = Vector3f(matrix.m[0][0], matrix.m[0][1], matrix.m[0][2]);
            Vector3f lb = Vector3f(matrix.m[1][0], matrix.m[1][1], matrix.m[1][2]);
            Vector3f lc = Vector3f(matrix.m[2][0], matrix.m[2][1], matrix.m[2][2]);

            float la2 = la.length_squared();
            float lb2 = lb.length_squared();
            float lc2 = lc.length_squared();

            // If the length of any of the transformed basis vectors is not 1,
            // then the Transform has scale.
            return (la2 < 0.999f || la2 > 1.001f ||
                lb2 < 0.999f || lb2 > 1.001f ||
                lc2 < 0.999f || lc2 > 1.001f);
        }

        /**
         * Checks if the transformation swaps the handedness of the coordinate system.
         * 
         * @return True if the transformation swaps handedness, false otherwise.
         */
        bool swaps_handedness() const {
            float det = matrix.m[0][0] * (matrix.m[1][1] * matrix.m[2][2] - matrix.m[1][2] * matrix.m[2][1]) -
                        matrix.m[0][1] * (matrix.m[1][0] * matrix.m[2][2] - matrix.m[1][2] * matrix.m[2][0]) +
                        matrix.m[0][2] * (matrix.m[1][0] * matrix.m[2][1] - matrix.m[1][1] * matrix.m[2][0]);
            return det < 0;
        }

        /**
         * Transforms a 3D point using the current transformation.
         * 
         * @tparam T The type of the point coordinates.
         * @param p The point to be transformed.
         * @return The transformed point.
         */
        template <typename T> inline Point3<T> transform_point(const Point3<T>& p) const {
            T x = p.x, y = p.y, z = p.z;
            T xp = m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z + m.m[0][3];
            T yp = m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z + m.m[1][3];
            T zp = m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z + m.m[2][3];
            T wp = m.m[3][0] * x + m.m[3][1] * y + m.m[3][2] * z + m.m[3][3];
            // "Normalize" the point's homogenous coordinates
            if (wp == 1) {
                return Point3<T>(xp, yp, zp);
            } else {
                return Point3<T>(xp, yp, zp) / wp;
            }
        }

        /**
         * Transforms a vector using the current transformation.
         * 
         * @param v The vector to be transformed.
         * @return The transformed vector.
         */
        template <typename T> inline Vector3<T> transform_vector(const Vector3<T>& v) const {
            T x = v.x, y = v.y, z = v.z;
            return Vector3<T>(m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z,
                              m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z,
                              m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z);
        }

        /**
         * Transforms a normal using the current transformation.
         * 
         * @param n The normal to be transformed.
         * @return The transformed normal.
         */
        template <typename T> inline Normal3<T> transform_normal(const Normal3<T>& n) const {
            T x = n.x, y = n.y, z = n.z;
            return Normal3<T>(inv_m.m[0][0] * x + inv_m.m[1][0] * y + inv_m.m[2][0] * z,
                              inv_m.m[0][1] * x + inv_m.m[1][1] * y + inv_m.m[2][1] * z,
                              inv_m.m[0][2] * x + inv_m.m[1][2] * y + inv_m.m[2][2] * z);
        }

        /**
         * Transforms a ray using the current transformation.
         * 
         * @param r The ray to be transformed.
         * @return The transformed ray.
         */
        Ray transform_ray(const Ray& r) const {
            Point3f new_origin = transform_point(r.origin);
            Vector3f new_direction = transform_vector(r.direction);
            return Ray(new_origin, new_direction, r.t_max, r.time);
        }

        /**
         * Transforms a bounding box using the current transformation.
         * 
         * @param bbox The bounding box to be transformed.
         * @return The transformed bounding box.
         */
        BBox3f transform_bbox(const BBox3f& bbox) const {
            Point3f p_min = bbox.p_min;
            Point3f p_max = bbox.p_max;
            BBox3f new_bbox = BBox3f(transform_point(p_min));
            new_bbox = union_bbox(new_bbox, transform_point(Point3f(p_max.x, p_min.y, p_min.z)));
            new_bbox = union_bbox(new_bbox, transform_point(Point3f(p_min.x, p_max.y, p_min.z)));
            new_bbox = union_bbox(new_bbox, transform_point(Point3f(p_min.x, p_min.y, p_max.z)));
            new_bbox = union_bbox(new_bbox, transform_point(Point3f(p_min.x, p_max.y, p_max.z)));
            new_bbox = union_bbox(new_bbox, transform_point(Point3f(p_max.x, p_max.y, p_min.z)));
            new_bbox = union_bbox(new_bbox, transform_point(Point3f(p_max.x, p_min.y, p_max.z)));
            new_bbox = union_bbox(new_bbox, transform_point(p_max));
            return new_bbox;
        }

};

/**
 * Creates a translation Transform.
 *
 * @param delta The translation vector.
 * @return The transformed vector.
 */
Transform translate(const Vector3f& delta);

/**
 * Creates a scaling Transform.
 *
 * @param x The x scaling factor.
 * @param y The y scaling factor.
 * @param z The z scaling factor.
 * @return The transformed vector.
 */
Transform scale(float x, float y, float z);

/**
 * Creates a rotation Transform about the x axis.
 *
 * @param theta The angle of rotation in degrees.
 * @return The transformed vector.
 */
Transform rotate_x(float theta);

/**
 * Creates a rotation Transform about the y axis.
 *
 * @param theta The angle of rotation in degrees.
 * @return The transformed vector.
 */
Transform rotate_y(float theta);

/**
 * Creates a rotation Transform about the z axis.
 *
 * @param theta The angle of rotation in degrees.
 * @return The transformed vector.
 */
Transform rotate_z(float theta);

/**
 * Creates a rotation Transform about an arbitrary axis.
 *
 * @param theta The angle of rotation in degrees.
 * @param axis The axis of rotation.
 * @return The transformed vector.
 */
Transform rotate(float theta, const Vector3f& axis);

/**
 * Creates the transformation matrix for "Look At" transformation.
 * 
 * @param pos The position of the camera.
 * @param look The point the camera should look at.
 * @param up The up vector of the camera.
 * @return The transformation matrix.
 */
Transform look_at(const Point3f& pos, const Point3f& look, const Vector3f& up);

#endif // !TRANSFORMS