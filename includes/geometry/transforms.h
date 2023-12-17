#ifndef TRANSFORMS
#define TRANSFORMS

#include <iostream>
#include "vectors.h"
#include "points.h"
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