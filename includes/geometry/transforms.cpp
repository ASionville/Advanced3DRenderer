#include <iostream>
#include "vectors.h"
#include "points.h"
#include "../utils/maths.h"
#include "transforms.h"

using namespace std;

/**
 * Transposes a 4x4 matrix.
 *
 * @param m The input matrix to be transposed.
 * @return The transposed matrix.
 */
Matrix4x4 transpose(const Matrix4x4& m) {
    return Matrix4x4(m.m[0][0], m.m[1][0], m.m[2][0], m.m[3][0],
                     m.m[0][1], m.m[1][1], m.m[2][1], m.m[3][1],
                     m.m[0][2], m.m[1][2], m.m[2][2], m.m[3][2],
                     m.m[0][3], m.m[1][3], m.m[2][3], m.m[3][3]);
}

/**
 * @brief Calculates the inverse of a 4x4 matrix.
 * @details Uses Gauss-Jordan with pivoting
 *
 * @param m The input matrix to invert.
 * @return The inverted matrix.
 */
Matrix4x4 invert(const Matrix4x4& m) {
    Matrix4x4 a = m;
    // b will contain the inverse of a, initialized to the identity
    Matrix4x4 b = Matrix4x4();
    int i, j, i1, j1, k, l;
    float max, t;

    // Gaussian elimination with pivoting
    for (i = 0; i < 4; i++) {
        // Find pivot
        max = a.m[i][i];
        i1 = i;
        for (j = i + 1; j < 4; j++) {
            if (fabs(a.m[j][i]) > fabs(max)) {
                max = a.m[j][i];
                i1 = j;
            }
        }

        // Swap rows
        if (i1 != i) {
            for (j = 0; j < 4; j++) {
                t = a.m[i][j];
                a.m[i][j] = a.m[i1][j];
                a.m[i1][j] = t;
                t = b.m[i][j];
                b.m[i][j] = b.m[i1][j];
                b.m[i1][j] = t;
            }
        }

        // Normalize row
        for (j = 0; j < 4; j++) {
            a.m[i][j] /= max;
            b.m[i][j] /= max;
        }

        // Eliminate column
        for (j = 0; j < 4; j++) {
            if (i != j) {
                t = a.m[j][i];
                for (k = 0; k < 4; k++) {
                    a.m[j][k] -= a.m[i][k] * t;
                    b.m[j][k] -= b.m[i][k] * t;
                }
            }
        }
    }

    return b;
}




/**
 * Creates a translation Transform.
 *
 * @param delta The translation vector.
 * @return The transformed vector.
 */
Transform translate(const Vector3f& delta) {
    Matrix4x4 mat = Matrix4x4(1, 0, 0, delta.x,
                              0, 1, 0, delta.y,
                              0, 0, 1, delta.z,
                              0, 0, 0, 1);
    Matrix4x4 inv_mat = Matrix4x4(1, 0, 0, -delta.x,
                                  0, 1, 0, -delta.y,
                                  0, 0, 1, -delta.z,
                                  0, 0, 0, 1);
    return Transform(mat, inv_mat);
}

/**
 * Creates a scaling Transform.
 *
 * @param x The x scaling factor.
 * @param y The y scaling factor.
 * @param z The z scaling factor.
 * @return The transformed vector.
 */
Transform scale(float x, float y, float z) {
    Matrix4x4 mat = Matrix4x4(x, 0, 0, 0,
                              0, y, 0, 0,
                              0, 0, z, 0,
                              0, 0, 0, 1);
    Matrix4x4 inv_mat = Matrix4x4(1 / x, 0, 0, 0,
                                  0, 1 / y, 0, 0,
                                  0, 0, 1 / z, 0,
                                  0, 0, 0, 1);
    return Transform(mat, inv_mat);
}

/**
 * Creates a rotation Transform about the x axis.
 *
 * @param theta The angle of rotation in degrees.
 * @return The transformed vector.
 */
Transform rotate_x(float theta) {
    float sin_theta = sin(radians(theta));
    float cos_theta = cos(radians(theta));
    Matrix4x4 mat = Matrix4x4(1, 0, 0, 0,
                              0, cos_theta, -sin_theta, 0,
                              0, sin_theta, cos_theta, 0,
                              0, 0, 0, 1);
    return Transform(mat, transpose(mat));
}

/**
 * Creates a rotation Transform about the y axis.
 *
 * @param theta The angle of rotation in degrees.
 * @return The transformed vector.
 */
Transform rotate_y(float theta) {
    float sin_theta = sin(radians(theta));
    float cos_theta = cos(radians(theta));
    Matrix4x4 mat = Matrix4x4(cos_theta, 0, sin_theta, 0,
                              0, 1, 0, 0,
                              -sin_theta, 0, cos_theta, 0,
                              0, 0, 0, 1);
    return Transform(mat, transpose(mat));
}

/**
 * Creates a rotation Transform about the z axis.
 *
 * @param theta The angle of rotation in degrees.
 * @return The transformed vector.
 */
Transform rotate_z(float theta) {
    float sin_theta = sin(radians(theta));
    float cos_theta = cos(radians(theta));
    Matrix4x4 mat = Matrix4x4(cos_theta, -sin_theta, 0, 0,
                              sin_theta, cos_theta, 0, 0,
                              0, 0, 1, 0,
                              0, 0, 0, 1);
    return Transform(mat, transpose(mat));
}

/**
 * Creates a rotation Transform about an arbitrary axis.
 *
 * @param theta The angle of rotation in degrees.
 * @param axis The axis of rotation.
 * @return The transformed vector.
 */
Transform rotate(float theta, const Vector3f& axis) {
    Vector3f a = axis.normalize();
    float sin_theta = sin(radians(theta));
    float cos_theta = cos(radians(theta));
    Matrix4x4 mat = Matrix4x4();
    // Compute the rotation matrix
    mat.m[0][0] = a.x * a.x + (1 - a.x * a.x) * cos_theta;
    mat.m[0][1] = a.x * a.y * (1 - cos_theta) - a.z * sin_theta;
    mat.m[0][2] = a.x * a.z * (1 - cos_theta) + a.y * sin_theta;
    mat.m[0][3] = 0;

    mat.m[1][0] = a.x * a.y * (1 - cos_theta) + a.z * sin_theta;
    mat.m[1][1] = a.y * a.y + (1 - a.y * a.y) * cos_theta;
    mat.m[1][2] = a.y * a.z * (1 - cos_theta) - a.x * sin_theta;
    mat.m[1][3] = 0;

    mat.m[2][0] = a.x * a.z * (1 - cos_theta) - a.y * sin_theta;
    mat.m[2][1] = a.y * a.z * (1 - cos_theta) + a.x * sin_theta;
    mat.m[2][2] = a.z * a.z + (1 - a.z * a.z) * cos_theta;
    mat.m[2][3] = 0;

    return Transform(mat, transpose(mat));
}

/**
 * Creates the transformation matrix for "Look At" transformation.
 * 
 * @param pos The position of the camera.
 * @param look The point the camera should look at.
 * @param up The up vector of the camera.
 * @return The transformation matrix.
 */
Transform look_at(const Point3f& pos, const Point3f& look, const Vector3f& up) {
    Matrix4x4 camera_to_world;

    // Initialize the camera basis vectors
    Vector3f direction = (look - pos).normalize();
    Vector3f left = cross(up, direction).normalize();
    Vector3f new_up = cross(direction, left);

    // Set the camera basis vectors
    camera_to_world.m[0][0] = left.x;
    camera_to_world.m[0][1] = left.y;
    camera_to_world.m[0][2] = left.z;
    camera_to_world.m[0][3] = 0; // 0 : vector

    camera_to_world.m[1][0] = new_up.x;
    camera_to_world.m[1][1] = new_up.y;
    camera_to_world.m[1][2] = new_up.z;
    camera_to_world.m[1][3] = 0;

    camera_to_world.m[2][0] = direction.x;
    camera_to_world.m[2][1] = direction.y;
    camera_to_world.m[2][2] = direction.z;
    camera_to_world.m[2][3] = 0;

    // Last row is the camera position
    camera_to_world.m[3][0] = pos.x;
    camera_to_world.m[3][1] = pos.y;
    camera_to_world.m[3][2] = pos.z;
    camera_to_world.m[3][3] = 1; // 1 : point

    // Return the inverse of the camera_to_world Transform
    // because we want to transform from world space to camera space
    return Transform(invert(camera_to_world), camera_to_world);
}