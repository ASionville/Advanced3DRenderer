#ifndef NORMALS
#define NORMALS

#include <iostream>
#include "vectors.h"

using namespace std;

/**
 * @brief Template class for 3D normals.
 * @tparam T The type of the components.
 */
template <typename T> class Normal3 {
    public :
        T x, y, z;

        /**
         * @brief Default constructor. Initializes x, y, and z to 0.
         */
        Normal3() : x(0), y(0), z(0) {}

        /**
         * @brief Constructor with components.
         * @param x The x component.
         * @param y The y component.
         * @param z The z component.
         */
        Normal3(T x, T y, T z) : x(x), y(y), z(z) {}

        /**
         * @brief Constructor from a vector.
         * @param v The vector to convert.
         */
        template <typename U> explicit Normal3(const Vector3<U>& v) : x(v.x), y(v.y), z(v.z) {
            Assert(!v.has_NaNs());
        }

        /**
         * @brief Check if the normal has any NaN components.
         * @return True if either x, y, or z is NaN, false otherwise.
         */
        bool has_NaNs() const {
            return isnan(x) || isnan(y) || isnan(z);
        }

        /**
         * @brief Access operator.
         * @param i The index (0 for x, 1 for y, 2 for z).
         * @return The component at the given index.
         */
        T operator[](int i) const {
            Assert(i >= 0 && i <= 2);
            if (i == 0) return x;
            if (i == 1) return y;
            return z;
        }

        /**
         * @brief Non-const access operator.
         * @param i The index (0 for x, 1 for y, 2 for z).
         * @return The component at the given index.
         */
        T& operator[](int i) {
            Assert(i >= 0 && i <= 2);
            if (i == 0) return x;
            if (i == 1) return y;
            return z;
        }

        /**
         * @brief Negation operator.
         * @return The negation of this normal.
         */
        Normal3<T> operator-() const {
            return Normal3<T>(-x, -y, -z);
        }

        /**
         * @brief Addition operator with a normal.
         * @param n The normal to add.
         * @return The result of the addition.
         */
        Normal3<T> operator+(const Normal3<T>& n) const {
            return Normal3<T>(x + n.x, y + n.y, z + n.z);
        }

        /**
         * @brief In-place addition operator with a normal.
         * @param n The normal to add.
         * @return This normal after the addition.
         */
        Normal3<T>& operator+=(const Normal3<T>& n) {
            x += n.x;
            y += n.y;
            z += n.z;
            return *this;
        }

        /**
         * @brief Subtraction operator with a normal.
         * @param n The normal to subtract.
         * @return The result of the subtraction.
         */
        Normal3<T> operator-(const Normal3<T>& n) const {
            return Normal3<T>(x - n.x, y - n.y, z - n.z);
        }

        /**
         * @brief In-place subtraction operator with a normal.
         * @param n The normal to subtract.
         * @return This normal after the subtraction.
         */
        Normal3<T>& operator-=(const Normal3<T>& n) {
            x -= n.x;
            y -= n.y;
            z -= n.z;
            return *this;
        }

        /**
         * @brief Multiplication operator with a scalar.
         * @param s The scalar to multiply.
         * @return The result of the multiplication.
         */
        Normal3<T> operator*(T s) const {
            return Normal3<T>(x * s, y * s, z * s);
        }

        /**
         * @brief In-place multiplication operator with a scalar.
         * @param s The scalar to multiply.
         * @return This normal after the multiplication.
         */
        Normal3<T>& operator*=(T s) {
            x *= s;
            y *= s;
            z *= s;
            return *this;
        }

        /**
         * @brief Division operator with a scalar.
         * @param s The scalar to divide.
         * @return The result of the division.
         */
        Normal3<T> operator/(T s) const {
            T inv = 1 / s;
            return Normal3<T>(x * inv, y * inv, z * inv);
        }

        /**
         * @brief In-place division operator with a scalar.
         * @param s The scalar to divide.
         * @return This normal after the division.
         */
        Normal3<T>& operator/=(T s) {
            T inv = 1 / s;
            x *= inv;
            y *= inv;
            z *= inv;
            return *this;
        }

        /**
         * Calculates the squared length of the vector.
         *
         * @return The squared length of the vector.
         */
        float length_squared() const {
            return x * x + y * y + z * z;
        }

        /**
         * Calculates the length of the vector.
         *
         * @return The length of the vector.
         */
        float length() const {
            return sqrt(length_squared());
        }

        /**
         * @brief Normalizes the current normal vector.
         * 
         * @return The normalized normal vector.
         */
        Normal3<T> normalize() const {
            return *this / length();
        }
};

/**
 * @brief Constructor that initializes a Vector3 from a Normal3.
 * @param n The Normal3 to convert.
 */
template <typename T> inline Vector3<T>::Vector3(const Normal3<T>& n) : x(n.x), y(n.y), z(n.z) {
    Assert(!n.has_NaNs());
}

/**
 * @brief Computes the dot product of a Normal3 and a Vector3.
 * @param n1 The Normal3.
 * @param v2 The Vector3.
 * @return The dot product of n1 and v2.
 */
template <typename T> inline float dot(const Normal3<T>& n1, const Vector3<T>& v2) {
    Assert(!n1.has_NaNs() && !v2.has_NaNs());
    return n1.x * v2.x + n1.y * v2.y + n1.z * v2.z;
}

/**
 * @brief Computes the dot product of a Vector3 and a Normal3.
 * @param v1 The Vector3.
 * @param n2 The Normal3.
 * @return The dot product of v1 and n2.
 */
template <typename T> inline float dot(const Vector3<T>& v1, const Normal3<T>& n2) {
    Assert(!v1.has_NaNs() && !n2.has_NaNs());
    return v1.x * n2.x + v1.y * n2.y + v1.z * n2.z;
}

/**
 * @brief Computes the dot product of two Normal3s.
 * @param n1 The first Normal3.
 * @param n2 The second Normal3.
 * @return The dot product of n1 and n2.
 */
template <typename T> inline float dot(const Normal3<T>& n1, const Normal3<T>& n2) {
    Assert(!n1.has_NaNs() && !n2.has_NaNs());
    return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
}

/**
 * @brief Adjusts the direction of a Normal3 to face towards a Vector3.
 * @param n The Normal3.
 * @param v The Vector3.
 * @return The adjusted Normal3.
 */
template <typename T> inline Normal3<T> face_forward(const Normal3<T>& n, const Vector3<T>& v) {
    return (dot(n, v) < 0.f) ? -n : n;
}

/**
 * @brief Adjusts the direction of a Normal3 to face towards another Normal3.
 * @param n The first Normal3.
 * @param n2 The second Normal3.
 * @return The adjusted Normal3.
 */
template <typename T> inline Normal3<T> face_forward(const Normal3<T>& n, const Normal3<T>& n2) {
    return (dot(n, n2) < 0.f) ? -n : n;
}

/**
 * @brief Adjusts the direction of a Vector3 to face towards a Normal3.
 * @param v The Vector3.
 * @param n The Normal3.
 * @return The adjusted Vector3.
 */
template <typename T> inline Vector3<T> face_forward(const Vector3<T>& v, const Normal3<T>& n) {
    return (dot(v, n) < 0.f) ? -v : v;
}

/**
 * @brief Adjusts the direction of a Vector3 to face towards another Vector3.
 * @param v The first Vector3.
 * @param v2 The second Vector3.
 * @return The adjusted Vector3.
 */
template <typename T> inline Vector3<T> face_forward(const Vector3<T>& v, const Vector3<T>& v2) {
    return (dot(v, v2) < 0.f) ? -v : v;
}

typedef Normal3<float> Normal3f;


#endif // !NORMALS