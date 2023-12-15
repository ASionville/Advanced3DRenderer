#ifndef VECTORS
#define VECTORS

#include <cmath>
#include <iostream>

using namespace std;

/**
 * @brief Multiplies a scalar value with a Vector2 object.
 * 
 * This function multiplies a scalar value with a Vector2 object and returns the result.
 * 
 * @tparam T The type of the elements in the Vector2 object.
 * @param s The scalar value to multiply.
 * @param v The Vector2 object to multiply.
 * @return The result of multiplying the scalar value with the Vector2 object.
 */
template <typename T> inline Vector2<T> operator* (T s, const Vector2<T> &v) {
    return v * s;
}

/**
 * Calculates the absolute value of each component of a 2D vector.
 * 
 * @tparam T The type of the vector components.
 * @param v The input vector.
 * @return A new vector with the absolute values of the input vector's components.
 */
template <typename T> inline Vector2<T> abs(const Vector2<T> &v) {
    return Vector2<T>(abs(v.x), abs(v.y));
}

/**
 * Calculates the dot product of two 2D vectors.
 * 
 * @tparam T The type of the vector components.
 * @param v1 The first vector.
 * @param v2 The second vector.
 * @return The dot product of the two vectors.
 */
template <typename T> inline T dot(const Vector2<T> &v1, const Vector2<T> &v2) {
    return v1.x * v2.x + v1.y * v2.y;
}

/**
 * Calculates the absolute value of the dot product between two Vector2 objects.
 * 
 * @param v1 The first Vector2 object.
 * @param v2 The second Vector2 object.
 * @return The absolute value of the dot product between v1 and v2.
 */
template <typename T> inline T abs_dot(const Vector2<T> &v1, const Vector2<T> &v2) {
    return abs(dot(v1, v2));
}

/**
 * Calculates the cross product of two 2D vectors.
 * 
 * @tparam T The type of the vector components.
 * @param v1 The first vector.
 * @param v2 The second vector.
 * @return The cross product of the two vectors.
 */
template <typename T> inline T cross(const Vector2<T> &v1, const Vector2<T> &v2) {
    return v1.x * v2.y - v1.y * v2.x;
}

/**
 * Returns the minimum component of a 2D vector.
 * 
 * @param v The 2D vector.
 * @return The minimum component of the vector.
 */
template <typename T> inline T min_component(const Vector2<T> &v) {
    return min(v.x, v.y);
}

/**
 * Returns the maximum component of a 2D vector.
 * 
 * @param v The 2D vector.
 * @return The maximum component of the vector.
 */
template <typename T> inline T max_component(const Vector2<T> &v) {
    return max(v.x, v.y);
}

/**
 * @brief Returns the index of the minimum dimension of the vector.
 * 
 * @tparam T The type of the vector's components.
 * @param v The vector.
 * @return int The index of the minimum dimension (0 for x, 1 for y).
 */
template <typename T> inline int min_dimension(const Vector2<T> &v) {
    return v.x < v.y ? 0 : 1;
}

/**
 * @brief Returns the index of the maximum dimension of the vector.
 * 
 * @tparam T The type of the vector's components.
 * @param v The vector.
 * @return int The index of the maximum dimension (0 for x, 1 for y).
 */
template <typename T> inline int max_dimension(const Vector2<T> &v) {
    return v.x > v.y ? 0 : 1;
}

/**
 * @brief Returns a new vector with the minimum components of the input vectors.
 * 
 * @tparam T The type of the vector's components.
 * @param v1 The first vector.
 * @param v2 The second vector.
 * @return Vector2<T> The new vector with the minimum components.
 */
template <typename T> inline Vector2<T> min(const Vector2<T> &v1, const Vector2<T> &v2) {
    return Vector2<T>(min(v1.x, v2.x), min(v1.y, v2.y));
}

/**
 * @brief Returns a new vector with the maximum components of the input vectors.
 * 
 * @tparam T The type of the vector's components.
 * @param v1 The first vector.
 * @param v2 The second vector.
 * @return Vector2<T> The new vector with the maximum components.
 */
template <typename T> inline Vector2<T> max(const Vector2<T> &v1, const Vector2<T> &v2) {
    return Vector2<T>(max(v1.x, v2.x), max(v1.y, v2.y));
}

/**
 * @brief Returns a new vector with the components permuted according to the input indices.
 * 
 * @tparam T The type of the vector's components.
 * @param v The vector.
 * @param x The index for the new vector's x component.
 * @param y The index for the new vector's y component.
 * @return Vector2<T> The new vector with permuted components.
 */
template <typename T> inline Vector2<T> permute(const Vector2<T> &v, int x, int y) {
    return Vector2<T>(v[x], v[y]);
}

/**
 * @brief Template class for a 3D vector.
 * 
 * @tparam T The type of the vector's components.
 */
template <typename T> class Vector3 {
public:
    T x, y, z;

    /**
     * @brief Default constructor. Initializes all components to zero.
     */
    Vector3() {x = y = z = 0;}

    /**
     * @brief Constructor with components.
     * @param _x The x component.
     * @param _y The y component.
     * @param _z The z component.
     */
    Vector3(T _x, T _y, T _z) : x(_x), y(_y), z(_z) {Assert(!has_NaNs());}

    /**
     * @brief Access operator.
     * @param i The index of the component to access.
     * @return The i-th component of the vector.
     */
    T operator[] (int i) const {
        if (i == 0) return x;
        else if (i == 1) return y;
        else if (i == 2) return z;
        else throw out_of_range("Index out of range");
    }

    /**
     * @brief Non-const access operator.
     * @param i The index of the component to access.
     * @return A reference to the i-th component of the vector.
     */
    T &operator[] (int i) {
        if (i == 0) return x;
        else if (i == 1) return y;
        else if (i == 2) return z;
        else throw out_of_range("Index out of range");
    }

    /**
     * @brief Checks if the vector has any NaN components.
     * @return true if any component is NaN, false otherwise.
     */
    bool has_NaNs() const {
        return isnan(x) || isnan(y) || isnan(z);
    }

    /**
     * @brief Addition operator.
     * @param other The vector to add.
     * @return The result of the addition.
     */
    Vector3<T> operator+ (const Vector3<T> &other) const {
        return Vector3<T>(x + other.x, y + other.y, z + other.z);
    }

    /**
     * @brief In-place addition operator.
     * @param other The vector to add.
     * @return A reference to this vector after the addition.
     */
    Vector3<T> &operator+= (const Vector3<T> &other) {
        x += other.x; y += other.y; z += other.z;
        return *this;
    }

    /**
     * @brief Subtraction operator.
     * @param other The vector to subtract.
     * @return The result of the subtraction.
     */
    Vector3<T> operator- (const Vector3<T> &other) const {
        return Vector3<T>(x - other.x, y - other.y, z - other.z);
    }

    /**
     * @brief In-place subtraction operator.
     * @param other The vector to subtract.
     * @return A reference to this vector after the subtraction.
     */
    Vector3<T> &operator-= (const Vector3<T> &other) {
        x -= other.x; y -= other.y; z -= other.z;
        return *this;
    }

    /**
     * @brief Negation operator.
     * @return The negation of this vector.
     */
    Vector3<T> operator- () const {
        return Vector3<T>(-x, -y, -z);
    }

    /**
     * @brief Multiplication operator.
     * @param s The scalar to multiply by.
     * @return The result of the multiplication.
     */
    Vector3<T> operator* (T s) const {
        return Vector3<T>(x * s, y * s, z * s);
    }

    /**
     * @brief In-place multiplication operator.
     * @param s The scalar to multiply by.
     * @return A reference to this vector after the multiplication.
     */
    Vector3<T> &operator*= (T s) {
        x *= s; y *= s; z *= s;
        return *this;
    }

    /**
     * @brief Division operator.
     * @param s The scalar to divide by.
     * @return The result of the division.
     */
    Vector3<T> operator/ (T s) const {
        Assert(s != 0);
        float inv = (float)1/s;
        return Vector3<T>(x * inv, y * inv, z * inv);
    }

    /**
     * @brief In-place division operator.
     * @param s The scalar to divide by.
     * @return A reference to this vector after the division.
     */
    Vector3<T> &operator/= (T s) {
        Assert(s != 0);
        float inv = (float)1/s;
        x *= inv; y *= inv; z *= inv;
        return *this;
    }

    /**
     * @brief Computes the squared length of the vector.
     * @return The squared length of the vector.
     */
    float length_squared() const {
        return x * x + y * y + z * z;
    }

    /**
     * @brief Computes the length of the vector.
     * @return The length of the vector.
     */
    float length() const {
        return sqrt(length_squared());
    }

    /**
     * @brief Normalizes the vector.
     * @return The normalized vector.
     */
    Vector3<T> normalize() const {
        return *this / length();
    }
};

/**
 * @brief Multiplies a scalar value with a vector.
 * 
 * This function multiplies a scalar value with a vector and returns the result.
 * 
 * @tparam T The type of the vector elements.
 * @param s The scalar value to multiply.
 * @param v The vector to multiply.
 * @return The result of multiplying the scalar value with the vector.
 */
template <typename T> inline Vector3<T> operator* (T s, const Vector3<T> &v) {
    return v * s;
}

/**
 * Calculates the absolute value of each component of a Vector3.
 * 
 * @param v The input Vector3.
 * @return A new Vector3 with the absolute values of the components.
 */
template <typename T> inline Vector3<T> abs(const Vector3<T> &v) {
    return Vector3<T>(abs(v.x), abs(v.y), abs(v.z));
}

/**
 * Calculates the dot product of two 3D vectors.
 * 
 * @tparam T The type of the vector components.
 * @param v1 The first vector.
 * @param v2 The second vector.
 * @return The dot product of the two vectors.
 */
template <typename T> inline T dot(const Vector3<T> &v1, const Vector3<T> &v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

/** @brief Computes the absolute value of the dot product of two vectors.
 * @param v1 The first vector.
 * @param v2 The second vector.
 * @return The absolute value of the dot product.
 */
template <typename T> inline T abs_dot(const Vector3<T> &v1, const Vector3<T> &v2) {
    return abs(dot(v1, v2));
}

/**
 * @brief Computes the cross product of two vectors.
 * @param v1 The first vector.
 * @param v2 The second vector.
 * @return The cross product of the two vectors.
 */
template <typename T> inline Vector3<T> cross(const Vector3<T> &v1, const Vector3<T> &v2) {
    return Vector3<T>(v1.y * v2.z - v1.z * v2.y,
                      v1.z * v2.x - v1.x * v2.z,
                      v1.x * v2.y - v1.y * v2.x);
}

/**
 * @brief Returns the minimum component of the vector.
 * @param v The vector.
 * @return The minimum component of the vector.
 */
template <typename T> inline T min_component(const Vector3<T> &v) {
    return min(v.x, min(v.y, v.z));
}

/**
 * @brief Returns the maximum component of the vector.
 * @param v The vector.
 * @return The maximum component of the vector.
 */
template <typename T> inline T max_component(const Vector3<T> &v) {
    return max(v.x, max(v.y, v.z));
}

/**
 * @brief Returns the index of the minimum dimension of the vector.
 * @param v The vector.
 * @return The index of the minimum dimension (0 for x, 1 for y, 2 for z).
 */
template <typename T> inline int min_dimension(const Vector3<T> &v) {
    return v.x < v.y ? (v.x < v.z ? 0 : 2) : (v.y < v.z ? 1 : 2);
}

/**
 * @brief Returns the index of the maximum dimension of the vector.
 * @param v The vector.
 * @return The index of the maximum dimension (0 for x, 1 for y, 2 for z).
 */
template <typename T> inline int max_dimension(const Vector3<T> &v) {
    return v.x > v.y ? (v.x > v.z ? 0 : 2) : (v.y > v.z ? 1 : 2);
}

/**
 * @brief Returns a new vector with the minimum components of the input vectors.
 * @param v1 The first vector.
 * @param v2 The second vector.
 * @return The new vector with the minimum components.
 */
template <typename T> inline Vector3<T> min(const Vector3<T> &v1, const Vector3<T> &v2) {
    return Vector3<T>(min(v1.x, v2.x), min(v1.y, v2.y), min(v1.z, v2.z));
}

/**
 * @brief Returns a new vector with the maximum components of the input vectors.
 * @param v1 The first vector.
 * @param v2 The second vector.
 * @return The new vector with the maximum components.
 */
template <typename T> inline Vector3<T> max(const Vector3<T> &v1, const Vector3<T> &v2) {
    return Vector3<T>(max(v1.x, v2.x), max(v1.y, v2.y), max(v1.z, v2.z));
}

/**
 * @brief Returns a new vector with the components permuted according to the input indices.
 * @param v The vector.
 * @param x The index for the new vector's x component.
 * @param y The index for the new vector's y component.
 * @param z The index for the new vector's z component.
 * @return The new vector with permuted components.
 */
template <typename T> inline Vector3<T> permute(const Vector3<T> &v, int x, int y, int z) {
    return Vector3<T>(v[x], v[y], v[z]);
}

// Coordinate system from vector
// Given a vector v, this function will return two vectors u and w such that u, v, w form an orthonormal basis
// This is useful for generating a coordinate system for shading
template <typename T> inline void coordinate_system(const Vector3<T> &v1, Vector3<T> *v2, Vector3<T> *v3) {
    // Create a normalized vector orthogonal to v1
    // We will use zero the largest component and swap the other two to create it
    if (abs(v1.x) > abs(v1.y)) {
        // If v1.x is the largest component, then the orthogonal vector is (-v1.z, 0, v1.x) / sqrt(v1.x * v1.x + v1.z * v1.z)
        *v2 = Vector3<T>(-v1.z, 0, v1.x) / sqrt(v1.x * v1.x + v1.z * v1.z);
    } else {
        // Otherwise, the orthogonal vector is (0, v1.z, -v1.y) / sqrt(v1.y * v1.y + v1.z * v1.z)
        *v2 = Vector3<T>(0, v1.z, -v1.y) / sqrt(v1.y * v1.y + v1.z * v1.z);
    }
    // The third vector is the cross product of the first two
    *v3 = cross(v1, *v2);
}

#endif // !VECTORS