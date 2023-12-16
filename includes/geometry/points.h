#ifndef POINTS
#define POINTS

#include <iostream>
#include "vectors.h"

using namespace std;

/**
 * @brief Template class for a 3D point.
 * 
 * @tparam T The type of the components.
 */
template <typename T> class Point2 {
    public:
        T x, y;

        /**
         * @brief Default constructor. Initializes x and y to 0.
         */
        Point2() : { x = y = 0; }

        /**
         * @brief Constructor with components.
         * @param _x The x component.
         * @param _y The y component.
         */
        Point2(T _x, T _y) : x(x), y(y) {Assert(!has_Nans);}

        /**
         * @brief Access operator.
         * @param i The index (0 for x, 1 for y).
         * @return The component at the given index.
         */
        template operator[] (int i) {
            Assert(i >= 0 && i <= 1);
            if (i == 0) return x;
            else if (i == 1) return y;
            else throw "Index out of bounds";
        }

        /**
         * @brief Non-const access operator.
         * @param i The index (0 for x, 1 for y).
         * @return The component at the given index.
         */
        template &operator[] (int i) {
            Assert(i >= 0 && i <= 1);
            if (i == 0) return x;
            else if (i == 1) return y;
            else throw "Index out of bounds";
        }

        /**
         * @brief Check if the point has any NaN components.
         * @return True if either x or y is NaN, false otherwise.
         */
        bool has_Nans() const {
            return isnan(x) || isnan(y);
        }

        /**
         * @brief Explicit conversion to another type.
         * @param p The point to convert.
         */
        template <typename U> explicit Point2(const Point2<U> &p) :
            x((T)p.x), y((T)p.y) { Assert(!has_Nans()); }

        /**
         * @brief Explicit conversion to a vector.
         * @return The vector representation of this point.
         */
        template <typename U> explicit operator Vector2<U>() const {
            return Vector2<U>(x, y);
        }

        /**
         * @brief Explicit conversion to a 3D point.
         * @return The point with z = 0.
         */
        template <typename U> explicit operator Point3<U>() const {
            return Vector3<U>(x, y, 0);
        }
        
        /**
         * @brief Addition operator with a vector.
         * @param v The vector to add.
         * @return The result of the addition.
         */
        Point2<T> operator+(const Vector2<T> &v) const {
            return Point2<T>(x + v.x, y + v.y);
        }

        /**
         * @brief In-place addition operator with a vector.
         * @param v The vector to add.
         * @return This point after the addition.
         */
        Point2<T> &operator+=(const Vector2<T> &v) {
            x += v.x; y += v.y;
            return *this;
        }

        /**
         * @brief Subtraction operator Point-Point -> Vector.
         * @param p The point to subtract.
         * @return The result of the subtraction.
         */
        Vector2<T> operator-(const Point2<T> &p) const {
            return Vector2<T>(x - p.x, y - p.y);
        }

        /**
         * @brief Subtraction operator Point-Vector -> Point.
         * @param v The vector to subtract.
         * @return The result of the subtraction.
         */
        Point2<T> operator-(const Vector2<T> &v) const {
            return Point2<T>(x - v.x, y - v.y);
        }

        /**
         * @brief In-place subtraction operator with a vector.
         * @param v The vector to subtract.
         * @return This point after the subtraction.
         */
        Point2<T> &operator-=(const Vector2<T> &v) {
            x -= v.x; y -= v.y;
            return *this;
        }

        /**
         * @brief In-place subtraction operator with a point.
         * @param p The point to subtract.
         * @return This point after the subtraction.
         */
        Point2<T> &operator-=(const Point2<T> &p) {
            x -= p.x; y -= p.y;
            return *this;
        }
        
        /**
         * @brief Negation operator.
         * @return The negation of this point.
         */
        Point2<T> operator-() const {
            return Point2<T>(-x, -y);
        }

        /**
         * @brief Multiplication operator with a scalar.
         * @param f The scalar to multiply.
         * @return The result of the multiplication.
         */
        Point2<T> operator*(T f) const {
            return Point2<T>(f * x, f * y);
        }

        /**
         * @brief In-place multiplication operator with a scalar.
         * @param f The scalar to multiply.
         * @return This point after the multiplication.
         */
        Point2<T> &operator*=(T f) {
            x *= f; y *= f;
            return *this;
        }

        /**
         * @brief Division operator with a scalar.
         * @param f The scalar to divide.
         * @return The result of the division.
         */
        Point2<T> operator/(T f) const {
            Assert(f != 0);
            float inv = 1.f / f;
            return Point2<T>(x * inv, y * inv);
        }

        /**
         * @brief In-place division operator with a scalar.
         * @param f The scalar to divide.
         * @return This point after the division.
         */
        Point2<T> &operator/=(T f) {
            Assert(f != 0);
            float inv = 1.f / f;
            x *= inv; y *= inv;
            return *this;
        }
};

/**
 * @brief Computes the squared distance between two points.
 * @param p1 The first point.
 * @param p2 The second point.
 * @return The squared distance between the points.
 */
template <typename T> inline float distance_squared(const Point2<T> &p1, const Point2<T> &p2) {
    Vector2<T> v = p1 - p2;
    return v.length_squared();
}

/**
 * @brief Computes the distance between two points.
 * @param p1 The first point.
 * @param p2 The second point.
 * @return The distance between the points.
 */
template <typename T> inline float distance(const Point2<T> &p1, const Point2<T> &p2) {
    return sqrt(distance_squared(p1, p2));
}

/**
 * @brief Performs linear interpolation between two points.
 * @param t The interpolation factor.
 * @param p0 The first point.
 * @param p1 The second point.
 * @return The interpolated point.
 */
template <typename T> inline Point2<T> lerp(float t, const Point2<T> &p0, const Point2<T> &p1) {
    return (1 - t) * p0 + t * p1;
}

/**
 * @brief Computes the point with the minimum components of two points.
 * @param p1 The first point.
 * @param p2 The second point.
 * @return The point with the minimum components.
 */
template <typename T> inline Point2<T> min(const Point2<T> &p1, const Point2<T> &p2) {
    return Point2<T>(min(p1.x, p2.x), min(p1.y, p2.y));
}

/**
 * @brief Computes the point with the maximum components of two points.
 * @param p1 The first point.
 * @param p2 The second point.
 * @return The point with the maximum components.
 */
template <typename T> inline Point2<T> max(const Point2<T> &p1, const Point2<T> &p2) {
    return Point2<T>(max(p1.x, p2.x), max(p1.y, p2.y));
}

/**
 * @brief Computes the floor of each component of a point.
 * @param p The point.
 * @return The point with each component being the floor of the original.
 */
template <typename T> inline Point2<T> floor(const Point2<T> &p) {
    return Point2<T>(floor(p.x), floor(p.y));
}

/**
 * @brief Computes the ceiling of each component of a point.
 * @param p The point.
 * @return The point with each component being the ceiling of the original.
 */
template <typename T> inline Point2<T> ceil(const Point2<T> &p) {
    return Point2<T>(ceil(p.x), ceil(p.y));
}

/**
 * @brief Computes the absolute value of each component of a point.
 * @param p The point.
 * @return The point with each component being the absolute value of the original.
 */
template <typename T> inline Point2<T> abs(const Point2<T> &p) {
    return Point2<T>(abs(p.x), abs(p.y));
}

/**
 * @brief Permutes the components of a point.
 * @param p The point.
 * @param x The index of the x component.
 * @param y The index of the y component.
 * @return The permuted point.
 */
template <typename T> Point2<T> permute(const Point2<T> &p, int x, int y) {
    return Point2<T>(p[x], p[y]);
}

/**
 * @brief Outputs a point to a stream.
 * @param os The output stream.
 * @param p The point.
 * @return The output stream.
 */
template <typename T> inline std::ostream &operator<<(std::ostream &os, const Point2<T> &p) {
    os << "[" << p.x << ", " << p.y << "]";
    return os;
}

/**
 * @brief Template class for a 3D point.
 * 
 * @tparam T The type of the components.
 */
template <typename T> class Point3 {
    public:
        T x, y, z;

        /**
         * @brief Default constructor. Initializes x, y and z to 0.
         */
        Point3() : { x = y = z = 0; }

        /**
         * @brief Constructor with components.
         * @param _x The x component.
         * @param _y The y component.
         * @param _z The z component.
         */
        Point3(T _x, T _y, T _z) : x(x), y(y), z(z) {Assert(!has_Nans);}

        /**
         * @brief Access operator.
         * @param i The index (0 for x, 1 for y, 2 for z).
         * @return The component at the given index.
         */
        template operator[] (int i) {
            Assert(i >= 0 && i <= 2);
            if (i == 0) return x;
            else if (i == 1) return y;
            else if (i == 2) return z;
            else throw "Index out of bounds";
        }

        /**
         * @brief Non-const access operator.
         * @param i The index (0 for x, 1 for y, 2 for z).
         * @return The component at the given index.
         */
        template &operator[] (int i) {
            Assert(i >= 0 && i <= 2);
            if (i == 0) return x;
            else if (i == 1) return y;
            else if (i == 2) return z;
            else throw "Index out of bounds";
        }

        /**
         * @brief Check if the point has any NaN components.
         * @return True if either x, y or z is NaN, false otherwise.
         */
        bool has_Nans() const {
            return isnan(x) || isnan(y) || isnan(z);
        }

        /**
         * @brief Explicit conversion to another type.
         * @param p The point to convert.
         */
        template <typename U> explicit Point3(const Point3<U> &p) :
            x((T)p.x), y((T)p.y), z((T)p.z) { Assert(!has_Nans()); }

        /**
         * @brief Explicit conversion to a vector.
         * @return The vector representation of this point.
         */
        template <typename U> explicit operator Vector3<U>() const {
            return Vector3<U>(x, y, z);
        }

        /**
         * @brief Explicit conversion to a 2D point.
         * @return The projection of this point on the xy plane.
         */
        template <typename U> explicit operator Point2<U>() const {
            return Vector2<U>(x, y);
        }

        /**
         * @brief Addition operator with a vector.
         * @param v The vector to add.
         * @return The result of the addition.
         */
        Point3<T> operator+(const Vector3<T> &v) const {
            return Point3<T>(x + v.x, y + v.y, z + v.z);
        }

        /**
         * @brief In-place addition operator with a vector.
         * @param v The vector to add.
         * @return This point after the addition.
         */
        Point3<T> &operator+=(const Vector3<T> &v) {
            x += v.x; y += v.y; z += v.z;
            return *this;
        }

        /**
         * @brief Subtraction operator Point-Point -> Vector.
         * @param p The point to subtract.
         * @return The result of the subtraction.
         */
        Vector3<T> operator-(const Point3<T> &p) const {
            return Vector3<T>(x - p.x, y - p.y, z - p.z);
        }

        /**
         * @brief Subtraction operator Point-Vector -> Point.
         * @param v The vector to subtract.
         * @return The result of the subtraction.
         */
        Point3<T> operator-(const Vector3<T> &v) const {
            return Point3<T>(x - v.x, y - v.y, z - v.z);
        }

        /**
         * @brief In-place subtraction operator with a vector.
         * @param v The vector to subtract.
         * @return This point after the subtraction.
         */
        Point3<T> &operator-=(const Vector3<T> &v) {
            x -= v.x; y -= v.y; z -= v.z;
            return *this;
        }

        /**
         * @brief In-place subtraction operator with a point.
         * @param p The point to subtract.
         * @return This point after the subtraction.
         */
        Point3<T> &operator-=(const Point3<T> &p) {
            x -= p.x; y -= p.y; z -= p.z;
            return *this;
        }

        /**
         * @brief Negation operator.
         * @return The negation of this point.
         */
        Point3<T> operator-() const {
            return Point3<T>(-x, -y, -z);
        }

        /**
         * @brief Multiplication operator with a scalar.
         * @param f The scalar to multiply.
         * @return The result of the multiplication.
         */
        Point3<T> operator*(T f) const {
            return Point3<T>(f * x, f * y, f * z);
        }

        /**
         * @brief In-place multiplication operator with a scalar.
         * @param f The scalar to multiply.
         * @return This point after the multiplication.
         */
        Point3<T> &operator*=(T f) {
            x *= f; y *= f; z *= f;
            return *this;
        }

        /**
         * @brief Division operator with a scalar.
         * @param f The scalar to divide.
         * @return The result of the division.
         */
        Point3<T> operator/(T f) const {
            Assert(f != 0);
            float inv = 1.f / f;
            return Point3<T>(x * inv, y * inv, z * inv);
        }

        /**
         * @brief In-place division operator with a scalar.
         * @param f The scalar to divide.
         * @return This point after the division.
         */
        Point3<T> &operator/=(T f) {
            Assert(f != 0);
            float inv = 1.f / f;
            x *= inv; y *= inv; z *= inv;
            return *this;
        }
};

/**
 * @brief Computes the squared distance between two points.
 * @param p1 The first point.
 * @param p2 The second point.
 * @return The squared distance between the points.
 */
template <typename T> inline float distance_squared(const Point3<T> &p1, const Point3<T> &p2) {
    Vector3<T> v = p1 - p2;
    return v.length_squared();
}

/**
 * @brief Computes the distance between two points.
 * @param p1 The first point.
 * @param p2 The second point.
 * @return The distance between the points.
 */
template <typename T> inline float distance(const Point3<T> &p1, const Point3<T> &p2) {
    return sqrt(distance_squared(p1, p2));
}

/**
 * @brief Performs linear interpolation between two points.
 * @param t The interpolation factor.
 * @param p0 The first point.
 * @param p1 The second point.
 * @return The interpolated point.
 */
template <typename T> inline Point3<T> lerp(float t, const Point3<T> &p0, const Point3<T> &p1) {
    return (1 - t) * p0 + t * p1;
}

/**
 * @brief Computes the point with the minimum components of two points.
 * @param p1 The first point.
 * @param p2 The second point.
 * @return The point with the minimum components.
 */
template <typename T> inline Point3<T> min(const Point3<T> &p1, const Point3<T> &p2) {
    return Point3<T>(min(p1.x, p2.x), min(p1.y, p2.y), min(p1.z, p2.z));
}

/**
 * @brief Computes the point with the maximum components of two points.
 * @param p1 The first point.
 * @param p2 The second point.
 * @return The point with the maximum components.
 */
template <typename T> inline Point3<T> max(const Point3<T> &p1, const Point3<T> &p2) {
    return Point3<T>(max(p1.x, p2.x), max(p1.y, p2.y), max(p1.z, p2.z));
}

/**
 * @brief Computes the floor of each component of a point.
 * @param p The point.
 * @return The point with each component being the floor of the original.
 */
template <typename T> inline Point3<T> floor(const Point3<T> &p) {
    return Point3<T>(floor(p.x), floor(p.y), floor(p.z));
}

/**
 * @brief Computes the ceiling of each component of a point.
 * @param p The point.
 * @return The point with each component being the ceiling of the original.
 */
template <typename T> inline Point3<T> ceil(const Point3<T> &p) {
    return Point3<T>(ceil(p.x), ceil(p.y), ceil(p.z));
}

/**
 * @brief Computes the absolute value of each component of a point.
 * @param p The point.
 * @return The point with each component being the absolute value of the original.
 */
template <typename T> inline Point3<T> abs(const Point3<T> &p) {
    return Point3<T>(abs(p.x), abs(p.y), abs(p.z));
}

/**
 * @brief Permutes the components of a point.
 * @param p The point.
 * @param x The index of the x component.
 * @param y The index of the y component.
 * @param z The index of the z component.
 * @return The permuted point.
 */
template <typename T> Point3<T> permute(const Point3<T> &p, int x, int y, int z) {
    return Point3<T>(p[x], p[y], p[z]);
}

/**
 * @brief Outputs a point to a stream.
 * @param os The output stream.
 * @param p The point.
 * @return The output stream.
 */
template <typename T> inline std::ostream &operator<<(std::ostream &os, const Point3<T> &p) {
    os << "[" << p.x << ", " << p.y << ", " << p.z << "]";
    return os;
}

typedef Point2<float> Point2f;
typedef Point2<int> Point2i;
typedef Point3<float> Point3f;
typedef Point3<int> Point3i;

#endif // !POINTS