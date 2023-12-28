#ifndef BBOXES
#define BBOXES

#include <iostream>
#include <limits>
#include "points.h"
#include "vectors.h"

using namespace std;

/**
 * @brief A template class for 2D bounding boxes.
 * 
 * @tparam T The type of the coordinates.
*/
template <typename T> class BBox2 {
    public:
        Point2<T> p_min, p_max;

        /**
         * @brief Default constructor. Initializes p_min to infinity and p_max to -infinity.
        */
        BBox2() {
            T min_num = numeric_limits<T>::lowest();
            T max_num = numeric_limits<T>::max();
            p_min = Point2<T>(max_num, max_num);
            p_max = Point2<T>(min_num, min_num);
        }

        /**
         * @brief Constructor with two points.
        */
        BBox2(const Point2<T>& p1, const Point2<T>& p2) :
            p_min(min(p1.x, p2.x), min(p1.y, p2.y)),
            p_max(max(p1.x, p2.x), max(p1.y, p2.y)) {}
        /**
         * @brief Constructor with a single point.
         * @param p The point.
        */
        BBox2(const Point2<T>& p) :
            p_min(p), p_max(p) {}

        /**
         * @brief Constructor with a single value.
         * @param v The value.
        */
        BBox2(T v) :
            p_min(v), p_max(v) {}

        /**
         * @brief Accesses the i-th point of the bounding box.
         * 
         * @param i The index of the point to access.
         * @return const Point2<T>& The i-th point of the bounding box.
         */
        const Point2<T>& operator[](int i) const {
            return (i == 0) ? p_min : p_max;
        }

        /**
         * @brief Accesses the i-th point of the bounding box.
         * 
         * @param i The index of the point to access.
         * @return Point2<T>& The i-th point of the bounding box.
         */
        Point2<T>& operator[](int i) {
            return (i == 0) ? p_min : p_max;
        }

        /**
         * @brief Returns the corner point of the bounding box at the specified index.
         * 
         * @param corner The index of the corner (0-3).
         * @return The corner point.
         */
        Point2<T> corner(int corner) const {
            return Point2<T>((*this)[(corner & 1)].x, (*this)[(corner & 2) ? 1 : 0].y);
        }

        /**
         * @brief Returns the surface area of the bounding box.
         * @return The surface area of the bounding box.
        */
        T surface_area() const {
            Vector2<T> d = p_max - p_min;
            return 2 * (d.x * d.y);
        }

        /**
         * @brief Returns the index of the maximum extent of the bounding box.
         * The maximum extent is determined by comparing the differences between the maximum and minimum coordinates in each dimension.
         * If the difference in the x dimension is greater than the difference in the y dimension, the maximum extent is 0.
         * Otherwise, the maximum extent is 1.
         *
         * @return The index of the maximum extent (0 for x, 1 for y).
         */
        int maximum_extent() const {
            Vector2<T> diag = p_max - p_min;
            if (diag.x > diag.y) {
                return 0;
            }
            else {
                return 1;
            }
        }

        /**
         * @brief Performs linear interpolation between the corners of the bounding box.
         * 
         * @param t The interpolation parameter.
         * @return The interpolated point.
         */
        Point2<T> lerp(const Point2f& t) const {
            return Point2<T>(::lerp(t.x, p_min.x, p_max.x), ::lerp(t.y, p_min.y, p_max.y));
        }

        /**
         * @brief Calculates the offset of a given point within the bounding box.
         * The offset is a normalized value between 0 and 1, representing the relative position of the point within the bounding box.
         * 
         * @param p The point for which to calculate the offset.
         * @return The offset of the point within the bounding box.
         */
        Point2<T> offset(const Point2<T>& p) const {
            Vector2<T> o = p - p_min;
            if (p_max.x > p_min.x) {
                o.x /= p_max.x - p_min.x;
            }
            if (p_max.y > p_min.y) {
                o.y /= p_max.y - p_min.y;
            }
            return o;
        }

        /**
         * @brief Calculates the bounding sphere of the bounding box.
         * 
         * @param center Pointer to a Point2 object to store the center of the bounding sphere.
         * @param radius Pointer to a float variable to store the radius of the bounding sphere.
         */
        void bounding_sphere(Point2<T>* center, float* radius) const {
            *center = (p_min + p_max) / 2;
            *radius = is_inside_bbox(*center, *this) ? distance(*center, p_max) : 0;
        }

        /**
         * @brief Outputs the bounding box to an ostream.
         * @param os The ostream.
         * @param r The bounding box.
         * @return The ostream.
         */
        friend ostream& operator<<(ostream& os, const BBox2<T>& b) {
            os << "BBox2(" << b.p_min << ", " << b.p_max << ")";
            return os;
        }
};

/**
 * @brief Iterator class for iterating over integer positions in a BBox2 object.
 * 
 * This class provides an iterator interface for iterating over the integer points within a BBox2 object.
 * It supports forward iteration and provides access to the x and y coordinates of each point.
 * 
 * @tparam T The type of the coordinates.
 */
template <typename T> class BBox2Iterator {
    public:
        using iterator_category = forward_iterator_tag;
        using value_type = Point2<T>;
        using reference = value_type&;
        using pointer = value_type*;
        using difference_type = ptrdiff_t;

    private:
        BBox2<T>* bounding_box;
        T x;
        T y;

    public:
        /**
         * @brief Constructs a BBox2Iterator object.
         * 
         * @param bounding_box Pointer to the BBox2 object to iterate over.
         */
        BBox2Iterator(BBox2<T>* bounding_box) : bounding_box(bounding_box), x(bounding_box->p_min.x), y(bounding_box->p_min.y) {}

        /**
         * @brief Advances the iterator to the next point.
         * 
         * If the current x coordinate is less than the maximum x coordinate of the bounding box,
         * the x coordinate is incremented. Otherwise, if the current y coordinate is less than the
         * maximum y coordinate of the bounding box, the x coordinate is reset to the minimum x coordinate
         * and the y coordinate is incremented.
         * 
         * @return Reference to the updated iterator.
         */
        BBox2Iterator& operator++() {
            if (x < (T)bounding_box->p_max.x) {
                x++;
                return *this;
            }

            if (y < (T)bounding_box->p_max.y) {
                x = (T)bounding_box->p_min.x;
                y++;
                return *this;
            }

            return *this;
        }

        /**
         * @brief Advances the iterator to the next point.
         * 
         * This postfix increment operator returns a copy of the iterator before it was incremented.
         * 
         * @return Copy of the iterator before incrementing.
         */
        BBox2Iterator operator++(int) {
            BBox2Iterator temp = *this;
            operator++();
            return temp;
        }

        /**
         * @brief Dereferences the iterator and returns a reference to the current point.
         * 
         * @return Reference to the current point.
         */
        reference operator*() const {
            return Point2<T>(x, y);
        }

        /**
         * @brief Returns a pointer to the current point.
         * 
         * @return Pointer to the current point.
         */
        pointer operator->() const {
            return &Point2<T>(x, y);
        }

        /**
         * @brief Checks if two iterators are equal.
         * 
         * Two iterators are considered equal if they are iterating over the same bounding box
         * and have the same current coordinates.
         * 
         * @param other The iterator to compare with.
         * @return True if the iterators are equal, false otherwise.
         */
        bool operator==(const BBox2Iterator& other) const {
            return bounding_box == other.bounding_box && x == other.x && y == other.y;
        }

        /**
         * @brief Checks if two iterators are not equal.
         * 
         * Two iterators are considered not equal if they are not iterating over the same bounding box
         * or have different current coordinates.
         * 
         * @param other The iterator to compare with.
         * @return True if the iterators are not equal, false otherwise.
         */
        bool operator!=(const BBox2Iterator& other) const {
            return !operator==(other);
        }
};


/**
 * @brief Computes the union of a bounding box and a point.
 * @details The union is the smallest bounding box that contains both the bounding box and the point.
 *
 * @tparam T The type of the coordinates of the bounding box and point.
 * @param b The original bounding box.
 * @param p The point to be included in the union.
 * @return The union of the bounding box and the point.
 */
template <typename T> BBox2<T> union_bbox(const BBox2<T>& b, const Point2<T>& p) {
    BBox2<T> ret;
    ret.p_min.x = min(b.p_min.x, p.x);
    ret.p_min.y = min(b.p_min.y, p.y);
    ret.p_max.x = max(b.p_max.x, p.x);
    ret.p_max.y = max(b.p_max.y, p.y);
    return ret;
}

/**
 * @brief Computes the union of two bounding boxes.
 * @details The union is the smallest bounding box that contains both bounding boxes.
 * 
 * @tparam T The type of the bounding box coordinates.
 * @param b1 The first bounding box.
 * @param b2 The second bounding box.
 * @return The union of the two bounding boxes.
 */
template <typename T> BBox2<T> union_bbox(const BBox2<T>& b1, const BBox2<T>& b2) {
    BBox2<T> ret;
    ret.p_min.x = min(b1.p_min.x, b2.p_min.x);
    ret.p_min.y = min(b1.p_min.y, b2.p_min.y);
    ret.p_max.x = max(b1.p_max.x, b2.p_max.x);
    ret.p_max.y = max(b1.p_max.y, b2.p_max.y);
    return ret;
}

/**
 * @brief Computes the intersection of two bounding boxes.
 * @details The intersection is the largest bounding box that is contained by both bounding boxes.
 *
 * @tparam T The type of the coordinates of the bounding boxes.
 * @param b1 The first bounding box.
 * @param b2 The second bounding box.
 * @return The intersection of the two bounding boxes.
 */
template <typename T> BBox2<T> intersect_bbox(const BBox2<T>& b1, const BBox2<T>& b2) {
    BBox2<T> ret;
    ret.p_min.x = max(b1.p_min.x, b2.p_min.x);
    ret.p_min.y = max(b1.p_min.y, b2.p_min.y);
    ret.p_max.x = min(b1.p_max.x, b2.p_max.x);
    ret.p_max.y = min(b1.p_max.y, b2.p_max.y);
    return ret;
}

/**
 * @brief Checks if two bounding boxes overlap.
 * 
 * @tparam T The type of the bounding box coordinates.
 * @param b1 The first bounding box.
 * @param b2 The second bounding box.
 * @return True if the bounding boxes overlap, false otherwise.
 */
template <typename T> bool overlaps_bbox(const BBox2<T>& b1, const BBox2<T>& b2) {
    bool x = (b1.p_max.x >= b2.p_min.x) && (b1.p_min.x <= b2.p_max.x);
    bool y = (b1.p_max.y >= b2.p_min.y) && (b1.p_min.y <= b2.p_max.y);
    return (x && y);
}

/**
 * @brief Checks if a point is inside a bounding box.
 * 
 * @param p The point to check.
 * @param b The bounding box to check against.
 * @return True if the point is inside the bounding box, false otherwise.
 */
template <typename T> bool is_inside_bbox(const Point2<T>& p, const BBox2<T>& b) {
    return (p.x >= b.p_min.x && p.x <= b.p_max.x &&
            p.y >= b.p_min.y && p.y <= b.p_max.y);
}

/**
 * @brief Checks if a point is inside a bounding box (exclusive).
 * \details Exclusive means that the point is not on the boundary of the bounding box.
 * Very useful for integer bounding boxes.
 * 
 * @tparam T The type of the coordinates.
 * @param p The point to check.
 * @param b The bounding box to check against.
 * @return True if the point is inside the bounding box, false otherwise.
 */
template <typename T> bool is_inside_bbox_exclusive(const Point2<T>& p, const BBox2<T>& b) {
    return (p.x > b.p_min.x && p.x < b.p_max.x &&
            p.y > b.p_min.y && p.y < b.p_max.y);
}

/**
 * @brief Expands the bounding box by a given delta value.
 * 
 * @tparam T The type of the bounding box coordinates.
 * @tparam U The type of the delta value.
 * @param b The original bounding box.
 * @param delta The delta value to expand the bounding box.
 * @return The expanded bounding box.
 */
template <typename T, typename U> inline BBox2<T> expand_bbox(const BBox2<T>& b, U delta) {
    return BBox2<T>(b.p_min - Vector2<T>(delta, delta), b.p_max + Vector2<T>(delta, delta));
}

/**
 * @brief A template class for 3D bounding boxes.
 * @tparam T The type of the coordinates.
 */
template <typename T> class BBox3 {
    public:
        Point3<T> p_min, p_max;

        /**
         * @brief Default constructor. Initializes p_min to infinity and p_max to -infinity.
        */
        BBox3() {
            T min_num = numeric_limits<T>::lowest();
            T max_num = numeric_limits<T>::max();
            p_min = Point3<T>(max_num, max_num, max_num);
            p_max = Point3<T>(min_num, min_num, min_num);
        }

        /**
         * @brief Constructor with two points.
        */
        BBox3(const Point3<T>& p1, const Point3<T>& p2) :
            p_min(min(p1.x, p2.x), min(p1.y, p2.y), min(p1.z, p2.z)),
            p_max(max(p1.x, p2.x), max(p1.y, p2.y), max(p1.z, p2.z)) {}
        /**
         * @brief Constructor with a single point.
         * @param p The point.
        */
        BBox3(const Point3<T>& p) :
            p_min(p), p_max(p) {}

        /**
         * @brief Constructor with a single value.
         * @param v The value.
        */
        BBox3(T v) :
            p_min(v), p_max(v) {}

        /**
         * @brief Accesses the i-th point of the bounding box.
         * 
         * @param i The index of the point to access.
         * @return const Point3<T>& The i-th point of the bounding box.
         */
        const Point3<T>& operator[](int i) const {
            return (i == 0) ? p_min : p_max;
        }

        /**
         * @brief Accesses the i-th point of the bounding box.
         * 
         * @param i The index of the point to access.
         * @return Point3<T>& The i-th point of the bounding box.
         */
        Point3<T>& operator[](int i) {
            return (i == 0) ? p_min : p_max;
        }

        /**
         * @brief Returns the corner point of the bounding box at the specified index.
         * 
         * @param corner The index of the corner (0-7).
         * @return The corner point.
         */
        Point3<T> corner(int corner) const {
            return Point3<T>((*this)[(corner & 1)].x, (*this)[(corner & 2) ? 1 : 0].y, (*this)[(corner & 4) ? 1 : 0].z);
        }

        /**
         * @brief Returns the surface area of the bounding box.
         * @return The surface area of the bounding box.
        */
        T surface_area() const {
            Vector3<T> d = p_max - p_min;
            return 2 * (d.x * d.y + d.x * d.z + d.y * d.z);
        }

        /**
         * Calculates the surface volume of the bounding box.
         * 
         * @tparam T The type of the bounding box coordinates.
         * @return The surface volume of the bounding box.
         */
        T surface_volume() const {
            Vector3<T> d = p_max - p_min;
            return d.x * d.y * d.z;
        }

        /**
         * @brief Returns the index of the maximum extent of the bounding box.
         * The maximum extent is determined by comparing the differences between the maximum and minimum coordinates in each dimension.
         * If the difference in the x dimension is greater than the difference in the y and z dimensions, the maximum extent is 0.
         * Otherwise, if the difference in the y dimension is greater than the difference in the z dimension, the maximum extent is 1.
         * Otherwise, the maximum extent is 2.
         *
         * @return The index of the maximum extent (0 for x, 1 for y, 2 for z).
         */
        int maximum_extent() const {
            Vector3<T> diag = p_max - p_min;
            if (diag.x > diag.y && diag.x > diag.z) {
                return 0;
            }
            else if (diag.y > diag.z) {
                return 1;
            }
            else {
                return 2;
            }
        }

        /**
         * @brief Performs linear interpolation between the corners of the bounding box.
         * 
         * @param t The interpolation parameter.
         * @return The interpolated point.
         */
        Point3<T> lerp(const Point3f& t) const {
            return Point3<T>(::lerp(t.x, p_min.x, p_max.x), ::lerp(t.y, p_min.y, p_max.y), ::lerp(t.z, p_min.z, p_max.z));
        }

        /**
         * @brief Calculates the offset of a given point within the bounding box.
         * The offset is a normalized value between 0 and 1, representing the relative position of the point within the bounding box.
         * 
         * @param p The point for which to calculate the offset.
         * @return The offset of the point within the bounding box.
         */
        Point3<T> offset(const Point3<T>& p) const {
            Vector3<T> o = p - p_min;
            if (p_max.x > p_min.x) {
                o.x /= p_max.x - p_min.x;
            }
            if (p_max.y > p_min.y) {
                o.y /= p_max.y - p_min.y;
            }
            if (p_max.z > p_min.z) {
                o.z /= p_max.z - p_min.z;
            }
            return o;
        }

        /**
         * @brief Calculates the bounding sphere of the bounding box.
         * 
         * @param center Pointer to a Point3 object to store the center of the bounding sphere.
         * @param radius Pointer to a float variable to store the radius of the bounding sphere.
         */
        void bounding_sphere(Point3<T>* center, float* radius) const {
            *center = (p_min + p_max) / 2;
            *radius = is_inside_bbox(*center, *this) ? distance(*center, p_max) : 0;
        }

        /**
         * @brief Outputs the bounding box to an ostream.
         * @param os The ostream.
         * @param r The bounding box.
         * @return The ostream.
         */
        friend ostream& operator<<(ostream& os, const BBox3<T>& b) {
            os << "BBox3(" << b.p_min << ", " << b.p_max << ")";
            return os;
        }
};

/**
 * @brief Iterator class for iterating over integer positions in a BBox3 object.
 * 
 * This class provides an iterator interface for iterating over the integer points within a BBox3 object.
 * It supports forward iteration and provides access to the x, y, and z coordinates of each point.
 * 
 * @tparam T The type of the coordinates.
 */
template <typename T> class BBox3Iterator {
    public:
        using iterator_category = forward_iterator_tag;
        using value_type = Point3<T>;
        using reference = value_type&;
        using pointer = value_type*;
        using difference_type = ptrdiff_t;

    private:
        BBox3<T>* bounding_box;
        T x;
        T y;
        T z;

    public:
        /**
         * @brief Constructs a BBox3Iterator object.
         * 
         * @param bounding_box Pointer to the BBox3 object to iterate over.
         */
        BBox3Iterator(BBox3<T>* bounding_box) : bounding_box(bounding_box), x(bounding_box->p_min.x), y(bounding_box->p_min.y), z(bounding_box->p_min.z) {}

        /**
         * @brief Advances the iterator to the next point.
         * 
         * If the current x coordinate is less than the maximum x coordinate of the bounding box,
         * the x coordinate is incremented. Otherwise, if the current y coordinate is less than the
         * maximum y coordinate of the bounding box, the x coordinate is reset to the minimum x coordinate
         * and the y coordinate is incremented. Otherwise, if the current z coordinate is less than the
         * maximum z coordinate of the bounding box, the x coordinate is reset to the minimum x coordinate,
         * the y coordinate is reset to the minimum y coordinate, and the z coordinate is incremented.
         * 
         * @return Reference to the updated iterator.
         */
        BBox3Iterator& operator++() {
            if (x < (T)bounding_box->p_max.x) {
                x++;
                return *this;
            }

            if (y < (T)bounding_box->p_max.y) {
                x = (T)bounding_box->p_min.x;
                y++;
                return *this;
            }

            if (z < (T)bounding_box->p_max.z) {
                x = (T)bounding_box->p_min.x;
                y = (T)bounding_box->p_min.y;
                z++;
                return *this;
            }

            return *this;
        }

        /**
         * @brief Advances the iterator to the next point.
         * 
         * This postfix increment operator returns a copy of the iterator before it was incremented.
         * 
         * @return Copy of the iterator before incrementing.
         */
        BBox3Iterator operator++(int) {
            BBox3Iterator temp = *this;
            operator++();
            return temp;
        }

        /**
         * @brief Dereferences the iterator and returns a reference to the current point.
         * 
         * @return Reference to the current point.
         */
        reference operator*() const {
            return Point3<T>(x, y, z);
        }

        /**
         * @brief Returns a pointer to the current point.
         * 
         * @return Pointer to the current point.
         */
        pointer operator->() const {
            return &Point3<T>(x, y, z);
        }

        /**
         * @brief Checks if two iterators are equal.
         * 
         * Two iterators are considered equal if they are iterating over the same bounding box
         * and have the same current coordinates.
         * 
         * @param other The iterator to compare with.
         * @return True if the iterators are equal, false otherwise.
         */
        bool operator==(const BBox3Iterator& other) const {
            return bounding_box == other.bounding_box && x == other.x && y == other.y && z == other.z;
        }

        /**
         * @brief Checks if two iterators are not equal.
         * 
         * Two iterators are considered not equal if they are not iterating over the same bounding box
         * or have different current coordinates.
         * 
         * @param other The iterator to compare with.
         * @return True if the iterators are not equal, false otherwise.
         */
        bool operator!=(const BBox3Iterator& other) const {
            return !operator==(other);
        }
};

/**
 * @brief Computes the union of a bounding box and a point.
 * @details The union is the smallest bounding box that contains both the bounding box and the point.
 *
 * @tparam T The type of the coordinates of the bounding box and point.
 * @param b The original bounding box.
 * @param p The point to be included in the union.
 * @return The union of the bounding box and the point.
 */
template <typename T> BBox3<T> union_bbox(const BBox3<T>& b, const Point3<T>& p) {
    BBox3<T> ret;
    ret.p_min.x = min(b.p_min.x, p.x);
    ret.p_min.y = min(b.p_min.y, p.y);
    ret.p_min.z = min(b.p_min.z, p.z);
    ret.p_max.x = max(b.p_max.x, p.x);
    ret.p_max.y = max(b.p_max.y, p.y);
    ret.p_max.z = max(b.p_max.z, p.z);
    return ret;
}

/**
 * @brief Computes the union of two bounding boxes.
 * @details The union is the smallest bounding box that contains both bounding boxes.
 * 
 * @tparam T The type of the bounding box coordinates.
 * @param b1 The first bounding box.
 * @param b2 The second bounding box.
 * @return The union of the two bounding boxes.
 */
template <typename T> BBox3<T> union_bbox(const BBox3<T>& b1, const BBox3<T>& b2) {
    BBox3<T> ret;
    ret.p_min.x = min(b1.p_min.x, b2.p_min.x);
    ret.p_min.y = min(b1.p_min.y, b2.p_min.y);
    ret.p_min.z = min(b1.p_min.z, b2.p_min.z);
    ret.p_max.x = max(b1.p_max.x, b2.p_max.x);
    ret.p_max.y = max(b1.p_max.y, b2.p_max.y);
    ret.p_max.z = max(b1.p_max.z, b2.p_max.z);
    return ret;
}

/**
 * @brief Computes the intersection of two bounding boxes.
 * @details The intersection is the largest bounding box that is contained by both bounding boxes.
 *
 * @tparam T The type of the coordinates of the bounding boxes.
 * @param b1 The first bounding box.
 * @param b2 The second bounding box.
 * @return The intersection of the two bounding boxes.
 */
template <typename T> BBox3<T> intersect_bbox(const BBox3<T>& b1, const BBox3<T>& b2) {
    BBox3<T> ret;
    ret.p_min.x = max(b1.p_min.x, b2.p_min.x);
    ret.p_min.y = max(b1.p_min.y, b2.p_min.y);
    ret.p_min.z = max(b1.p_min.z, b2.p_min.z);
    ret.p_max.x = min(b1.p_max.x, b2.p_max.x);
    ret.p_max.y = min(b1.p_max.y, b2.p_max.y);
    ret.p_max.z = min(b1.p_max.z, b2.p_max.z);
    return ret;
}

/**
 * @brief Checks if two bounding boxes overlap.
 * 
 * @tparam T The type of the bounding box coordinates.
 * @param b1 The first bounding box.
 * @param b2 The second bounding box.
 * @return True if the bounding boxes overlap, false otherwise.
 */
template <typename T> bool overlaps_bbox(const BBox3<T>& b1, const BBox3<T>& b2) {
    bool x = (b1.p_max.x >= b2.p_min.x) && (b1.p_min.x <= b2.p_max.x);
    bool y = (b1.p_max.y >= b2.p_min.y) && (b1.p_min.y <= b2.p_max.y);
    bool z = (b1.p_max.z >= b2.p_min.z) && (b1.p_min.z <= b2.p_max.z);
    return (x && y && z);
}

/**
 * @brief Checks if a point is inside a bounding box.
 * 
 * @param p The point to check.
 * @param b The bounding box to check against.
 * @return True if the point is inside the bounding box, false otherwise.
 */
template <typename T> bool is_inside_bbox(const Point3<T>& p, const BBox3<T>& b) {
    return (p.x >= b.p_min.x && p.x <= b.p_max.x &&
            p.y >= b.p_min.y && p.y <= b.p_max.y &&
            p.z >= b.p_min.z && p.z <= b.p_max.z);
}

/**
 * @brief Checks if a point is inside a bounding box (exclusive).
 * \details Exclusive means that the point is not on the boundary of the bounding box.
 * Very useful for integer bounding boxes.
 * 
 * @tparam T The type of the coordinates.
 * @param p The point to check.
 * @param b The bounding box to check against.
 * @return True if the point is inside the bounding box, false otherwise.
 */
template <typename T> bool is_inside_bbox_exclusive(const Point3<T>& p, const BBox3<T>& b) {
    return (p.x > b.p_min.x && p.x < b.p_max.x &&
            p.y > b.p_min.y && p.y < b.p_max.y &&
            p.z > b.p_min.z && p.z < b.p_max.z);
}

/**
 * @brief Expands the bounding box by a given delta value.
 * 
 * @tparam T The type of the bounding box coordinates.
 * @tparam U The type of the delta value.
 * @param b The original bounding box.
 * @param delta The delta value to expand the bounding box.
 * @return The expanded bounding box.
 */
template <typename T, typename U> inline BBox3<T> expand_bbox(const BBox3<T>& b, U delta) {
    return BBox3<T>(b.p_min - Vector3<T>(delta, delta, delta), b.p_max + Vector3<T>(delta, delta, delta));
}

typedef BBox2<int> BBox2i;
typedef BBox2<float> BBox2f;
typedef BBox3<int> BBox3i;
typedef BBox3<float> BBox3f;

#endif // !BBOXES