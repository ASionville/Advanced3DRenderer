#include "maths.h"

using namespace std;

/**
 * Clamps a value between a specified range.
 *
 * @param val The value to be clamped.
 * @param low The lower bound of the range.
 * @param high The upper bound of the range.
 * @return The clamped value.
 */
template <typename T, typename U, typename V> inline T clamp(T val, U low, V high) {
    if (val < low) return low;
    if (val > high) return high;
    return val;
}

/**
 * Performs linear interpolation between two values.
 * 
 * @tparam T The type of the values to interpolate.
 * @param t The interpolation factor, ranging from 0 to 1.
 * @param v1 The starting value.
 * @param v2 The ending value.
 * @return The interpolated value between v1 and v2.
 */
template <typename T> inline T lerp(float t, T v1, T v2) {
    return v1 * (1 - t) + v2 * t;
}


/**
 * Converts degrees to radians.
 * 
 * @param degrees The angle in degrees.
 * @return The angle in radians.
 */
inline float radians(float degrees) {
    return degrees * (PI / 180);
}

/**
 * Converts an angle from radians to degrees.
 *
 * @param radians The angle in radians.
 * @return The angle in degrees.
 */
inline float degrees(float radians) {
    return radians * (180 / PI);
}