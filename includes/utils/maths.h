#ifndef UTILS_MATHS
#define UTILS_MATHS

// Constants
static const float PI = 3.14159265358979323846f;
static const float INV_PI = 0.31830988618379067154f;
static const float INV_TWO_PI = 0.15915494309189533577f;
static const float INV_FOUR_PI = 0.07957747154594766788f;
static const float PI_OVER_TWO = 1.57079632679489661923f;
static const float PI_OVER_FOUR = 0.78539816339744830961f;
static const float SQRT_TWO = 1.41421356237309504880f;

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
 * Converts degrees to radians.
 * 
 * @param degrees The angle in degrees.
 * @return The angle in radians.
 */
inline float radians(float degrees);

/**
 * Converts an angle from radians to degrees.
 *
 * @param radians The angle in radians.
 * @return The angle in degrees.
 */
inline float degrees(float radians);

#endif // !UTILS_MATHS