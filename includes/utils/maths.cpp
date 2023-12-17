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