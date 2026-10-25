// ============================================================================
// Copyright (c) 2026, Anhui Jianzhu University, China
// Author: Jia-Chun Guo
// Contact: jchguo@ahjzu.edu.cn; jchguo@whu.edu.cn
// License: MIT (https://opensource.org/licenses/MIT)
// ============================================================================
#pragma once

// GeodesicDefinitions.hpp
// Shared multiprecision type, WGS84 parameters, angle conventions,
// and elementary geodesic transformations.
// Angles are in radians unless explicitly indicated otherwise.
// Requires C++17 and Boost.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

namespace Geodesic_Defs {

using real = boost::multiprecision::cpp_dec_float_100;

using boost::multiprecision::abs;
using boost::multiprecision::atan2;
using boost::multiprecision::cos;
using boost::multiprecision::floor;
using boost::multiprecision::hypot;
using boost::multiprecision::sin;
using boost::multiprecision::sqrt;

// Construct decimal constants from strings to avoid binary64 rounding.
// These parameters are immutable; cached series tables use this ellipsoid.
inline const real PI = boost::math::constants::pi<real>();
inline const real TWO_PI = real(2) * PI;
inline const real HALF_PI = PI / real(2);
inline const real DEG_TO_RAD = PI / real(180);
inline const real RAD_TO_DEG = real(180) / PI;
inline const real a("6378137.0");
inline const real f = real(1) / real("298.257223563");
inline const real one_minus_f = real(1) - f;
inline const real e2 = f * (real(2) - f);
inline const real e = sqrt(e2);
inline const real n = f / (real(2) - f);
inline const real n2 = n * n;

inline void require_finite(const real& value, const char* message)
{
    if (!boost::math::isfinite(value))
        throw std::invalid_argument(message);
}

inline real deg_to_rad(real degrees) { return degrees * DEG_TO_RAD; }
inline real rad_to_deg(real radians) { return radians * RAD_TO_DEG; }

// Longitudes and signed azimuths use (-PI, PI] or (-180, 180].
// In particular, both signs of an exact half-turn map to the positive endpoint.
inline real normalize_lon_rad(real angle)
{
    require_finite(angle, "Non-finite angle in normalize_lon_rad.");
    angle -= floor(angle / TWO_PI) * TWO_PI;
    if (angle > PI) angle -= TWO_PI;
    return angle;
}

inline real normalize_lon_deg(real angle)
{
    require_finite(angle, "Non-finite angle in normalize_lon_deg.");
    angle -= floor(angle / real(360)) * real(360);
    if (angle > real(180)) angle -= real(360);
    return angle;
}

// Unsigned azimuths use [0, 2*PI) or [0, 360).
inline real normalize_azi_rad(real angle)
{
    require_finite(angle, "Non-finite angle in normalize_azi_rad.");
    angle -= floor(angle / TWO_PI) * TWO_PI;
    // A rounded quotient can fall just below an exact integer.
    if (angle >= TWO_PI) angle -= TWO_PI;
    if (angle < 0) angle += TWO_PI;
    return angle;
}

inline real normalize_azi_deg(real angle)
{
    require_finite(angle, "Non-finite angle in normalize_azi_deg.");
    angle -= floor(angle / real(360)) * real(360);
    if (angle >= real(360)) angle -= real(360);
    if (angle < 0) angle += real(360);
    return angle;
}

template<class T>
inline int sign(const T& value)
{
    return (value > 0) - (value < 0);
}

// Direction sign used by the endpoint-azimuth branch rules.
// At exact cardinal azimuths, 0 and 90 degrees give +1,
// while 180 and 270 degrees give -1.
// Elsewhere, return sign(sin(2*azimuth)).
inline int Sign_sin2_deg(real azimuth)
{
    azimuth = normalize_azi_deg(azimuth);
    if (azimuth == 0 || azimuth == 90) return +1;
    if (azimuth == 180 || azimuth == 270) return -1;
    return (azimuth < 90 || (azimuth > 180 && azimuth < 270)) ? +1 : -1;
}

// Convert numeric DMS components to decimal degrees.
// Any negative component makes the complete angle negative.
inline real dms_to_deg(real degrees, real minutes, real seconds)
{
    const int direction = (degrees < 0 || minutes < 0 || seconds < 0) ? -1 : 1;
    return real(direction) * (abs(degrees) + abs(minutes) / real(60)
        + abs(seconds) / real(3600));
}

// Read decimal degrees, D M, or D M S from one input line. The leading token
// supplies the sign (including "-0"); minutes and seconds must lie in [0, 60).
// Reject trailing text rather than silently accepting a partial number.
inline bool read_dms_space(real& degrees)
{
    std::string line;
    if (!std::getline(std::cin >> std::ws, line)) return false;
    std::istringstream input(line);
    std::vector<std::string> tokens;
    for (std::string token; input >> token;) tokens.push_back(token);
    if (tokens.empty() || tokens.size() > 3) return false;
    real values[3] = {real(0), real(0), real(0)};
    for (std::size_t i = 0; i < tokens.size(); ++i)
    {
        std::istringstream component(tokens[i]);
        if (!(component >> values[i]) || !boost::math::isfinite(values[i]))
            return false;
        component >> std::ws;
        if (!component.eof()) return false;
    }
    if (values[1] < 0 || values[1] >= 60 || values[2] < 0 || values[2] >= 60)
        return false;
    const int direction = tokens.front().front() == '-' ? -1 : 1;
    degrees = real(direction) * (abs(values[0]) + values[1] / real(60)
        + values[2] / real(3600));
    return true;
}

// Reduced latitude beta from geodetic latitude, and its inverse.
// Exact poles remain exact poles under either transformation.
inline real beta_of_lat(real latitude)
{
    if (latitude == HALF_PI || latitude == -HALF_PI || latitude == 0)
        return latitude;
    return atan2(one_minus_f * sin(latitude), cos(latitude));
}

inline real lat_of_beta(real beta)
{
    if (beta == HALF_PI || beta == -HALF_PI || beta == 0) return beta;
    return atan2(sin(beta), one_minus_f * cos(beta));
}

namespace detail {

// q^2 = (cos(beta)-c)*(cos(beta)+c).
// Clamp a negative result to zero only when c exceeds cos(beta)
// by at most 64 times the arithmetic machine epsilon.
// Larger violations raise a domain error.
inline real vertex_q2(const real& cosine, const real& c)
{
    require_finite(c, "Non-finite Clairaut constant.");
    if (c < 0) throw std::domain_error("The Clairaut constant must be nonnegative.");
    real q2 = (cosine - c) * (cosine + c);
    if (q2 < 0)
    {
        const real allowance = real(64) * std::numeric_limits<real>::epsilon();
        if (c - cosine > allowance)
            throw std::domain_error("Clairaut constant outside the admissible interval.");
        q2 = 0;
    }
    return q2;
}

} // namespace detail

// Direct problem, Eq.(22): phi = atan2(|cos(beta)*cos(alpha)|, sin(beta)).
// This evaluates the principal vertex-based amplitude without recovering c.
inline real phi1_of_beta1(real beta, real alpha)
{
    return atan2(abs(cos(beta) * cos(alpha)), sin(beta));
}

// Vertex-based amplitude phi in [0, PI].
// At beta = 0, assign PI/2, including the degenerate case c = 1.
inline real phi_of_beta(real beta, real c)
{
    const real cosine = cos(beta);
    const real q2 = detail::vertex_q2(cosine, c);
    if (beta == 0) return HALF_PI;
    if (beta == HALF_PI) return real(0);
    if (beta == -HALF_PI) return PI;
    return atan2(sqrt(q2), sin(beta));
}

// Direct problem, Eqs.(56)-(57). hypot avoids cancellation in sqrt(1-c^2).
inline real beta2_of_phi(real alpha1, real beta1,
    real phi2, real c)
{
    const real p = hypot(cos(alpha1), sin(alpha1) * sin(beta1));
    const real sine = sin(phi2);
    const real cosine = cos(phi2);
    return atan2(p * cosine, hypot(sine, c * cosine));
}

// Principal terminal azimuth before the directional adjustments
// of the direct-problem branch rules.
inline real Azi2_of_phi(real alpha1, real beta1,
    real phi2, real c)
{
    const real p = hypot(cos(alpha1), sin(alpha1) * sin(beta1));
    return atan2(c, p * sin(phi2));
}

// Principal azimuth in [0, PI/2] from the reduced latitude
// and the nonnegative Clairaut constant.
inline real Azi(real c, real beta)
{
    return atan2(c, sqrt(detail::vertex_q2(cos(beta), c)));
}

// Sum the supplied values using Kahan compensated summation.
template<class T>
inline T KahanSum(std::initializer_list<T> values)
{
    T total = 0;
    T correction = 0;
    for (const T& value : values)
    {
        const T adjusted = value - correction;
        const T updated = total + adjusted;
        correction = (updated - total) - adjusted;
        total = updated;
    }
    return total;
}

} // namespace Geodesic_Defs
