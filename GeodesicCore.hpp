// ============================================================================
// Copyright (c) 2026, Anhui Jianzhu University, China
// Author: Jia-Chun Guo
// Contact: jchguo@ahjzu.edu.cn
// License: MIT (https://opensource.org/licenses/MIT)
// ============================================================================
#pragma once

// GeodesicCore.hpp
// Shared numerical kernels for direct and inverse geodesic calculations
// on an ellipsoid of revolution. Angles are expressed in radians.
// Provides arc lengths, longitude corrections, amplitude inversion,
// and solvers for the Clairaut constant.
// Arc lengths use elliptic integrals of the second kind by default.
// Polynomial and sine series use Horner's method and Clenshaw recurrence.

#include "GeodesicDefinitions.hpp"
#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_3.hpp>

namespace GeodesicCore {

using Geodesic_Defs::real;
using Geodesic_Defs::PI;
using Geodesic_Defs::a;
using Geodesic_Defs::e;
using Geodesic_Defs::e2;
using Geodesic_Defs::n;  // Third flattening.
using Geodesic_Defs::n2;
using Geodesic_Defs::deg_to_rad;
using Geodesic_Defs::beta_of_lat;
using Geodesic_Defs::lat_of_beta;
using Geodesic_Defs::phi1_of_beta1;
using Geodesic_Defs::phi_of_beta;
using Geodesic_Defs::beta2_of_phi;
using Geodesic_Defs::Azi;
using Geodesic_Defs::Azi2_of_phi;
using Geodesic_Defs::KahanSum;
using boost::math::ellint_1;
using boost::math::ellint_2;
using boost::math::ellint_3;
using boost::multiprecision::floor;
using boost::multiprecision::abs;
using boost::multiprecision::sqrt;
using boost::multiprecision::sign;
using boost::multiprecision::sin;
using boost::multiprecision::cos;
using boost::multiprecision::atan2;
using boost::multiprecision::hypot;

// Constants representing powers of the third flattening n.
inline const real n3 = n2 * n;
inline const real n4 = n2 * n2;
inline const real n5 = n4 * n;
inline const real n6 = n3 * n3;
inline const real n7 = n6 * n;
inline const real n8 = n4 * n4;

namespace detail {

inline real ratio(long long numerator, long long denominator)
{
    return real(numerator) / real(denominator);
}

struct PolynomialValue
{
    real value;
    real derivative; // Derivative with respect to the polynomial argument x.
};

// Evaluate a polynomial using Horner's method with the first count coefficients,
// ordered by ascending powers of x.
// Requires 1 <= count <= N.
template<std::size_t N>
inline real polynomial_value(const std::array<real, N>& coefficients,
    unsigned count, const real& x)
{
    real value = coefficients[count - 1];
    for (unsigned i = count - 1; i > 0; --i)
        value = value * x + coefficients[i - 1];
    return value;
}

// Evaluate the polynomial and its derivative with respect to x
// in a single Horner pass, using the same coefficient ordering.
// Requires 1 <= count <= N.
template<std::size_t N>
inline PolynomialValue polynomial(const std::array<real, N>& coefficients,
    unsigned count, const real& x)
{
    real value = coefficients[count - 1];
    real derivative = 0;
    for (unsigned i = count - 1; i > 0; --i)
    {
        derivative = derivative * x + value;
        value = value * x + coefficients[i - 1];
    }
    return {value, derivative};
}

// Shared series coefficient tables, initialized once per program.
// Longitude contraction epsilon(c): through n^10.
// Longitude harmonics f_j(c) and inversion of epsilon(c): through n^8.
// Arc-length evaluation and inversion: expansions in k1.
// Normalization factor F_E(c): numerator polynomial in k1^2.
struct SeriesTable
{
    std::array<real, 10> epsilon{};
    std::array<std::array<real, 8>, 8> harmonic{};
    std::array<real, 8> inverse_epsilon{};
    std::array<std::array<real, 6>, 8> forward_arc{};
    std::array<std::array<real, 6>, 8> inverse_arc{};
    std::array<real, 7> fe_numerator{};

    SeriesTable()
    {
        std::array<real, 11> np{};
        np[0] = 1;
        for (std::size_t i = 1; i < np.size(); ++i) np[i] = np[i - 1] * n;
        epsilon[0] = PI * (np[10] * ratio(-82591, 16384)
            + np[9] * ratio(40193, 8192)
            + np[8] * ratio(-4871, 1024)
            + np[7] * ratio(587, 128)
            + np[6] * ratio(-281, 64)
            + np[5] * ratio(133, 32)
            + np[4] * ratio(-31, 8)
            + np[3] * ratio(7, 2)
            + np[2] * real(-3)
            + np[1] * real(2));
        epsilon[1] = PI * (np[10] * ratio(675581, 16384)
            + np[9] * ratio(-34001, 1024)
            + np[8] * ratio(26651, 1024)
            + np[7] * ratio(-1261, 64)
            + np[6] * ratio(911, 64)
            + np[5] * ratio(-77, 8)
            + np[4] * ratio(47, 8)
            + np[3] * real(-3)
            + np[2] * real(1));
        epsilon[2] = PI * (np[10] * ratio(-1632483, 4096)
            + np[9] * ratio(527607, 2048)
            + np[8] * ratio(-161271, 1024)
            + np[7] * ratio(11457, 128)
            + np[6] * ratio(-1473, 32)
            + np[5] * ratio(327, 16)
            + np[4] * ratio(-57, 8)
            + np[3] * ratio(3, 2));
        epsilon[3] = PI * (np[10] * ratio(9422165, 4096)
            + np[9] * ratio(-1204295, 1024)
            + np[8] * ratio(562115, 1024)
            + np[7] * ratio(-7255, 32)
            + np[6] * ratio(2515, 32)
            + np[5] * ratio(-165, 8)
            + np[4] * ratio(25, 8));
        epsilon[4] = PI * (np[10] * ratio(-65821385, 8192)
            + np[9] * ratio(13013315, 4096)
            + np[8] * ratio(-1120525, 1024)
            + np[7] * ratio(39725, 128)
            + np[6] * ratio(-4165, 64)
            + np[5] * ratio(245, 32));
        epsilon[5] = PI * (np[10] * ratio(143275923, 8192)
            + np[9] * ratio(-5284503, 1024)
            + np[8] * ratio(1266489, 1024)
            + np[7] * ratio(-13797, 64)
            + np[6] * ratio(1323, 64));
        epsilon[6] = PI * (np[10] * ratio(-97724319, 4096)
            + np[9] * ratio(10130967, 2048)
            + np[8] * ratio(-754677, 1024)
            + np[7] * ratio(7623, 128));
        epsilon[7] = PI * (np[10] * ratio(81162081, 4096)
            + np[9] * ratio(-2637921, 1024)
            + np[8] * ratio(184041, 1024));
        epsilon[8] = PI * (np[10] * ratio(-149993415, 16384)
            + np[9] * ratio(4601025, 8192));
        epsilon[9] = PI * (np[10] * ratio(29548805, 16384));
        harmonic[0][0] = np[8] * ratio(435, 32768)
            + np[7] * ratio(-243, 16384)
            + np[6] * ratio(21, 1024)
            + np[5] * ratio(-3, 128)
            + np[4] * ratio(5, 128)
            + np[3] * ratio(-3, 64)
            + np[2] * ratio(1, 8)
            + np[1] * ratio(-1, 4);
        harmonic[0][1] = np[8] * ratio(-33331, 32768)
            + np[7] * ratio(7785, 8192)
            + np[6] * ratio(-895, 1024)
            + np[5] * ratio(51, 64)
            + np[4] * ratio(-89, 128)
            + np[3] * ratio(19, 32)
            + np[2] * ratio(-3, 8);
        harmonic[0][2] = np[8] * ratio(563683, 32768)
            + np[7] * ratio(-198013, 16384)
            + np[6] * ratio(4077, 512)
            + np[5] * ratio(-19, 4)
            + np[4] * ratio(307, 128)
            + np[3] * ratio(-51, 64);
        harmonic[0][3] = np[8] * ratio(-3742507, 32768)
            + np[7] * ratio(239599, 4096)
            + np[6] * ratio(-13331, 512)
            + np[5] * ratio(589, 64)
            + np[4] * ratio(-255, 128);
        harmonic[0][4] = np[8] * ratio(11971921, 32768)
            + np[7] * ratio(-2109245, 16384)
            + np[6] * ratio(35537, 1024)
            + np[5] * ratio(-701, 128);
        harmonic[0][5] = np[8] * ratio(-19678337, 32768)
            + np[7] * ratio(1066345, 8192)
            + np[6] * ratio(-16411, 1024);
        harmonic[0][6] = np[8] * ratio(15990169, 32768)
            + np[7] * ratio(-803251, 16384);
        harmonic[0][7] = np[8] * ratio(-5080225, 32768);
        harmonic[1][0] = np[8] * ratio(287, 32768)
            + np[7] * ratio(-187, 16384)
            + np[6] * ratio(27, 2048)
            + np[5] * ratio(-5, 256)
            + np[4] * ratio(3, 128)
            + np[3] * ratio(-3, 64)
            + np[2] * ratio(1, 16);
        harmonic[1][1] = np[8] * ratio(-10329, 16384)
            + np[7] * ratio(9425, 16384)
            + np[6] * ratio(-265, 512)
            + np[5] * ratio(113, 256)
            + np[4] * ratio(-23, 64)
            + np[3] * ratio(13, 64);
        harmonic[1][2] = np[8] * ratio(304541, 32768)
            + np[7] * ratio(-50651, 8192)
            + np[6] * ratio(3809, 1024)
            + np[5] * ratio(-483, 256)
            + np[4] * ratio(79, 128);
        harmonic[1][3] = np[8] * ratio(-425467, 8192)
            + np[7] * ratio(194053, 8192)
            + np[6] * ratio(-4401, 512)
            + np[5] * ratio(487, 256);
        harmonic[1][4] = np[8] * ratio(4291513, 32768)
            + np[7] * ratio(-599951, 16384)
            + np[6] * ratio(12299, 2048);
        harmonic[1][5] = np[8] * ratio(-2465457, 16384)
            + np[7] * ratio(317733, 16384);
        harmonic[1][6] = np[8] * ratio(2092939, 32768);
        harmonic[2][0] = np[8] * ratio(243, 32768)
            + np[7] * ratio(-139, 16384)
            + np[6] * ratio(3, 256)
            + np[5] * ratio(-7, 512)
            + np[4] * ratio(3, 128)
            + np[3] * ratio(-5, 192);
        harmonic[2][1] = np[8] * ratio(-14141, 32768)
            + np[7] * ratio(4709, 12288)
            + np[6] * ratio(-41, 128)
            + np[5] * ratio(193, 768)
            + np[4] * ratio(-17, 128);
        harmonic[2][2] = np[8] * ratio(264185, 49152)
            + np[7] * ratio(-26561, 8192)
            + np[6] * ratio(1261, 768)
            + np[5] * ratio(-271, 512);
        harmonic[2][3] = np[8] * ratio(-388881, 16384)
            + np[7] * ratio(107509, 12288)
            + np[6] * ratio(-47, 24);
        harmonic[2][4] = np[8] * ratio(1374599, 32768)
            + np[7] * ratio(-115563, 16384);
        harmonic[2][5] = np[8] * ratio(-824769, 32768);
        harmonic[3][0] = np[8] * ratio(193, 32768)
            + np[7] * ratio(-127, 16384)
            + np[6] * ratio(9, 1024)
            + np[5] * ratio(-7, 512)
            + np[4] * ratio(7, 512);
        harmonic[3][1] = np[8] * ratio(-2477, 8192)
            + np[7] * ratio(4071, 16384)
            + np[6] * ratio(-97, 512)
            + np[5] * ratio(49, 512);
        harmonic[3][2] = np[8] * ratio(48321, 16384)
            + np[7] * ratio(-24377, 16384)
            + np[6] * ratio(485, 1024);
        harmonic[3][3] = np[8] * ratio(-74691, 8192)
            + np[7] * ratio(33633, 16384);
        harmonic[3][4] = np[8] * ratio(273437, 32768);
        harmonic[4][0] = np[8] * ratio(179, 32768)
            + np[7] * ratio(-99, 16384)
            + np[6] * ratio(9, 1024)
            + np[5] * ratio(-21, 2560);
        harmonic[4][1] = np[8] * ratio(-6583, 32768)
            + np[7] * ratio(1227, 8192)
            + np[6] * ratio(-75, 1024);
        harmonic[4][2] = np[8] * ratio(45125, 32768)
            + np[7] * ratio(-7107, 16384);
        harmonic[4][3] = np[8] * ratio(-70753, 32768);
        harmonic[5][0] = np[8] * ratio(143, 32768)
            + np[7] * ratio(-99, 16384)
            + np[6] * ratio(11, 2048);
        harmonic[5][1] = np[8] * ratio(-6017, 49152)
            + np[7] * ratio(957, 16384);
        harmonic[5][2] = np[8] * ratio(13211, 32768);
        harmonic[6][0] = np[8] * ratio(143, 32768)
            + np[7] * ratio(-429, 114688);
        harmonic[6][1] = np[8] * ratio(-1573, 32768);
        harmonic[7][0] = np[8] * ratio(715, 262144);
        inverse_epsilon[0] = np[8] * ratio(-3651, 2048)
            + np[7] * ratio(3329, 2048)
            + np[6] * ratio(-375, 256)
            + np[5] * ratio(83, 64)
            + np[4] * ratio(-9, 8)
            + np[3] * ratio(15, 16)
            + np[2] * ratio(-3, 4)
            + np[1] * ratio(1, 2);
        inverse_epsilon[1] = np[8] * ratio(9759, 1024)
            + np[7] * ratio(-12465, 2048)
            + np[6] * ratio(57, 16)
            + np[5] * ratio(-117, 64)
            + np[4] * ratio(3, 4)
            + np[3] * ratio(-3, 16);
        inverse_epsilon[2] = np[8] * ratio(-123771, 4096)
            + np[7] * ratio(28929, 2048)
            + np[6] * ratio(-1473, 256)
            + np[5] * ratio(123, 64)
            + np[4] * ratio(-15, 32)
            + np[3] * ratio(1, 16);
        inverse_epsilon[3] = np[8] * ratio(395355, 8192)
            + np[7] * ratio(-32375, 2048)
            + np[6] * ratio(1065, 256)
            + np[5] * ratio(-25, 32)
            + np[4] * ratio(5, 64);
        inverse_epsilon[4] = np[8] * ratio(-176025, 4096)
            + np[7] * ratio(19185, 2048)
            + np[6] * ratio(-369, 256)
            + np[5] * ratio(15, 128);
        inverse_epsilon[5] = np[8] * ratio(43995, 2048)
            + np[7] * ratio(-5733, 2048)
            + np[6] * ratio(49, 256);
        inverse_epsilon[6] = np[8] * ratio(-23095, 4096)
            + np[7] * ratio(681, 2048);
        inverse_epsilon[7] = np[8] * ratio(9909, 16384);
        forward_arc[0][0] = ratio(1, 2);
        forward_arc[0][1] = ratio(-3, 16);
        forward_arc[0][2] = ratio(1, 32);
        forward_arc[0][3] = ratio(-19, 2048);
        forward_arc[0][4] = ratio(3, 4096);
        forward_arc[0][5] = ratio(-53, 65536);
        forward_arc[1][0] = ratio(-1, 16);
        forward_arc[1][1] = ratio(1, 32);
        forward_arc[1][2] = ratio(-9, 2048);
        forward_arc[1][3] = ratio(7, 4096);
        forward_arc[1][4] = ratio(1, 65536);
        forward_arc[1][5] = ratio(27, 131072);
        forward_arc[2][0] = ratio(1, 48);
        forward_arc[2][1] = ratio(-3, 256);
        forward_arc[2][2] = ratio(3, 2048);
        forward_arc[2][3] = ratio(-17, 24576);
        forward_arc[2][4] = ratio(-3, 65536);
        forward_arc[3][0] = ratio(-5, 512);
        forward_arc[3][1] = ratio(3, 512);
        forward_arc[3][2] = ratio(-11, 16384);
        forward_arc[3][3] = ratio(3, 8192);
        forward_arc[3][4] = ratio(651, 16777216);
        forward_arc[4][0] = ratio(7, 1280);
        forward_arc[4][1] = ratio(-7, 2048);
        forward_arc[4][2] = ratio(3, 8192);
        forward_arc[4][3] = ratio(-117, 524288);
        forward_arc[5][0] = ratio(-7, 2048);
        forward_arc[5][1] = ratio(9, 4096);
        forward_arc[5][2] = ratio(-117, 524288);
        forward_arc[5][3] = ratio(467, 3145728);
        forward_arc[6][0] = ratio(33, 14336);
        forward_arc[6][1] = ratio(-99, 65536);
        forward_arc[6][2] = ratio(77, 524288);
        forward_arc[7][0] = ratio(-429, 262144);
        forward_arc[7][1] = ratio(143, 131072);
        forward_arc[7][2] = ratio(-429, 4194304);
        inverse_arc[0][0] = ratio(-1, 2);
        inverse_arc[0][1] = ratio(9, 32);
        inverse_arc[0][2] = ratio(-205, 1536);
        inverse_arc[0][3] = ratio(4879, 73728);
        inverse_arc[0][4] = ratio(-9039, 327680);
        inverse_arc[0][5] = ratio(1050467, 88473600);
        inverse_arc[1][0] = ratio(5, 16);
        inverse_arc[1][1] = ratio(-37, 96);
        inverse_arc[1][2] = ratio(1335, 4096);
        inverse_arc[1][3] = ratio(-86171, 368640);
        inverse_arc[1][4] = ratio(4119073, 28311552);
        inverse_arc[1][5] = ratio(-18357853, 220200960);
        inverse_arc[2][0] = ratio(-29, 96);
        inverse_arc[2][1] = ratio(75, 128);
        inverse_arc[2][2] = ratio(-2901, 4096);
        inverse_arc[2][3] = ratio(443327, 655360);
        inverse_arc[2][4] = ratio(-1152507, 2097152);
        inverse_arc[3][0] = ratio(539, 1536);
        inverse_arc[3][1] = ratio(-2391, 2560);
        inverse_arc[3][2] = ratio(1082857, 737280);
        inverse_arc[3][3] = ratio(-2722891, 1548288);
        inverse_arc[3][4] = ratio(6190623251LL, 3523215360LL);
        inverse_arc[4][0] = ratio(-3467, 7680);
        inverse_arc[4][1] = ratio(28223, 18432);
        inverse_arc[4][2] = ratio(-1361343, 458752);
        inverse_arc[4][3] = ratio(211942939, 49545216);
        inverse_arc[5][0] = ratio(38081, 61440);
        inverse_arc[5][1] = ratio(-733437, 286720);
        inverse_arc[5][2] = ratio(10820079, 1835008);
        inverse_arc[5][3] = ratio(-547525831, 55050240);
        inverse_arc[6][0] = ratio(-459485, 516096);
        inverse_arc[6][1] = ratio(709743, 163840);
        inverse_arc[6][2] = ratio(-983638957, 84934656);
        inverse_arc[7][0] = ratio(109167851, 82575360);
        inverse_arc[7][1] = ratio(-550835669, 74317824);
        inverse_arc[7][2] = ratio(29797006823LL, 1321205760);

        // Numerator coefficients for the k1-series approximation to F_E(c).
        fe_numerator = {real(1), ratio(1, 4), ratio(1, 64), ratio(1, 256),
            ratio(25, 16384), ratio(49, 65536), ratio(441, 1048576)};
    }
};

inline const SeriesTable& series_table()
{
    static const SeriesTable table;
    return table;
}

// Evaluate sum_{j=1}^N coefficients[j-1]*sin(2*j*phi)
// using Clenshaw recurrence.
// Return exact zero at phi = 0, PI/2, and PI.
// The caller handles any linear term in the unreduced amplitude phi.
template<std::size_t N>
inline real sine_series(const std::array<real, N>& coefficients, const real& phi)
{
    if (phi == 0 || phi == Geodesic_Defs::HALF_PI || phi == PI) return real(0);
    const real twice_phi = real(2) * phi;
    const real multiplier = real(2) * cos(twice_phi);
    real b1 = 0, b2 = 0;
    for (std::size_t i = N; i > 0; --i)
    {
        const real b0 = coefficients[i - 1] + multiplier * b1 - b2;
        b2 = b1;
        b1 = b0;
    }
    return sin(twice_phi) * b1;
}

// Squared elliptic modulus k^2 = e2*(1-c^2)/(1-e2*c^2).
// Requires 0 <= c <= 1.
inline real modulus_squared(const real& c)
{
    Geodesic_Defs::require_finite(c, "Non-finite Clairaut constant.");
    if (c < 0 || c > 1)
        throw std::domain_error("The Clairaut constant must lie in [0, 1].");
    return e2 * (real(1) - c) * (real(1) + c) / (real(1) - e2 * c * c);
}

// Transformed modulus k1 = (1-sqrt(1-k^2))/(1+sqrt(1-k^2)).
inline real landen_parameter(const real& c)
{
    const real k2 = modulus_squared(c);
    const real denominator = real(1) + sqrt(real(1) - k2);
    // Rationalized form avoids subtracting nearly equal numbers near c = 1.
    return k2 / (denominator * denominator);
}

// Evaluate a_j(k1) = k1^j*P_j(k1^2), j = 1,...,8,
// from the selected forward or inverse arc-length coefficient table.
inline std::array<real, 8> arc_coefficients(const real& k1,
    const std::array<std::array<real, 6>, 8>& table)
{
    constexpr unsigned counts[8] = {6, 6, 5, 5, 4, 4, 3, 3};
    const real k12 = k1 * k1;
    real power = k1;
    std::array<real, 8> coefficients{};
    for (unsigned j = 0; j < 8; ++j)
    {
        coefficients[j] = power * polynomial_value(table[j], counts[j], k12);
        power *= k1;
    }
    return coefficients;
}

} // namespace detail

// Normalization factor F_E(c) = 2*E(k)/PI, Eq.(33).
inline real FE(real c)
{
    return real(2) * ellint_2(sqrt(detail::modulus_squared(c))) / PI;
}

// Arc-length evaluation for a fixed Clairaut constant c.
// Stores a' = a*sqrt(1-e2*c*c) and the elliptic modulus k.
// Arc lengths are returned in metres.
// half_period: s_pi(c) = 2*a'*E(k), Eq.(24).
struct ArcLength
{
    real scale;
    real modulus;

    explicit ArcLength(const real& c)
        : scale(a * sqrt(real(1) - e2 * c * c)),
          modulus(sqrt(detail::modulus_squared(c))) {}

    real operator()(const real& phi) const { return scale * ellint_2(modulus, phi); }
    real half_period() const { return real(2) * scale * ellint_2(modulus); }
};

// Arc length s(phi,c) = a'*E(phi,k), Eq.(18).
// The amplitude phi is in radians and may include complete periods.
inline real s_of_phi(real phi, real c) { return ArcLength(c)(phi); }

// Optional series approximation to F_E(c), Eq.(33).
inline real FE_series(real c)
{
    const real k1 = detail::landen_parameter(c);
    return detail::polynomial_value(detail::series_table().fe_numerator, 7, real(k1 * k1))
        / (real(1) + k1);
}

// Eq.(31), optional series approximation to s(phi,c), Eq.(18).
// The linear phi term retains accumulated complete periods.
inline real s_of_phi_series(real phi, real c)
{
    const real k1 = detail::landen_parameter(c);
    const auto coefficients = detail::arc_coefficients(k1, detail::series_table().forward_arc);
    const real factor = detail::polynomial_value(detail::series_table().fe_numerator, 7,
        real(k1 * k1)) / (real(1) + k1);
    return a * sqrt(real(1) - e2 * c * c) * factor
        * (phi + detail::sine_series(coefficients, phi));
}

// Recover phi from the normalized arc length phi_prime, Eq.(32).
// The linear term retains accumulated complete periods.
inline real phi_from_phi_prime(real phi_prime, real c)
{
    const real k1 = detail::landen_parameter(c);
    const auto coefficients = detail::arc_coefficients(k1, detail::series_table().inverse_arc);
    return phi_prime + detail::sine_series(coefficients, phi_prime);
}

// DL(phi,c) = atan2(sin(phi), c*cos(phi)), Eq.(36).
// Returns the principal atan2 value.
inline real DL(real phi, real c) { return atan2(sin(phi), c * cos(phi)); }

// Longitude contraction epsilon(c), Eqs.(41) and (A-4).
inline real Epsilon(real c)
{
    return c * detail::polynomial_value(detail::series_table().epsilon, 10, real(c * c));
}

// Evaluate epsilon(c) and its derivative with respect to c
// from the finite series used for longitude evaluation.
struct EpsilonEvaluation
{
    real value;
    real derivative; // With respect to c.
};

inline EpsilonEvaluation Epsilon_with_derivative(const real& c)
{
    const real x = c * c;
    const auto p = detail::polynomial(detail::series_table().epsilon, 10, x);
    return {c * p.value, p.value + real(2) * x * p.derivative};
}

// Longitude contraction and harmonic coefficients at a fixed c,
// shared by evaluations at both endpoints.
struct LongitudeCoefficients
{
    real epsilon;
    std::array<real, 8> harmonic{};

    explicit LongitudeCoefficients(const real& c)
    {
        const auto& table = detail::series_table();
        const real x = c * c;
        const real g = (real(1) - c) * (real(1) + c);
        epsilon = c * detail::polynomial_value(table.epsilon, 10, x);
        real power = g;
        for (unsigned j = 0; j < 8; ++j)
        {
            harmonic[j] = power * detail::polynomial_value(table.harmonic[j], 8 - j, x);
            power *= g;
        }
    }
};

// Longitude contraction, harmonic coefficients, and their derivatives
// with respect to c, evaluated from the same finite expansions.
struct LongitudeCoefficientsWithDerivative
{
    real epsilon;
    real epsilon_derivative;
    std::array<real, 8> harmonic{};
    std::array<real, 8> harmonic_derivative{};

    explicit LongitudeCoefficientsWithDerivative(const real& c)
    {
        const auto& table = detail::series_table();
        const real x = c * c;
        const real g = (real(1) - c) * (real(1) + c);
        const auto p = detail::polynomial(table.epsilon, 10, x);
        epsilon = c * p.value;
        epsilon_derivative = p.value + real(2) * x * p.derivative;
        real previous_power = 1;
        for (unsigned j = 0; j < 8; ++j)
        {
            const real power = previous_power * g;
            const auto h = detail::polynomial(table.harmonic[j], 8 - j, x);
            harmonic[j] = power * h.value;
            harmonic_derivative[j] = real(2) * c
                * (power * h.derivative - real(j + 1) * previous_power * h.value);
            previous_power = power;
        }
    }
};

// Longitude correction dL(phi,c), Eq.(40), using precomputed coefficients.
inline real dL_with_coefficients(const real& phi, const LongitudeCoefficients& coefficients)
{
    return coefficients.epsilon / PI
        * (phi + detail::sine_series(coefficients.harmonic, phi));
}

// Longitude correction dL(phi,c), Eq.(40).
inline real dL(real phi, real c) { return dL_with_coefficients(phi, LongitudeCoefficients(c)); }

// Recover c for a symmetric half geodesic from the nonnegative
// longitude difference dlon12, expressed in radians.
// Invert epsilon(c) = PI - dlon12 using the finite inverse n-series.
inline real c_from_Epsilon(real dlon12)
{
    static const real normalization = real(2)
        * (ellint_2(e) - (real(1) - e2) * ellint_1(e));
    const real delta = abs(PI - dlon12) / normalization;
    const real delta2 = delta * delta;
    const real correction = detail::polynomial_value(
        detail::series_table().inverse_epsilon, 8, delta2);
    return abs(delta - delta * delta2 * correction);
}

// Geometric branches of the longitude-matching equation.
enum class ClairautCase
{
    A_same_side = 1,
    B_north_vertex = 2,
    C_south_vertex = 3
};

// Clairaut constant and selected geometric branch.
struct ClairautResult
{
    real c;
    ClairautCase which;
};

// Best available Clairaut estimate and its convergence information.
// resid is the smallest evaluated absolute longitude residual, in radians.
// Convergence requires resid <= tolerance.
struct ClairautIterationReport
{
    real c;
    ClairautCase which;
    int iters = 0;
    real resid = std::numeric_limits<real>::infinity(); // Radians.
    bool converged = false;
};

namespace detail {

// Sine and cosine of an endpoint's reduced latitude.
struct ReducedLatitude
{
    real sine;
    real cosine;

    explicit ReducedLatitude(const real& beta)
        : sine(beta == Geodesic_Defs::HALF_PI ? real(1)
            : beta == -Geodesic_Defs::HALF_PI ? real(-1) : real(sin(beta))),
          cosine(abs(beta) == Geodesic_Defs::HALF_PI ? real(0) : real(cos(beta))) {}

    real amplitude(const real& c) const
    {
        const real q2 = Geodesic_Defs::detail::vertex_q2(cosine, c);
        if (sine == 0) return Geodesic_Defs::HALF_PI;
        return atan2(sqrt(q2), sine);
    }
};

// Endpoint geometry shared by the Clairaut updates.
struct ClairautGeometry
{
    ReducedLatitude first;
    ReducedLatitude second;
    real sine_sum;
    real sine_difference;
    real c_max;

    ClairautGeometry(const real& beta1, const real& beta2)
        : first(beta1), second(beta2), sine_sum(sin(beta2 + beta1)),
          sine_difference(sin(beta2 - beta1)),
          c_max((std::min)(first.cosine, second.cosine)) {}

    real from_D(const real& D) const
    {
        const real sine_D = sin(D);
        const real u = second.cosine * sine_D;
        real v;
        if (abs(D) <= Geodesic_Defs::HALF_PI)
        {
            const real sine_half = sin(D / real(2));
            v = sine_difference + real(2) * first.sine * second.cosine
                * sine_half * sine_half;
        }
        else
        {
            const real cosine_half = cos(D / real(2));
            v = sine_sum - real(2) * first.sine * second.cosine
                * cosine_half * cosine_half;
        }
        // Norm form of Eq.(80), using half-angle identities to reduce
        // cancellation near D = 0 and |D| = PI.
        const real denominator = hypot(u, v);
        if (denominator == 0) return real(0);
        const real c = abs(first.cosine * u) / denominator;
        return (std::min)(c, c_max);
    }
};

struct LongitudeState
{
    real epsilon;
    real dL1;
    real dL2;
    std::array<real, 3> longitude;
};

inline LongitudeState longitude_state(const ClairautGeometry& geometry, const real& c)
{
    const LongitudeCoefficients coefficients(c);
    const real phi1 = geometry.first.amplitude(c);
    const real phi2 = geometry.second.amplitude(c);
    const real DL1 = DL(phi1, c);
    const real DL2 = DL(phi2, c);
    const real dL1 = dL_with_coefficients(phi1, coefficients);
    const real dL2 = dL_with_coefficients(phi2, coefficients);
    return {coefficients.epsilon, dL1, dL2,
        {abs(DL1 - DL2) - abs(dL1 - dL2),
         abs(DL1 + DL2) - abs(dL1 + dL2),
         real(2) * (PI - coefficients.epsilon) - abs(DL1 + DL2) + abs(dL1 + dL2)}};
}

inline ClairautIterationReport iterate_clairaut(real initial, const real& beta1,
    const real& beta2, const real& target, int maximum, const real& tolerance,
    bool same_side_only)
{
    if (!boost::math::isfinite(initial) || !boost::math::isfinite(target)
        || !boost::math::isfinite(beta1) || !boost::math::isfinite(beta2)
        || !boost::math::isfinite(tolerance) || tolerance <= 0 || maximum < 0
        || target < 0 || target > PI
        || abs(beta1) >= Geodesic_Defs::HALF_PI
        || abs(beta2) >= Geodesic_Defs::HALF_PI)
        throw std::invalid_argument("Invalid Clairaut fixed-point iteration input.");
    const ClairautGeometry geometry(beta1, beta2);
    real c = (std::max)(real(0), (std::min)(initial, geometry.c_max));
    ClairautIterationReport best{c, ClairautCase::A_same_side};
    int performed = 0;
    for (int iteration = 0; iteration < maximum; ++iteration)
    {
        performed = iteration + 1;
        const auto state = longitude_state(geometry, c);
        unsigned selected = 0;
        if (!same_side_only)
        {
            const real r1 = abs(state.longitude[0] - target);
            const real r2 = abs(state.longitude[1] - target);
            const real r3 = abs(state.longitude[2] - target);
            if (r2 <= r1 && r2 <= r3) selected = 1;
            else if (r3 <= r1 && r3 < r2) selected = 2;
        }
        const real residual = abs(state.longitude[selected] - target);
        if (residual < best.resid)
            best = {c, static_cast<ClairautCase>(selected + 1), iteration + 1,
                residual, residual <= tolerance};
        if (best.converged) return best;
        real D;
        if (selected == 0) D = target + abs(state.dL1 - state.dL2);
        else if (selected == 1) D = target + abs(state.dL1 + state.dL2);
        else D = target + real(2) * state.epsilon - abs(state.dL1 + state.dL2);
        const real next = geometry.from_D(D);
        if (next == c) break; // Stagnation is not convergence.
        c = next;
    }
    best.iters = performed;
    return best; // Return the best trial and its convergence status.
}

} // namespace detail

// Compute the Clairaut constant from D and the endpoint reduced latitudes,
// using the norm form of Eq.(80).
inline real c_from_D(real D, real beta1, real beta2)
{
    return detail::ClairautGeometry(beta1, beta2).from_D(D);
}

// Fixed-point iteration with branch selection and convergence reporting.
// The tolerance applies to the absolute longitude residual, in radians.
inline ClairautIterationReport solve_clairaut_c_iter_D_checked(real c0,
    real beta1, real beta2, real dlon12, int max_iter = 80,
    real tol = real("1e-30"))
{
    return detail::iterate_clairaut(c0, beta1, beta2, dlon12, max_iter, tol, false);
}

// Return the Clairaut constant and selected branch.
// Throw std::runtime_error if the longitude-residual tolerance is not met.
inline ClairautResult solve_clairaut_c_iter_D(real c0, real beta1, real beta2,
    real dlon12, int max_iter = 80, real tol = real("1e-30"))
{
    const auto result = solve_clairaut_c_iter_D_checked(c0, beta1, beta2,
        dlon12, max_iter, tol);
    if (!result.converged)
        throw std::runtime_error("Clairaut iteration did not satisfy the longitude residual tolerance.");
    return {result.c, result.which};
}

// Fixed-point iteration restricted to Case A, with convergence reporting.
inline ClairautIterationReport solve_clairaut_c_iter_same_side_checked(real c0,
    real beta1, real beta2, real dlon12, int max_iter = 80,
    real tol = real("1e-30"))
{
    return detail::iterate_clairaut(c0, beta1, beta2, dlon12, max_iter, tol, true);
}

// Return the Clairaut constant for Case A.
// Throw std::runtime_error if the longitude-residual tolerance is not met.
inline real solve_clairaut_c_iter_same_side(real c0, real beta1, real beta2,
    real dlon12, int max_iter = 80, real tol = real("1e-30"))
{
    const auto result = solve_clairaut_c_iter_same_side_checked(c0, beta1,
        beta2, dlon12, max_iter, tol);
    if (!result.converged)
        throw std::runtime_error("Same-side Clairaut iteration did not satisfy the longitude residual tolerance.");
    return result.c;
}

} // namespace GeodesicCore
