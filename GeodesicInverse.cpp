// ============================================================================
// Copyright (c) 2026, Anhui Jianzhu University, China
// Author: Jia-Chun Guo
// Contact: jchguo@ahjzu.edu.cn
// License: MIT (https://opensource.org/licenses/MIT)
// ============================================================================

// GeodesicInverse.cpp
// Solve the shortest inverse geodesic between two given endpoints.
// Input and printed angles are in degrees; distances are in metres.
// The terminal azimuth is reported in both backward and continuing forward directions.
// No arguments: interactive decimal-degree/DMS input.
// Command line: GeodesicInverse lat1_deg lon1_deg lat2_deg lon2_deg.
// Build as a separate C++17 console executable.

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include "GeodesicDefinitions.hpp"
#include "GeodesicCore.hpp"
#include "NearAntipodalClairaut.hpp"

namespace GeodesicInverseSolver {

using Geodesic_Defs::real;
using Geodesic_Defs::PI;
using Geodesic_Defs::HALF_PI;
using Geodesic_Defs::a;
using Geodesic_Defs::f;
using Geodesic_Defs::abs;
using Geodesic_Defs::atan2;
using Geodesic_Defs::cos;
using Geodesic_Defs::sqrt;
using Geodesic_Defs::beta_of_lat;
using Geodesic_Defs::deg_to_rad;
using Geodesic_Defs::rad_to_deg;
using Geodesic_Defs::normalize_lon_deg;
using Geodesic_Defs::normalize_lon_rad;
using NearAntipodal::LCBranch;

// Absolute longitude-residual tolerance for iterative convergence, in radians.
inline const real longitude_tolerance("1e-30");
// Tolerance for identifying the symmetric vertex-limit configuration, in radians.
// This geometric classification tolerance is independent of iterative convergence.
inline const real vertex_longitude_tolerance("1e-30");
inline constexpr int fixed_point_iterations = 80;
inline constexpr int safeguarded_iterations = 512;

// Internal inverse result: distance in metres and azimuths in radians.
// azi1 is forward; azi2_backward is backward.
struct InverseSolution
{
    real s12;
    real azi1;
    real azi2_backward;
    real c = 0; // Nonnegative Clairaut constant of the selected geodesic.
};

// Evaluate distance and azimuths once after c and the geometric branch are known.
inline InverseSolution solution_from_branch(const real& c, LCBranch branch,
    const real& beta1, const real& beta2, bool increasing_latitude, int direction)
{
    const real phi1 = Geodesic_Defs::phi_of_beta(beta1, c);
    const real phi2 = Geodesic_Defs::phi_of_beta(beta2, c);
    const GeodesicCore::ArcLength arc(c);
    const real s1 = arc(phi1);
    const real s2 = arc(phi2);
    const real A1 = Geodesic_Defs::Azi(c, beta1);
    const real A2 = Geodesic_Defs::Azi(c, beta2);
    const real signed_A1 = real(direction) * A1;
    const real signed_A2 = real(direction) * A2;
    switch (branch)
    {
    case LCBranch::A:
        if (increasing_latitude)
            return {abs(s1 - s2), signed_A1, PI + signed_A2, c};
        return {abs(s1 - s2), PI - signed_A1, -signed_A2, c};
    case LCBranch::B:
        return {abs(s1 + s2), signed_A1, -signed_A2, c};
    case LCBranch::C:
        return {abs(real(2) * arc.half_period() - abs(s1 + s2)),
            PI - signed_A1, PI + signed_A2, c};
    }
    throw std::logic_error("Invalid inverse-geodesic branch.");
}

// Solve epsilon(c) = PI - longitude for a symmetric half geodesic.
// Initialize with the inverse n-series, then apply safeguarded Newton
// to the finite epsilon(c) series on [0, c_max].
inline real symmetric_half_constant(const real& longitude, const real& c_max)
{
    const real target = PI - longitude;
    real lower = 0, upper = c_max;
    real c = std::max(lower, std::min(GeodesicCore::c_from_Epsilon(longitude), upper));
    for (int iteration = 0; iteration < fixed_point_iterations; ++iteration)
    {
        const auto epsilon = GeodesicCore::Epsilon_with_derivative(c);
        const real residual = epsilon.value - target;
        if (abs(residual) <= longitude_tolerance) return c;
        if (residual < 0) lower = c;
        else upper = c;
        const real derivative = epsilon.derivative;
        real next = (lower + upper) / real(2);
        if (derivative != 0)
        {
            const real candidate = c - residual / derivative;
            if (candidate > lower && candidate < upper) next = candidate;
        }
        if (next == c) break;
        c = next;
    }
    throw std::runtime_error("Symmetric half-period inversion did not converge.");
}

// Solve Case A by fixed-point iteration,
// with safeguarded Newton iteration as a fallback.
inline real same_side_constant(const real& initial, const real& beta1,
    const real& beta2, const real& longitude)
{
    const auto ordinary = GeodesicCore::solve_clairaut_c_iter_same_side_checked(initial,
        beta1, beta2, longitude, fixed_point_iterations, longitude_tolerance);
    if (ordinary.converged) return ordinary.c;
    const auto safeguarded = NearAntipodal::solve_c_near_antipodal_one_branch(
        LCBranch::A, beta1, beta2, longitude, ordinary.c,
        safeguarded_iterations, longitude_tolerance);
    if (!safeguarded.converged)
        throw std::runtime_error("Same-side inverse solution did not converge.");
    return safeguarded.c;
}

// Solve a shortest inverse geodesic.
// When two symmetric shortest solutions exist, select the initially northward one.
inline InverseSolution inverse(const real& lat1_deg, const real& lon1_deg,
    const real& lat2_deg, const real& lon2_deg)
{
    for (const real& value : {lat1_deg, lon1_deg, lat2_deg, lon2_deg})
        Geodesic_Defs::require_finite(value, "Non-finite endpoint coordinate.");
    if (abs(lat1_deg) > 90 || abs(lat2_deg) > 90)
        throw std::invalid_argument("Latitude must lie in [-90, 90] degrees.");
    if (abs(lon1_deg) > 180 || abs(lon2_deg) > 180)
        throw std::invalid_argument("Longitude must lie in [-180, 180] degrees.");

    const real dlon_deg = normalize_lon_deg(lon2_deg - lon1_deg);
    const real signed_longitude = deg_to_rad(dlon_deg);
    const real longitude = abs(signed_longitude);
    const int direction = Geodesic_Defs::sign(dlon_deg);
    const real beta1 = beta_of_lat(deg_to_rad(lat1_deg));
    const real beta2 = beta_of_lat(deg_to_rad(lat2_deg));
    const bool first_pole = abs(lat1_deg) == 90;
    const bool second_pole = abs(lat2_deg) == 90;

    // Coincident endpoints have no unique azimuth. Longitude labels at a pole
    // identify the same physical point, and -180 and +180 are equivalent.
    if (lat1_deg == lat2_deg && (dlon_deg == 0 || first_pole))
        throw std::invalid_argument("Coincident endpoints have no unique inverse azimuth.");

    // Case I1: poles and meridians. Polar azimuths use the meridional
    // convention (0 or pi), independent of the supplied polar longitude.
    if (first_pole || second_pole || dlon_deg == 0 || abs(dlon_deg) == 180)
    {
        const real phi1 = Geodesic_Defs::phi_of_beta(beta1, real(0));
        const real phi2 = Geodesic_Defs::phi_of_beta(beta2, real(0));
        const GeodesicCore::ArcLength arc(real(0));
        const real s1 = arc(phi1), s2 = arc(phi2);
        if (first_pole || second_pole)
        {
            // Assign forward and reverse azimuths along the selected meridian.
            if (lat1_deg < lat2_deg)
                return { abs(s1 - s2), real(0), PI };
            return { abs(s1 - s2), PI, real(0) };
        }
        if (dlon_deg == 0)
        {
            if (lat1_deg < lat2_deg) return {abs(s1 - s2), real(0), PI};
            return {abs(s1 - s2), PI, real(0)};
        }
        if (lat2_deg >= -lat1_deg)
            return {abs(s1 + s2), real(0), real(0)};
        return {abs(real(2) * arc.half_period() - s1 - s2), PI, PI};
    }

    // Case I2: two equatorial endpoints.
    // Beyond the equatorial shortest-arc limit,
    // select the initially northward member of the symmetric pair.
    const real equatorial_limit = PI * (real(1) - f);
    if (lat1_deg == 0 && lat2_deg == 0)
    {
        if (longitude <= equatorial_limit)
            return {a * longitude, real(direction) * HALF_PI,
                -real(direction) * HALF_PI, real(1)};
        const real c = symmetric_half_constant(longitude, real(1));
        const real A = atan2(c, sqrt((real(1) - c) * (real(1) + c)));
        return {GeodesicCore::ArcLength(c).half_period(),
            real(direction) * A, -real(direction) * A, c};
    }

    const real initial = GeodesicCore::c_from_D(longitude, beta1, beta2);

    // Case I3: symmetric latitudes, lat2 = -lat1.
    // Compare the longitude difference with PI - epsilon(c_max) to distinguish
    // same-side, vertex-limit, and symmetric half-geodesic configurations.
    if (lat2_deg == -lat1_deg)
    {
        if (longitude <= equatorial_limit)
        {
            const real c = same_side_constant(initial, beta1, beta2, longitude);
            return solution_from_branch(c, LCBranch::A, beta1, beta2,
                lat1_deg < lat2_deg, direction);
        }
        const real c_max = abs(cos(beta1));
        const real vertex_longitude = PI - GeodesicCore::Epsilon(c_max);
        if (abs(longitude - vertex_longitude) <= vertex_longitude_tolerance)
            return {GeodesicCore::ArcLength(c_max).half_period(),
                real(direction) * HALF_PI, -real(direction) * HALF_PI, c_max};
        if (longitude < vertex_longitude)
        {
            const real c = same_side_constant(initial, beta1, beta2, longitude);
            return solution_from_branch(c, LCBranch::A, beta1, beta2,
                lat1_deg < lat2_deg, direction);
        }
        const real c = symmetric_half_constant(longitude, c_max);
        const real A = Geodesic_Defs::Azi(c, beta1);
        return {GeodesicCore::ArcLength(c).half_period(),
            real(direction) * A, -real(direction) * A, c};
    }

    // Case I5: operational near-antipodal domain, Eq.(89).
    // The safeguarded solver selects the branch from the endpoint geometry.
    const bool near_antipodal = abs(lat1_deg + lat2_deg) <= 1
        && abs(dlon_deg) >= 179 && abs(dlon_deg) < 180;
    if (near_antipodal)
    {
        const auto result = NearAntipodal::solve_c_near_antipodal_all_branches(
            beta1, beta2, longitude, initial, safeguarded_iterations, longitude_tolerance);
        return solution_from_branch(result.c, result.br, beta1, beta2,
            lat1_deg < lat2_deg, direction);
    }

    // Case I4: fixed-point iteration with a longitude-residual check.
    // If needed, initialize the safeguarded solver with the best fixed-point estimate.
    const auto ordinary = GeodesicCore::solve_clairaut_c_iter_D_checked(initial, beta1,
        beta2, longitude, fixed_point_iterations, longitude_tolerance);
    if (ordinary.converged)
        return solution_from_branch(ordinary.c,
            static_cast<LCBranch>(static_cast<int>(ordinary.which)), beta1, beta2,
            lat1_deg < lat2_deg, direction);
    const auto safeguarded = NearAntipodal::solve_c_near_antipodal_all_branches(
        beta1, beta2, longitude, ordinary.c, safeguarded_iterations, longitude_tolerance);
    return solution_from_branch(safeguarded.c, safeguarded.br, beta1, beta2,
        lat1_deg < lat2_deg, direction);
}

} // namespace GeodesicInverseSolver

#ifndef GEODESIC_CONSOLE_INPUT_DEFINED
#define GEODESIC_CONSOLE_INPUT_DEFINED
namespace ConsoleInput {
using Geodesic_Defs::real;

inline std::string trim(const std::string& text)
{
    const auto first = text.find_first_not_of(" \t\r");
    if (first == std::string::npos) return {};
    const auto last = text.find_last_not_of(" \t\r");
    return text.substr(first, last - first + 1);
}

inline real decimal(const std::string& token)
{
    std::istringstream input(token);
    real value;
    if (!(input >> value) || !boost::math::isfinite(value))
        throw std::invalid_argument("Expected a finite decimal number.");
    input >> std::ws;
    if (!input.eof()) throw std::invalid_argument("Unexpected text after the number.");
    return value;
}

inline bool read_angle(const char* prompt, bool dms, real& value)
{
    std::cout << prompt;
    if (dms)
    {
        if (!Geodesic_Defs::read_dms_space(value))
        {
            if (std::cin.eof()) return false;
            throw std::invalid_argument("Invalid DMS input; minutes and seconds must lie in [0,60).");
        }
        return true;
    }
    std::string line;
    if (!std::getline(std::cin, line)) return false;
    value = decimal(line);
    return true;
}

inline bool read_decimal(const char* prompt, real& value)
{
    return read_angle(prompt, false, value);
}
} // namespace ConsoleInput
#endif

namespace InverseConsole {
using Geodesic_Defs::real;
inline void display(const GeodesicInverseSolver::InverseSolution& solution)
{
    const real r = Geodesic_Defs::sqrt((real(1) - solution.c) * (real(1) + solution.c));
    const real lat_vertex = Geodesic_Defs::rad_to_deg(
        Geodesic_Defs::atan2(r, Geodesic_Defs::one_minus_f * solution.c));
    const real beta_vertex = Geodesic_Defs::rad_to_deg(Geodesic_Defs::atan2(r, solution.c));
    const real forward_terminal = solution.azi2_backward - Geodesic_Defs::PI;
    std::cout << "\n======== INVERSE RESULTS ========\n"
        << std::fixed << std::setprecision(30)
        << "s12 = " << solution.s12 << " m\n"
        << "azi1 = " << Geodesic_Defs::rad_to_deg(Geodesic_Defs::normalize_azi_rad(solution.azi1))
        << " deg  (forward direction, [0,360))\n"
        << "azi2_backward = " << Geodesic_Defs::rad_to_deg(
            Geodesic_Defs::normalize_azi_rad(solution.azi2_backward))
        << " deg  (backward direction, [0,360))\n"
        << "azi2_forward = " << Geodesic_Defs::rad_to_deg(
            Geodesic_Defs::normalize_azi_rad(forward_terminal))
        << " deg  (Karney forward direction, [0,360))\n"
        << "azi1_signed = " << Geodesic_Defs::rad_to_deg(
            Geodesic_Defs::normalize_lon_rad(solution.azi1)) << " deg\n"
        << "azi2_backward_signed = " << Geodesic_Defs::rad_to_deg(
            Geodesic_Defs::normalize_lon_rad(solution.azi2_backward)) << " deg\n"
        << "azi2_forward_signed = " << Geodesic_Defs::rad_to_deg(
            Geodesic_Defs::normalize_lon_rad(forward_terminal)) << " deg\n"
        << "lat_V = " << lat_vertex << " deg\n"
        << "beta_V = " << beta_vertex << " deg\n"
        << "c = " << solution.c << '\n';
}
} // namespace InverseConsole

#ifndef GEODESIC_INVERSE_NO_MAIN
int main(int argc, char* argv[])
{
    using Geodesic_Defs::real;
    if (argc != 1)
    {
        if (argc != 5)
        {
            std::cerr << "Usage: GeodesicInverse [lat1_deg lon1_deg lat2_deg lon2_deg]\n";
            return 1;
        }
        try
        {
            InverseConsole::display(GeodesicInverseSolver::inverse(
                ConsoleInput::decimal(argv[1]), ConsoleInput::decimal(argv[2]),
                ConsoleInput::decimal(argv[3]), ConsoleInput::decimal(argv[4])));
            return 0;
        }
        catch (const std::exception& exception)
        {
            std::cerr << "ERROR: " << exception.what() << '\n';
            return 1;
        }
    }
    for (std::string command;;)
    {
        std::cout << "\n===== Inverse Geodesic Calculator =====\n"
            << "1 = decimal degrees\n2 = degrees-minutes-seconds (DMS)\nq = quit\n> ";
        if (!std::getline(std::cin, command)) break;
        command = ConsoleInput::trim(command);
        if (command.empty()) continue;
        if (command == "q" || command == "Q") break;
        if (command != "1" && command != "2")
        {
            std::cerr << "ERROR: choose 1, 2, or q.\n";
            continue;
        }
        try
        {
            real lat1, lon1, lat2, lon2;
            const bool dms = command == "2";
            if (!ConsoleInput::read_angle("Enter lat1: ", dms, lat1)
                || !ConsoleInput::read_angle("Enter lon1: ", dms, lon1)
                || !ConsoleInput::read_angle("Enter lat2: ", dms, lat2)
                || !ConsoleInput::read_angle("Enter lon2: ", dms, lon2)) break;
            InverseConsole::display(GeodesicInverseSolver::inverse(lat1, lon1, lat2, lon2));
        }
        catch (const std::exception& exception)
        {
            std::cerr << "ERROR: " << exception.what() << '\n';
        }
    }
    std::cout << "Program exited.\n";
    return 0;
}
#endif
