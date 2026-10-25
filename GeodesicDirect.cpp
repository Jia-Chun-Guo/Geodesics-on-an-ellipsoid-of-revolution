// ============================================================================
// Copyright (c) 2026, Anhui Jianzhu University, China
// Author: Jia-Chun Guo
// Contact: jchguo@ahjzu.edu.cn
// License: MIT (https://opensource.org/licenses/MIT)
// ============================================================================

// GeodesicDirect.cpp
// Single-point direct geodesics, including multiple complete latitude cycles.
// No arguments: interactive decimal-degree/DMS input. Four arguments:
// GeodesicDirect lat1_deg lon1_deg azi1_deg s12_m
// Angles in the public input/output are degrees; core computations use radians.
// Build this source as its own C++17 console executable.

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include "GeodesicDefinitions.hpp"
#include "GeodesicCore.hpp"

namespace GeodesicDirectSolver {

using Geodesic_Defs::real;
using Geodesic_Defs::PI;
using Geodesic_Defs::HALF_PI;
using Geodesic_Defs::TWO_PI;
using Geodesic_Defs::a;
using Geodesic_Defs::abs;
using Geodesic_Defs::atan2;
using Geodesic_Defs::cos;
using Geodesic_Defs::sin;
using Geodesic_Defs::hypot;
using Geodesic_Defs::floor;
using Geodesic_Defs::deg_to_rad;
using Geodesic_Defs::rad_to_deg;
using Geodesic_Defs::normalize_lon_deg;
using Geodesic_Defs::normalize_lon_rad;
using Geodesic_Defs::normalize_azi_deg;

struct DirectSolution
{
    real lat2_deg;
    real lon2_deg;
    real azi2_deg; // Continuing forward direction, normalized to (-180, 180].
    real c;
    // Half of the distance-reduction cycle, in metres:
    // s_pi for non-equatorial geodesics; PI*a for the Equator.
    real half_period_m;
    real complete_cycles; // Integer-valued count of complete distance-reduction cycles.
    real reduced_distance_m;
    const char* geometric_case;
};

struct DistanceReduction
{
    real cycles;
    real remainder;
};

inline DistanceReduction reduce_distance(const real& distance, const real& period)
{
    real cycles = floor(distance / period);
    real remainder = distance - cycles * period;
    // Correct a quotient rounded across an integer boundary.
    if (remainder < 0) { remainder += period; cycles -= 1; }
    if (remainder >= period) { remainder -= period; cycles += 1; }
    // Products involving a large integer cycle count can leave a remainder
    // at the working-precision level even for a represented complete cycle.
    const real allowance = real(64) * std::numeric_limits<real>::epsilon()
        * (real(1) + abs(distance) + period);
    if (cycles != 0 && remainder <= allowance) remainder = 0;
    else if (period - remainder <= allowance) { remainder = 0; cycles += 1; }
    if (remainder < 0 || remainder >= period)
        throw std::runtime_error("Distance reduction exceeded the available precision.");
    return {cycles, remainder};
}

inline real amplitude_from_distance(const real& distance, const real& half_period,
    const real& mean_scale, const real& c)
{
    if (distance == 0) return real(0);
    if (distance == half_period) return PI;
    if (distance == half_period / real(2)) return HALF_PI;
    return GeodesicCore::phi_from_phi_prime(distance / mean_scale, c);
}

// Solve from the initial point, forward azimuth, and nonnegative distance.
// Return the endpoint and its continuing forward azimuth.
inline DirectSolution direct(real lat1_deg, real lon1_deg, real azi1_deg, real s12)
{
    for (const real& value : {lat1_deg, lon1_deg, azi1_deg, s12})
        Geodesic_Defs::require_finite(value, "Non-finite direct-problem input.");
    if (abs(lat1_deg) > 90) throw std::invalid_argument("Latitude must lie in [-90, 90] degrees.");
    if (abs(lon1_deg) > 180) throw std::invalid_argument("Longitude must lie in [-180, 180] degrees.");
    if (azi1_deg < -180 || azi1_deg > 360)
        throw std::invalid_argument("Azimuth must lie in [-180, 360] degrees.");
    if (s12 < 0) throw std::invalid_argument("Distance must be nonnegative.");
    azi1_deg = normalize_azi_deg(azi1_deg);
    const real azi1 = deg_to_rad(azi1_deg);
    const real lon1 = deg_to_rad(lon1_deg);
    real sine_alpha, cosine_alpha;
    if (azi1_deg == 0) { sine_alpha = 0; cosine_alpha = 1; }
    else if (azi1_deg == 90) { sine_alpha = 1; cosine_alpha = 0; }
    else if (azi1_deg == 180) { sine_alpha = 0; cosine_alpha = -1; }
    else if (azi1_deg == 270) { sine_alpha = -1; cosine_alpha = 0; }
    else { sine_alpha = sin(azi1); cosine_alpha = cos(azi1); }

    const bool at_pole = abs(lat1_deg) == 90;
    const real beta1 = lat1_deg == 90 ? HALF_PI : lat1_deg == -90 ? real(-HALF_PI)
        : Geodesic_Defs::beta_of_lat(deg_to_rad(lat1_deg));
    const real sine_beta = at_pole ? (lat1_deg > 0 ? real(1) : real(-1)) : real(sin(beta1));
    const real cosine_beta = at_pole ? real(0) : real(cos(beta1));
    const real c = (std::min)(real(1), real(abs(sine_alpha * cosine_beta)));

    // Case D1: an equatorial geodesic has the equatorial circumference as its
    // physical cycle. It has no nontrivial latitude oscillation.
    if (lat1_deg == 0 && (azi1_deg == 90 || azi1_deg == 270))
    {
        const real half = a * PI;
        const auto reduced = reduce_distance(s12, real(2) * half);
        const int east_west = azi1_deg == 90 ? +1 : -1;
        return {real(0), normalize_lon_deg(lon1_deg
                + real(east_west) * rad_to_deg(reduced.remainder / a)),
            normalize_lon_deg(azi1_deg), c, half, reduced.cycles,
            reduced.remainder, "D1: Equator"};
    }

    const GeodesicCore::ArcLength arc(c);
    const real half = arc.half_period();
    const real period = real(2) * half;
    const real mean_scale = half / PI;  // Mean arc-length scale: s_pi/PI = a'*F_E(c).
    const auto reduced = reduce_distance(s12, period);
    const real S12 = reduced.remainder;
    if (s12 == 0)
        return {lat1_deg, normalize_lon_deg(lon1_deg), normalize_lon_deg(azi1_deg),
            c, half, real(0), real(0), "D0: Zero distance"};

    // Case D3: meridians, including arbitrary azimuths specified at a pole.
    // At c=0 the longitude jumps by a half-turn at a pole; evaluating atan2(0,0)
    // in the generic longitude formula would lose that geometric information.
    if (c == 0)
    {
        if (S12 == 0)
            return {lat1_deg, normalize_lon_deg(lon1_deg), normalize_lon_deg(azi1_deg),
                c, half, reduced.cycles, S12, "D3: Complete meridian cycles"};
        real plane_longitude = lon1_deg;
        int southward;
        if (lat1_deg == 90) { plane_longitude += real(180) - azi1_deg; southward = +1; }
        else if (lat1_deg == -90) { plane_longitude += azi1_deg; southward = -1; }
        else southward = azi1_deg == 0 ? -1 : +1;
        const real phi1 = HALF_PI - beta1;
        const real s1 = lat1_deg == 90 ? real(0) : lat1_deg == -90 ? half : arc(phi1);
        auto position = reduce_distance(s1 + real(southward) * S12, period);
        const real allowance = real(64) * std::numeric_limits<real>::epsilon() * period;
        if (position.remainder <= allowance || period - position.remainder <= allowance)
            position.remainder = 0;
        else if (abs(position.remainder - half) <= allowance) position.remainder = half;
        const real folded = position.remainder <= half ? position.remainder
            : real(period - position.remainder);
        // At an exact pole, use the incoming meridian's longitude frame.
        const bool opposite_half = position.remainder > half
            || (position.remainder == half && southward < 0)
            || (position.remainder == 0 && southward > 0);
        const real phi2 = amplitude_from_distance(folded, half, mean_scale, c);
        const real latitude = folded == 0 ? real(90) : folded == half ? real(-90)
            : rad_to_deg(Geodesic_Defs::lat_of_beta(HALF_PI - phi2));
        const bool northward = (southward < 0) != opposite_half;
        return {latitude, normalize_lon_deg(plane_longitude + (opposite_half ? real(180) : real(0))),
            northward ? real(0) : real(180), c, half, reduced.cycles, S12, "D3: Meridian"};
    }

    // Case D2: non-equatorial, non-meridional geodesics.
    // Select one of six endpoint branches from the initial direction
    // and the reduced arc-length increment S12.
    const GeodesicCore::LongitudeCoefficients coefficients(c);
    const real contraction = coefficients.epsilon;
    const int east_west = azi1_deg < 180 ? +1 : -1;
    // Signed longitude contraction accumulated over complete latitude cycles.
    const real cycle_correction = real(2) * reduced.cycles * contraction * real(east_west); 
    if (S12 == 0)
        return {lat1_deg, rad_to_deg(normalize_lon_rad(lon1 - cycle_correction)),
            normalize_lon_deg(azi1_deg), c, half, reduced.cycles, S12, "D2: Complete latitude cycles"};

    const real numerator = abs(cosine_beta * cosine_alpha);
    const real phi1 = atan2(numerator, sine_beta);
    const real s1 = arc(phi1);
    const real longitude1 = atan2(numerator, c * sine_beta)
        - GeodesicCore::dL_with_coefficients(phi1, coefficients);
    const int direction = Geodesic_Defs::Sign_sin2_deg(azi1_deg);
    const real half_longitude = PI - contraction;
    const bool initially_northward = azi1_deg <= 90 || azi1_deg >= 270;
    real endpoint_distance, longitude_offset = 0, azimuth_offset = 0;
    int endpoint_sign = -1, azimuth_sign = direction;
    const char* geometric_case;
    if (initially_northward)
    {
        if (S12 <= s1)
        {
            endpoint_distance = s1 - S12;
            geometric_case = "D2.1.1";
        }
        else if (S12 <= s1 + half)
        {
            endpoint_distance = S12 - s1;
            endpoint_sign = +1;
            azimuth_offset = PI;
            azimuth_sign = -direction;
            geometric_case = "D2.1.2";
        }
        else
        {
            endpoint_distance = period + s1 - S12;
            longitude_offset = real(2) * half_longitude;
            geometric_case = "D2.1.3";
        }
    }
    else
    {
        if (S12 <= half - s1)
        {
            endpoint_distance = s1 + S12;
            azimuth_offset = PI;
            geometric_case = "D2.2.1";
        }
        else if (S12 <= period - s1)
        {
            endpoint_distance = period - s1 - S12;
            endpoint_sign = +1;
            longitude_offset = -real(2) * half_longitude;
            azimuth_sign = -direction;
            geometric_case = "D2.2.2";
        }
        else
        {
            endpoint_distance = s1 + S12 - period;
            longitude_offset = -real(2) * half_longitude;
            azimuth_offset = PI;
            geometric_case = "D2.2.3";
        }
    }
    // Enforce the theoretical interval [0, s_pi] for endpoint_distance
    // before arc-length inversion.
    endpoint_distance = (std::max)(real(0), (std::min)(endpoint_distance, half));
    const real phi2 = amplitude_from_distance(endpoint_distance, half, mean_scale, c);
    const real sine_phi = sin(phi2), cosine_phi = cos(phi2);
    const real p = hypot(cosine_alpha, sine_alpha * sine_beta);
    const real beta2 = atan2(p * cosine_phi, hypot(sine_phi, c * cosine_phi));
    const real longitude2 = atan2(sine_phi, c * cosine_phi)
        - GeodesicCore::dL_with_coefficients(phi2, coefficients);
    const real longitude_increment = longitude1 + real(endpoint_sign) * longitude2 + longitude_offset;
    const real principal_azimuth = atan2(c, p * sine_phi);
    const real azimuth = azimuth_offset + real(azimuth_sign) * principal_azimuth;
    return {rad_to_deg(Geodesic_Defs::lat_of_beta(beta2)),
        rad_to_deg(normalize_lon_rad(lon1 + real(direction) * longitude_increment - cycle_correction)),
        normalize_lon_deg(rad_to_deg(azimuth)), c, half, reduced.cycles, S12, geometric_case};
}

} // namespace GeodesicDirectSolver

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

namespace DirectConsole {
using Geodesic_Defs::real;

inline void display(const real& lat1, const real& lon1, const real& azi1,
    const real& distance, const GeodesicDirectSolver::DirectSolution& solution)
{
    const real backward = Geodesic_Defs::normalize_azi_deg(solution.azi2_deg + real(180));
    std::cout << "\n======== DIRECT RESULTS ========\n"
        << std::fixed << std::setprecision(30)
        << "lat1 = " << lat1 << " deg\nlon1 = " << lon1 << " deg\n"
        << "azi1 = " << Geodesic_Defs::normalize_azi_deg(azi1) << " deg\n"
        << "s12 = " << distance << " m\n"
        << "lat2 = " << solution.lat2_deg << " deg\n"
        << "lon2 = " << solution.lon2_deg << " deg\n"
        << "azi2_backward = " << backward << " deg  (backward direction, [0,360))\n"
        << "azi2_forward = " << Geodesic_Defs::normalize_azi_deg(solution.azi2_deg)
        << " deg  (Karney forward direction, [0,360))\n"
        << "s_PI = " << solution.half_period_m << " m  (length of a half geodesic)\n"
        << "c = " << solution.c << '\n'
        << std::setprecision(0) << "complete_cycles = " << solution.complete_cycles << '\n'
        << std::setprecision(30) << "remaining_distance = " << solution.reduced_distance_m << " m\n"
        << "geometric_case = " << solution.geometric_case << "\n";
}
} // namespace DirectConsole

#ifndef GEODESIC_DIRECT_NO_MAIN
int main(int argc, char* argv[])
{
    using Geodesic_Defs::real;
    if (argc != 1)
    {
        if (argc != 5)
        {
            std::cerr << "Usage: GeodesicDirect [lat1_deg lon1_deg azi1_deg s12_m]\n";
            return 1;
        }
        try
        {
            const real lat1 = ConsoleInput::decimal(argv[1]), lon1 = ConsoleInput::decimal(argv[2]);
            const real azi1 = ConsoleInput::decimal(argv[3]), distance = ConsoleInput::decimal(argv[4]);
            DirectConsole::display(lat1, lon1, azi1, distance,
                GeodesicDirectSolver::direct(lat1, lon1, azi1, distance));
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
        std::cout << "\n===== Direct Geodesic Calculator =====\n"
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
            real lat1, lon1, azi1, distance;
            const bool dms = command == "2";
            if (!ConsoleInput::read_angle("Enter lat1: ", dms, lat1)
                || !ConsoleInput::read_angle("Enter lon1: ", dms, lon1)
                || !ConsoleInput::read_angle("Enter azi1: ", dms, azi1)
                || !ConsoleInput::read_decimal("Enter s12 (metres, nonnegative): ", distance)) break;
            DirectConsole::display(lat1, lon1, azi1, distance,
                GeodesicDirectSolver::direct(lat1, lon1, azi1, distance));
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
