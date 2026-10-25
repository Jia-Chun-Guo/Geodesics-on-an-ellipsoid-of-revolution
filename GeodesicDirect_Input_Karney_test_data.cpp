// ============================================================================
// Copyright (c) 2026, Anhui Jianzhu University, China
// Author: Jia-Chun Guo
// Contact: jchguo@ahjzu.edu.cn
// License: MIT (https://opensource.org/licenses/MIT)
// ============================================================================

// GeodesicDirect_Input_Karney_test_data.cpp
// Batch direct solutions for GeodTest's ten-column format:
// lat1 lon1 azi1 lat2 lon2 azi2 s12 a12 m12 S12.
// Only lat1, lon1, azi1, and s12 are used as direct-problem input.
// Endpoint azimuth comparisons use Karney's continuing forward convention.
// No arguments: Windows file dialog or terminal prompt. Optional arguments:
// GeodesicDirectTest input.dat [output.csv]
// Build this source as its own C++17 console executable.
// Solve from lat1, lon1, azi1, and s12.
// Compare the computed lat2, lon2, and continuing forward azi2 with the reference values.

#include <algorithm>
#include <array>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <commdlg.h>
#ifdef _MSC_VER
#pragma comment(lib, "Comdlg32.lib")
#endif
#endif
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
    real half_period_m;   // s_pi for non-equatorial geodesics; PI*a for the Equator.
    real complete_cycles; // Integer-valued count of complete distance-reduction cycles..
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
    const real mean_scale = half / PI; // Mean arc-length scale: s_pi/PI = a'*F_E(c).
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

namespace GeodesicDirectTest {
using Geodesic_Defs::real;

// Calculated minus reference, wrapped to [-180, 180) degrees.
inline real angular_difference_deg(const real& calculated, const real& reference)
{
    real difference = calculated - reference;
    difference -= real(360) * Geodesic_Defs::floor((difference + real(180)) / real(360));
    return difference;
}

std::filesystem::path select_karney_test_file()
{
#ifdef _WIN32
    MessageBoxExW(nullptr, L"Please select the Karney GeodTest data file.",
        L"GeodesicDirect", MB_OK | MB_ICONINFORMATION,
        MAKELANGID(LANG_ENGLISH, SUBLANG_ENGLISH_US));
    wchar_t file_name[32768] = L"";
    OPENFILENAMEW dialog{};
    dialog.lStructSize = sizeof(dialog);
    dialog.lpstrFile = file_name;
    dialog.nMaxFile = static_cast<DWORD>(sizeof(file_name) / sizeof(file_name[0]));
    dialog.lpstrFilter = L"Karney geodesic test files (*.dat)\0*.dat\0All files (*.*)\0*.*\0\0";
    dialog.nFilterIndex = 1;
    dialog.lpstrTitle = L"Select a Karney GeodTest data file";
    dialog.lpstrDefExt = L"dat";
    dialog.Flags = OFN_EXPLORER | OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST | OFN_NOCHANGEDIR;
    if (GetOpenFileNameW(&dialog)) return std::filesystem::path(file_name);
    const DWORD error = CommDlgExtendedError();
    if (error != 0)
        throw std::runtime_error("File-selection dialog failed; error code " + std::to_string(error));
    return {};
#else
    std::cout << "Enter the path of the Karney GeodTest data file: ";
    std::string name;
    std::getline(std::cin >> std::ws, name);
    return std::filesystem::path(name);
#endif
}

// Latitude difference is calculated minus reference.
// Longitude and azimuth differences are wrapped to [-180, 180) degrees.
// Write values in fixed notation and differences in scientific notation.
int run_file(const std::filesystem::path& input_path, std::filesystem::path output_path)
{
    if (input_path.empty())
    {
        std::cout << "No input file selected. Program exited.\n";
        return 0;
    }
    if (output_path.empty())
    {
        output_path = input_path.parent_path() / input_path.stem();
        output_path += "_GeodesicDirect_results.CSV";
    }
    std::error_code error;
    if (std::filesystem::equivalent(input_path, output_path, error)
        || std::filesystem::absolute(input_path).lexically_normal()
            == std::filesystem::absolute(output_path).lexically_normal())
        throw std::invalid_argument("Input and output paths must differ.");
    std::ifstream input(input_path);
    if (!input) throw std::runtime_error("Cannot open the selected input file.");
    std::ofstream output(output_path);
    if (!output) throw std::runtime_error("Cannot create the CSV output file.");
    output << "lat1_deg,lon1_deg,azi1_deg,s12_m,"
        "lat2_calc_deg,lon2_calc_deg,azi2_calc_deg,"
        "lat2_error_deg,lon2_error_deg,azi2_error_deg\n";

    std::string line;
    std::size_t line_number = 0, records = 0;
    while (std::getline(input, line))
    {
        ++line_number;
        if (line_number == 1 && line.compare(0, 3, "\xEF\xBB\xBF") == 0) line.erase(0, 3);
        const std::size_t first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos || line[first] == '#') continue;
        std::array<real, 10> fields{};
        try
        {
            std::istringstream row(line);
            for (real& value : fields)
            {
                if (!(row >> value) || !boost::math::isfinite(value))
                    throw std::invalid_argument("Expected ten finite numeric fields.");
            }
            row >> std::ws;
            if (!row.eof() && row.peek() != '#')
                throw std::invalid_argument("Unexpected text after the tenth field.");
            if (fields[6] < 0) throw std::invalid_argument("Reference distance must be nonnegative.");
            const auto result = GeodesicDirectSolver::direct(fields[0], fields[1], fields[2], fields[6]);
            if (!boost::math::isfinite(result.lat2_deg) || !boost::math::isfinite(result.lon2_deg)
                || !boost::math::isfinite(result.azi2_deg))
                throw std::runtime_error("The direct solution is not finite.");
            const real latitude_error = result.lat2_deg - fields[3];
            const real longitude_error = angular_difference_deg(result.lon2_deg, fields[4]);
            const real azimuth_error = angular_difference_deg(result.azi2_deg, fields[5]);
            // Preserve the original CSV column order, precision, and error signs.
            output << std::fixed << std::setprecision(20)
                << fields[0] << ',' << fields[1] << ',' << fields[2] << ',' << fields[6] << ','
                << result.lat2_deg << ',' << result.lon2_deg << ',' << result.azi2_deg << ','
                << std::scientific << std::setprecision(20)
                << latitude_error << ',' << longitude_error << ',' << azimuth_error << '\n';
            ++records;
        }
        catch (const std::exception& exception)
        {
            std::ostringstream message;
            message << std::setprecision(40) << "GeodTest direct failure at line " << line_number
                << ": " << exception.what() << '\n'
                << "lat1=" << fields[0] << ", lon1=" << fields[1]
                << ", azi1=" << fields[2] << ", s12=" << fields[6];
            throw std::runtime_error(message.str());
        }
    }
    if (input.bad()) throw std::runtime_error("Failed while reading the input file.");
    output.flush();
    if (!output) throw std::runtime_error("Failed while writing the CSV output file.");
    std::cout << "Processed " << records << " records.\n";
#ifdef _WIN32
    std::wcout << L"CSV output: " << output_path.wstring() << L'\n';
#else
    std::cout << "CSV output: " << output_path.string() << '\n';
#endif
    return 0;
}

} // namespace GeodesicDirectTest

#ifndef GEODESIC_DIRECT_TEST_NO_MAIN
// MSVC's wide-character entry point preserves Unicode command-line paths.
#if defined(_WIN32) && defined(_MSC_VER)
int wmain(int argc, wchar_t* argv[])
#else
int main(int argc, char* argv[])
#endif
{
    try
    {
        if (argc > 3)
        {
            std::cerr << "Usage: GeodesicDirectTest [input.dat [output.csv]]\n";
            return 1;
        }
        const std::filesystem::path input = argc >= 2
            ? std::filesystem::path(argv[1]) : GeodesicDirectTest::select_karney_test_file();
        const std::filesystem::path output = argc == 3
            ? std::filesystem::path(argv[2]) : std::filesystem::path();
        return GeodesicDirectTest::run_file(input, output);
    }
    catch (const std::exception& exception)
    {
        std::cerr << "ERROR: " << exception.what() << '\n';
        return 1;
    }
}
#endif
