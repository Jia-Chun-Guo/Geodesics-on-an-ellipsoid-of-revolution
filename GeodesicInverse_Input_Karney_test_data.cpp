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
// Solve from lat1, lon1, lat2, and lon2.
// Compare the computed s12 and forward azimuths with the reference values.

#include <algorithm>
#include <array>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
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
#include "NearAntipodalClairaut.hpp"

namespace GeodesicInverseTest {

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
inline const real vertex_longitude_tolerance("2e-20");
inline constexpr int fixed_point_iterations = 80;
inline constexpr int safeguarded_iterations = 256;

struct InverseSolution
{
    real s12;
    real azi1;
    real azi2_backward;
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
            return {abs(s1 - s2), signed_A1, PI + signed_A2};
        return {abs(s1 - s2), PI - signed_A1, -signed_A2};
    case LCBranch::B:
        return {abs(s1 + s2), signed_A1, -signed_A2};
    case LCBranch::C:
        return {abs(real(2) * arc.half_period() - abs(s1 + s2)),
            PI - signed_A1, PI + signed_A2};
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
                -real(direction) * HALF_PI};
        const real c = symmetric_half_constant(longitude, real(1));
        const real A = atan2(c, sqrt((real(1) - c) * (real(1) + c)));
        return {GeodesicCore::ArcLength(c).half_period(),
            real(direction) * A, -real(direction) * A};
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
                real(direction) * HALF_PI, -real(direction) * HALF_PI};
        if (longitude < vertex_longitude)
        {
            const real c = same_side_constant(initial, beta1, beta2, longitude);
            return solution_from_branch(c, LCBranch::A, beta1, beta2,
                lat1_deg < lat2_deg, direction);
        }
        const real c = symmetric_half_constant(longitude, c_max);
        const real A = Geodesic_Defs::Azi(c, beta1);
        return {GeodesicCore::ArcLength(c).half_period(),
            real(direction) * A, -real(direction) * A};
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

// Errors use the signed periodic difference in [-180, 180) degrees.
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
        L"GeodesicInverse", MB_OK | MB_ICONINFORMATION,
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

// Distance difference is calculated minus reference.
// Azimuth differences are wrapped to [-180, 180) degrees.
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
        output_path += "_GeodesicInverse_results.CSV";
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
    output << "lat1_deg,lon1_deg,lat2_deg,lon2_deg,"
        "s12_calc_m,azi1_calc_deg,azi2_calc_deg,"
        "s12_error_m,azi1_error_deg,azi2_error_deg\n";

    std::string line;
    std::size_t line_number = 0, records = 0;
    while (std::getline(input, line))
    {
        ++line_number;
        if (line_number == 1 && line.compare(0, 3, "\xEF\xBB\xBF") == 0) line.erase(0, 3);
        const std::size_t first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos || line[first] == '#') continue;
        std::array<real, 10> fields{};
        // std::array<real, 7> fields{}; 10 
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
            const auto result = inverse(fields[0], fields[1], fields[3], fields[4]);
            const real azi1 = rad_to_deg(normalize_lon_rad(result.azi1));
            // Convert the terminal azimuth from backward to continuing forward direction
            // before comparison with the GeodTest reference.
            const real azi2 = rad_to_deg(normalize_lon_rad(result.azi2_backward - PI));
            if (!boost::math::isfinite(result.s12) || !boost::math::isfinite(azi1)
                || !boost::math::isfinite(azi2) || result.s12 < 0)
                throw std::runtime_error("The inverse solution is not finite or has negative distance.");
            const real distance_error = result.s12 - fields[6];
            const real azi1_error = angular_difference_deg(azi1, fields[2]);
            const real azi2_error = angular_difference_deg(azi2, fields[5]);
            // Preserve the original DAT column order, precision, and error signs.
            // output << std::fixed << std::setprecision(30)
            output << std::fixed << std::setprecision(20)
                << fields[0] << ',' << fields[1] << ',' << fields[3] << ',' << fields[4] << ','
                << result.s12 << ',' << azi1 << ',' << azi2 << ','
                << std::scientific << std::setprecision(20)
                << distance_error << ',' << azi1_error << ',' << azi2_error << '\n';
            ++records;
        }
        catch (const std::exception& exception)
        {
            std::ostringstream message;
            message << std::setprecision(40) << "GeodTest failure at line " << line_number
                << ": " << exception.what() << '\n'
                << "lat1=" << fields[0] << ", lon1=" << fields[1]
                << ", lat2=" << fields[3] << ", lon2=" << fields[4];
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

} // namespace GeodesicInverseTest

#ifndef GEODESIC_INVERSE_TEST_NO_MAIN
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
            std::cerr << "Usage: GeodesicInverse [input.dat [output.CSV]]\n";
            return 1;
        }
        const std::filesystem::path input = argc >= 2
            ? std::filesystem::path(argv[1]) : GeodesicInverseTest::select_karney_test_file();
        const std::filesystem::path output = argc == 3
            ? std::filesystem::path(argv[2]) : std::filesystem::path();
        return GeodesicInverseTest::run_file(input, output);
    }
    catch (const std::exception& exception)
    {
        std::cerr << "ERROR: " << exception.what() << '\n';
        return 1;
    }
}
#endif
