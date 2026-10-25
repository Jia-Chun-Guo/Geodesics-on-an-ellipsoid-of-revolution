// ============================================================================
//
// Copyright (c) 2026, Anhui Jianzhu University, China
// Author: Jia-Chun Guo
// Contact: jchguo@ahjzu.edu.cn
//
// License: This source code is released under the MIT License.
//          https://opensource.org/licenses/MIT
//
// ============================================================================
#pragma once
// ============================================================
// NearAntipodalClairaut.hpp
// Determine the Clairaut constant and geometric branch
// for near-antipodal inverse geodesics.
// Candidate branches follow Eq.(91); longitude matching uses Eq.(92).
// Safeguarded Newton iteration uses u = sqrt(1-c/c_max)
// and analytic derivatives of the finite longitude series.
// Angles and longitude residuals are expressed in radians.
// Mathematical definitions and series coefficients are provided by
// GeodesicDefinitions.hpp and GeodesicCore.hpp.
// ============================================================

// =================== includes ===================
#include <algorithm>
#include <array>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <boost/math/special_functions/fpclassify.hpp>

#include "GeodesicDefinitions.hpp"
#include "GeodesicCore.hpp"

namespace NearAntipodal {

using Geodesic_Defs::real;
using Geodesic_Defs::PI;
using Geodesic_Defs::phi_of_beta;

using boost::multiprecision::abs;
using boost::multiprecision::cos;

using GeodesicCore::DL;
using GeodesicCore::dL;
using GeodesicCore::Epsilon;

// =====================================================
// Near-antipodal branch identifiers, Eq.(90).
// =====================================================
enum class LCBranch
{
    A = 1,     // Same side of a vertex
    B = 2,     // Across the northern vertex
    C = 3      // Across the southern vertex
};

// =====================================================
// Result of the Clairaut constant search.
// =====================================================
struct ClairautRootResult
{
    real c;
    LCBranch br;
    real resid;       // Absolute longitude residual (radians)
    int iters;
    bool converged;
};

// =====================================================
// Restrict a trial value to the given interval.
// =====================================================
inline real clamp_real(real x, real lo, real hi)
{
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

namespace detail {

// Longitude-series coefficients and their derivatives with respect to c.
using Coefficients = GeodesicCore::LongitudeCoefficientsWithDerivative;

struct Endpoint
{
    real sine;     // abs(sin(beta)).
    real cosine;   // cos(beta).
    bool southern; // True for beta < 0.

    explicit Endpoint(const real& beta)
        : sine(abs(boost::multiprecision::sin(beta))),
          cosine(boost::multiprecision::cos(beta)), southern(beta < 0)
    {}
};

struct EndpointValue
{
    real longitude;
    real derivative;
    bool derivative_valid;
};

inline EndpointValue north_longitude(
    const Endpoint& endpoint, const real& c, const Coefficients& coefficients)
{
    const real half_period = PI - coefficients.epsilon;
    // The equatorial-node limit is also regular at c = 1.
    if (endpoint.sine == 0)
        return {half_period / real(2), -coefficients.epsilon_derivative / real(2), true};

    if (c == endpoint.cosine)
        return {real(0), real(0), false};

    const real q2 = (endpoint.cosine - c) * (endpoint.cosine + c);
    if (q2 < 0)
        throw std::domain_error("Clairaut constant outside the admissible interval.");
    const real q = boost::multiprecision::sqrt(q2);
    const real g = (real(1) - c) * (real(1) + c);
    const real sine2 = real(2) * endpoint.sine * q / g;
    const real cosine2 = (endpoint.sine * endpoint.sine - q2) / g;

    // Evaluate sum f_j*sin(2*j*phi), sum f_j'(c)*sin(2*j*phi),
    // and sum 2*j*f_j*cos(2*j*phi) by Clenshaw recurrences.
    real b1 = 0, b2 = 0, d1 = 0, d2 = 0, a1 = 0, a2 = 0;
    const real twice_cosine = real(2) * cosine2;
    for (int j = 7; j >= 0; --j)
    {
        const real b0 = coefficients.harmonic[j] + twice_cosine * b1 - b2;
        const real d0 = coefficients.harmonic_derivative[j] + twice_cosine * d1 - d2;
        const real a0 = real(2 * (j + 1)) * coefficients.harmonic[j]
            + twice_cosine * a1 - a2;
        b2 = b1; b1 = b0;
        d2 = d1; d1 = d0;
        a2 = a1; a1 = a0;
    }

    const real phi = boost::multiprecision::atan2(q, endpoint.sine);
    const real spherical = (c == 0) ? PI / real(2)
        : real(boost::multiprecision::atan2(q, c * endpoint.sine));
    const real amplitude_sum = phi + sine2 * b1;
    const real scaled_epsilon = coefficients.epsilon / PI;
    const real longitude = spherical - scaled_epsilon * amplitude_sum;
    const real base_derivative = -endpoint.sine / (q * g);
    const real amplitude_derivative = c * base_derivative;
    const real derivative = base_derivative
        - coefficients.epsilon_derivative * amplitude_sum / PI
        - scaled_epsilon * (amplitude_derivative * (real(1) + cosine2 * a1 - a2)
            + sine2 * d1);
    return {longitude, derivative, true};
}

struct Trial
{
    real c;
    std::array<real, 3> longitude;
    std::array<real, 3> derivative;
    bool derivative_valid;
};

inline unsigned branch_index(LCBranch branch)
{
    const unsigned value = static_cast<unsigned>(branch);
    if (value < 1 || value > 3)
        throw std::invalid_argument("Invalid near-antipodal branch identifier.");
    return value - 1;
}

struct Context
{
    Endpoint first;
    Endpoint second;
    real c_max;
    LCBranch vertex_branch;
    std::array<Trial, 3> common;
    unsigned common_count = 0;
    mutable Trial last;
    mutable bool last_valid = false;

    Context(const real& beta1, const real& beta2)
        : first(beta1), second(beta2), c_max((std::min)(first.cosine, second.cosine)),
          vertex_branch(beta2 >= -beta1 ? LCBranch::B : LCBranch::C)
    {}

    Trial compute(const real& c) const
    {
        if (!boost::math::isfinite(c) || c < 0 || c > c_max)
            throw std::domain_error("Clairaut constant outside the endpoint interval.");
        const Coefficients coefficients(c);
        const real half_period = PI - coefficients.epsilon;
        auto x = north_longitude(first, c, coefficients);
        auto y = north_longitude(second, c, coefficients);
        if (first.southern)
        {
            x.longitude = half_period - x.longitude;
            x.derivative = -coefficients.epsilon_derivative - x.derivative;
        }
        if (second.southern)
        {
            y.longitude = half_period - y.longitude;
            y.derivative = -coefficients.epsilon_derivative - y.derivative;
        }
        const real difference = x.longitude - y.longitude;
        const real sum = x.longitude + y.longitude;
        const real sum_derivative = sum >= 0
            ? real(x.derivative + y.derivative) : real(-x.derivative - y.derivative);
        const real complement = real(2) * half_period - abs(sum);
        const real complement_derivative = -real(2) * coefficients.epsilon_derivative - sum_derivative;
        return {c,
            {abs(difference), abs(sum), abs(complement)},
            {difference >= 0 ? real(x.derivative - y.derivative) : real(y.derivative - x.derivative),
             sum_derivative,
             complement >= 0 ? complement_derivative : real(-complement_derivative)},
            x.derivative_valid && y.derivative_valid};
    }

    const Trial& evaluate(const real& c) const
    {
        for (unsigned i = 0; i < common_count; ++i)
            if (c == common[i].c) return common[i];
        if (last_valid && c == last.c) return last;
        last = compute(c);
        last_valid = true;
        return last;
    }

    void prepare_common(const real& initial)
    {
        common_count = 0;
        last_valid = false;
        for (const real& c : {real(0), c_max, clamp_real(initial, real(0), c_max)})
        {
            bool present = false;
            for (unsigned i = 0; i < common_count; ++i)
                present = present || common[i].c == c;
            if (!present) common[common_count++] = compute(c);
        }
    }
};

} // namespace detail

inline real dlon_of_c(LCBranch branch, real c, real beta1, real beta2)
{
    const detail::Context context(beta1, beta2);
    return context.compute(c).longitude[detail::branch_index(branch)];
}

// =====================================================
// Evaluate the signed longitude residual, Eq.(92), and
// retain the trial with the smallest absolute residual.
// =====================================================
inline real longitude_residual(
    LCBranch br, real c,
    real beta1, real beta2, real dlon12,
    ClairautRootResult& result
)
{
    real R = dlon_of_c(br, c, beta1, beta2) - dlon12;
    real r = abs(R);

    if (r < result.resid)
    {
        result.c = c;
        result.resid = r;
    }

    return R;
}

// =====================================================
// Solve the longitude-matching equation by safeguarded Newton iteration
// in u, using analytic derivatives and an interval-midpoint fallback.
// Convergence requires abs(R_j(c)) <= tolerance.
// =====================================================
namespace detail {

// c = c_max*(1-u^2), evaluated in factored form.
// u = 0 is the vertex limit; u = 1 corresponds to c = 0.
inline real c_of_u(const real& u, const real& c_max)
{
    return c_max * (real(1) - u) * (real(1) + u);
}

inline real u_of_c(const real& c, const real& c_max)
{
    return boost::multiprecision::sqrt((c_max - c) / c_max);
}

inline ClairautRootResult solve_branch(
    const Context& context, LCBranch branch,
    const real& target, const real& initial,
    int max_iterations, const real& tolerance)
{
    ClairautRootResult best{real(0), branch,
        std::numeric_limits<real>::infinity(), 0, false};
    if (context.c_max <= 0) return best;
    const unsigned index = branch_index(branch);
    auto residual = [&](const real& c) -> real {
        const real value = context.evaluate(c).longitude[index] - target;
        const real magnitude = abs(value);
        if (magnitude < best.resid)
        {
            best.c = c;
            best.resid = magnitude;
            best.converged = magnitude <= tolerance;
        }
        return value;
    };
    // Test both physical endpoints and the initial estimate.
    // Accept a trial when its absolute longitude residual meets the tolerance.
    real lower = 0, upper = 1;
    real r_lower = residual(context.c_max);
    if (best.converged) return best;
    real r_upper = residual(real(0));
    if (best.converged) return best;
    const real c0 = clamp_real(initial, real(0), context.c_max);
    residual(c0);
    if (best.converged) return best;

    // Monotonicity excludes a root on an eligible shortest branch
    // when its endpoint residuals have the same nonzero sign.
    if (r_lower * r_upper > 0
        && (branch == LCBranch::A || branch == context.vertex_branch))
        return best;

    // For an explicitly requested vertex branch outside the candidate set,
    // scan the interval for a sign change.
    // geometrically excluded vertex branch.
    if (r_lower * r_upper > 0)
    {
        constexpr int count = 8;
        real previous = lower, previous_residual = r_lower;
        for (int i = 1; i <= count; ++i)
        {
            const real u = real(i) / real(count);
            const real value = (i == count) ? r_upper : residual(c_of_u(u, context.c_max));
            if (best.converged) return best;
            if (previous_residual * value <= 0)
            {
                lower = previous; r_lower = previous_residual;
                upper = u; r_upper = value;
                break;
            }
            previous = u;
            previous_residual = value;
        }
    }

    real u = u_of_c(best.c, context.c_max);
    if (!(u > lower && u < upper))
    {
        // Use an interior endpoint-secant estimate when available;
        // otherwise start from the interval midpoint.
        u = (lower + upper) / real(2);
        if (r_lower * r_upper < 0)
        {
            const real secant = lower - r_lower * (upper - lower) / (r_upper - r_lower);
            if (secant > lower && secant < upper) u = secant;
        }
    }
    real c = c_of_u(u, context.c_max);
    real value = residual(c);
    real previous_width = upper - lower;
    real previous_residual = abs(value);
    bool previous_newton = false;

    for (int iteration = 0; iteration < max_iterations && !best.converged; ++iteration)
    {
        best.iters = iteration + 1;
        if (r_lower * value <= 0)
        {
            upper = u; r_upper = value;
        }
        else if (r_upper * value <= 0)
        {
            lower = u; r_lower = value;
        }
        const real width = upper - lower;
        const bool force_midpoint = previous_newton
            && r_lower * r_upper <= 0
            && width > previous_width / real(2)
            && abs(value) > previous_residual / real(2);
        previous_width = width;

        real next = lower + width / real(2);
        real next_value = 0;
        bool accepted_newton = false;
        if (!force_midpoint)
        {
            const auto& trial = context.evaluate(c);
            const real derivative = -real(2) * context.c_max * u * trial.derivative[index];
            if (trial.derivative_valid && derivative != 0 && boost::math::isfinite(derivative))
            {
                const real candidate = u - value / derivative;
                if (boost::math::isfinite(candidate) && candidate > lower && candidate < upper)
                {
                    const real candidate_value = residual(c_of_u(candidate, context.c_max));
                    if (best.converged) break;
                    if (abs(candidate_value) <= abs(value))
                    {
                        next = candidate;
                        next_value = candidate_value;
                        accepted_newton = true;
                    }
                }
            }
        }
        if (!accepted_newton)
            next_value = residual(c_of_u(next, context.c_max));
        if (best.converged) break;
        if (next == u || next == lower || next == upper) break;
        previous_newton = accepted_newton;
        previous_residual = abs(value);
        u = next;
        c = c_of_u(u, context.c_max);
        value = next_value;
    }
    return best;
}

inline void validate(
    const real& beta1, const real& beta2, const real& target,
    const real& initial, int maximum, const real& tolerance)
{
    if (!boost::math::isfinite(beta1) || !boost::math::isfinite(beta2)
        || !boost::math::isfinite(target) || !boost::math::isfinite(initial)
        || !boost::math::isfinite(tolerance)
        || abs(beta1) >= PI / real(2) || abs(beta2) >= PI / real(2)
        || target < 0 || target > PI || maximum < 0 || tolerance <= 0)
        throw std::invalid_argument("Invalid near-antipodal solver input.");
}

} // namespace detail

// Search the specified branch and return its best trial,
// absolute longitude residual, and convergence status.
inline ClairautRootResult solve_c_near_antipodal_one_branch(
    LCBranch branch, real beta1, real beta2, real dlon12, real c0,
    int max_iter = 256, real tol_dlon = real("1e-30"))
{
    detail::validate(beta1, beta2, dlon12, c0, max_iter, tol_dlon);
    detail::Context context(beta1, beta2);
    context.prepare_common(c0);
    return detail::solve_branch(context, branch, dlon12, c0, max_iter, tol_dlon);
}

// Select an eligible branch geometrically and solve the longitude equation.
// Throw std::runtime_error if neither eligible branch meets the residual tolerance.
inline ClairautRootResult solve_c_near_antipodal_all_branches(
    real beta1, real beta2, real dlon12, real c0,
    int max_iter = 256, real tol_dlon = real("1e-30"))
{
    detail::validate(beta1, beta2, dlon12, c0, max_iter, tol_dlon);
    detail::Context context(beta1, beta2);
    context.prepare_common(c0);
    const auto vertex_branch = context.vertex_branch;

    // For nonsymmetric endpoints, c_max is the endpoint-at-vertex
    // limit of Case A. The retained vertex-branch expression gives
    // the same longitude at this boundary. Accept the Case A limit
    // if its longitude residual is within tolerance; otherwise,
    // use the residual sign to select the branch.
    const real vertex_residual = context.evaluate(context.c_max).longitude[0] - dlon12;

    if (abs(vertex_residual) <= tol_dlon)
        return {context.c_max, LCBranch::A, abs(vertex_residual), 0, true};

    const LCBranch selected = vertex_residual > 0 ? LCBranch::A : vertex_branch;

    auto best = detail::solve_branch(context, selected, dlon12, c0, max_iter, tol_dlon);

    // Try the other eligible branch if the selected search does not converge.
    if (!best.converged)
    {
        const LCBranch alternative = selected == LCBranch::A ? vertex_branch : LCBranch::A;
        auto candidate = detail::solve_branch(context, alternative, dlon12, c0, max_iter, tol_dlon);
        if (candidate.resid < best.resid) best = candidate;
    }
    if (!best.converged)
    {
        std::ostringstream message;
        message << "Near-antipodal inverse solution did not converge: "
            << "best longitude residual = " << best.resid
            << " rad, tolerance = " << tol_dlon << " rad.";
        throw std::runtime_error(message.str());
    }
    return best;
}

} // namespace NearAntipodal
