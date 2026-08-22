#!/usr/bin/env python3
"""Exact controls for THM-3210's factorial exterior double-cancellation ray.

Extends THM-3186's exit-time amplitude E_L.  Exact symbolic and rational
arithmetic only; every gate is an explicit ``require`` so that ordinary and
``-O`` replay are byte-identical.
"""

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


n_sym, d_sym, v_sym = sp.symbols("n d v")


def amplitude(length, index, scalar_d, weight_v):
    """THM-3186 (21): E_L = sum_j c_(n+j) C_(j-1) prod_(h>j) u_(n+h)."""

    delta = 1 - 4 * scalar_d * weight_v
    a_of = lambda i: 2 * (i + 1) * (2 * i + 1) * weight_v
    b_of = lambda i: i * (i + 1) * delta
    c_of = lambda i: scalar_d - i - 1
    u_of = lambda i: -b_of(i)
    alpha_of = lambda i: a_of(i) * scalar_d
    beta_of = lambda i: b_of(i) * scalar_d
    continuant = {0: sp.Integer(1), 1: alpha_of(index + 1)}
    for r in range(2, length + 1):
        continuant[r] = sp.expand(alpha_of(index + r) * continuant[r - 1]
                                  + scalar_d * beta_of(index + r)
                                  * continuant[r - 2])
    total = sp.Integer(0)
    for j in range(1, length):
        tail = sp.Integer(1)
        for h in range(j + 1, length):
            tail *= u_of(index + h)
        total += c_of(index + j) * continuant[j - 1] * tail
    return sp.expand(total)


# ------------------------------------------- 1.  degree and the length-3 locus

DEGREE_ROWS = []
for index in (1, 2, 3, 4):
    degrees = tuple(sp.Poly(amplitude(length, index, d_sym, v_sym),
                            v_sym).degree()
                    for length in (2, 3, 4, 5, 6))
    require(degrees == (0, 1, 2, 3, 4),
            "deg_v E_L = L - 2")
    DEGREE_ROWS.append((index, degrees))

LENGTH3_LOCUS = ((n_sym + 3) * (d_sym - n_sym - 2)
                 / (2 * d_sym * ((4 * n_sym + 9) * d_sym
                                 - (n_sym + 3) * (4 * n_sym + 7))))

LOCUS_ROWS = []
for index in range(1, 10):
    solutions = sp.solve(amplitude(3, index, d_sym, v_sym), v_sym)
    require(len(solutions) == 1, "E_3 is affine in v with one root")
    require(sp.simplify(solutions[0] - LENGTH3_LOCUS.subs(n_sym, index)) == 0,
            "closed form of the length-three cancellation locus")
    LOCUS_ROWS.append((index, sp.simplify(solutions[0])))

# THM-3186's published hostile is the n=1 point of that locus at d=5.
require(sp.nsimplify(LENGTH3_LOCUS.subs({n_sym: 1, d_sym: 5}))
        == sp.Rational(4, 105),
        "THM-3186 hostile lies on the length-three locus")


# ------------------------------------------------- 2.  the double-cancellation ray

RAY_D = n_sym + 4
RAY_V = (n_sym + 3) / (3 * (n_sym + 4) * (2 * n_sym + 5))

RAY_E2 = sp.simplify(amplitude(2, n_sym, RAY_D, RAY_V))
RAY_E3 = sp.simplify(amplitude(3, n_sym, RAY_D, RAY_V))
RAY_E4 = sp.simplify(amplitude(4, n_sym, RAY_D, RAY_V))
RAY_E5 = sp.factor(sp.simplify(amplitude(5, n_sym, RAY_D, RAY_V)))

require(RAY_E2 == 2, "E_2 = 2 on the ray")
require(RAY_E3 == 0, "E_3 vanishes identically on the ray")
require(RAY_E4 == 0, "E_4 vanishes identically on the ray")
require(RAY_E5 != 0, "E_5 does not vanish on the ray")

RAY_DELTA = sp.factor(sp.simplify(1 - 4 * RAY_D * RAY_V))
RAY_BETA = sp.factor(sp.simplify(n_sym * (n_sym + 1) * RAY_DELTA * RAY_D))
RAY_C1 = sp.simplify(RAY_D - n_sym - 2)
require(sp.simplify(RAY_DELTA - (2 * n_sym + 3) / (3 * (2 * n_sym + 5))) == 0,
        "Delta on the ray")
require(RAY_C1 == 2, "c_(n+1) = 2 on the ray")

RAY_E5_COFACTOR = sp.expand(10 * n_sym ** 3 + 101 * n_sym ** 2
                            + 336 * n_sym + 366)
require(sp.simplify(RAY_E5
                    + 4 * (n_sym + 2) * (n_sym + 3) ** 2 * (n_sym + 4)
                    * (2 * n_sym + 3) * RAY_E5_COFACTOR
                    / (27 * (2 * n_sym + 5) ** 2)) == 0,
        "closed form of E_5 on the ray")

# Every factor of Delta, beta_n, v and E_5 is strictly positive for n >= 1,
# so the ray is admissible under THM-3183's standing hypothesis and the
# invisible window is genuinely two steps wide.
POSITIVITY_ROWS = []
for index in range(1, 41):
    delta_value = sp.Rational(2 * index + 3, 3 * (2 * index + 5))
    v_value = sp.Rational(index + 3, 3 * (index + 4) * (2 * index + 5))
    beta_value = sp.nsimplify(RAY_BETA.subs(n_sym, index))
    e5_value = sp.nsimplify(RAY_E5.subs(n_sym, index))
    require(delta_value > 0 and v_value > 0 and beta_value > 0,
            "standing hypothesis holds on the ray")
    require(e5_value != 0, "E_5 nonzero on the ray")
    require(sp.simplify(amplitude(3, index, index + 4, v_value)) == 0
            and sp.simplify(amplitude(4, index, index + 4, v_value)) == 0,
            "numeric replay of the double cancellation")
    POSITIVITY_ROWS.append(index)
require(len(POSITIVITY_ROWS) == 40, "positivity bank is complete")

# Visibility profile V_L = -beta_n E_L: visible, invisible, invisible, visible.
VISIBILITY = tuple("visible" if value != 0 else "invisible"
                   for value in (RAY_E2, RAY_E3, RAY_E4, RAY_E5))
require(VISIBILITY == ("visible", "invisible", "invisible", "visible"),
        "non-monotone visibility profile at lengths 2,3,4,5")


# ------------------------------- 3.  no triple cancellation, and uniqueness

TRIPLE_ROWS = []
UNIQUENESS_ROWS = []
for index in (1, 2, 3):
    poly3 = sp.Poly(amplitude(3, index, d_sym, v_sym), v_sym)
    poly4 = sp.Poly(amplitude(4, index, d_sym, v_sym), v_sym)
    poly5 = sp.Poly(amplitude(5, index, d_sym, v_sym), v_sym)
    res34 = sp.Poly(sp.resultant(poly3, poly4, v_sym), d_sym)
    res35 = sp.Poly(sp.resultant(poly3, poly5, v_sym), d_sym)
    shared = sp.factor(sp.gcd(res34, res35).as_expr())
    require(sp.simplify(shared / d_sym ** 2).is_number,
            "the only common resultant root is the excluded d = 0")
    TRIPLE_ROWS.append((index, str(shared)))

    factored = sp.factor(res34.as_expr())
    rational_roots = [root for root in sp.solve(sp.Eq(factored, 0), d_sym)
                      if root.is_rational and root != 0]
    require(rational_roots == [sp.Integer(index + 4)],
            "d = n+4 is the unique nonzero rational double-cancellation ray")
    UNIQUENESS_ROWS.append((index, [int(root) for root in rational_roots]))

# Positive control: off the ray, length four is visible.
CONTROL_INDEX = 1
CONTROL_V = sp.Rational(4, 105)
require(sp.simplify(amplitude(3, CONTROL_INDEX, 5, CONTROL_V)) == 0,
        "control: on-ray length three cancels")
require(sp.simplify(amplitude(4, CONTROL_INDEX, 6, CONTROL_V)) != 0,
        "control: moving d off the ray restores length-four visibility")


print("THM-3210 FACTORIAL EXTERIOR DOUBLE-CANCELLATION RAY EXACT CONTROL")
print("amplitude_degree_rows=" + repr(DEGREE_ROWS))
print("length3_locus=v=(n+3)(d-n-2)/(2d[(4n+9)d-(n+3)(4n+7)])")
print("length3_locus_rows=" + repr([(i, str(x)) for i, x in LOCUS_ROWS]))
print("thm3186_hostile_is_locus_point=(n=1,d=5,v=4/105)")
print("ray=(d=n+4, v=(n+3)/(3(n+4)(2n+5)))")
print("ray_Delta=" + str(RAY_DELTA))
print("ray_beta_n=" + str(RAY_BETA))
print("ray_c_(n+1)=2")
print("ray_E2=2  ray_E3=0  ray_E4=0")
print("ray_E5=" + str(RAY_E5))
print("visibility_profile_lengths_2_3_4_5=" + repr(VISIBILITY))
print("positivity_bank_n=1..40")
print("triple_cancellation_shared_resultant=" + repr(TRIPLE_ROWS))
print("rational_double_ray_uniqueness=" + repr(UNIQUENESS_ROWS))
print("scope=exit_time_amplitude_only_no_PRS_depth_or_GMC2_consequence")
print("ALL EXACT CHECKS PASSED")
