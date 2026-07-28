#!/usr/bin/env python3
"""Exact referee for THM-2641.

Dependency-free, deterministic, and optimized-mode safe.
"""

from math import gcd, lcm


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def mmul(a, b):
    return tuple(
        tuple(sum(a[i][k] * b[k][j] for k in range(2)) for j in range(2))
        for i in range(2)
    )


def mneg(a):
    return tuple(tuple(-x for x in row) for row in a)


def minv(a):
    det = a[0][0] * a[1][1] - a[0][1] * a[1][0]
    require(det == 1, "matrix is not in SL2(Z)")
    return ((a[1][1], -a[0][1]), (-a[1][0], a[0][0]))


def mpow(a, exponent):
    require(exponent >= 0, "negative exponent")
    out = ((1, 0), (0, 1))
    base = a
    n = exponent
    while n:
        if n & 1:
            out = mmul(out, base)
        base = mmul(base, base)
        n //= 2
    return out


def mmod(a, modulus):
    return tuple(tuple(x % modulus for x in row) for row in a)


def pcanon2(a):
    a = mmod(a, 2)
    return min(a, mmod(mneg(a), 2))


def mvec(a, v, modulus=2):
    return tuple(sum(a[i][j] * v[j] for j in range(2)) % modulus for i in range(2))


def cycle_conjugate(c, t, exponent):
    cp = mpow(c, exponent)
    return mmul(mmul(cp, t), minv(cp))


T = ((1, 1), (0, 1))
S = ((0, -1), (1, 0))
C = ((0, -1), (1, 1))
I = ((1, 0), (0, 1))
require(mpow(S, 2) == mneg(I), "S is not projectively order two")
require(mpow(C, 3) == mneg(I), "C is not projectively order three")
require(mmul(S, C) == mneg(T), "T is not the product of the free-factor generators")

# With chi(T)=epsilon, the other sign convention sends every displayed value
# to its negative.  The unordered fibres and all conclusions are unchanged.
epsilon = 1
require({(3 + 4) % 6, (3 + 2) % 6} == {1, 5}, "free-factor orientations do not make T a C6 generator")

M = lcm(49, 91)
STEP = 2 * M
require(M == 637 and STEP == 1274, "wrong 49/91 invisibility level")
require(STEP == lcm(2, 49, 91), "1274 is not the least simultaneous congruence step")
require(STEP % 6 == 2, "invisible step does not move the ternary class")

parents = tuple(mpow(T, 1 + STEP * j) for j in range(3))
for modulus in (2, 49, 91):
    require(len({mmod(g, modulus) for g in parents}) == 1, f"parents differ modulo {modulus}")
chi_values = tuple((1 + STEP * j) % 6 for j in range(3))
require(set(chi_values) == {1, 3, 5}, "parents do not exhaust one C3 fibre")
for sign in (1, -1):
    require({sign * x % 6 for x in chi_values} == {1, 3, 5}, "sign convention changed the fibre")

# General level-M statement, checked over a representative finite range.
general_checks = 0
for odd_m in range(1, 250, 2):
    if odd_m % 3 == 0:
        continue
    values = {(2 * odd_m * j) % 6 for j in range(3)}
    require(values == {0, 2, 4}, "level-M lifts miss the ternary kernel")
    for j in range(3):
        lift = mpow(T, 2 * odd_m * j)
        require(mmod(lift, 2) == mmod(I, 2), "general lift moved modulo two")
        require(mmod(lift, odd_m) == mmod(I, odd_m), "general lift moved modulo M")
    general_checks += 1

# Conjugacy preserves abelianization but rotates through all three mod-two
# transvections/theta directions.
conjugates = tuple(cycle_conjugate(C, T, j) for j in range(3))
expected = (
    ((1, 1), (0, 1)),
    ((1, 0), (-1, 1)),
    ((2, 1), (-1, 0)),
)
require(conjugates == expected, "wrong integral conjugate triple")
mod2_conjugates = tuple(pcanon2(g) for g in conjugates)
require(len(set(mod2_conjugates)) == 3, "conjugates do not exhaust three mod-two transvections")
directions = ((1, 0), (0, 1), (1, 1))
fixed_directions = []
for g in mod2_conjugates:
    fixed = tuple(v for v in directions if mvec(g, v) == v)
    require(len(fixed) == 1, "a transvection does not select one theta direction")
    fixed_directions.append(fixed[0])
require(set(fixed_directions) == set(directions), "conjugacy does not exhaust theta directions")

# Exact modulus bookkeeping for the proposed spectral refinement.
require(gcd(6, 91) == 1 and lcm(6, 91) == 546, "wrong post-91 C6 modulus")
require(gcd(49, 91) == 7 and lcm(49, 91) == 637, "wrong full ancestry fibre product")
require(lcm(6, 637) == 3822, "wrong full C6/ancestry modulus")
require(lcm(3, 91) == 273 and lcm(3, 637) == 1911, "wrong C3-only refinements")

# Boolean spectral hostile F(x)=1_{{6x}<1/2}.  The analytic coefficient
# formula is zero unless n=6k with k odd (apart from the mean).  Check that
# every class modulo 91 contains such a frequency, while all have residue
# zero modulo six.
inv6 = pow(6, -1, 91)
spectral_witnesses = []
for residue in range(91):
    k = (inv6 * residue) % 91
    if k % 2 == 0:
        k += 91
    n = 6 * k
    require(k % 2 == 1 and n % 91 == residue and n % 6 == 0, "bad spectral hostile witness")
    spectral_witnesses.append(n)
require(len(set(n % 91 for n in spectral_witnesses)) == 91, "hostile misses a mod-91 class")
require({n % 6 for n in spectral_witnesses} == {0}, "hostile leaks into another C6 class")

print("THM-2641 MODULAR CONNECTOR EXACT REFEREE")
print(f"free_factors=S_order2,C_order3 T_class=+-1_in_C6")
print(f"invisibility_level M={M} step={STEP} residues=(2,49,91)")
print(f"positive_parent_exponents={(1, 1 + STEP, 1 + 2 * STEP)} chi_C6={chi_values}")
print(f"general_odd_M_not3_checks={general_checks} ternary_kernel=surjective")
print(f"theta_conjugates={conjugates}")
print(f"theta_fixed_directions={tuple(fixed_directions)} common_abelian_class=+-1")
print("spectral_moduli post91=546 full637=3822 C3_only=(273,1911)")
print(f"boolean_hostile_mod91_classes={len(spectral_witnesses)} C6_support={{0}} max_frequency={max(spectral_witnesses)}")
print("ALL CHECKS PASSED")
