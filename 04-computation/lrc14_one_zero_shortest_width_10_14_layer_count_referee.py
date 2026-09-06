#!/usr/bin/env python3
"""Exact referee for one-zero shortest relations of l1 width 10, 12, or 14.

The analytic tail is a rational, finite signed-pattern calculation.  The
finite H101 base is then rebuilt directly from the carrier inequalities.  No
repository implementation is imported, and all gates survive python -O.
"""

from fractions import Fraction
from itertools import combinations, product
from math import gcd
import sys


sys.stdout.reconfigure(newline="\n")
CHECKS = 0


def gate(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def gcd3(v):
    return gcd(gcd(abs(v[0]), abs(v[1])), abs(v[2]))


def norm1(v):
    return sum(abs(x) for x in v)


def dot(u, v):
    return sum(x * y for x, y in zip(u, v))


def cross(u, v):
    return (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )


def canonical_sign(v):
    for x in v:
        if x:
            return v if x > 0 else tuple(-y for y in v)
    raise ValueError("zero vector")


def one_zero(v):
    return sum(x % 3 == 0 for x in v) == 1


def one_zero_patterns(weight):
    result = set()
    for v in product(range(-weight, weight + 1), repeat=3):
        if norm1(v) == weight and gcd3(v) == 1 and one_zero(v):
            result.add(canonical_sign(v))
    return result


def relations_up_to(bound):
    result = set()
    for v in product(range(-bound, bound + 1), repeat=3):
        weight = norm1(v)
        if weight and weight <= bound and gcd3(v) == 1:
            result.add(canonical_sign(v))
    return sorted(result, key=lambda v: (norm1(v), v))


def bezout_vector(n):
    """Deterministic small q with n.q=1."""
    for radius in range(1, max(abs(x) for x in n) + 1):
        for q in product(range(-radius, radius + 1), repeat=3):
            if dot(n, q) == 1:
                return q
    raise RuntimeError(f"Bezout search failed for {n}")


def feasible_packet(n):
    """Normalize c=1 and return w(t)=base+t*direction, t in (lo,hi)."""
    if n[0]:
        A = Fraction(-n[1], n[0])
        B = Fraction(-n[2], n[0])
        lower, upper = Fraction(0), Fraction(1)
        # a=A*b+B; impose a>0 and b-a>0 in addition to 0<b<1.
        for slope, intercept in ((A, B), (1 - A, -B)):
            if slope > 0:
                lower = max(lower, -intercept / slope)
            elif slope < 0:
                upper = min(upper, -intercept / slope)
            elif intercept <= 0:
                return None
        if lower < upper:
            return (
                (B, Fraction(0), Fraction(1)),
                (A, Fraction(1), Fraction(0)),
                lower,
                upper,
            )
        return None
    if not n[1]:
        return None
    b = Fraction(-n[2], n[1])
    if Fraction(0) < b < Fraction(1):
        return (
            (Fraction(0), b, Fraction(1)),
            (Fraction(1), Fraction(0), Fraction(0)),
            Fraction(0),
            b,
        )
    return None


def fadd(x, y):
    return tuple(x[i] + y[i] for i in range(3))


def fsub(x, y):
    return tuple(x[i] - y[i] for i in range(3))


def fscale(q, x):
    return tuple(q * x[i] for i in range(3))


def feval(x, w):
    return sum(x[i] * w[i] for i in range(3))


E = (
    (Fraction(1), Fraction(0), Fraction(0)),
    (Fraction(0), Fraction(1), Fraction(0)),
    (Fraction(0), Fraction(0), Fraction(1)),
)


def affine_endpoint_forms(n, beta):
    """Relaxed k-roof endpoints for beta*m+k*n.

    If n_i=0, that coordinate can only delete a line, so omitting its constant
    roof gives a valid upper bound for the chord length.
    """
    q = bezout_vector(n)
    m = (
        fsub(fscale(q[2], E[1]), fscale(q[1], E[2])),
        fsub(fscale(q[0], E[2]), fscale(q[2], E[0])),
        fsub(fscale(q[1], E[0]), fscale(q[0], E[1])),
    )
    lower, upper = [], []
    for i in range(3):
        if not n[i]:
            continue
        B = fscale(Fraction(3, 14), fadd(E[(i + 1) % 3], E[(i + 2) % 3]))
        beta_m = fscale(beta, m[i])
        left = fscale(Fraction(1, n[i]), fsub(fscale(-1, B), beta_m))
        right = fscale(Fraction(1, n[i]), fsub(B, beta_m))
        if n[i] > 0:
            lower.append(left)
            upper.append(right)
        else:
            lower.append(right)
            upper.append(left)
    return lower, upper


def as_affine(form, base, direction):
    return feval(form, direction), feval(form, base)


def affine_value(function, t):
    return function[0] * t + function[1]


def relaxed_chord_at(groups, t):
    result = Fraction(0)
    for lower, upper in groups:
        left = max(affine_value(f, t) for f in lower)
        right = min(affine_value(f, t) for f in upper)
        result += max(Fraction(0), right - left)
    return result


def exact_two_chord_supremum(n):
    """Exact sup of ell_1+ell_2 on the closed normalized speed cone."""
    base, direction, interval_left, interval_right = feasible_packet(n)
    groups = []
    breakpoints = {interval_left, interval_right}
    for beta in (1, 2):
        lower_forms, upper_forms = affine_endpoint_forms(n, beta)
        lower = [as_affine(f, base, direction) for f in lower_forms]
        upper = [as_affine(f, base, direction) for f in upper_forms]
        groups.append((lower, upper))
        # Active minima/maxima and the zero/nonzero chord status can change
        # only at an equality among these affine endpoint functions.
        functions = lower + upper
        for i, first in enumerate(functions):
            for second in functions[i + 1 :]:
                if first[0] == second[0]:
                    continue
                t = Fraction(second[1] - first[1], first[0] - second[0])
                if interval_left <= t <= interval_right:
                    breakpoints.add(t)
    values = [(relaxed_chord_at(groups, t), t) for t in breakpoints]
    return max(values)


PATTERNS = {weight: one_zero_patterns(weight) for weight in (10, 12, 14)}
EXPECTED_SHAPES = {
    10: {(2, 3, 5)},
    12: {(0, 1, 11), (0, 5, 7), (1, 2, 9), (1, 3, 8), (1, 5, 6), (2, 3, 7), (3, 4, 5)},
    14: {(0, 1, 13), (1, 1, 12), (1, 3, 10), (1, 4, 9), (1, 6, 7), (3, 4, 7)},
}
EXPECTED_ORIENTED = {10: 24, 12: 144, 14: 120}
EXPECTED_FEASIBLE = {10: 7, 12: 49, 14: 39}
for weight in (10, 12, 14):
    shapes = {tuple(sorted(abs(x) for x in n)) for n in PATTERNS[weight]}
    gate(shapes == EXPECTED_SHAPES[weight], f"shape classification at {weight}")
    gate(len(PATTERNS[weight]) == EXPECTED_ORIENTED[weight], f"oriented count at {weight}")

FEASIBLE = {
    weight: {n for n in PATTERNS[weight] if feasible_packet(n) is not None}
    for weight in (10, 12, 14)
}
for weight in (10, 12, 14):
    gate(len(FEASIBLE[weight]) == EXPECTED_FEASIBLE[weight], f"feasible count at {weight}")

# For M>=8, each chord is shorter than 6c/(7M), so their sum is at most
# 3c/14.  Only the following small-M shapes need the exact arrangement audit.
SPECIAL_SHAPES = {
    (10, (2, 3, 5)): Fraction(13, 70),
    (12, (0, 5, 7)): Fraction(9, 49),
    (12, (1, 5, 6)): Fraction(3, 14),
    (12, (2, 3, 7)): Fraction(4, 21),
    (12, (3, 4, 5)): Fraction(71, 420),
    (14, (1, 6, 7)): Fraction(19, 98),
    (14, (3, 4, 7)): Fraction(121, 588),
}
special_pattern_count = 0
computed_shape_suprema = {}
for weight in (10, 12, 14):
    for n in sorted(FEASIBLE[weight]):
        shape = tuple(sorted(abs(x) for x in n))
        if max(shape) >= 8:
            continue
        special_pattern_count += 1
        value, location = exact_two_chord_supremum(n)
        gate(value <= Fraction(3, 14), f"two-chord bound at n={n}, t={location}")
        key = (weight, shape)
        old = computed_shape_suprema.get(key, Fraction(0))
        computed_shape_suprema[key] = max(old, value)

gate(special_pattern_count == 49, "special signed-pattern count")
gate(computed_shape_suprema == SPECIAL_SHAPES, "special shape supremum table")
gate(Fraction(1, 14) + 2 <= Fraction(1, 11) * 103, "c>=103 tail arithmetic")


def ceil_div(x, y):
    return -((-x) // y)


def integer_k_range(w, n, base):
    lower, upper = -10**30, 10**30
    total = sum(w)
    for i in range(3):
        A = 3 * (total - w[i])
        coefficient, offset = n[i], base[i]
        if coefficient < 0:
            coefficient, offset = -coefficient, -offset
        if not coefficient:
            if 14 * abs(offset) >= A:
                return None
            continue
        lower = max(lower, ceil_div(-A + 1 - 14 * offset, 14 * coefficient))
        upper = min(upper, (A - 1 - 14 * offset) // (14 * coefficient))
    return None if lower > upper else (lower, upper)


def layer_points(w, n, base):
    bounds = integer_k_range(w, n, base)
    if bounds is None:
        return set()
    lower, upper = bounds
    return {
        tuple(base[i] + k * n[i] for i in range(3))
        for k in range(lower, upper + 1)
        if all((base[i] + k * n[i]) % 3 for i in range(3))
    }


def reconstruct_four_layers(w, n):
    q = bezout_vector(n)
    m = cross(w, q)
    gate(dot(n, q) == 1 and cross(n, m) == w, f"basis identity at w={w}, n={n}")
    layers = {}
    for beta in (-2, -1, 1, 2):
        base = tuple(beta * x for x in m)
        residues = [r for r in range(3) if all((base[i] + r * n[i]) % 3 for i in range(3))]
        gate(len(residues) == 1, f"unique live residue at w={w}, n={n}, beta={beta}")
        layers[beta] = layer_points(w, n, base)
    gate(len(layers[1]) == len(layers[-1]), f"beta 1 symmetry at {w}")
    gate(len(layers[2]) == len(layers[-2]), f"beta 2 symmetry at {w}")
    union = set().union(*layers.values())
    half_count = len(layers[1]) + len(layers[2])
    return union, half_count


def direct_carriers(w):
    total = sum(w)
    B = tuple((3 * (total - w[i]) - 1) // 14 for i in range(3))
    result = set()
    for x in range(-B[0], B[0] + 1):
        for y in range(-B[1], B[1] + 1):
            numerator = -w[0] * x - w[1] * y
            if numerator % w[2]:
                continue
            z = numerator // w[2]
            if abs(z) <= B[2] and x % 3 and y % 3 and z % 3:
                result.add((x, y, z))
    return result


def projection_sums(w, carriers):
    cap = Fraction(3, 7 * w[2])
    result = []
    for i in range(3):
        j, k = [r for r in range(3) if r != i]
        result.append(
            sum(
                min(
                    cap,
                    Fraction(3 * (w[j] + w[k]) - 14 * abs(C[i]), 14 * w[j] * w[k]),
                )
                for C in carriers
            )
        )
    return tuple(result)


HEIGHT = 101
VALUES = [x for x in range(1, HEIGHT + 1, 2) if x % 3]
RELATIONS = relations_up_to(14)
eligible = 0
shortest_rows = {10: [], 12: [], 14: []}
line_reconstructions = 0
for w in combinations(VALUES, 3):
    if gcd3(w) != 1:
        continue
    eligible += 1
    relations = [n for n in RELATIONS if dot(n, w) == 0]
    for n in relations:
        gate(norm1(n) % 2 == 0, f"odd relation norm for odd w={w}")
    if not relations:
        continue
    mu = min(norm1(n) for n in relations)
    if mu not in (10, 12, 14):
        continue
    minimizers = [n for n in relations if norm1(n) == mu]
    one_zero_minimizers = [n for n in minimizers if one_zero(n)]
    if not one_zero_minimizers:
        continue
    carriers = direct_carriers(w)
    for n in one_zero_minimizers:
        rebuilt, half_count = reconstruct_four_layers(w, n)
        gate(rebuilt == carriers, f"four-layer/direct mismatch at w={w}, n={n}")
        gate(2 * half_count == len(carriers), f"four-layer count mismatch at w={w}, n={n}")
        line_reconstructions += 1
    automatic = 11 * len(carriers) <= 2 * w[2]
    shortest_rows[mu].append((w, len(carriers), automatic))

gate(len(VALUES) == 34 and eligible == 5937, "H101 finite universe")
gate(len(shortest_rows[10]) == 373, "mu10 base count")
gate(len(shortest_rows[12]) == 1288, "mu12 base count")
gate(len(shortest_rows[14]) == 578, "mu14 base count")
nonautomatic = [
    (mu, w, N)
    for mu in (10, 12, 14)
    for w, N, automatic in shortest_rows[mu]
    if not automatic
]
gate(nonautomatic == [(10, (19, 23, 29), 6)], "unique H101 nonautomatic row")

A2 = (19, 23, 29)
A2_CARRIERS = {
    (-11, -1, 8),
    (-10, 7, 1),
    (-1, -8, 7),
    (1, 8, -7),
    (10, -7, -1),
    (11, 1, -8),
}
a2_carriers = direct_carriers(A2)
a2_sums = projection_sums(A2, a2_carriers)
gate(a2_carriers == A2_CARRIERS, "A2 carrier circuit")
gate(
    a2_sums == (Fraction(156, 4669), Fraction(192, 3857), Fraction(3840, 88711)),
    "A2 projection sums",
)
gate(min(a2_sums) < Fraction(6, 77), "A2 projection target")

max_density = {}
for mu in (10, 12, 14):
    max_density[mu] = max(
        shortest_rows[mu],
        key=lambda row: (Fraction(row[1], row[0][2]), row[0]),
    )

print("LRC14 ONE-ZERO SHORTEST WIDTH 10--14 LAYER-COUNT REFEREE")
print("arithmetic=exact integers/Fraction; standalone_no_repo_import=true; optimizable_assertions=0")
print("scope=sorted primitive distinct positive odd ternary-unit triples w=(a,b,c)")
print("layer_equality=carriers are exactly beta*m+k*n for beta in {-2,-1,1,2}, one k class mod3")
print("analytic=ell_1+ell_2<=3c/14 for every feasible signed pattern of norm 10,12,14")
print("large_coefficient=M>=8 uses two single-roof bounds; special_M<8 signed_patterns=49")
print(
    "special_shape_suprema="
    + ",".join(f"{weight}:{shape}:{value}" for (weight, shape), value in sorted(SPECIAL_SHAPES.items()))
)
print("count_tail=half_count<c/14+2<=c/11 for c>=103; complete remaining base H101")
print(f"finite_base=H{HEIGHT}; values={len(VALUES)}; eligible={eligible}; line_reconstructions={line_reconstructions}")
for mu in (10, 12, 14):
    row = max_density[mu]
    failures = sum(not automatic for _, _, automatic in shortest_rows[mu])
    print(
        f"mu{mu}_shortest_onezero={len(shortest_rows[mu])}; nonautomatic={failures}; "
        f"max_base_density={row[1]}/{row[0][2]}@{row[0]}"
    )
print("unique_nonautomatic=mu10@(19,23,29):N=6; A2=true")
print("A2_projection=(156/4669,192/3857,3840/88711); min<6/77")
print("consequence=every shortest one-zero row with mu<=14 satisfies the THM-4422 projection target")
print(f"checks={CHECKS}")
print("verdict=PASS")
