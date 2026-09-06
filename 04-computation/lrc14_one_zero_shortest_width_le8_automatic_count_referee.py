#!/usr/bin/env python3
"""Independent exact referee for the one-zero shortest-width <= 8 closure.

This file is deliberately standalone: it imports no repository implementation.
All comparisons use integer or Fraction arithmetic, and every verification gate
remains active under ``python -O``.
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
    raise ValueError("zero vector has no canonical sign")


def one_zero_mod_three(v):
    return sum(x % 3 == 0 for x in v) == 1


def primitive_one_zero_patterns(weight):
    out = set()
    for v in product(range(-weight, weight + 1), repeat=3):
        if norm1(v) != weight or gcd3(v) != 1:
            continue
        if not one_zero_mod_three(v):
            continue
        out.add(canonical_sign(v))
    return out


def all_primitive_relations_up_to(bound):
    out = set()
    for v in product(range(-bound, bound + 1), repeat=3):
        weight = norm1(v)
        if not weight or weight > bound or gcd3(v) != 1:
            continue
        out.add(canonical_sign(v))
    return sorted(out, key=lambda v: (norm1(v), v))


def sorted_positive_feasible_interval(n):
    """Exact open parameter interval for n.(a,b,c)=0, 0<a<b<c, c=1.

    Returns w(t)=base+t*direction and the open interval (lower, upper), or
    None.  This is a one-dimensional rational feasibility calculation, not a
    sampled search.
    """
    if n[0]:
        # t=b and a=A*t+B.
        A = Fraction(-n[1], n[0])
        B = Fraction(-n[2], n[0])
        lower, upper = Fraction(0), Fraction(1)
        # a>0 and b-a>0.
        for slope, intercept in ((A, B), (1 - A, -B)):
            if slope > 0:
                lower = max(lower, -intercept / slope)
            elif slope < 0:
                upper = min(upper, -intercept / slope)
            elif intercept <= 0:
                return None
        if not lower < upper:
            return None
        return (
            (B, Fraction(0), Fraction(1)),
            (A, Fraction(1), Fraction(0)),
            lower,
            upper,
        )
    if not n[1]:
        return None
    b = Fraction(-n[2], n[1])
    if not Fraction(0) < b < Fraction(1):
        return None
    # t=a in (0,b).
    return (
        (Fraction(0), b, Fraction(1)),
        (Fraction(1), Fraction(0), Fraction(0)),
        Fraction(0),
        b,
    )


def form_add(x, y):
    return tuple(x[i] + y[i] for i in range(3))


def form_sub(x, y):
    return tuple(x[i] - y[i] for i in range(3))


def form_scale(q, x):
    return tuple(q * x[i] for i in range(3))


def form_eval(x, w):
    return sum(x[i] * w[i] for i in range(3))


E = (
    (Fraction(1), Fraction(0), Fraction(0)),
    (Fraction(0), Fraction(1), Fraction(0)),
    (Fraction(0), Fraction(0), Fraction(1)),
)


def forms_equivalent_on_relation(x, y, n):
    """Whether x-y is a rational multiple of the relation n."""
    difference = form_sub(x, y)
    multiplier = None
    for i in range(3):
        if n[i]:
            q = difference[i] / n[i]
            if multiplier is None:
                multiplier = q
            elif multiplier != q:
                return False
        elif difference[i]:
            return False
    return True


def roof_endpoint_forms(n):
    """Return q, m=w cross q, and lower/upper k-roof forms."""
    unit_index = next(i for i, x in enumerate(n) if abs(x) == 1)
    q = [0, 0, 0]
    q[unit_index] = n[unit_index]
    # Linear-form version of m=w cross q.
    m = (
        form_sub(form_scale(q[2], E[1]), form_scale(q[1], E[2])),
        form_sub(form_scale(q[0], E[2]), form_scale(q[2], E[0])),
        form_sub(form_scale(q[1], E[0]), form_scale(q[0], E[1])),
    )
    lower, upper = [None] * 3, [None] * 3
    for i in range(3):
        if not n[i]:
            continue
        B = form_scale(Fraction(3, 14), form_add(E[(i + 1) % 3], E[(i + 2) % 3]))
        minus_endpoint = form_sub(form_scale(-1, B), m[i])
        plus_endpoint = form_sub(B, m[i])
        if n[i] > 0:
            lower[i] = form_scale(Fraction(1, n[i]), minus_endpoint)
            upper[i] = form_scale(Fraction(1, n[i]), plus_endpoint)
        else:
            lower[i] = form_scale(Fraction(1, n[i]), plus_endpoint)
            upper[i] = form_scale(Fraction(1, n[i]), minus_endpoint)
    return tuple(q), m, lower, upper


def linear_form_nonnegative_on_feasible_cone(form, n):
    packet = sorted_positive_feasible_interval(n)
    if packet is None:
        return False
    base, direction, lower, upper = packet
    intercept = form_eval(form, base)
    slope = form_eval(form, direction)
    return intercept + slope * lower >= 0 and intercept + slope * upper >= 0


def ceil_div(x, y):
    return -((-x) // y)


def integer_k_range(w, n, m):
    """Exact integer k range cut out by every strict carrier roof."""
    lower, upper = -10**30, 10**30
    total = sum(w)
    for i in range(3):
        A = 3 * (total - w[i])
        coefficient, offset = n[i], m[i]
        if coefficient < 0:
            coefficient, offset = -coefficient, -offset
        if not coefficient:
            if 14 * abs(offset) >= A:
                return None
            continue
        # -A < 14(offset+k*coefficient) < A, with integral sides.
        lower = max(lower, ceil_div(-A + 1 - 14 * offset, 14 * coefficient))
        upper = min(upper, (A - 1 - 14 * offset) // (14 * coefficient))
    return None if lower > upper else (lower, upper)


def layer_points(w, n, m):
    bounds = integer_k_range(w, n, m)
    if bounds is None:
        return set()
    lower, upper = bounds
    out = set()
    for k in range(lower, upper + 1):
        C = tuple(m[i] + k * n[i] for i in range(3))
        if all(x % 3 for x in C):
            out.add(C)
    return out


def line_reconstruction(w, n):
    unit_index = next(i for i, x in enumerate(n) if abs(x) == 1)
    q = [0, 0, 0]
    q[unit_index] = n[unit_index]
    m = cross(w, q)
    gate(dot(n, q) == 1, f"Bezout failure for {n}")
    gate(cross(n, m) == w, f"basis failure for w={w}, n={n}")
    residues = [r for r in range(3) if all((m[i] + r * n[i]) % 3 for i in range(3))]
    gate(len(residues) == 1, f"live residue is not unique for w={w}, n={n}")
    positive = layer_points(w, n, m)
    negative = layer_points(w, n, tuple(-x for x in m))
    gate(len(positive) == len(negative), f"opposite layer imbalance for w={w}, n={n}")
    return positive | negative, len(positive), residues[0]


def direct_carriers(w):
    total = sum(w)
    B = tuple((3 * (total - w[i]) - 1) // 14 for i in range(3))
    out = set()
    for x in range(-B[0], B[0] + 1):
        for y in range(-B[1], B[1] + 1):
            numerator = -w[0] * x - w[1] * y
            if numerator % w[2]:
                continue
            z = numerator // w[2]
            if abs(z) > B[2] or x % 3 == 0 or y % 3 == 0 or z % 3 == 0:
                continue
            out.add((x, y, z))
    return out


def projection_sums(w, carriers):
    q = Fraction(3, 7 * w[2])
    result = []
    for i in range(3):
        j, k = [r for r in range(3) if r != i]
        result.append(
            sum(
                min(
                    q,
                    Fraction(3 * (w[j] + w[k]) - 14 * abs(C[i]), 14 * w[j] * w[k]),
                )
                for C in carriers
            )
        )
    return tuple(result)


# Exact coefficient-pattern classification.
P2 = primitive_one_zero_patterns(2)
P4 = primitive_one_zero_patterns(4)
P6 = primitive_one_zero_patterns(6)
P8 = primitive_one_zero_patterns(8)
gate({tuple(sorted(abs(x) for x in n)) for n in P2} == {(0, 1, 1)}, "norm-two shape")
gate(not P4, "a primitive one-zero norm-four pattern unexpectedly exists")
gate(
    {tuple(sorted(abs(x) for x in n)) for n in P6} == {(0, 1, 5), (1, 2, 3)},
    "norm-six shapes",
)
gate(
    {tuple(sorted(abs(x) for x in n)) for n in P8}
    == {(0, 1, 7), (1, 1, 6), (1, 3, 4)},
    "norm-eight shapes",
)
gate(len(P2) == 6 and len(P6) == 36 and len(P8) == 48, "oriented pattern counts")

FEASIBLE6 = {n for n in P6 if sorted_positive_feasible_interval(n) is not None}
EXPECTED_FEASIBLE6 = {
    (0, 5, -1),
    (1, -3, 2),
    (1, 3, -2),
    (2, -3, 1),
    (2, 3, -1),
    (3, -2, 1),
    (3, 1, -2),
    (3, 2, -1),
    (5, -1, 0),
    (5, 0, -1),
}
gate(FEASIBLE6 == EXPECTED_FEASIBLE6, "sorted norm-six signed-placement classification")
gate(not any(sorted_positive_feasible_interval(n) for n in P2), "norm-two would repeat a speed")

# Each table row chooses one upper roof U_i and one lower roof L_j.  Their
# difference is an upper bound for the whole k-interval.  The displayed exact
# form is understood modulo n.(a,b,c)=0.
NORM6_ROOF_TABLE = {
    (0, 5, -1): (1, 2, (Fraction(2, 35), Fraction(3, 7), Fraction(0))),
    (1, -3, 2): (0, 2, (Fraction(0), Fraction(1, 7), Fraction(0))),
    (1, 3, -2): (1, 0, (Fraction(0), Fraction(0), Fraction(2, 21))),
    (2, -3, 1): (0, 2, (Fraction(0), Fraction(1, 7), Fraction(0))),
    (2, 3, -1): (1, 0, (Fraction(0), Fraction(0), Fraction(1, 21))),
    (3, -2, 1): (0, 2, (Fraction(0), Fraction(2, 21), Fraction(0))),
    (3, 1, -2): (1, 0, (Fraction(0), Fraction(0), Fraction(2, 21))),
    (3, 2, -1): (1, 0, (Fraction(0), Fraction(0), Fraction(1, 21))),
    (5, -1, 0): (0, 1, (Fraction(3, 7), Fraction(0), Fraction(2, 35))),
    (5, 0, -1): (2, 0, (Fraction(3, 7), Fraction(2, 35), Fraction(0))),
}
gate(set(NORM6_ROOF_TABLE) == FEASIBLE6, "norm-six roof table coverage")
for n in sorted(FEASIBLE6):
    upper_index, lower_index, exact_width = NORM6_ROOF_TABLE[n]
    q, _, lower, upper = roof_endpoint_forms(n)
    gate(dot(n, q) == 1, f"symbolic Bezout failure for {n}")
    selected_width = form_sub(upper[upper_index], lower[lower_index])
    gate(
        forms_equivalent_on_relation(selected_width, exact_width, n),
        f"roof-width identity failure for {n}",
    )
    gap = form_sub((Fraction(0), Fraction(0), Fraction(1, 7)), exact_width)
    gate(linear_form_nonnegative_on_feasible_cone(gap, n), f"roof-width bound failure for {n}")

# The norm-eight tail needs only its largest coordinate.  Each exact shape has
# M>=4, so one roof gives a k-interval of length
# 3(w_j+w_k)/(7M)<3c/14; one live residue modulo three gives m<c/14+1.
gate(min(max(abs(x) for x in n) for n in P8) == 4, "norm-eight maximum coefficient")
gate(3 * 53 >= 154, "norm-eight c>=53 tail arithmetic")
gate(all(x % 2 == 0 for x in (6, 8)), "parity bookkeeping")


# Complete finite base for the norm-eight tail, plus independent controls for
# the stronger shortest-one-zero conclusion.
HEIGHT = 49
VALUES = [x for x in range(1, HEIGHT + 1, 2) if x % 3]
RELATIONS_LE8 = all_primitive_relations_up_to(8)
eligible = 0
admit = {6: [], 8: []}
shortest = {6: [], 8: []}
line_checks = 0
for w in combinations(VALUES, 3):
    if gcd3(w) != 1:
        continue
    eligible += 1
    relations = [n for n in RELATIONS_LE8 if dot(n, w) == 0]
    for n in relations:
        gate(norm1(n) % 2 == 0, f"odd relation norm for odd w={w}")
    carriers = None
    admitting = {}
    for weight, patterns in ((6, P6), (8, P8)):
        admitting[weight] = [n for n in patterns if dot(n, w) == 0]
        if admitting[weight]:
            if carriers is None:
                carriers = direct_carriers(w)
            admit[weight].append((w, len(carriers)))
            for n in admitting[weight]:
                rebuilt, half_count, _ = line_reconstruction(w, n)
                gate(rebuilt == carriers, f"line/direct mismatch for w={w}, n={n}")
                gate(2 * half_count == len(carriers), f"line count mismatch for w={w}, n={n}")
                line_checks += 1
    if relations:
        mu = min(norm1(n) for n in relations)
        minimizers = [n for n in relations if norm1(n) == mu]
        if mu in (6, 8) and any(one_zero_mod_three(n) for n in minimizers):
            if carriers is None:
                carriers = direct_carriers(w)
            gate(11 * len(carriers) <= 2 * w[2], f"shortest one-zero count failure at {w}")
            shortest[mu].append((w, len(carriers)))

gate(len(VALUES) == 17 and eligible == 678, "finite-base universe")
gate(len(admit[6]) == 204 and len(admit[8]) == 190, "finite-base admitting counts")
gate(len(shortest[6]) == 202 and len(shortest[8]) == 163, "finite-base shortest counts")

for weight in (6, 8):
    failures = [(w, N) for w, N in admit[weight] if 11 * N > 2 * w[2]]
    gate(failures == [((1, 5, 7), 2)], f"unexpected admitting norm-{weight} hostile")
    gate(all(11 * N <= 2 * w[2] for w, N in shortest[weight]), f"shortest norm-{weight}")

EXCEPTION = (1, 5, 7)
exception_carriers = direct_carriers(EXCEPTION)
exception_relations = [n for n in RELATIONS_LE8 if dot(n, EXCEPTION) == 0]
exception_mu = min(norm1(n) for n in exception_relations)
exception_S = projection_sums(EXCEPTION, exception_carriers)
gate(exception_mu == 4, "exception must have a shorter unit relation")
gate(exception_carriers == {(-2, -1, 1), (2, 1, -1)}, "exception carrier set")
gate(exception_S == (Fraction(8, 245), Fraction(6, 49), Fraction(4, 35)), "exception projections")
gate(min(exception_S) == Fraction(8, 245) < Fraction(6, 77), "exception projection closure")

CONTROL6 = (1, 5, 61)
CONTROL8 = (17, 23, 25)
control6_carriers = direct_carriers(CONTROL6)
control8_carriers = direct_carriers(CONTROL8)
gate(len(control6_carriers) == 4 and 11 * 4 <= 2 * 61, "norm-six rank-two control")
gate(len(control8_carriers) == 4 and 11 * 4 <= 2 * 25, "norm-eight rank-two control")

max6 = max(shortest[6], key=lambda row: Fraction(row[1], row[0][2]))
max8 = max(shortest[8], key=lambda row: Fraction(row[1], row[0][2]))

print("LRC14 ONE-ZERO SHORTEST WIDTH <=8 AUTOMATIC COUNT REFEREE")
print("arithmetic=exact integers/Fraction; standalone_no_repo_import=true; optimizable_assertions=0")
print("scope=sorted primitive distinct positive odd ternary-unit triples w=(a,b,c)")
print("pattern_classification=mu parity even; norm2 repeats a speed; norm4 one-zero empty")
print("norm6_absolute_shapes=(0,1,5),(1,2,3); oriented=36; sorted_feasible=10")
print("norm6_sorted_patterns=" + ",".join(str(n) for n in sorted(FEASIBLE6)))
print("norm6_analytic=each live half lies in one k residue mod3 inside interval length <=c/7")
print("norm6_count=m<=ceil(c/21); automatic for c>=11; sole lower eligible row=(1,5,7)")
print("norm8_absolute_shapes=(0,1,7),(1,1,6),(1,3,4); oriented=48; max_coefficient>=4")
print("norm8_analytic=m<c/14+1<=c/11 for c>=53; complete remaining base c<=49")
print(f"finite_base=H{HEIGHT}; values={len(VALUES)}; eligible={eligible}; line_reconstructions={line_checks}")
print(f"admit_norm6={len(admit[6])}; admit_norm8={len(admit[8])}; unique_count_hostile=(1,5,7):N=2")
print(f"shortest_onezero_mu6={len(shortest[6])}; all_automatic=true; max_base_density={max6[1]}/{max6[0][2]}@{max6[0]}")
print(f"shortest_onezero_mu8={len(shortest[8])}; all_automatic=true; max_base_density={max8[1]}/{max8[0][2]}@{max8[0]}")
print("hostile_scope=(1,5,7) admits norm6/norm8 one-zero relations but mu1=4 unit; excluded by shortest hypothesis")
print("hostile_projection=S=(8/245,6/49,4/35); min=8/245<6/77")
print(f"control_norm6={CONTROL6}:N={len(control6_carriers)}; control_norm8={CONTROL8}:N={len(control8_carriers)}")
print(f"checks={CHECKS}")
print("verdict=PASS")
