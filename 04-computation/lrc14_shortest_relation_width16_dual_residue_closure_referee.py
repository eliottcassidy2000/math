#!/usr/bin/env python3
"""Independent referee for the complete shortest-relation l1=16 shell.

The script proves rational coefficient-cone chord bounds for both possible
mod-three residue types and rebuilds the finite base directly.  It imports no
repository implementation; all gates remain active under python -O.
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


def residue_type(n):
    zeros = sum(x % 3 == 0 for x in n)
    if zeros == 0:
        return "unit"
    if zeros == 1:
        return "onezero"
    return "forbidden"


def weight16_patterns(kind):
    result = set()
    for v in product(range(-16, 17), repeat=3):
        if norm1(v) != 16 or gcd3(v) != 1:
            continue
        v = canonical_sign(v)
        if residue_type(v) == kind:
            result.add(v)
    return result


BEZOUT_CACHE = {}


def bezout_vector(n):
    if n in BEZOUT_CACHE:
        return BEZOUT_CACHE[n]
    for radius in range(1, 17):
        for q in product(range(-radius, radius + 1), repeat=3):
            if dot(n, q) == 1:
                BEZOUT_CACHE[n] = q
                return q
    raise RuntimeError(f"Bezout search failed for {n}")


def feasible_packet(n):
    """Normalize c=1 and solve the exact open cone 0<a<b<1."""
    if n[0]:
        A = Fraction(-n[1], n[0])
        B = Fraction(-n[2], n[0])
        lower, upper = Fraction(0), Fraction(1)
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


def feval(form, w):
    return sum(form[i] * w[i] for i in range(3))


E = (
    (Fraction(1), Fraction(0), Fraction(0)),
    (Fraction(0), Fraction(1), Fraction(0)),
    (Fraction(0), Fraction(0), Fraction(1)),
)


def endpoint_forms(n, beta):
    """Return a relaxed collection of lower/upper k-roof linear forms."""
    q = bezout_vector(n)
    m = (
        fsub(fscale(q[2], E[1]), fscale(q[1], E[2])),
        fsub(fscale(q[0], E[2]), fscale(q[2], E[0])),
        fsub(fscale(q[1], E[0]), fscale(q[0], E[1])),
    )
    lower, upper = [], []
    for i in range(3):
        if not n[i]:
            # A constant-coordinate roof can delete a line but cannot enlarge
            # its chord, so omission is a valid analytic upper relaxation.
            continue
        B = fscale(Fraction(3, 14), fadd(E[(i + 1) % 3], E[(i + 2) % 3]))
        beta_m = fscale(beta, m[i])
        first = fscale(Fraction(1, n[i]), fsub(fscale(-1, B), beta_m))
        second = fscale(Fraction(1, n[i]), fsub(B, beta_m))
        if n[i] > 0:
            lower.append(first)
            upper.append(second)
        else:
            lower.append(second)
            upper.append(first)
    return lower, upper


def as_affine(form, base, direction):
    return feval(form, direction), feval(form, base)


def affine_value(function, t):
    return function[0] * t + function[1]


def chord_sum_at(groups, t):
    total = Fraction(0)
    for lower, upper in groups:
        left = max(affine_value(f, t) for f in lower)
        right = min(affine_value(f, t) for f in upper)
        total += max(Fraction(0), right - left)
    return total


def exact_chord_sum_supremum(n, betas):
    """Exact maximum on the closed normalized speed cone.

    Each chord is piecewise affine.  All changes of active roof or of empty
    status occur where two endpoint forms cross.
    """
    base, direction, interval_left, interval_right = feasible_packet(n)
    groups = []
    breakpoints = {interval_left, interval_right}
    for beta in betas:
        lower_forms, upper_forms = endpoint_forms(n, beta)
        lower = [as_affine(f, base, direction) for f in lower_forms]
        upper = [as_affine(f, base, direction) for f in upper_forms]
        groups.append((lower, upper))
        functions = lower + upper
        for i, first in enumerate(functions):
            for second in functions[i + 1 :]:
                if first[0] == second[0]:
                    continue
                t = Fraction(second[1] - first[1], first[0] - second[0])
                if interval_left <= t <= interval_right:
                    breakpoints.add(t)
    return max((chord_sum_at(groups, t), t) for t in breakpoints)


UNIT = weight16_patterns("unit")
ONEZERO = weight16_patterns("onezero")
gate(len(UNIT) == 144 and len(ONEZERO) == 96, "oriented residue-type counts")
gate(
    {tuple(sorted(abs(x) for x in n)) for n in UNIT}
    == {(1, 1, 14), (1, 2, 13), (1, 4, 11), (1, 5, 10), (1, 7, 8), (2, 7, 7), (4, 5, 7)},
    "unit absolute-shape classification",
)
gate(
    {tuple(sorted(abs(x) for x in n)) for n in ONEZERO}
    == {(0, 5, 11), (2, 3, 11), (2, 5, 9), (3, 5, 8), (5, 5, 6)},
    "one-zero absolute-shape classification",
)

FEASIBLE_UNIT = {n for n in UNIT if feasible_packet(n) is not None}
FEASIBLE_ONEZERO = {n for n in ONEZERO if feasible_packet(n) is not None}
gate(len(FEASIBLE_UNIT) == 51, "feasible unit signed patterns")
gate(len(FEASIBLE_ONEZERO) == 32, "feasible one-zero signed patterns")

EXPECTED_SUPREMA = {
    "unit": {
        (1, 1, 14): Fraction(9, 49),
        (1, 2, 13): Fraction(15, 91),
        (1, 4, 11): Fraction(27, 154),
        (1, 5, 10): Fraction(33, 175),
        (1, 7, 8): Fraction(45, 196),
        (2, 7, 7): Fraction(57, 343),
        (4, 5, 7): Fraction(69, 490),
    },
    "onezero": {
        (0, 5, 11): Fraction(1, 7),
        (2, 3, 11): Fraction(12, 77),
        (2, 5, 9): Fraction(11, 63),
        (3, 5, 8): Fraction(11, 56),
        (5, 5, 6): Fraction(1, 7),
    },
}

COMPUTED_SUPREMA = {"unit": {}, "onezero": {}}
for kind, patterns, betas in (
    ("unit", FEASIBLE_UNIT, (-3, 0, 3)),
    ("onezero", FEASIBLE_ONEZERO, (1, 2)),
):
    for n in sorted(patterns):
        value, location = exact_chord_sum_supremum(n, betas)
        shape = tuple(sorted(abs(x) for x in n))
        old = COMPUTED_SUPREMA[kind].get(shape, Fraction(0))
        COMPUTED_SUPREMA[kind][shape] = max(old, value)
        ceiling = Fraction(45, 196) if kind == "unit" else Fraction(11, 56)
        gate(value <= ceiling, f"{kind} chord ceiling at n={n}, t={location}")
gate(COMPUTED_SUPREMA == EXPECTED_SUPREMA, "shape chord-supremum tables")

# One live congruence class per coset and geometric line gives the half-count
# bounds H<15c/196+3 (unit) and H<11c/168+2 (one-zero).
gate(Fraction(15, 196) * 209 + 3 <= Fraction(209, 11), "unit c>=209 tail")
gate(Fraction(11, 168) * 79 + 2 <= Fraction(79, 11), "one-zero c>=79 tail")


def has_relation_below_16(w):
    """Exhaust every nonzero relation with l1 norm at most 15."""
    for x in range(-15, 16):
        remainder = 15 - abs(x)
        for y in range(-remainder, remainder + 1):
            numerator = -w[0] * x - w[1] * y
            if numerator % w[2]:
                continue
            z = numerator // w[2]
            weight = abs(x) + abs(y) + abs(z)
            if 0 < weight < 16:
                return True
    return False


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


def real_k_interval(w, n, base):
    lower, upper = None, None
    total = sum(w)
    for i in range(3):
        B = Fraction(3 * (total - w[i]), 14)
        if not n[i]:
            if abs(base[i]) >= B:
                return None
            continue
        first = Fraction(-B - base[i], n[i])
        second = Fraction(B - base[i], n[i])
        coordinate_lower, coordinate_upper = min(first, second), max(first, second)
        lower = coordinate_lower if lower is None else max(lower, coordinate_lower)
        upper = coordinate_upper if upper is None else min(upper, coordinate_upper)
    return None if lower is None or lower >= upper else (lower, upper)


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


def reconstruct_layers(w, n, kind):
    q = bezout_vector(n)
    m = cross(w, q)
    gate(dot(n, q) == 1 and cross(n, m) == w, f"basis identity at w={w}, n={n}")
    betas = (-3, 0, 3) if kind == "unit" else (-2, -1, 1, 2)
    layers = {}
    for beta in betas:
        base = tuple(beta * x for x in m)
        expected_residues = 2 if kind == "unit" else 1
        residues = [r for r in range(3) if all((base[i] + r * n[i]) % 3 for i in range(3))]
        gate(len(residues) == expected_residues, f"live residues at w={w}, n={n}, beta={beta}")
        layers[beta] = layer_points(w, n, base)
    return set().union(*layers.values()), layers, m


def ceil_fraction(x):
    return (x.numerator + x.denominator - 1) // x.denominator


def phase_blind_full_bound(w, n, kind):
    """Chord lengths plus ceilings, but without endpoint residue phases."""
    q = bezout_vector(n)
    m = cross(w, q)
    betas = (-3, 0, 3) if kind == "unit" else (1, 2)
    half_bound = 0
    lengths = []
    for beta in betas:
        base = tuple(beta * x for x in m)
        interval = real_k_interval(w, n, base)
        length = Fraction(0) if interval is None else interval[1] - interval[0]
        lengths.append((beta, length))
        half_bound += ceil_fraction(length / 3)
    return 2 * half_bound, tuple(lengths)


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


def canonical_direction(v):
    divisor = gcd3(v)
    u = tuple(x // divisor for x in v)
    return canonical_sign(u)


HEIGHT = 205
VALUES = [x for x in range(1, HEIGHT + 1, 2) if x % 3]
FEASIBLE16 = sorted(FEASIBLE_UNIT | FEASIBLE_ONEZERO)
eligible = 0
admitting = 0
rows = []
type_rows = {"unit": [], "onezero": []}
mixed_rows = 0
line_reconstructions = 0
phase_blind_hostiles = []
for w in combinations(VALUES, 3):
    if gcd3(w) != 1:
        continue
    eligible += 1
    minimizers = [n for n in FEASIBLE16 if dot(n, w) == 0]
    if not minimizers:
        continue
    admitting += 1
    if has_relation_below_16(w):
        continue
    kinds = {residue_type(n) for n in minimizers}
    if len(kinds) == 2:
        mixed_rows += 1
    carriers = direct_carriers(w)
    candidate_bounds = []
    for kind in sorted(kinds):
        type_rows[kind].append((w, len(carriers)))
        for n in minimizers:
            if residue_type(n) != kind:
                continue
            rebuilt, layers, _ = reconstruct_layers(w, n, kind)
            gate(rebuilt == carriers, f"layer/direct mismatch at w={w}, n={n}")
            line_reconstructions += 1
            bound, lengths = phase_blind_full_bound(w, n, kind)
            candidate_bounds.append((bound, kind, n, lengths, layers))
    best = min(candidate_bounds, key=lambda packet: (packet[0], packet[1], packet[2]))
    gate(11 * len(carriers) < 2 * w[2], f"nonstrict shortest width-16 row at {w}")
    rows.append((w, len(carriers), kinds, best))
    if 11 * best[0] > 2 * w[2]:
        phase_blind_hostiles.append((w, len(carriers), best[0], best[1], best[2], best[3]))

gate(len(VALUES) == 69 and eligible == 51873, "H205 finite universe")
gate(admitting == 6228, "H205 norm-16 admitting count")
gate(len(rows) == 5230, "H205 shortest width-16 row count")
gate(len(type_rows["unit"]) == 2210, "unit row count")
gate(len(type_rows["onezero"]) == 3106, "one-zero row count")
gate(
    mixed_rows == 86 and line_reconstructions == 5358,
    f"mixed/reconstruction counts: mixed={mixed_rows}, reconstructions={line_reconstructions}",
)

EXPECTED_PHASE_HOSTILES = [
    ((53, 71, 73), 10, 14),
    ((55, 71, 73), 8, 14),
    ((61, 103, 109), 14, 20),
]
gate(
    [(w, N, bound) for w, N, bound, _, _, _ in phase_blind_hostiles]
    == EXPECTED_PHASE_HOSTILES,
    "endpoint-phase hostile classification",
)

first_row = min(rows, key=lambda row: (row[0][2], row[0]))
gate(first_row[0] == (43, 53, 61) and first_row[1] == 6, "first width-16 row")
first_directions = {canonical_direction(C) for C in direct_carriers(first_row[0])}
gate(len(first_directions) == 2, "first width-16 row direction count")

unit_leader = max(type_rows["unit"], key=lambda row: (Fraction(row[1], row[0][2]), row[0]))
onezero_leader = max(type_rows["onezero"], key=lambda row: (Fraction(row[1], row[0][2]), row[0]))
gate(unit_leader == ((71, 95, 97), 16), "unit density leader")
gate(onezero_leader == ((109, 119, 125), 18), "one-zero density leader")


def phase_word(w, n):
    q = bezout_vector(n)
    m = cross(w, q)
    word = []
    for beta in (-3, 0, 3):
        base = tuple(beta * x for x in m)
        interval = real_k_interval(w, n, base)
        counts = []
        for residue in (1, 2):
            bounds = integer_k_range(w, n, base)
            if bounds is None:
                count = 0
            else:
                lower, upper = bounds
                first = lower + (residue - lower) % 3
                count = 0 if first > upper else 1 + (upper - first) // 3
            counts.append(count)
        word.append((beta, interval, tuple(counts)))
    return tuple(word)


EXPECTED_PHASE_WORD_COUNTS = {
    (53, 71, 73): ((1, 2), (2, 2), (2, 1)),
    (55, 71, 73): ((1, 1), (2, 2), (1, 1)),
    (61, 103, 109): ((2, 2), (3, 3), (2, 2)),
}
phase_words = {}
for w, N, _, kind, n, _ in phase_blind_hostiles:
    gate(kind == "unit", f"unexpected nonunit phase hostile at {w}")
    word = phase_word(w, n)
    counts = tuple(packet[2] for packet in word)
    gate(counts == EXPECTED_PHASE_WORD_COUNTS[w], f"phase word at {w}")
    gate(sum(sum(pair) for pair in counts) == N, f"phase total at {w}")
    phase_words[w] = counts

print("LRC14 SHORTEST-RELATION WIDTH-16 DUAL-RESIDUE CLOSURE REFEREE")
print("arithmetic=exact integers/Fraction; standalone_no_repo_import=true; optimizable_assertions=0")
print("scope=sorted primitive distinct positive odd ternary-unit triples with mu1=16")
print("residue_types=unit:three geometric layers/six coset-layers; onezero:four geometric/coset-layers")
print("unit_shapes=" + ",".join(str(shape) for shape in sorted(EXPECTED_SUPREMA["unit"])))
print("onezero_shapes=" + ",".join(str(shape) for shape in sorted(EXPECTED_SUPREMA["onezero"])))
print("signed_feasible=unit:51/144; onezero:32/96")
print("analytic_unit=ell_-3+ell_0+ell_3<=45c/196; automatic tail c>=209")
print("analytic_onezero=ell_1+ell_2<=11c/56; automatic tail c>=79")
print(f"finite_base=H{HEIGHT}; values={len(VALUES)}; eligible={eligible}; admitting={admitting}")
print(
    f"shortest_rows={len(rows)}; unit_rows={len(type_rows['unit'])}; onezero_rows={len(type_rows['onezero'])}; "
    f"mixed_rows={mixed_rows}; line_reconstructions={line_reconstructions}"
)
print("dense_rows=0; verdict_count=ALL_AUTOMATIC")
print(f"first_shell_row={first_row[0]}:N={first_row[1]}:directions={len(first_directions)}")
print(f"unit_density_leader={unit_leader[0]}:N={unit_leader[1]}:ratio={unit_leader[1]}/{unit_leader[0][2]}")
print(
    f"onezero_density_leader={onezero_leader[0]}:N={onezero_leader[1]}:"
    f"ratio={onezero_leader[1]}/{onezero_leader[0][2]}"
)
print("phase_blind_count_hostiles=" + ",".join(f"{w}:N={N}:bound={bound}" for w, N, bound in EXPECTED_PHASE_HOSTILES))
print("endpoint_phase_words=" + ",".join(f"{w}:{phase_words[w]}" for w in sorted(phase_words)))
print("cheapest_extra_invariant=(beta,live k-residue set,open-endpoint floor pairs)")
print(f"checks={CHECKS}")
print("verdict=PASS")
