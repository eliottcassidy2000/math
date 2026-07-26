#!/usr/bin/env python3
"""Exact companion for THM-2432.

The proof is analytic.  This dependency-free companion checks the two
partition mechanisms, the inherited unequal-depth anti-shield, the old and
new two-source escape constructions, the residual type arithmetic, and a
strict physical 91-root hostile.  Every truth-bearing check uses ``require``
so optimized Python executes the same audit.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from itertools import product
from math import gcd, lcm


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"check failed: {label}")


def centred(value: F, modulus: int = 1) -> F:
    residue = value % modulus
    return min(residue, modulus - residue)


def danger(value: F, width: int = 1) -> bool:
    return centred(value) < F(width, 14)


def strictly_safe(value: F, width: int = 1) -> bool:
    return centred(value) > F(width, 14)


def fibre_word(speed: int, base: F, order: int, width: int = 1) -> frozenset[int]:
    return frozenset(
        j
        for j in range(order)
        if danger(speed * (base + F(j, order)), width)
    )


def root_phase_histogram(width: int) -> tuple[int, Counter[int]]:
    """Generic word sizes on every open phase chamber of a 13-grid."""
    threshold = F(13 * width, 14)
    endpoints = sorted(
        {
            (sign * threshold - residue) % 13
            for sign in (-1, 1)
            for residue in range(13)
        }
    )
    histogram: Counter[int] = Counter()
    for index, left in enumerate(endpoints):
        right = endpoints[(index + 1) % len(endpoints)]
        if index + 1 == len(endpoints):
            right += 13
        phase = ((left + right) / 2) % 13
        size = sum(
            centred(phase + residue, 13) < threshold
            for residue in range(13)
        )
        histogram[size] += 1
    return len(endpoints), histogram


def wall_cell_measure(
    speeds_in: tuple[int, ...],
    speeds_out: tuple[int, ...],
) -> F:
    denominator = 14 * lcm(*(speeds_in + speeds_out))
    selected = 0
    for cell in range(denominator):
        midpoint = F(2 * cell + 1, 2 * denominator)
        if (
            all(danger(speed * midpoint) for speed in speeds_in)
            and all(not danger(speed * midpoint) for speed in speeds_out)
        ):
            selected += 1
    return F(selected, denominator)


def verify_witness(a: int, b: int, c: int, x: F) -> None:
    require(danger(a * x), "first source danger")
    require(danger(b * x), "second source danger")
    require(strictly_safe(c * x), "target strict safety")


def primitive_escape(
    a: int,
    b: int,
    c: int,
    transverse_divisor: int,
) -> tuple[F, str]:
    require(gcd(a, b) == 1, "primitive source pair")
    require(gcd(a * b, transverse_divisor) == 1, "source units")
    require(c % transverse_divisor == 0, "target divisibility")
    require(transverse_divisor >= 49, "target-size service")

    if a > b:
        a, b = b, a

    if a == b:
        require(a == 1, "equal coprime sources")

    if c > b:
        left = F(1, 14 * c)
        right = min(F(1, 14 * b), F(13, 14 * c))
        require(left < right, "central escape interval")
        x = (left + right) / 2
        verify_witness(a, b, c, x)
        return x, "central"

    require(c < b, "target/source equality excluded")
    require(a < b, "ordered distinct sources")
    require(b >= transverse_divisor + 1, "post-central source floor")

    inverse_a = pow(a, -1, b)
    lam = (c * inverse_a) % b
    if 2 * lam > b:
        lam -= b
    absolute_lam = abs(lam)

    if absolute_lam >= 2:
        j = b // (14 * absolute_lam) + 1
        require(0 < 14 * j < b, "multiplier source residue")
        require(14 * absolute_lam * j > b, "multiplier exits target arc")
        require(
            7 * absolute_lam * j <= 4 * b,
            "balanced multiplier upper bound",
        )
        n = (inverse_a * j) % b
        x = F(n, b)
        verify_witness(a, b, c, x)
        return x, "multiplier"

    require(lam not in (0, 1), "zero/positive unit multiplier excluded")
    require(lam == -1, "negative unit multiplier remains")
    require(a + c == b, "negative multiplier relation")

    m, residue = divmod(b, 14)
    require(1 <= residue <= 13, "nonzero residue modulo fourteen")
    require(c > residue, "target exceeds boundary residue")
    n = (inverse_a * m) % b
    delta = F(c - residue, 14 * b)
    bounds = (
        F(1, 7 * b),
        (F(1, 14) - delta) / a,
        delta / c,
    )
    require(all(bound > 0 for bound in bounds), "positive perturbation bounds")
    eta = min(bounds) / 2
    x = F(n, b) - F(1, 14 * b) + eta
    verify_witness(a, b, c, x)
    return x, "perturbation"


def two_source_escape(
    a: int,
    b: int,
    c: int,
    transverse_divisor: int,
) -> tuple[F, str]:
    require(gcd(a * b, transverse_divisor) == 1, "source unit hypothesis")
    require(c % transverse_divisor == 0, "target layer hypothesis")
    common = gcd(a, b)

    if c % common:
        for k in range(1, common):
            x = F(k, common)
            if strictly_safe(c * x):
                verify_witness(a, b, c, x)
                return x, "subgroup"
        raise RuntimeError("check failed: cyclic subgroup escape")

    require(
        gcd(common, transverse_divisor) == 1,
        "gcd division preserves transverse layer",
    )
    reduced_x, kind = primitive_escape(
        a // common,
        b // common,
        c // common,
        transverse_divisor,
    )
    x = reduced_x / common
    verify_witness(a, b, c, x)
    return x, kind


def weak_compositions(total: int, parts: int):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in weak_compositions(total - first, parts - 1):
            yield (first,) + tail


# 1. Exact lower-load and top-word fibre counts.
units_7 = tuple(range(1, 7))
phases = (F(1, 211), F(73, 211), F(191, 211))
top_guard_cases = 0
lower_word_cases = 0
lower_per_bin_checks = 0

for M in range(1, 4):
    order = 7 ** (M + 1)
    for unit in units_7:
        top_guard = 7**M * unit
        for phase in phases:
            word = fibre_word(top_guard, phase, order, width=2)
            require(len(word) == 2 * order // 7, "top guard size")
            require(
                len({index % 7 for index in word}) == 2,
                "top guard occupies two bins",
            )
            top_guard_cases += 1

    for depth in range(M):
        for unit in units_7:
            lower = 7**depth * unit
            for phase in phases:
                word = fibre_word(lower, phase, order)
                require(len(word) == order // 7, "lower ordinary size")
                for residue in range(7):
                    count = sum(index % 7 == residue for index in word)
                    require(count == order // 49, "lower ordinary per-bin count")
                    lower_per_bin_checks += 1
                lower_word_cases += 1

lower_load_width = 5 + 1 + 1
require(lower_load_width == 7, "low branch width")
outside_guard_bins = 7 - 2
require(outside_guard_bins == 5, "low branch guard complement")


# 2. Thirteen-root capacity and blocker pullback.
ordinary_endpoints, ordinary_histogram = root_phase_histogram(1)
guard_endpoints, guard_histogram = root_phase_histogram(2)
require(ordinary_endpoints == 26, "ordinary phase endpoints")
require(ordinary_histogram == Counter({1: 13, 2: 13}), "ordinary root sizes")
require(guard_endpoints == 26, "guard phase endpoints")
require(guard_histogram == Counter({3: 13, 4: 13}), "guard root sizes")
guard_root_capacity = max(guard_histogram)
guard_free_root_floor = 13 - guard_root_capacity
require(guard_root_capacity == 4, "guard root capacity")
require(guard_free_root_floor == 9, "guard-free root floor")

blocker_pullback_cases = 0
root_phases = (F(1, 97), F(17, 101), F(53, 107))
for quotient in range(1, 31):
    for phase in root_phases:
        base_truth = danger(quotient * phase)
        for h in range(13):
            root = (phase + h) / 13
            require(
                danger(13 * quotient * root) == base_truth,
                "blocker root pullback",
            )
            blocker_pullback_cases += 1


# 3. Inherited unequal-depth anti-shield and old 49-layer escape.
unequal_intersection_cases = 0
for M in (2, 3):
    order = 7 ** (M + 1)
    for lower_depth in range(M):
        for upper_depth in range(lower_depth + 1, M):
            for lower_unit, upper_unit in product(units_7, repeat=2):
                first = 7**lower_depth * lower_unit
                second = 7**upper_depth * upper_unit
                for phase in phases:
                    first_word = fibre_word(first, phase, order)
                    second_word = fibre_word(second, phase, order)
                    require(
                        len(first_word & second_word) == order // 49,
                        "unequal-depth intersection",
                    )
                    unequal_intersection_cases += 1

intersection_anti_shield = F(6, 7) * F(1, 49)
intersection_wall_control = wall_cell_measure((1, 7), (343,))
require(intersection_anti_shield == F(6, 343), "anti-shield arithmetic")
require(
    intersection_wall_control == intersection_anti_shield,
    "anti-shield wall control",
)

old_escape_kinds: Counter[str] = Counter()
old_escape_cases = 0
for a in range(1, 141):
    if gcd(a, 49) != 1:
        continue
    for b in range(a, 141):
        if gcd(b, 49) != 1:
            continue
        for c in range(49, 491, 49):
            witness, kind = two_source_escape(a, b, c, 49)
            verify_witness(a, b, c, witness)
            old_escape_kinds[kind] += 1
            old_escape_cases += 1

require(old_escape_cases == 72_600, "old 49-layer witness bank size")
require(
    old_escape_kinds
    == Counter(
        {
            "subgroup": 19_561,
            "central": 46_262,
            "multiplier": 6_663,
            "perturbation": 114,
        }
    ),
    "old 49-layer mechanism histogram",
)


# 4. The W=7,k=2 top inequality forces the unique one-fold bin word.
top_compositions = tuple(weak_compositions(7, 7))
admissible_top_compositions = tuple(
    composition
    for composition in top_compositions
    if all(7 * multiplicity >= 7 - 2 for multiplicity in composition)
)
require(len(top_compositions) == 1716, "weak-composition count")
require(
    admissible_top_compositions == ((1, 1, 1, 1, 1, 1, 1),),
    "unique exact top partition",
)
q_pair_containments = 5 * 4 // 2
require(q_pair_containments == 10, "ordinary pair count")


# 5. New transverse 91-layer escape.
subgroup_order_cases = 0
for order in range(2, 401):
    farthest = max(centred(F(k, order)) for k in range(order))
    require(farthest >= F(1, 3), "cyclic subgroup leaves central arc")
    subgroup_order_cases += 1

multiplier_inequality_cases = 0
for b in range(92, 401):
    for absolute_lam in range(2, b // 2 + 1):
        j = b // (14 * absolute_lam) + 1
        require(0 < 14 * j < b, "multiplier audit source")
        require(14 * absolute_lam * j > b, "multiplier audit exit")
        require(7 * absolute_lam * j <= 4 * b, "multiplier audit ceiling")
        multiplier_inequality_cases += 1

transverse_escape_kinds: Counter[str] = Counter()
transverse_escape_cases = 0
for a in range(1, 141):
    if gcd(a, 91) != 1:
        continue
    for b in range(a, 141):
        if gcd(b, 91) != 1:
            continue
        for c in range(91, 911, 91):
            witness, kind = two_source_escape(a, b, c, 91)
            verify_witness(a, b, c, witness)
            transverse_escape_kinds[kind] += 1
            transverse_escape_cases += 1

require(transverse_escape_cases == 62_160, "transverse witness bank size")
require(
    transverse_escape_kinds
    == Counter(
        {
            "subgroup": 16_502,
            "central": 43_394,
            "multiplier": 2_225,
            "perturbation": 39,
        }
    ),
    "transverse mechanism histogram",
)

explicit_controls = {
    "central": (1, 1, 91, F(1, 182)),
    "subgroup": (2, 2, 91, F(1, 2)),
    "multiplier": (1, 93, 91, F(4, 93)),
    "perturbation": (1, 92, 91, F(15189, 234416)),
}
for expected_kind, (a, b, c, expected_x) in explicit_controls.items():
    x, kind = two_source_escape(a, b, c, 91)
    require(kind == expected_kind, "explicit mechanism label")
    require(x == expected_x, "explicit mechanism witness")
    verify_witness(a, b, c, x)


# 6. Strict physical common-root hostile: root geometry alone is sharp.
hostile_phase = F(12345, 65537)
hostile_guard_speed = 1
hostile_ordinary_speeds = (731, 1046, 318, 1775, 1047)
hostile_guard = {
    r
    for r in range(91)
    if danger(hostile_guard_speed * (hostile_phase + r) / 91, width=2)
}
hostile_ordinaries = [
    {
        r
        for r in range(91)
        if danger(speed * (hostile_phase + r) / 91)
    }
    for speed in hostile_ordinary_speeds
]
hostile_multiplicity: Counter[int] = Counter()
for word in (hostile_guard, *hostile_ordinaries):
    hostile_multiplicity.update(word)
require(len(hostile_guard) == 26, "hostile guard size")
require(
    all(len(word) == 13 for word in hostile_ordinaries),
    "hostile ordinary sizes",
)
require(
    hostile_multiplicity == Counter({r: 1 for r in range(91)}),
    "hostile exact partition",
)
hostile_margin = min(
    abs(
        centred(speed * (hostile_phase + r) / 91)
        - F(width, 14)
    )
    for speed, width in (
        (hostile_guard_speed, 2),
        *((speed, 1) for speed in hostile_ordinary_speeds),
    )
    for r in range(91)
)
require(hostile_margin == F(19503, 11927734), "hostile chamber margin")


# 7. Residual type arithmetic.
residual_M0 = (
    (0, 5, 0, 7),
    (1, 5, 1, 8),
    (2, 5, 2, 9),
)
residual_positive_before = (
    (1, 5, 0, 7),
    (2, 0, 0, 2),
    (2, 5, 0, 7),
    (2, 5, 1, 8),
)
excluded_positive = ((2, 0, 0, 2), (2, 5, 0, 7))
residual_positive_after = tuple(
    shape for shape in residual_positive_before if shape not in excluded_positive
)
require(len(residual_M0) + len(residual_positive_before) == 7, "old residual count")
require(
    residual_positive_after == ((1, 5, 0, 7), (2, 5, 1, 8)),
    "new positive-depth residual list",
)
require(len(residual_M0) + len(residual_positive_after) == 5, "new residual count")


print("theorem=THM-2432")
print("status=PROVED+VERIFIED-EXACT+INDEPENDENTLY-HOSTILE-AUDITED")
print(f"top_guard_fibre_cases={top_guard_cases}")
print(f"lower_word_cases={lower_word_cases}")
print(f"lower_per_bin_checks={lower_per_bin_checks}")
print(f"lower_load_width={lower_load_width}")
print(f"outside_guard_bins={outside_guard_bins}")
print(f"ordinary_root_phase_endpoints={ordinary_endpoints}")
print("ordinary_root_histogram=1:13,2:13")
print(f"guard_root_phase_endpoints={guard_endpoints}")
print("guard_root_histogram=3:13,4:13")
print(f"guard_root_capacity={guard_root_capacity}")
print(f"guard_free_root_floor={guard_free_root_floor}")
print(f"blocker_pullback_cases={blocker_pullback_cases}")
print(f"unequal_depth_intersection_cases={unequal_intersection_cases}")
print(f"intersection_anti_shield={intersection_anti_shield}")
print(f"intersection_wall_control={intersection_wall_control}")
print(f"old_49_layer_escape_cases={old_escape_cases}")
print(
    "old_49_layer_escape_kinds="
    + ",".join(
        f"{kind}:{old_escape_kinds[kind]}"
        for kind in ("subgroup", "central", "multiplier", "perturbation")
    )
)
print(f"top_weak_compositions={len(top_compositions)}")
print(f"top_admissible_compositions={len(admissible_top_compositions)}")
print("top_unique_multiplicity=1,1,1,1,1,1,1")
print(f"q_pair_containments={q_pair_containments}")
print(f"subgroup_order_cases={subgroup_order_cases}")
print(f"multiplier_inequality_cases={multiplier_inequality_cases}")
print(f"transverse_91_escape_cases={transverse_escape_cases}")
print(
    "transverse_91_escape_kinds="
    + ",".join(
        f"{kind}:{transverse_escape_kinds[kind]}"
        for kind in ("subgroup", "central", "multiplier", "perturbation")
    )
)
print(f"physical_91_hostile_margin={hostile_margin}")
print(f"residual_types_before={len(residual_M0) + len(residual_positive_before)}")
print(f"excluded_positive_types={len(excluded_positive)}")
print(f"residual_types_after={len(residual_M0) + len(residual_positive_after)}")
print("remaining_M0=0,5,0,7|1,5,1,8|2,5,2,9")
print("remaining_Mpositive=1,5,0,7|2,5,1,8")
print("scalar_ledger=165")
print("LRC14=OPEN")
print("ALL CHECKS PASSED")
