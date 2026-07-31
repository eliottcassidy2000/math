#!/usr/bin/env python3
r"""Exact referee for the reflected ``k=1`` pair and conic-tail theorems.

This promotion-ready scratch artifact audits the following addendum to
THM-2941 without changing the canonical proof surface.

Let ``E`` be a six-element subset of ``{1,...,14}``, let

    L = 14*lcm(E),

and give the six reflected drifts independently varying levels

    z_e = q_e*L-e,       q_e >= 1.

The strongest statement audited here is a labelled-pair Bonferroni theorem.
Put

    epsilon(E,q)=sum_e e/[7(q_e*L-e)].

Choose any positive integer base ``b`` and two labels ``e,f`` with unequal
levels.  Write

    lambda_e=q_e-b-e/L,   lambda_f=q_f-b-f/L,
    S=|lambda_e|+|lambda_f|,   v=lambda_e-lambda_f,

and define the exact translation-free tent floor

    G(ell)=[floor(ell)/49+max(0,{ell}-5/7)^2/4]/ell.

If

    epsilon(E,q) + [4S+|v|/2]/b < G(|v|),               (P)

then every positive body-safe cell has six-clause union mass below ``6/7``.
Indeed the exact singleton sum is ``6/7+epsilon``.  Toothpick subdivision
compares the selected pair with its limiting tent: the two dynamic arcs cost
``4S/b``, and the 1-Lipschitz tent's left-grid error is ``|v|/(2b)``.  Hence
the exact pair overlap exceeds ``epsilon``, and Bonferroni closes the union.

This produces a strong uniform nearest-level corollary.  Let ``m=min q_e``
and let

    Delta=min{q_e-m:q_e>m}.

Choose a minimum-level label and a nearest higher-level label, and take
``b=m``.  Their offsets have opposite signs, so ``S=|v|<=Delta+1/14``.
Moreover

    epsilon<=1/[42(m-1/14)]<=1/(39m).

Since ``G(|v|)>=gamma=35/2976``, every nonconstant packet closes whenever

    m > (13392/35)Delta + 93992/3185.                  (P')

In particular ``m>=383Delta+30`` is sufficient, so a packet not closed by
this theorem must have ``m<=383Delta+29``.  The other four levels are
arbitrary: this is pair/Kakeya-strip transversality, not a bounded-box claim.

The earlier all-six conic estimate is a useful secondary sidecar.  More
intrinsically, choose any positive integer base ``m`` and put

    d_e=q_e-m,       Lambda_m=sum_e |d_e-e/L|.

If the level vector is nonconstant and there exists such a base with

    m > (2976/7)*Lambda_m,                              (C)

then every positive body-safe cell closes.  Since

    Lambda_m <= sum_e |q_e-m|+1/6,

a body-independent sufficient condition is

    m > (496/7)*(6*sum_e |q_e-m|+1).                   (C')

An integer median minimizes the L1 spread in (C'), but need not maximize the
ratio of its two sides.  Between consecutive level values that ratio is
fractional-linear, so the best coarse certificate is attained at one of the
six values ``m in {q_e}``; checking those six is exact.

For the exact condition (C), ``Lambda_m`` is affine between its real kinks
``q_e-e/L``.  The ratio ``m/Lambda_m`` is therefore fractional-linear and
monotone on each such segment; on the positive integer lattice a maximizer is
adjacent to a kink, hence belongs to ``{q_e-1,q_e}`` (discarding zero), and the
two exterior tails are monotone.  Thus the exact existential test also has a
twelve-candidate finite optimizer ledger.

For the range corollary, put
``m=min_e q_e``, ``d_e=q_e-m``, and ``D=max_e d_e``.  If ``D=0``, this is the
already-proved common-level diagonal.  If ``D>=1`` and

    m >= floor(496*(30*D+1)/7)+1,

then on *every* positive body-safe ``1/L`` cell the six reflected clauses have
union mass strictly below ``6/7``.  Hence the THM-2941 projected residual has
mass greater than ``1/7`` and cannot fit in the one remaining aligned comb.
Equivalently, a reflected packet not closed by this argument must satisfy

    min(q_e) <= floor(496*(30*D+1)/7).                   (BS)

The proof of (C) is an exact toothpick/Riemann argument.  On a safe cell
``t=(j+u)/L``, split ``u=(r+x)/m`` and put

    lambda_e = d_e-e/L,       s_r=r/m.

Modulo an integer, the exact phase is

    z_e*(j+u)/L = x + lambda_e*(j+s_r+x/m).              (1)

Its limiting arc is ``||x+lambda_e*(j+s)||<1/14``.  If two offsets differ,
their relative velocity

    v = (d_e-d_f)-(e-f)/L

has ``|v|>=1-13/168=155/168`` because ``L>=168``.  The intersection of two
limiting arcs is the one-periodic tent

    phi(y)=max(0,1/7-||y||).

For an interval of length ``ell`` the least possible integral of ``phi`` is

    floor(ell)/49 + max(0, frac(ell)-5/7)^2/4.           (2)

Thus its average is at least ``35/2976`` whenever ``ell>=155/168``.  Below
one, (2) is increasing and starts at that value.  For ``ell=k+w>=1``: if
``k>=2``, retain the complete periods to get at least ``2/147``; if ``k=1``
and ``w<=5/7``, get at least ``1/84``; and if ``k=1,w>=5/7``, subtracting
``gamma(1+w)`` from the numerator in (2), for ``gamma=35/2976``, gives a
quadratic whose exact minimum is ``47111/433972224>0``.  Consequently the
limiting six-arc union has average at most ``6/7-35/2976``.

For (1), the phase perturbation relative to the limiting arc at ``s_r`` is
bounded by ``|lambda_e|/m``.  Membership can change only within that distance
of either endpoint, a set of mass at most ``4|lambda_e|/m``.  Unioning the six
clauses gives the exact-to-sampled-limit error

    <= 4*sum_e |lambda_e|/m.                             (3)

Translation of one limiting arc by ``lambda_e*(s-t)`` changes it in mass by
at most ``2|lambda_e||s-t|``.  Integrating on every left Riemann subinterval
gives sampled-limit-to-integral error

    <= sum_e |lambda_e|/m.                               (4)

For arbitrary signed offsets, the exact error is

    5*Lambda_m/m.

This gives (C).  The triangle inequality and the exact body census
``sum(E)/L<=1/6`` give (C').  Under minimum normalization the ``d_e`` are
nonnegative and

    sum_e |lambda_e| <= sum_e d_e + sum_e e/L <= 5D+1/6,

because at least one of the six offsets is zero.

Thus the general finite union mass is at most

    6/7 - 35/2976 + 5*Lambda_m/m,

which is below ``6/7`` precisely under (C).  With ``m=min q``, the range
corollary follows from the displayed ``5D+1/6`` bound as soon as

    m > (2976/35)(25D+5/6) = 496(30D+1)/7.

The displayed integer threshold is therefore strict.  A constant level vector
has no differing pair and is routed to the common diagonal.  After minimum
normalization, ``D=0`` means exactly that constant case; there is no omitted
"nonconstant D=0" lane.

The referee below checks all body constants, the exact branch identity (1)
on every body in both the constant and nonconstant lanes (plus two deeper
patterns on the hostile bodies), the exact tent primitive and
translation lower bound, all permitted pair velocities with offset difference
through eight, and a complete moving-cell box ``q_e in {1,2,3}`` on two hostile
bodies.  Direct multiplier pullback and the independent affine branch route
must agree exactly.  Assertions are ``RuntimeError`` gates under ``python -O``.
"""

from __future__ import annotations

import argparse
import hashlib
from fractions import Fraction as Q
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, combinations_with_replacement, product
from math import lcm
from pathlib import Path


BASE_RELATIVE = Path("04-computation") / "lrc14_j7_reflected_levels_all_q_mass_closure_thm2941.py"
ROOT = next(parent for parent in Path(__file__).resolve().parents if (parent / BASE_RELATIVE).is_file())
BASE = ROOT / BASE_RELATIVE
OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_reflected_nearest_level_matching_tail_referee_thm2941.out"
)
EXPECTED_BASE_SHA256 = "2cf0866932f775cc493f97093333e81e65ac3aa76a8e439de969aa700c993f31"

H = Q(1, 14)
THRESHOLD = Q(6, 7)
TENT_FLOOR = Q(35, 2976)
MIN_COARSE_LEVEL = lambda D: (496 * (30 * D + 1)) // 7 + 1
PAIR_CLEAN_LEVEL = lambda delta: 383 * delta + 30
ALL_BODY_PATTERNS = (
    (1, (0, 0, 0, 0, 0, 0)),
    (2, (-1, 0, 0, 0, 0, 0)),
)
EXTRA_HOSTILE_PATTERNS = (
    (3, (0, 1, 0, 1, 0, 1)),
    (2, (0, 2, 1, 0, 2, 1)),
)
MOVING_CELL_BODIES = (
    (1, 2, 3, 4, 6, 12),
    (1, 2, 4, 6, 9, 12),
)

EXPECTED_BRANCH_DIGEST = "a5243d543ec62accc4d775262d2cb4314a8e090e44d7439732585b33d72ae92c"
EXPECTED_PAIR_DIGEST = "58783e1e6926d6e9109693e477b6eb7250c485ddcb97b6ebdd809334f50ffdda"
EXPECTED_POSITIVE_PAIR_DIGEST = "13625540cdf995f3c982829fab9f394e7d6357460483d25078f95ec8bd6cd35a"
EXPECTED_TENT_DIGEST = "80de5ccde3cf0946a98aaafcad5c2c27c32cb61c9d22e035eb2909894bf4ee59"
EXPECTED_MOVING_CELL_ROWS = (
    (
        (1, 2, 3, 4, 6, 12),
        (Q(522473986, 843926985), (1, 1, 2, 1, 2, 2), 18, 88),
    ),
    (
        (1, 2, 4, 6, 9, 12),
        (Q(35402863, 57635875), (1, 1, 1, 2, 1, 2), 161, 240),
    ),
)
EXPECTED_SEMANTIC_SHA256 = "d386be566a00055a955ab63b8cd30b5ae5d8e749782b2d83fbafaf5a422ce3e5"


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    """Hash the repository's declared LF-normalized evidence image."""
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


require(file_sha256(BASE) == EXPECTED_BASE_SHA256, "all-q reflected engine changed")
SPEC = spec_from_file_location("bounded_spread_base", BASE)
require(SPEC is not None and SPEC.loader is not None, "cannot load all-q engine")
R = module_from_spec(SPEC)
SPEC.loader.exec_module(R)


def qtext(value: Q) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def strict_integer_above(value: Q) -> int:
    return value.numerator // value.denominator + 1


def l1_coarse_level(spread: int) -> int:
    require(spread >= 0, ("negative L1 spread", spread))
    return (496 * (6 * spread + 1)) // 7 + 1


def coarse_certificate_ratio(levels: tuple[int, ...], base: int) -> Q:
    require(base >= 1 and all(q >= 1 for q in levels), ("bad coarse ratio input", levels, base))
    spread = sum(abs(q - base) for q in levels)
    return Q(base, 6 * spread + 1)


def exact_certificate_ratio(
    E: tuple[int, ...], L: int, levels: tuple[int, ...], base: int
) -> Q:
    require(base >= 1 and all(q >= 1 for q in levels), ("bad exact ratio input", levels, base))
    amplitude = sum(
        (abs(Q(q - base) - Q(e, L)) for e, q in zip(E, levels)),
        Q(0),
    )
    require(amplitude > 0, ("zero exact amplitude", E, L, levels, base))
    return Q(base) / amplitude


def clip(left: Q, right: Q) -> tuple[Q, Q] | None:
    left = max(Q(0), left)
    right = min(Q(1), right)
    return (left, right) if left < right else None


def affine_circle_arcs(coefficient: Q, phase: Q) -> tuple[tuple[Q, Q], ...]:
    """Pull back ``||coefficient*x+phase||<1/14`` on ``0<=x<=1``."""
    require(coefficient > 0, ("nonpositive affine coefficient", coefficient, phase))
    low_value = phase
    high_value = phase + coefficient
    lo = low_value.numerator // low_value.denominator - 1
    hi = high_value.numerator // high_value.denominator + 1
    arcs = []
    for integer in range(lo, hi + 1):
        piece = clip(
            (Q(integer) - H - phase) / coefficient,
            (Q(integer) + H - phase) / coefficient,
        )
        if piece is not None:
            arcs.append(piece)
    return R.merge_intervals(tuple(arcs))


def mapped_direct_branch_arcs(
    L: int, e: int, m: int, d: int, j: int, r: int
) -> tuple[tuple[Q, Q], ...]:
    """Restrict direct coarse ``u`` arcs to branch ``r`` and map to ``x``."""
    z = (m + d) * L - e
    coarse = R.direct_multiplier_arcs(L, z, j)
    left_branch = Q(r, m)
    right_branch = Q(r + 1, m)
    mapped = []
    for left, right in coarse:
        piece_left = max(left, left_branch)
        piece_right = min(right, right_branch)
        if piece_left < piece_right:
            mapped.append((m * piece_left - r, m * piece_right - r))
    return R.merge_intervals(tuple(mapped))


def affine_exact_branch_arcs(
    L: int, e: int, m: int, d: int, j: int, r: int
) -> tuple[tuple[Q, Q], ...]:
    lam = Q(d * L - e, L)
    s = Q(r, m)
    return affine_circle_arcs(Q(1) + lam / m, lam * (j + s))


def limiting_branch_arcs(
    L: int, e: int, d: int, j: int, s: Q
) -> tuple[tuple[Q, Q], ...]:
    lam = Q(d * L - e, L)
    return affine_circle_arcs(Q(1), lam * (j + s))


def interval_union_mass(rows: tuple[tuple[Q, Q], ...]) -> Q:
    return R.interval_mass(R.merge_intervals(rows))


def intersection_mass(
    left: tuple[tuple[Q, Q], ...], right: tuple[tuple[Q, Q], ...]
) -> Q:
    left_mass = R.interval_mass(R.merge_intervals(left))
    right_mass = R.interval_mass(R.merge_intervals(right))
    union = interval_union_mass(tuple(left) + tuple(right))
    return left_mass + right_mass - union


def singleton_excess(E: tuple[int, ...], L: int, levels: tuple[int, ...]) -> Q:
    return sum(
        (Q(e, 7 * (q * L - e)) for e, q in zip(E, levels)),
        Q(0),
    )


def direct_cell_mass(L: int, E: tuple[int, ...], levels: tuple[int, ...], j: int) -> Q:
    arcs = tuple(
        arc
        for e, q in zip(E, levels)
        for arc in R.direct_multiplier_arcs(L, q * L - e, j)
    )
    return interval_union_mass(arcs)


def exact_branch_profile(
    E: tuple[int, ...], L: int, j: int, m: int, offsets: tuple[int, ...]
) -> tuple[Q, Q, Q, int]:
    """Audit (1), its union form, and the pointwise perturbation budget."""
    exact_masses = []
    limit_masses = []
    clause_checks = 0
    lam_sum = sum((abs(Q(d * L - e, L)) for e, d in zip(E, offsets)), Q(0))
    for r in range(m):
        exact_arcs = []
        mapped_arcs = []
        limit_arcs = []
        for e, d in zip(E, offsets):
            route_a = mapped_direct_branch_arcs(L, e, m, d, j, r)
            route_b = affine_exact_branch_arcs(L, e, m, d, j, r)
            require(route_a == route_b, ("branch identity failed", E, L, j, m, offsets, r, e))
            exact_arcs.extend(route_b)
            mapped_arcs.extend(route_a)
            limit_arcs.extend(limiting_branch_arcs(L, e, d, j, Q(r, m)))
            clause_checks += 1
        exact_mass = interval_union_mass(tuple(exact_arcs))
        mapped_mass = interval_union_mass(tuple(mapped_arcs))
        limit_mass = interval_union_mass(tuple(limit_arcs))
        require(exact_mass == mapped_mass, ("branch union routes differ", E, m, offsets, r))
        require(
            abs(exact_mass - limit_mass) <= 4 * lam_sum / m,
            ("phase perturbation budget failed", E, m, offsets, r, exact_mass, limit_mass),
        )
        exact_masses.append(exact_mass)
        limit_masses.append(limit_mass)
    branch_average = sum(exact_masses, Q(0)) / m
    sampled_limit = sum(limit_masses, Q(0)) / m
    levels = tuple(m + d for d in offsets)
    direct_mass = direct_cell_mass(L, E, levels, j)
    require(branch_average == direct_mass, ("toothpick average failed", E, m, offsets))
    require(
        abs(branch_average - sampled_limit) <= 4 * lam_sum / m,
        ("averaged perturbation budget failed", E, m, offsets),
    )
    return direct_mass, sampled_limit, lam_sum, clause_checks


def pair_branch_profile(
    E: tuple[int, ...],
    L: int,
    j: int,
    m: int,
    offsets: tuple[int, ...],
    a: int,
    b: int,
) -> tuple[Q, Q, Q, Q, Q, Q]:
    """Exact pair overlap versus its sampled and integrated tent floors."""
    require(offsets[a] != offsets[b], ("pair levels are equal", offsets, a, b))
    lam_a = Q(offsets[a] * L - E[a], L)
    lam_b = Q(offsets[b] * L - E[b], L)
    velocity = lam_a - lam_b
    pair_size = abs(lam_a) + abs(lam_b)
    exact_rows = []
    limit_rows = []
    for r in range(m):
        exact_a = affine_exact_branch_arcs(L, E[a], m, offsets[a], j, r)
        exact_b = affine_exact_branch_arcs(L, E[b], m, offsets[b], j, r)
        limit_a = limiting_branch_arcs(L, E[a], offsets[a], j, Q(r, m))
        limit_b = limiting_branch_arcs(L, E[b], offsets[b], j, Q(r, m))
        exact_rows.append(intersection_mass(exact_a, exact_b))
        limit_rows.append(intersection_mass(limit_a, limit_b))
    exact_average = sum(exact_rows, Q(0)) / m
    sampled_average = sum(limit_rows, Q(0)) / m
    integrated_average = translated_tent_average(velocity, j)
    abstract_floor = minimum_tent_average(abs(velocity))
    require(
        abs(exact_average - sampled_average) <= 4 * pair_size / m,
        ("pair dynamic error failed", E, m, offsets, a, b),
    )
    require(
        abs(sampled_average - integrated_average) <= abs(velocity) / (2 * m),
        ("pair Riemann error failed", E, m, offsets, a, b),
    )
    require(integrated_average >= abstract_floor, ("pair tent floor failed", E, offsets, a, b))
    levels = tuple(m + d for d in offsets)
    direct_a = R.direct_multiplier_arcs(L, levels[a] * L - E[a], j)
    direct_b = R.direct_multiplier_arcs(L, levels[b] * L - E[b], j)
    direct_overlap = intersection_mass(direct_a, direct_b)
    require(exact_average == direct_overlap, ("direct pair route failed", E, m, offsets, a, b))
    certified_floor = abstract_floor - (4 * pair_size + abs(velocity) / 2) / m
    require(direct_overlap >= certified_floor, ("pair certified floor failed", E, m, offsets, a, b))
    return (
        direct_overlap,
        sampled_average,
        integrated_average,
        abstract_floor,
        pair_size,
        velocity,
    )


def tent(value: Q) -> Q:
    frac = value - value.numerator // value.denominator
    distance = min(frac, 1 - frac)
    return max(Q(0), Q(1, 7) - distance)


def minimum_tent_average(length: Q) -> Q:
    """Exact minimum over translations of an interval of this length."""
    require(length > 0, ("nonpositive tent interval", length))
    whole = length.numerator // length.denominator
    remainder = length - whole
    extra = Q(0) if remainder <= Q(5, 7) else (remainder - Q(5, 7)) ** 2 / 4
    return (Q(whole, 49) + extra) / length


def universal_tent_lower_bound(length: Q) -> Q:
    """Analytic branches proving the all-velocity ``35/2976`` floor."""
    coarse_velocity = Q(155, 168)
    require(length >= coarse_velocity, ("velocity below theorem range", length))
    if length < 1:
        # ((ell-5/7)^2)/(4 ell) has derivative
        # (1-(5/(7 ell))^2)/4 > 0 on this interval.
        return Q(35, 2976)
    whole = length.numerator // length.denominator
    remainder = length - whole
    if whole >= 2:
        # Keep complete periods and use whole/(whole+1)>=2/3.
        return Q(2, 147)
    if remainder <= Q(5, 7):
        return Q(1, 84)
    # In the final branch, with gamma=35/2976,
    #   1/49+(w-5/7)^2/4-gamma(1+w)
    # has minimum 47111/433972224 at w=5/7+2gamma.
    quadratic_margin = Q(47111, 433972224)
    require(quadratic_margin > 0, "tent quadratic margin changed")
    return TENT_FLOOR


def translated_tent_average(velocity: Q, j: int) -> Q:
    """Integrate ``tent(velocity*(j+s))`` exactly for ``0<=s<=1``."""
    require(velocity != 0, ("zero relative velocity", velocity))
    y0 = velocity * j
    y1 = velocity * (j + 1)
    low = min(y0, y1)
    high = max(y0, y1)
    lo = low.numerator // low.denominator - 2
    hi = high.numerator // high.denominator + 2
    events = {Q(0), Q(1)}
    for integer in range(lo, hi + 1):
        for offset in (Q(0), Q(1, 7), Q(6, 7)):
            s = (Q(integer) + offset) / velocity - j
            if 0 < s < 1:
                events.add(s)
    ordered = sorted(events)
    total = Q(0)
    for left, right in zip(ordered, ordered[1:]):
        total += (right - left) * (
            tent(velocity * (j + left)) + tent(velocity * (j + right))
        ) / 2
    return total


def first_different_pair(offsets: tuple[int, ...]) -> tuple[int, int]:
    for a, b in combinations(range(6), 2):
        if offsets[a] != offsets[b]:
            return a, b
    raise RuntimeError(("offset vector is constant", offsets))


def moving_cell_box(E: tuple[int, ...]) -> tuple[Q, tuple[int, ...], int, int]:
    """Complete exact q=1..3 box: maximize the best available safe-cell mass."""
    L, safe = R.safe_cell_ranges(E)
    cells = tuple(j for left, right in safe for j in range(left, right))
    require(cells, ("no body-safe cell", E))
    table = {
        (j, index, q): R.direct_multiplier_arcs(L, q * L - e, j)
        for j in cells
        for index, e in enumerate(E)
        for q in range(1, 4)
    }
    maximum = None
    for levels in product(range(1, 4), repeat=6):
        rows = tuple(
            (
                interval_union_mass(
                    tuple(
                        arc
                        for index, q in enumerate(levels)
                        for arc in table[j, index, q]
                    )
                ),
                j,
            )
            for j in cells
        )
        best = min(rows)
        row = (best[0], levels, best[1], len(cells))
        if maximum is None or row > maximum:
            maximum = row
    require(maximum is not None and maximum[0] < THRESHOLD, ("moving-cell survivor", E, maximum))
    return maximum


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    bodies = tuple(combinations(range(1, 15), 6))
    require(len(bodies) == 3003, "body universe changed")
    body_rows = tuple((14 * lcm(*E), E) for E in bodies)
    min_ruler = min(body_rows)
    sum_ratio_max = max((Q(sum(E), L), E, L) for L, E in body_rows)
    pair_ratio_max = max(
        (Q(abs(e - f), L), E, L, e, f)
        for L, E in body_rows
        for e, f in combinations(E, 2)
    )
    label_ratio_max = max((Q(e, L), E, L, e) for L, E in body_rows for e in E)
    require(min_ruler == (168, (1, 2, 3, 4, 6, 12)), ("minimum ruler changed", min_ruler))
    require(sum_ratio_max == (Q(1, 6), (1, 2, 3, 4, 6, 12), 168), ("sum ratio changed", sum_ratio_max))
    require(
        pair_ratio_max == (Q(11, 168), (1, 2, 3, 4, 6, 12), 168, 1, 12),
        ("pair ratio changed", pair_ratio_max),
    )
    require(label_ratio_max[0] == Q(1, 14), ("label ratio changed", label_ratio_max))

    coarse_velocity = Q(155, 168)
    coarse_subunit_tent = minimum_tent_average(coarse_velocity)
    require(coarse_subunit_tent == Q(35, 2976), ("coarse tent constant changed", coarse_subunit_tent))
    require(coarse_subunit_tent >= TENT_FLOOR, "subunit tent floor failed")
    require(MIN_COARSE_LEVEL(1) == 2197, "integer threshold changed")
    require(
        Q(MIN_COARSE_LEVEL(1)) * TENT_FLOOR > 25 + Q(5, 6),
        "strict threshold lost",
    )
    require(Q(9, 2) / TENT_FLOOR == Q(13392, 35), "pair slope constant changed")
    require(
        (Q(9, 28) + Q(1, 39)) / TENT_FLOOR == Q(93992, 3185),
        "pair intercept constant changed",
    )
    for delta in range(1, 1001):
        m = PAIR_CLEAN_LEVEL(delta)
        require(
            Q(9, 2) * (delta + Q(1, 14)) / m + Q(1, 39 * m) < TENT_FLOOR,
            ("clean nearest-level threshold failed", delta, m),
        )

    branch_digest = hashlib.sha256()
    pair_digest = hashlib.sha256()
    branch_rows = 0
    branch_clause_checks = 0
    translated_tent_checks = 0
    finite_bound_checks = 0
    pair_rows = 0
    pair_certificate_fires = 0
    for L, E in body_rows:
        safe = R.safe_cell_ranges(E)[1]
        j = safe[0][0]
        require(R.body_cell_is_safe(L, E, j), ("hostile cell unsafe", E, L, j))
        patterns = ALL_BODY_PATTERNS + (
            EXTRA_HOSTILE_PATTERNS if E in MOVING_CELL_BODIES else ()
        )
        for m, offsets in patterns:
            direct_mass, sampled_limit, lam_sum, clause_checks = exact_branch_profile(
                E, L, j, m, offsets
            )
            require(all(m + d >= 1 for d in offsets), ("nonpositive level", m, offsets))
            require(
                all(Q(m + d, m) - Q(e, m * L) > 0 for e, d in zip(E, offsets)),
                ("signed affine coefficient lost positivity", E, L, m, offsets),
            )
            l1_spread = sum(abs(d) for d in offsets)
            require(
                lam_sum <= l1_spread + Q(1, 6),
                ("signed lambda budget failed", E, offsets, lam_sum),
            )
            if len(set(offsets)) == 1:
                require(len(set(m + d for d in offsets)) == 1, ("constant route changed", offsets))
            else:
                a, b = first_different_pair(offsets)
                velocity = Q(offsets[a] - offsets[b]) - Q(E[a] - E[b], L)
                require(abs(velocity) >= coarse_velocity, ("velocity floor failed", E, offsets, velocity))
                translated = translated_tent_average(velocity, j)
                abstract_minimum = minimum_tent_average(abs(velocity))
                require(translated >= abstract_minimum >= TENT_FLOOR, ("tent floor failed", E, offsets))
                translated_tent_checks += 1
                pair_row = pair_branch_profile(E, L, j, m, offsets, a, b)
                levels = tuple(m + d for d in offsets)
                epsilon = singleton_excess(E, L, levels)
                direct_singleton_sum = sum(
                    (
                        R.interval_mass(
                            R.merge_intervals(
                                R.direct_multiplier_arcs(L, q * L - e, j)
                            )
                        )
                        for e, q in zip(E, levels)
                    ),
                    Q(0),
                )
                require(
                    direct_singleton_sum == THRESHOLD + epsilon,
                    ("singleton excess identity failed", E, levels, direct_singleton_sum, epsilon),
                )
                require(
                    direct_mass <= THRESHOLD + epsilon - pair_row[0],
                    ("Bonferroni direction failed", E, levels, direct_mass, pair_row[0]),
                )
                pair_certificate_left = epsilon + (
                    4 * pair_row[4] + abs(pair_row[5]) / 2
                ) / m
                if pair_certificate_left < pair_row[3]:
                    require(direct_mass < THRESHOLD, ("pair certificate fired falsely", E, levels))
                    pair_certificate_fires += 1
                pair_digest.update(
                    f"{E}|{L}|{j}|{m}|{offsets}|{a}|{b}|{pair_row}|{epsilon}|{pair_certificate_left}\n".encode()
                )
                pair_rows += 1
                exact_level = strict_integer_above(Q(2976, 7) * lam_sum)
                require(
                    Q(exact_level) * TENT_FLOOR > 5 * lam_sum,
                    ("exact conic strict bound failed", E, offsets, exact_level),
                )
                coarse_level = l1_coarse_level(l1_spread)
                require(
                    Q(coarse_level) * TENT_FLOOR
                    > 5 * (l1_spread + Q(1, 6)),
                    ("L1 conic strict bound failed", offsets, coarse_level),
                )
                if min(offsets) == 0:
                    D = max(offsets)
                    range_level = MIN_COARSE_LEVEL(D)
                    require(
                        sum(offsets) <= 5 * D,
                        ("minimum-normalized sum bound failed", offsets, D),
                    )
                    require(
                        Q(range_level) * TENT_FLOOR > 25 * D + Q(5, 6),
                        ("range strict bound failed", D, range_level),
                    )
                finite_bound_checks += 1
            row = (E, L, j, m, offsets, direct_mass, sampled_limit, lam_sum)
            branch_digest.update(f"{row}\n".encode())
            branch_rows += 1
            branch_clause_checks += clause_checks

    velocity_rows = set()
    tent_digest = hashlib.sha256()
    velocity_minimum = None
    for L, E in body_rows:
        for e, f in combinations(E, 2):
            for offset_difference in range(-8, 9):
                if offset_difference == 0:
                    continue
                key = (L, e, f, offset_difference)
                if key in velocity_rows:
                    continue
                velocity_rows.add(key)
                velocity = Q(offset_difference) - Q(e - f, L)
                minimum = minimum_tent_average(abs(velocity))
                analytic_floor = universal_tent_lower_bound(abs(velocity))
                require(abs(velocity) >= coarse_velocity, ("velocity census floor failed", key, velocity))
                require(
                    minimum >= analytic_floor >= TENT_FLOOR,
                    ("velocity census tent failed", key, minimum, analytic_floor),
                )
                row = (minimum, abs(velocity), offset_difference, L, e, f)
                if velocity_minimum is None or row < velocity_minimum:
                    velocity_minimum = row
                tent_digest.update(f"{row}\n".encode())
    require(
        velocity_minimum
        == (Q(1369, 105504), Q(157, 168), -1, 168, 1, 12),
        ("hostile velocity minimum changed", velocity_minimum),
    )

    breakpoint_digest = hashlib.sha256()
    breakpoint_rows = 0
    for levels in combinations_with_replacement(range(1, 9), 6):
        if len(set(levels)) == 1:
            continue
        all_integer_rows = tuple(
            (coarse_certificate_ratio(levels, base), base)
            for base in range(1, max(levels) + 9)
        )
        level_rows = tuple(
            (coarse_certificate_ratio(levels, base), base)
            for base in sorted(set(levels))
        )
        require(max(all_integer_rows) == max(level_rows), ("breakpoint certificate failed", levels))
        breakpoint_digest.update(f"{levels}|{max(level_rows)}\n".encode())
        breakpoint_rows += 1
    nonmedian_example = (1, 1, 1, 2, 1000, 1000)
    nonmedian_rows = tuple(
        sorted(
            (
                coarse_certificate_ratio(nonmedian_example, base),
                base,
                sum(abs(q - base) for q in nonmedian_example),
            )
            for base in sorted(set(nonmedian_example))
        )
    )
    require(nonmedian_rows[-1][1] == 1000, ("nonmedian control changed", nonmedian_rows))
    positive_example = (10000, 10000, 10001, 10001, 10002, 10002)
    positive_rows = tuple(
        sorted(
            (coarse_certificate_ratio(positive_example, base), base)
            for base in sorted(set(positive_example))
        )
    )
    require(positive_rows[-1] == (Q(10001, 25), 10001), ("positive L1 control changed", positive_rows))
    require(positive_rows[-1][0] > Q(496, 7), ("positive L1 cone did not fire", positive_rows))

    exact_breakpoint_digest = hashlib.sha256()
    exact_breakpoint_rows = 0
    for E in MOVING_CELL_BODIES:
        L = 14 * lcm(*E)
        for levels in product(range(1, 5), repeat=6):
            if len(set(levels)) == 1:
                continue
            all_integer_rows = tuple(
                (exact_certificate_ratio(E, L, levels, base), base)
                for base in range(1, max(levels) + 9)
            )
            candidate_bases = tuple(
                sorted(
                    {
                        candidate
                        for q in levels
                        for candidate in (q - 1, q)
                        if candidate >= 1
                    }
                )
            )
            candidate_rows = tuple(
                (exact_certificate_ratio(E, L, levels, base), base)
                for base in candidate_bases
            )
            require(
                max(all_integer_rows) == max(candidate_rows),
                ("exact breakpoint certificate failed", E, levels, candidate_bases),
            )
            exact_breakpoint_digest.update(f"{E}|{levels}|{max(candidate_rows)}\n".encode())
            exact_breakpoint_rows += 1

    positive_pair_digest = hashlib.sha256()
    positive_pair_rows = 0
    positive_levels = (413, 414, 413, 413, 413, 413)
    for L, E in body_rows:
        base = 413
        a, b = 0, 1
        lam_a = Q(positive_levels[a] - base) - Q(E[a], L)
        lam_b = Q(positive_levels[b] - base) - Q(E[b], L)
        pair_size = abs(lam_a) + abs(lam_b)
        velocity = lam_a - lam_b
        epsilon = singleton_excess(E, L, positive_levels)
        left = epsilon + (4 * pair_size + abs(velocity) / 2) / base
        floor = minimum_tent_average(abs(velocity))
        require(left < floor, ("positive pair certificate did not fire", E, left, floor))
        positive_pair_digest.update(f"{E}|{L}|{left}|{floor}|{pair_size}|{velocity}\n".encode())
        positive_pair_rows += 1
    direct_positive_E = (1, 2, 3, 4, 6, 12)
    direct_positive_L, direct_positive_safe = R.safe_cell_ranges(direct_positive_E)
    direct_positive_j = direct_positive_safe[0][0]
    direct_positive_mass = direct_cell_mass(
        direct_positive_L, direct_positive_E, positive_levels, direct_positive_j
    )
    require(direct_positive_mass < THRESHOLD, ("positive pair direct control failed", direct_positive_mass))
    require(
        positive_pair_digest.hexdigest() == EXPECTED_POSITIVE_PAIR_DIGEST,
        "positive pair digest changed",
    )

    moving_rows = tuple((E, moving_cell_box(E)) for E in MOVING_CELL_BODIES)
    if EXPECTED_MOVING_CELL_ROWS is not None:
        require(moving_rows == EXPECTED_MOVING_CELL_ROWS, ("moving-cell rows changed", moving_rows))
    if EXPECTED_BRANCH_DIGEST is not None:
        require(branch_digest.hexdigest() == EXPECTED_BRANCH_DIGEST, "branch digest changed")
    if EXPECTED_PAIR_DIGEST is not None:
        require(pair_digest.hexdigest() == EXPECTED_PAIR_DIGEST, "pair digest changed")
    if EXPECTED_TENT_DIGEST is not None:
        require(tent_digest.hexdigest() == EXPECTED_TENT_DIGEST, "tent digest changed")

    semantic_payload = (
        min_ruler,
        sum_ratio_max,
        pair_ratio_max,
        label_ratio_max,
        coarse_velocity,
        coarse_subunit_tent,
        TENT_FLOOR,
        ALL_BODY_PATTERNS,
        EXTRA_HOSTILE_PATTERNS,
        branch_rows,
        branch_clause_checks,
        pair_rows,
        pair_certificate_fires,
        translated_tent_checks,
        finite_bound_checks,
        len(velocity_rows),
        velocity_minimum,
        breakpoint_rows,
        breakpoint_digest.hexdigest(),
        nonmedian_rows,
        positive_rows,
        exact_breakpoint_rows,
        exact_breakpoint_digest.hexdigest(),
        positive_pair_rows,
        positive_pair_digest.hexdigest(),
        direct_positive_mass,
        moving_rows,
        branch_digest.hexdigest(),
        pair_digest.hexdigest(),
        tent_digest.hexdigest(),
    )
    semantic_hash = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    lines = [
        "LRC14 reflected k=1 nearest-level pair and conic tail exact referee",
        f"all_q_engine_lf_sha256={file_sha256(BASE)}",
        "scope=all E subset {1,...,14}, |E|=6; z_e=q_e*L-e; q_e positive integers",
        "primary_pair=epsilon=sum e/[7(q_eL-e)];lambda_e=q_e-b-e/L;S=|lambda_e|+|lambda_f|;v=lambda_e-lambda_f",
        "primary_pair_certificate=epsilon+[4S+|v|/2]/b<G(|v|),where G(ell)=[floor(ell)/49+max(0,{ell}-5/7)^2/4]/ell",
        "nearest_level_corollary=m=min q;Delta=least positive q_e-m;m>(13392/35)Delta+93992/3185;clean m>=383Delta+30",
        "secondary_cone=nonconstant q closes if there exists integer base m>=1 with Lambda=sum|q_e-m-e/L| and m>(2976/7)Lambda",
        "body_independent_L1_cone=m>(496/7)(6*sum|q_e-m|+1);best coarse certificate occurs at one of the six level values",
        "median_note=an integer median minimizes L1 spread but need not maximize the certificate ratio",
        "range_corollary=m=min(q_e);d_e=q_e-m;D=max(d_e)=max(q_e)-min(q_e)",
        "D=0=common-level diagonal already closed for every m>=1",
        "exact_base_optimizer=it suffices to test positive integers in {q_e-1,q_e};coarse_L1_optimizer needs only {q_e}",
        "D>=1 conclusion=every body-safe cell has six-clause union mass<6/7 when m>=floor(496(30D+1)/7)+1",
        "contrapositive_wedge=min(q_e)<=floor(496(30D+1)/7)",
        "geometry=conic tail reduction at bounded projective spread;not an all-k1 closure",
        "exact_subdivision=z_e(j+(r+x)/m)/L=x+(d_e-e/L)(j+r/m+x/m) mod1",
        "limit_pair_overlap=phi(v(j+s));phi(y)=max(0,1/7-||y||);v=(d_e-d_f)-(e-f)/L",
        f"minimum_ruler={min_ruler};max_sumE_over_L={sum_ratio_max};max_pair_ratio={pair_ratio_max}",
        f"max_label_over_L={label_ratio_max};nearest_level_exact_bound=m>(13392/35)Delta+93992/3185;clean=m>=383Delta+30",
        f"coarse_velocity_floor={qtext(coarse_velocity)};subunit_tent_average={qtext(coarse_subunit_tent)}",
        f"uniform_tent_floor={qtext(TENT_FLOOR)};limiting_union_gap={qtext(TENT_FLOOR)}",
        "all_velocity_tent_gate=ell<1 monotone;ell=k+w,k>=2 gives 2/147;k=1,w<=5/7 gives 1/84;k=1,w>=5/7 has quadratic margin 47111/433972224",
        "exact_to_sampled_error<=4*sum|d_e-e/L|/m;left_Riemann_error<=sum|d_e-e/L|/m",
        "min_normalization=sum d_e<=5D;sum|d_e-e/L|<=5D+1/6;total_error<35/2976 at integer m>=floor(496(30D+1)/7)+1",
        f"branch_rows={branch_rows};branch_clause_checks={branch_clause_checks};translated_tent_checks={translated_tent_checks}",
        f"pair_rows={pair_rows};pair_certificate_fires={pair_certificate_fires};pair_digest_sha256={pair_digest.hexdigest()}",
        f"velocity_hostile_rows={len(velocity_rows)};velocity_hostile_minimum={velocity_minimum}",
        f"coarse_breakpoint_rows={breakpoint_rows};breakpoint_digest_sha256={breakpoint_digest.hexdigest()}",
        f"median_not_optimal_control={nonmedian_example};rows={nonmedian_rows}",
        f"positive_L1_cone_control={positive_example};rows={positive_rows};conclusion=CLOSES",
        f"exact_breakpoint_q1_to4_rows={exact_breakpoint_rows};exact_breakpoint_digest_sha256={exact_breakpoint_digest.hexdigest()}",
        f"positive_nearest_pair_rows={positive_pair_rows};levels={positive_levels};positive_pair_digest_sha256={positive_pair_digest.hexdigest()};direct_hostile_mass={qtext(direct_positive_mass)}",
    ]
    for E, row in moving_rows:
        lines.append(
            f"moving_cell_q1_to3;E={E};max_over_vectors_min_over_safe_cells="
            f"{qtext(row[0])};levels={row[1]};best_j={row[2]};safe_cells={row[3]};gap={qtext(THRESHOLD-row[0])}"
        )
    lines.extend(
        (
            f"branch_digest_sha256={branch_digest.hexdigest()}",
            f"tent_digest_sha256={tent_digest.hexdigest()}",
            f"semantic_sha256={semantic_hash}",
            "normal_vs_python_O=BYTE_IDENTICAL",
            f"source_lf_sha256={file_sha256(Path(__file__))}",
            "status=FINITE-EXACT referee for the displayed analytic proof;not arbitrary reflected k=1 closure",
            "all_exact_controls=PASS",
        )
    )
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
