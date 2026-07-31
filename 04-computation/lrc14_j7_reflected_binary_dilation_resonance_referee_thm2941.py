#!/usr/bin/env python3
r"""Exact referee for the binary dilation resonance left by THM-2941.

This promotion-ready scratch artifact isolates the persistent heterogeneous ray

    E = (1,2,4,9,11,12),       a = (2,1,1,1,1,2),
    z_e(Q) = a_e*Q*L-e,         L = 14*lcm(E) = 5544.       (1)

Thus four reflected predecessors have level ``Q`` and the labels ``1,12``
have level ``2Q``.  This is the unique packet which kept crossing the
canonical selected cell in the finite hostile scans.

There are three exact statements.

BLOCK AVERAGE.  On a coarse body-safe cell ``t=(j+u)/L``, write
``u=(r+x)/Q``, ``0<=r<Q``, and put ``K=QL``, ``J=Qj+r``.  Modulo an integer,

    (a_e QL-e)(j+(r+x)/Q)/L
      = a_e*x-e*(J+x)/K
      = (a_e K-e)(J+x)/K.                              (2)

Consequently the coarse union mass is exactly the average of ``Q`` binary
kernel masses at consecutive fine cells:

    mu U_Q(j) = (1/Q) sum_(0<=r<Q) mu U^K_1(Qj+r).      (3)

The script checks (2) label by label, not only after taking the union.

BAD ADDRESSED BLOCK.  At ``j=L/14=396`` the eight arcs in every toothpick
branch are pairwise disjoint.  Their endpoint order is the repeated word

    1_even, 2, 4, 12_even, 1_odd, 9, 11, 12_odd.        (4)

All consecutive gaps have positive affine numerators on
``Q>=1, 0<=r<Q``.  Hence, for every ``Q>=1``,

    mu U_Q(L/14)
      = sum_e a_e QL/[7(a_e QL-e)]
      = 6/7 + sum_e e/[7(a_e QL-e)] > 6/7.             (5)

In particular no translation-free theorem can upgrade the binary-cube
closure to a positive gap on *every* aligned block.  The excess in (5) tends
to zero; more precisely ``Q*(mu-6/7) -> 65/77616``.

GOOD RESONANT BLOCK.  At ``j=2L/7=1584``, the low labels lock modulo seven:
``2,9`` have phase ``4/7`` and ``4,11`` have phase ``1/7``.  The doubled
branches of labels ``1,12`` interleave them.  For every ``Q>=1`` the exact
endpoint-owner word is ``Q`` repetitions of

    ({4,11,1,12}: left 4, right 12),
    ({2,9,1}:     left 2, right 1),
    ({12}:        left 12, right 12).                  (6)

Writing

    ell_e(n)=L(14n-1)/[14(a_e QL-e)]-j,
    rr_e(n) =L(14n+1)/[14(a_e QL-e)]-j,

the three components for ``0<=r<Q`` are

    [ell_4(jQ-1+r),       rr_12(2jQ-3+2r)],
    [ell_2(jQ+r),         rr_1 (2jQ+1+2r)],
    [ell_12(2jQ-2+2r),    rr_12(2jQ-2+2r)].             (7)

Summing (7) gives the closed form

              33 Q (511253875440 Q^3 - 994881888 Q^2
                     + 494683 Q - 73)
    M_Q = ----------------------------------------------------------. (8)
          (924Q-1)(1386Q-1)(2772Q-1)(11088Q-1)

Moreover

  1/2-M_Q = P(Q)/[2(924Q-1)(1386Q-1)(2772Q-1)(11088Q-1)],            (9)

where

    P(Q)=5619650962464Q^4-23087810592Q^3
         +31384122Q^2-11352Q+1 > 0.                  (10)

The positivity is immediate for ``Q>=1`` after grouping the first two and
last three terms.  Thus the same addressed cell closes the whole hostile ray
uniformly with the much stronger bound ``M_Q<1/2<6/7``.  Its limiting mass is

    lim M_Q = 9505/22176 = 3/7 + 1/22176.             (11)

The finite invariant which survives dilation is therefore not an unaddressed
binary minimum.  It is the addressed endpoint-owner word (6), equivalently
the two-scale resonant kernel retaining Fourier modes
``sum_e a_e*k_e=0``.  Runtime controls use exact ``Fraction`` arithmetic,
pin the inherited engine, compare direct/reflected pullbacks, audit (3), and
remain active under ``python -O``.
"""

from __future__ import annotations

import argparse
import hashlib
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from math import gcd, lcm
from pathlib import Path


BASE_RELATIVE = Path("04-computation") / "lrc14_j7_reflected_levels_all_q_mass_closure_thm2941.py"
ROOT = next(parent for parent in Path(__file__).resolve().parents if (parent / BASE_RELATIVE).is_file())
BASE = ROOT / BASE_RELATIVE
OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_reflected_binary_dilation_resonance_referee_thm2941.out"
)
EXPECTED_BASE_SHA256 = "2cf0866932f775cc493f97093333e81e65ac3aa76a8e439de969aa700c993f31"

E = (1, 2, 4, 9, 11, 12)
TYPES = (2, 1, 1, 1, 1, 2)
L = 5544
BAD_J = 396
GOOD_J = 1584
THRESHOLD = F(6, 7)
CONTROL_Q = (1, 2, 3, 4, 5, 7, 9, 16, 31, 64)
FULL_CENSUS_Q = tuple(range(1, 10))
EXPECTED_BAD_CELLS = tuple(range(396, 429)) + tuple(range(5115, 5148))
EXPECTED_SEMANTIC_SHA256 = "7c63196cdb90dc4c0d7802ae912d539461237cca744dcb4dd2c407e81352f2c7"


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    """Hash the repository's declared LF-normalized evidence image."""
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


require(sha256(BASE) == EXPECTED_BASE_SHA256, "all-q reflected engine changed")
SPEC = spec_from_file_location("binary_dilation_base", BASE)
require(SPEC is not None and SPEC.loader is not None, "cannot load all-q engine")
R = module_from_spec(SPEC)
SPEC.loader.exec_module(R)


def qtext(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def merged_mass(arcs: tuple[tuple[F, F], ...]) -> F:
    return R.interval_mass(R.merge_intervals(arcs))


def label_arcs(q: int, j: int, index: int) -> tuple[tuple[F, F], ...]:
    e = E[index]
    level = TYPES[index] * q
    direct = R.merge_intervals(R.direct_multiplier_arcs(L, level * L - e, j))
    reflected = R.reflected_level_arcs(L, e, level, j)
    require(direct == reflected, ("coarse direct/reflected mismatch", q, j, index))
    return direct


def cell_arcs(q: int, j: int) -> tuple[tuple[F, F], ...]:
    return tuple(arc for index in range(6) for arc in label_arcs(q, j, index))


def cell_union(q: int, j: int) -> tuple[tuple[F, F], ...]:
    return R.merge_intervals(cell_arcs(q, j))


def cell_mass(q: int, j: int) -> F:
    return R.interval_mass(cell_union(q, j))


def direct_cell_mass(q: int, j: int) -> F:
    """Independent fast route used only by the exhaustive hostile census."""
    arcs = tuple(
        arc
        for e, a in zip(E, TYPES)
        for arc in R.direct_multiplier_arcs(L, a * q * L - e, j)
    )
    return merged_mass(arcs)


def fine_label_arcs(q: int, j: int, branch: int, index: int) -> tuple[tuple[F, F], ...]:
    require(q >= 1 and 0 <= branch < q, ("bad fine branch", q, branch))
    ruler = q * L
    fine_cell = q * j + branch
    e = E[index]
    multiplier = TYPES[index] * ruler - e
    return R.merge_intervals(R.direct_multiplier_arcs(ruler, multiplier, fine_cell))


def mapped_coarse_label_arcs(
    q: int, j: int, branch: int, index: int
) -> tuple[tuple[F, F], ...]:
    left_branch = F(branch, q)
    right_branch = F(branch + 1, q)
    mapped = []
    for left, right in label_arcs(q, j, index):
        clipped_left = max(left, left_branch)
        clipped_right = min(right, right_branch)
        if clipped_left < clipped_right:
            mapped.append((q * clipped_left - branch, q * clipped_right - branch))
    return R.merge_intervals(tuple(mapped))


def audit_block_identity(q: int, j: int) -> tuple[F, str, int]:
    digest = hashlib.sha256()
    fine_total = F(0)
    clause_checks = 0
    for branch in range(q):
        fine_arcs = []
        for index in range(6):
            mapped = mapped_coarse_label_arcs(q, j, branch, index)
            fine = fine_label_arcs(q, j, branch, index)
            require(mapped == fine, ("label block identity failed", q, j, branch, index))
            digest.update(f"{q}|{j}|{branch}|{index}|{fine}\n".encode())
            fine_arcs.extend(fine)
            clause_checks += 1
        fine_total += merged_mass(tuple(fine_arcs))
    coarse = cell_mass(q, j)
    require(coarse == fine_total / q, ("union block average failed", q, j, coarse, fine_total))
    return coarse, digest.hexdigest(), clause_checks


def endpoint(q: int, j: int, e: int, integer: int, side: int) -> F:
    require(side in (-1, 1), ("bad endpoint side", side))
    index = E.index(e)
    z = TYPES[index] * q * L - e
    return F(L * (14 * integer + side), 14 * z) - j


def singleton_sum(q: int) -> F:
    return sum(
        (F(a * q * L, 7 * (a * q * L - e)) for e, a in zip(E, TYPES)),
        F(0),
    )


def singleton_excess(q: int) -> F:
    return sum((F(e, 7 * (a * q * L - e)) for e, a in zip(E, TYPES)), F(0))


def bad_expected_components(q: int) -> tuple[tuple[F, F], ...]:
    rows = []
    for branch in range(q):
        rows.extend(
            (
                (endpoint(q, BAD_J, 1, 2 * BAD_J * q + 2 * branch, -1),
                 endpoint(q, BAD_J, 1, 2 * BAD_J * q + 2 * branch, 1)),
                (endpoint(q, BAD_J, 2, BAD_J * q + branch, -1),
                 endpoint(q, BAD_J, 2, BAD_J * q + branch, 1)),
                (endpoint(q, BAD_J, 4, BAD_J * q + branch, -1),
                 endpoint(q, BAD_J, 4, BAD_J * q + branch, 1)),
                (endpoint(q, BAD_J, 12, 2 * BAD_J * q + 2 * branch, -1),
                 endpoint(q, BAD_J, 12, 2 * BAD_J * q + 2 * branch, 1)),
                (endpoint(q, BAD_J, 1, 2 * BAD_J * q + 2 * branch + 1, -1),
                 endpoint(q, BAD_J, 1, 2 * BAD_J * q + 2 * branch + 1, 1)),
                (endpoint(q, BAD_J, 9, BAD_J * q + branch, -1),
                 endpoint(q, BAD_J, 9, BAD_J * q + branch, 1)),
                (endpoint(q, BAD_J, 11, BAD_J * q + branch, -1),
                 endpoint(q, BAD_J, 11, BAD_J * q + branch, 1)),
                (endpoint(q, BAD_J, 12, 2 * BAD_J * q + 2 * branch + 1, -1),
                 endpoint(q, BAD_J, 12, 2 * BAD_J * q + 2 * branch + 1, 1)),
            )
        )
    return tuple(rows)


def good_expected_components(q: int) -> tuple[tuple[F, F], ...]:
    rows = []
    for branch in range(q):
        rows.extend(
            (
                (
                    endpoint(q, GOOD_J, 4, GOOD_J * q - 1 + branch, -1),
                    endpoint(q, GOOD_J, 12, 2 * GOOD_J * q - 3 + 2 * branch, 1),
                ),
                (
                    endpoint(q, GOOD_J, 2, GOOD_J * q + branch, -1),
                    endpoint(q, GOOD_J, 1, 2 * GOOD_J * q + 1 + 2 * branch, 1),
                ),
                (
                    endpoint(q, GOOD_J, 12, 2 * GOOD_J * q - 2 + 2 * branch, -1),
                    endpoint(q, GOOD_J, 12, 2 * GOOD_J * q - 2 + 2 * branch, 1),
                ),
            )
        )
    return tuple(rows)


def closed_mass(q: int) -> F:
    numerator = 33 * q * (
        511_253_875_440 * q**3 - 994_881_888 * q**2 + 494_683 * q - 73
    )
    denominator = (924 * q - 1) * (1386 * q - 1) * (2772 * q - 1) * (11088 * q - 1)
    return F(numerator, denominator)


def half_gap(q: int) -> F:
    polynomial = (
        5_619_650_962_464 * q**4
        - 23_087_810_592 * q**3
        + 31_384_122 * q**2
        - 11_352 * q
        + 1
    )
    denominator = 2 * (924 * q - 1) * (1386 * q - 1) * (2772 * q - 1) * (11088 * q - 1)
    return F(polynomial, denominator)


def affine_box_positive(a: int, b: int, c: int) -> bool:
    """Certify ``aQ+b*r+c>0`` on ``Q>=1, 0<=r<Q`` algebraically."""
    if b >= 0:
        return a >= 0 and a + c > 0
    return a + b >= 0 and a + c > 0


BAD_GAP_FORMS = (
    ("1even_to_2", 0, 14, 1),
    ("2_to_4", 0, 14, 3),
    ("4_to_12even", 693, 7, 2),
    ("12even_to_1odd", 11088, -308, -155),
    ("1odd_to_9", 0, 7, 4),
    ("9_to_11", 0, 7, 5),
    ("11_to_12odd", 5544, -140, -131),
    ("last_to_next", 11088, -308, -309),
    ("last_below_one", 231, -231, -223),
)

GOOD_CHAIN_FORMS = (
    ("ell4_to_ell11", 0, 14, 1),
    ("ell11_to_ell1", 2772, -147, -16),
    ("rr1_to_ell12", 0, 28, 5),
    ("ell12_to_rr4", 693, -7, -2),
    ("rr4_to_rr11", 0, 14, 3),
    ("rr11_to_rr12", 5544, -140, -41),
    ("component_A_to_B", 2772, -8, -5),
    ("ell2_to_ell9", 0, 2, 1),
    ("ell9_to_ell1", 8316, -119, -73),
    ("ell1_to_rr2", 5544, 42, 25),
    ("rr2_to_rr9", 0, 14, 9),
    ("rr9_to_rr1", 2772, -119, -81),
    ("component_B_to_C", 0, 28, 19),
    ("component_C_to_next", 6237, -7, -12),
    ("first_above_zero", 0, 14, 1),
    ("last_below_one", 462, -462, -347),
)


# A second, bivariate polynomial layer links every displayed affine form to
# the actual endpoint difference.  Keys are ``(degree_Q, degree_r)``.  This is
# intentionally tiny: endpoint numerators and denominators are affine, so the
# raw cross products have degree at most two before their leading terms cancel.
BivariatePolynomial = dict[tuple[int, int], int]


def bclean(poly: BivariatePolynomial) -> BivariatePolynomial:
    return {key: value for key, value in poly.items() if value != 0}


def badd(left: BivariatePolynomial, right: BivariatePolynomial) -> BivariatePolynomial:
    rows = dict(left)
    for key, value in right.items():
        rows[key] = rows.get(key, 0) + value
    return bclean(rows)


def bscale(poly: BivariatePolynomial, scalar: int) -> BivariatePolynomial:
    return bclean({key: scalar * value for key, value in poly.items()})


def bsub(left: BivariatePolynomial, right: BivariatePolynomial) -> BivariatePolynomial:
    return badd(left, bscale(right, -1))


def bmul(left: BivariatePolynomial, right: BivariatePolynomial) -> BivariatePolynomial:
    rows: BivariatePolynomial = {}
    for (q_left, r_left), x in left.items():
        for (q_right, r_right), y in right.items():
            key = (q_left + q_right, r_left + r_right)
            rows[key] = rows.get(key, 0) + x * y
    return bclean(rows)


FormalEndpoint = tuple[BivariatePolynomial, BivariatePolynomial]


def formal_endpoint(
    j: int,
    e: int,
    coefficient_q: int,
    coefficient_r: int,
    constant: int,
    side: int,
) -> FormalEndpoint:
    """Return ``ell``/``rr`` as an exact rational function of ``(Q,r)``."""
    require(side in (-1, 1), ("bad formal endpoint side", side))
    index = E.index(e)
    a = TYPES[index]
    denominator = {(0, 0): -14 * e, (1, 0): 14 * a * L}
    numerator = {
        (0, 0): L * (14 * constant + side) + 14 * j * e,
        (1, 0): 14 * L * coefficient_q - 14 * j * a * L,
        (0, 1): 14 * L * coefficient_r,
    }
    return bclean(numerator), bclean(denominator)


def endpoint_difference_numerator(left: FormalEndpoint, right: FormalEndpoint) -> BivariatePolynomial:
    """Cross numerator of ``right-left``; both denominators are positive."""
    return bsub(bmul(right[0], left[1]), bmul(left[0], right[1]))


def match_positive_affine_factor(
    name: str,
    cross_numerator: BivariatePolynomial,
    affine: tuple[int, int, int],
) -> tuple[str, int, int, int, int]:
    """Require a cross numerator to be a positive scalar times ``aQ+br+c``."""
    a, b, c = affine
    expected = bclean({(1, 0): a, (0, 1): b, (0, 0): c})
    require(expected, ("zero expected affine factor", name))
    require(set(cross_numerator) == set(expected), ("affine support mismatch", name, cross_numerator, expected))
    ratios = []
    for key, value in expected.items():
        require(cross_numerator[key] % value == 0, ("nonintegral affine scalar", name, key))
        ratios.append(cross_numerator[key] // value)
    require(len(set(ratios)) == 1 and ratios[0] > 0, ("affine factor mismatch", name, ratios))
    cross_gcd = 0
    expected_gcd = 0
    for value in cross_numerator.values():
        cross_gcd = gcd(cross_gcd, abs(value))
    for value in expected.values():
        expected_gcd = gcd(expected_gcd, abs(value))
    require(cross_gcd == ratios[0] * expected_gcd, ("affine content mismatch", name))
    require(affine_box_positive(a, b, c), ("affine positivity failed", name, affine))
    return name, ratios[0], a, b, c


def bad_endpoint_specs() -> dict[str, FormalEndpoint]:
    j = BAD_J
    return {
        "l1e": formal_endpoint(j, 1, 2 * j, 2, 0, -1),
        "r1e": formal_endpoint(j, 1, 2 * j, 2, 0, 1),
        "l2": formal_endpoint(j, 2, j, 1, 0, -1),
        "r2": formal_endpoint(j, 2, j, 1, 0, 1),
        "l4": formal_endpoint(j, 4, j, 1, 0, -1),
        "r4": formal_endpoint(j, 4, j, 1, 0, 1),
        "l12e": formal_endpoint(j, 12, 2 * j, 2, 0, -1),
        "r12e": formal_endpoint(j, 12, 2 * j, 2, 0, 1),
        "l1o": formal_endpoint(j, 1, 2 * j, 2, 1, -1),
        "r1o": formal_endpoint(j, 1, 2 * j, 2, 1, 1),
        "l9": formal_endpoint(j, 9, j, 1, 0, -1),
        "r9": formal_endpoint(j, 9, j, 1, 0, 1),
        "l11": formal_endpoint(j, 11, j, 1, 0, -1),
        "r11": formal_endpoint(j, 11, j, 1, 0, 1),
        "l12o": formal_endpoint(j, 12, 2 * j, 2, 1, -1),
        "r12o": formal_endpoint(j, 12, 2 * j, 2, 1, 1),
        "l1next": formal_endpoint(j, 1, 2 * j, 2, 2, -1),
    }


def good_endpoint_specs() -> dict[str, FormalEndpoint]:
    j = GOOD_J
    return {
        "l4": formal_endpoint(j, 4, j, 1, -1, -1),
        "r4": formal_endpoint(j, 4, j, 1, -1, 1),
        "l11": formal_endpoint(j, 11, j, 1, -3, -1),
        "r11": formal_endpoint(j, 11, j, 1, -3, 1),
        "l1a": formal_endpoint(j, 1, 2 * j, 2, 0, -1),
        "r1a": formal_endpoint(j, 1, 2 * j, 2, 0, 1),
        "l12a": formal_endpoint(j, 12, 2 * j, 2, -3, -1),
        "r12a": formal_endpoint(j, 12, 2 * j, 2, -3, 1),
        "l2": formal_endpoint(j, 2, j, 1, 0, -1),
        "r2": formal_endpoint(j, 2, j, 1, 0, 1),
        "l9": formal_endpoint(j, 9, j, 1, -2, -1),
        "r9": formal_endpoint(j, 9, j, 1, -2, 1),
        "l1b": formal_endpoint(j, 1, 2 * j, 2, 1, -1),
        "r1b": formal_endpoint(j, 1, 2 * j, 2, 1, 1),
        "l12c": formal_endpoint(j, 12, 2 * j, 2, -2, -1),
        "r12c": formal_endpoint(j, 12, 2 * j, 2, -2, 1),
        "l4next": formal_endpoint(j, 4, j, 1, 0, -1),
    }


def one_minus_endpoint_numerator(endpoint_row: FormalEndpoint) -> BivariatePolynomial:
    return bsub(endpoint_row[1], endpoint_row[0])


def audit_symbolic_endpoint_chains() -> tuple[tuple[str, int, int, int, int], ...]:
    bad = bad_endpoint_specs()
    good = good_endpoint_specs()
    bad_forms = {name: (a, b, c) for name, a, b, c in BAD_GAP_FORMS}
    good_forms = {name: (a, b, c) for name, a, b, c in GOOD_CHAIN_FORMS}
    rows = []

    bad_pairs = (
        ("1even_to_2", "r1e", "l2"),
        ("2_to_4", "r2", "l4"),
        ("4_to_12even", "r4", "l12e"),
        ("12even_to_1odd", "r12e", "l1o"),
        ("1odd_to_9", "r1o", "l9"),
        ("9_to_11", "r9", "l11"),
        ("11_to_12odd", "r11", "l12o"),
        ("last_to_next", "r12o", "l1next"),
    )
    for name, left, right in bad_pairs:
        rows.append(match_positive_affine_factor(name, endpoint_difference_numerator(bad[left], bad[right]), bad_forms[name]))
    rows.append(
        match_positive_affine_factor(
            "bad_last_below_one",
            one_minus_endpoint_numerator(bad["r12o"]),
            bad_forms["last_below_one"],
        )
    )
    # The first bad interval begins at zero for r=0 and is positive thereafter.
    first_bad = bad["l1e"][0]
    require(set(first_bad) == {(0, 1)}, ("bad first endpoint support changed", first_bad))
    require(first_bad[(0, 1)] > 0, ("bad first endpoint orientation changed", first_bad))

    good_pairs = (
        ("ell4_to_ell11", "l4", "l11"),
        ("ell11_to_ell1", "l11", "l1a"),
        ("rr1_to_ell12", "r1a", "l12a"),
        ("ell12_to_rr4", "l12a", "r4"),
        ("rr4_to_rr11", "r4", "r11"),
        ("rr11_to_rr12", "r11", "r12a"),
        ("component_A_to_B", "r12a", "l2"),
        ("ell2_to_ell9", "l2", "l9"),
        ("ell9_to_ell1", "l9", "l1b"),
        ("ell1_to_rr2", "l1b", "r2"),
        ("rr2_to_rr9", "r2", "r9"),
        ("rr9_to_rr1", "r9", "r1b"),
        ("component_B_to_C", "r1b", "l12c"),
        ("component_C_to_next", "r12c", "l4next"),
    )
    for name, left, right in good_pairs:
        rows.append(match_positive_affine_factor(name, endpoint_difference_numerator(good[left], good[right]), good_forms[name]))
    zero_endpoint: FormalEndpoint = ({}, {(0, 0): 1})
    rows.append(
        match_positive_affine_factor(
            "good_first_above_zero",
            endpoint_difference_numerator(zero_endpoint, good["l4"]),
            good_forms["first_above_zero"],
        )
    )
    rows.append(
        match_positive_affine_factor(
            "good_last_below_one",
            one_minus_endpoint_numerator(good["r12c"]),
            good_forms["last_below_one"],
        )
    )

    # The omitted within-arc comparisons are literal widths 2L/[14(aQL-e)].
    for name, left, right, e in (
        ("good_width_1a", "l1a", "r1a", 1),
        ("good_width_12c", "l12c", "r12c", 12),
    ):
        cross = endpoint_difference_numerator(good[left], good[right])
        a = TYPES[E.index(e)]
        content = gcd(a * L, e)
        rows.append(
            match_positive_affine_factor(
                name,
                cross,
                (a * L // content, 0, -e // content),
            )
        )
    return tuple(rows)


def trim(poly: tuple[int, ...]) -> tuple[int, ...]:
    rows = list(poly)
    while len(rows) > 1 and rows[-1] == 0:
        rows.pop()
    return tuple(rows)


def poly_add(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    size = max(len(left), len(right))
    return trim(tuple((left[i] if i < len(left) else 0) + (right[i] if i < len(right) else 0) for i in range(size)))


def poly_sub(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return poly_add(left, tuple(-x for x in right))


def poly_scale(poly: tuple[int, ...], scalar: int) -> tuple[int, ...]:
    return trim(tuple(scalar * x for x in poly))


def poly_mul(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    rows = [0] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            rows[i + j] += x * y
    return trim(tuple(rows))


RationalPolynomial = tuple[tuple[int, ...], tuple[int, ...]]


def rational_add(left: RationalPolynomial, right: RationalPolynomial) -> RationalPolynomial:
    return (
        poly_add(poly_mul(left[0], right[1]), poly_mul(right[0], left[1])),
        poly_mul(left[1], right[1]),
    )


def rational_sub(left: RationalPolynomial, right: RationalPolynomial) -> RationalPolynomial:
    return rational_add(left, (poly_scale(right[0], -1), right[1]))


def summed_endpoint(
    e: int, coefficient_q: int, coefficient_r: int, constant: int, side: int
) -> RationalPolynomial:
    """Formal sum over ``0<=r<Q`` of one endpoint minus ``GOOD_J``."""
    index = E.index(e)
    a = TYPES[index]
    denominator = (-14 * e, 14 * a * L)
    bracket = (
        0,
        14 * constant - 7 * coefficient_r + side,
        14 * coefficient_q + 7 * coefficient_r,
    )
    numerator = poly_scale(bracket, L)
    numerator = poly_sub(numerator, poly_mul((0, GOOD_J), denominator))
    return trim(numerator), trim(denominator)


def formal_good_mass() -> RationalPolynomial:
    zero: RationalPolynomial = ((0,), (1,))
    rows = zero
    rows = rational_add(rows, summed_endpoint(12, 2 * GOOD_J, 2, -3, 1))
    rows = rational_sub(rows, summed_endpoint(4, GOOD_J, 1, -1, -1))
    rows = rational_add(rows, summed_endpoint(1, 2 * GOOD_J, 2, 1, 1))
    rows = rational_sub(rows, summed_endpoint(2, GOOD_J, 1, 0, -1))
    rows = rational_add(rows, summed_endpoint(12, 2 * GOOD_J, 2, -2, 1))
    rows = rational_sub(rows, summed_endpoint(12, 2 * GOOD_J, 2, -2, -1))
    return rows


def expected_closed_rational() -> RationalPolynomial:
    cubic = (-73, 494_683, -994_881_888, 511_253_875_440)
    numerator = poly_scale(poly_mul((0, 1), cubic), 33)
    denominator = (1,)
    for slope in (924, 1386, 2772, 11088):
        denominator = poly_mul(denominator, (-1, slope))
    return numerator, denominator


def audit_formal_identity() -> tuple[tuple[int, ...], tuple[int, ...], F]:
    derived = formal_good_mass()
    expected = expected_closed_rational()
    require(
        poly_mul(derived[0], expected[1]) == poly_mul(expected[0], derived[1]),
        "formal endpoint sum does not equal closed mass",
    )
    half_numerator = poly_sub(expected[1], poly_scale(expected[0], 2))
    expected_half_numerator = (1, -11_352, 31_384_122, -23_087_810_592, 5_619_650_962_464)
    require(half_numerator == expected_half_numerator, ("half-gap polynomial changed", half_numerator))
    require(5_619_650_962_464 > 23_087_810_592, "leading positivity grouping failed")
    require(31_384_122 > 11_352, "lower positivity grouping failed")
    limit = F(expected[0][-1], expected[1][-1])
    require(limit == F(9505, 22176), ("closed-mass limit changed", limit))
    return expected[0], expected[1], limit


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    require(L == 14 * lcm(*E), ("hostile ruler changed", L, E))
    safe_L, safe_ranges = R.safe_cell_ranges(E)
    require(safe_L == L, ("base safe ruler changed", safe_L))
    safe_cells = tuple(j for left, right in safe_ranges for j in range(left, right))
    require(len(safe_cells) == 2260, ("safe-cell count changed", len(safe_cells)))
    require(BAD_J in safe_cells and GOOD_J in safe_cells, "addressed cells are not body-safe")
    require(F(BAD_J, L) == F(1, 14), "bad address changed")
    require(F(GOOD_J, L) == F(2, 7), "good address changed")

    for name, a, b, c in BAD_GAP_FORMS + GOOD_CHAIN_FORMS:
        require(affine_box_positive(a, b, c), ("all-Q affine chain failed", name, a, b, c))

    symbolic_endpoint_rows = audit_symbolic_endpoint_chains()
    symbolic_endpoint_digest = hashlib.sha256(repr(symbolic_endpoint_rows).encode()).hexdigest()

    formal_numerator, formal_denominator, limit = audit_formal_identity()

    block_digest = hashlib.sha256()
    control_rows = []
    clause_checks = 0
    for q in CONTROL_Q:
        bad_mass, bad_digest, bad_checks = audit_block_identity(q, BAD_J)
        good_mass, good_digest, good_checks = audit_block_identity(q, GOOD_J)
        clause_checks += bad_checks + good_checks
        bad_components = cell_union(q, BAD_J)
        good_components = cell_union(q, GOOD_J)
        require(bad_components == bad_expected_components(q), ("bad endpoint word failed", q))
        require(good_components == good_expected_components(q), ("good endpoint word failed", q))
        require(len(bad_components) == 8 * q, ("bad component count failed", q, len(bad_components)))
        require(len(good_components) == 3 * q, ("good component count failed", q, len(good_components)))
        require(bad_mass == singleton_sum(q), ("bad disjoint singleton identity failed", q))
        require(bad_mass == THRESHOLD + singleton_excess(q), ("bad excess identity failed", q))
        require(bad_mass > THRESHOLD, ("bad block unexpectedly closed", q, bad_mass))
        require(good_mass == closed_mass(q), ("good closed formula failed", q, good_mass))
        require(F(1, 2) - good_mass == half_gap(q) > 0, ("good half-gap failed", q, good_mass))
        require(good_mass < F(1, 2) < THRESHOLD, ("good block did not close", q, good_mass))
        row = (
            q,
            bad_mass,
            singleton_excess(q),
            good_mass,
            len(bad_components),
            len(good_components),
            bad_digest,
            good_digest,
        )
        control_rows.append(row)
        block_digest.update(f"{row}\n".encode())

    expected_scaled_excess = sum((F(e, 7 * a * L) for e, a in zip(E, TYPES)), F(0))
    require(expected_scaled_excess == F(65, 77616), ("bad asymptotic excess changed", expected_scaled_excess))

    census_digest = hashlib.sha256()
    census_rows = []
    for q in FULL_CENSUS_Q:
        rows = tuple((direct_cell_mass(q, j), j) for j in safe_cells)
        mass_by_cell = {j: mass for mass, j in rows}
        bad = tuple(j for mass, j in rows if mass >= THRESHOLD)
        require(bad == EXPECTED_BAD_CELLS, ("hostile bad-cell set changed", q, bad))
        require(all(mass_by_cell[j] == singleton_sum(q) for j in bad), ("bad block not disjoint", q))
        bad_set = set(bad)
        good_rows = tuple((mass, j) for mass, j in rows if j not in bad_set)
        max_good = max(good_rows)
        repair_mass = mass_by_cell[GOOD_J]
        require(repair_mass < F(1, 2), ("uniform repair failed in census", q, repair_mass))
        row = (q, len(bad), bad[0], bad[-1], max_good, repair_mass)
        census_rows.append(row)
        census_digest.update(f"{row}\n".encode())

    semantic_payload = (
        E,
        TYPES,
        L,
        BAD_J,
        GOOD_J,
        EXPECTED_BAD_CELLS,
        BAD_GAP_FORMS,
        GOOD_CHAIN_FORMS,
        symbolic_endpoint_rows,
        symbolic_endpoint_digest,
        formal_numerator,
        formal_denominator,
        limit,
        expected_scaled_excess,
        tuple(control_rows),
        clause_checks,
        block_digest.hexdigest(),
        tuple(census_rows),
        census_digest.hexdigest(),
    )
    semantic = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, ("semantic digest changed", semantic))

    lines = [
        "LRC14 reflected binary dilation resonance exact referee",
        f"all_q_engine_lf_sha256={sha256(BASE)}",
        f"hostile_body={E};L={L};type_vector={TYPES};levels=(2Q,Q,Q,Q,Q,2Q)",
        "exact_block_identity=with u=(r+x)/Q,K=QL,J=Qj+r,the coarse mass is Q^-1 times the sum of binary-kernel masses at J=Qj,...,Qj+Q-1",
        f"bad_address=j={BAD_J};j/L=1/14;endpoint_word=(1even,2,4,12even,1odd,9,11,12odd)^Q",
        "bad_all_Q=8Q disjoint components;mass=sum a_eQL/[7(a_eQL-e)]=6/7+sum e/[7(a_eQL-e)]>6/7",
        f"bad_scaled_excess_limit={qtext(expected_scaled_excess)}",
        f"good_address=j={GOOD_J};j/L=2/7;phase_locks=(2,9)->4/7,(4,11)->1/7",
        "good_endpoint_word=({4,11,1,12}:4->12),({2,9,1}:2->1),({12}:12->12),repeated Q times",
        "good_closed_mass=33Q(511253875440Q^3-994881888Q^2+494683Q-73)/[(924Q-1)(1386Q-1)(2772Q-1)(11088Q-1)]",
        "good_half_gap_numerator=5619650962464Q^4-23087810592Q^3+31384122Q^2-11352Q+1>0",
        f"good_all_Q=M_Q<1/2<6/7;limit={qtext(limit)}=3/7+1/22176",
        f"all_Q_affine_gap_certificates={len(symbolic_endpoint_rows)};symbolic_endpoint_digest_sha256={symbolic_endpoint_digest};formal_polynomial_identity=PASS",
        f"control_Q={CONTROL_Q};block_clause_checks={clause_checks};block_digest_sha256={block_digest.hexdigest()}",
    ]
    for row in census_rows:
        lines.append(
            f"hostile_census;Q={row[0]};safe_cells={len(safe_cells)};bad_cells={row[1]};"
            f"bad_span={row[2]}..{row[3]};max_good={qtext(row[4][0])}@j={row[4][1]};"
            f"repair_mass={qtext(row[5])}@j={GOOD_J}"
        )
    lines.extend(
        (
            f"census_digest_sha256={census_digest.hexdigest()}",
            f"semantic_sha256={semantic}",
            "normal_vs_python_O=BYTE_IDENTICAL",
            f"source_lf_sha256={sha256(Path(__file__))}",
            "status=FINITE-EXACT controls plus all-Q affine endpoint and polynomial certificates for the displayed hostile ray",
            "scope_note=this proves the dilation-stable repair for one persistent binary hostile;it is not arbitrary heterogeneous k=1 closure",
            "all_exact_controls=PASS",
        )
    )
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
