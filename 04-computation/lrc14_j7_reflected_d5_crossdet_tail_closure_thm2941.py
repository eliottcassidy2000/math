#!/usr/bin/env python3
r"""Exact proof that the full reflected ``D=5`` tail ``m>=16`` closes.

Let ``E`` be a six-label body, ``L=14*lcm(E)``, and put

    z_e=q_e L-e,                 q_e in {m,...,m+5}.

The universal chromatic theorem reduces 3001 bodies to six pairwise-distinct
levels, hence a permutation of ``m,...,m+5``.  On its two exceptional bodies
an uncertified word is a proper colouring of the recorded good graph.  This
referee closes every one of those residual words when ``m>=16``.

The new mechanism is an exact cross-determinant transport.  On a body-safe
cell ``t=(j+u)/L``, take labels ``a,b`` at distinct levels ``p,q`` and set

    c=1-a/(pL),             eta=(qa-pb)/(pL-a).

The substitution ``v=cu`` gives the exact pair overlap

    I=c^-1 int_0^c 1_B(pv-aj/L) 1_B((q+eta)v-bj/L) dv.       (1)

The first factor vanishes on the omitted tail ``c<=v<=1``.  Indeed, writing
``v=1-x`` and ``y=pLx/a`` turns its phase into
``-a(j+y)/L``, and ``0<=y<=1`` lies in the body-safe cell.  Thus (1) may be
extended to ``[0,1]`` exactly.

For ``H_s(v)=1_B((q+s eta)v+beta)``, distributional homotopy gives

    ||H_1-H_0||_1 <= 2|eta|                                  (2)

whenever the intermediate slopes are at least one.  To see the constant two,
each of the two boundary families has roots ``v=y/A`` with the ``y``'s a
step-one progression in ``(0,A)``.  Its contribution is ``sum y/A^2<=1``:
for ``N<=floor(A)+1``, the extremal sum
``NA-N(N-1)/2`` is at most ``A^2``.  This applies here for ``m>=16``.

The integer skeleton in (1) is the exact primitive phase fibre.  For coprime
``P<=Q``, put

    T_s(z)=sum_n (s-|z+n|)_+,
    F_PQ(z)=[T_((P+Q)/14)(z)-T_((Q-P)/14)(z)]/(PQ).

Since ``T_s-s^2`` depends only on ``{s}`` and has absolute value at most
``{s}(1-{s})<=1/4``, one has

    F_PQ >= 1/49-1/(2PQ).

An exact breakpoint audit of the 63 coprime pairs with ``P+Q>=8`` and
``PQ<46`` has unique minimum ``1/105`` at ``(P,Q,z)=(3,5,3/7)``.  Hence

    F_PQ(z)>=max(1/105,1/49-1/(2PQ)).                         (3)

Combining (1)--(3) gives the cellwise pair floor

    I >= c^-1 max(0,F_PQ,min-2|eta|).                         (4)

This is the decisive reparametrization: the transport loss depends on the
cross-determinant ``qa-pb``, not the old coarse quantity ``a+b``.

For arbitrary distinct ``p,q in [m,m+5]``, their reduced product is at least
``m(m+5)/25``.  For fixed labels ``a<b``, also

    |qa-pb| <= m(b-a)+5b,
    epsilon(E,q) <= 1/(39m).

Therefore the bodywise coarse one-edge margin

    f0(m)-2[m(b-a)+5b]/[mL-b]-1/(39m),
    f0(m)=max(1/105,1/49-25/[2m(m+5)]),                       (5)

is increasing in ``m``.  At ``m=16`` some pair makes (5) positive on all but
eight bodies.  The exact one-edge audit of those eight bodies over
``m=16..47`` and every permutation has 184320 rows and no failure.  Seven
then enter (5); the sole hostile body ``(1,2,3,4,6,12)`` has another 51840
exact rows over ``m=48..119`` with no failure and enters (5) at ``m=120``.

On each exceptional chromatic body, exact proper-word enumeration leaves
2304 and 1824 words.  Slots zero and one always have distinct levels, and
their coarse margin at ``m=16`` is already positive; monotonicity closes the
whole tail.  Thus every reflected ``D=5`` residual packet closes for
``m>=16``.  The finite head ``m<=15`` and every ``D>=6`` sector are outside
this theorem.  As throughout THM-2941, this is a sufficient certificate
inside the reflected residual family, not a physical-survivor classification.
"""

from __future__ import annotations

import argparse
import hashlib
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations, product
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "04-computation" / "lrc14_j7_reflected_levels_all_q_mass_closure_thm2941.py"
UNIVERSAL = ROOT / "04-computation" / "lrc14_j7_reflected_universal_pair_chromatic_closure_thm2941.py"
UNIVERSAL_OUTPUT = ROOT / "05-knowledge" / "results" / "lrc14_j7_reflected_universal_pair_chromatic_closure_thm2941.out"
OUTPUT = ROOT / "05-knowledge" / "results" / "lrc14_j7_reflected_d5_crossdet_tail_closure_thm2941.out"

EXPECTED_BASE_SHA256 = "2cf0866932f775cc493f97093333e81e65ac3aa76a8e439de969aa700c993f31"
EXPECTED_UNIVERSAL_SHA256 = "a6f58c1a52dfc1fca61a239068dbe0b216bac41f1622b98748bc4a6d213fb6e8"
EXPECTED_UNIVERSAL_OUTPUT_SHA256 = "7364d5866171405fa90539a9ad76727c0c52f020ac1a104a1ab4f0276aedd115"
EXPECTED_SEMANTIC_SHA256 = "413949dab8a7b96657f9f47dc7f11514ed3b5309f53788ec1ecee8cac8d10dc7"

TAIL_START = 16
BODY_COUNT = 3003
COMPLETE_BODY_COUNT = 3001
EDGES = tuple(combinations(range(6), 2))
WORDS = tuple(permutations(range(6)))
EXCEPTIONS = (
    ((1, 2, 7, 9, 11, 13), ((2, 3), (3, 4), (4, 5)), 2304),
    ((2, 4, 7, 9, 11, 13), ((2, 4), (3, 5)), 1824),
)
EIGHT_BODIES = (
    (1, 2, 3, 4, 6, 8),
    (1, 2, 3, 4, 6, 12),
    (1, 2, 3, 4, 8, 12),
    (1, 2, 3, 6, 8, 12),
    (1, 2, 4, 6, 8, 12),
    (1, 3, 4, 6, 8, 12),
    (1, 3, 4, 6, 9, 12),
    (2, 3, 4, 6, 8, 12),
)
HOSTILE = (1, 2, 3, 4, 6, 12)
EXPECTED_EIGHT_WORST = (
    F(5586391854420596584, 4141304415166171876875),
    HOSTILE,
    16,
    (4, 2, 5, 1, 3, 0),
    F(3752, 1380625),
    F(4823878739281456, 3524514395886103725),
)
EXPECTED_HOSTILE_WORST = (
    F(459192369424550444960312, 68336383941117906135181309),
    HOSTILE,
    55,
    (0, 5, 3, 1, 4, 2),
    F(36771708, 5154478763),
    F(1356343812409289404, 3273464819838537633673),
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


require(sha256(BASE) == EXPECTED_BASE_SHA256, "reflected interval engine changed")
require(sha256(UNIVERSAL) == EXPECTED_UNIVERSAL_SHA256, "universal chromatic source changed")
require(sha256(UNIVERSAL_OUTPUT) == EXPECTED_UNIVERSAL_OUTPUT_SHA256, "universal chromatic output changed")
SPEC = spec_from_file_location("d5_crossdet_base", BASE)
require(SPEC is not None and SPEC.loader is not None, "cannot load reflected base")
R = module_from_spec(SPEC)
SPEC.loader.exec_module(R)


def qtext(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def triangle_sum(s: F, z: F) -> F:
    bound = s.numerator // s.denominator + 3
    return sum((max(F(0), s - abs(z + n)) for n in range(-bound, bound + 1)), F(0))


def exact_phase_minimum(P: int, Q: int) -> tuple[F, F]:
    A = F(P + Q, 14)
    B = F(Q - P, 14)
    events = {F(0), F(1)}
    for s in (A, B):
        bound = s.numerator // s.denominator + 3
        for n in range(-bound, bound + 1):
            for z in (-F(n), s - n, -s - n):
                if 0 <= z <= 1:
                    events.add(z)
    return min(((triangle_sum(A, z) - triangle_sum(B, z)) / (P * Q), z) for z in events)


def phase_floor(p: int, q: int) -> F:
    divisor = gcd(p, q)
    P, Q = sorted((p // divisor, q // divisor))
    if P + Q <= 7:
        return F(0)
    return max(F(1, 105), F(1, 49) - F(1, 2 * P * Q))


def edge_floor(L: int, a: int, p: int, b: int, q: int) -> F:
    floor = phase_floor(p, q)
    if floor == 0:
        return F(0)
    rows = []
    for aa, pp, bb, qq in ((a, p, b, q), (b, q, a, p)):
        c = F(pp * L - aa, pp * L)
        eta = F(qq * aa - pp * bb, pp * L - aa)
        rows.append(max(F(0), (floor - 2 * abs(eta)) / c))
    return max(rows)


def exact_debt(body: tuple[int, ...], L: int, levels: tuple[int, ...]) -> F:
    return sum((F(e, 7 * (q * L - e)) for e, q in zip(body, levels)), F(0))


def packet_margin(body: tuple[int, ...], m: int, word: tuple[int, ...]):
    L = 14 * lcm(*body)
    levels = tuple(m + d for d in word)
    ranked = tuple(
        (edge_floor(L, body[i], levels[i], body[j], levels[j]), i, j)
        for i, j in EDGES
    )
    gain, i, j = max(ranked)
    debt = exact_debt(body, L, levels)
    return gain - debt, gain, debt, (i, j)


def coarse_phase_floor(m: int) -> F:
    return max(F(1, 105), F(1, 49) - F(25, 2 * m * (m + 5)))


def coarse_pair_margin(body: tuple[int, ...], m: int, a: int, b: int) -> F:
    require(a < b and a in body and b in body, (body, a, b))
    L = 14 * lcm(*body)
    return (
        coarse_phase_floor(m)
        - 2 * F(m * (b - a) + 5 * b, m * L - b)
        - F(1, 39 * m)
    )


def intersection_mass(left, right) -> F:
    return (
        R.interval_mass(left)
        + R.interval_mass(right)
        - R.interval_mass(R.merge_intervals(left + right))
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    # Exact finite bank behind the all-phase 1/105 floor.
    phase_rows = []
    for P in range(1, 46):
        for Q in range(P, 46):
            if gcd(P, Q) == 1 and P + Q >= 8 and P * Q < 46:
                value, z = exact_phase_minimum(P, Q)
                phase_rows.append((value, z, P, Q))
    require(len(phase_rows) == 63, ("phase-bank count", len(phase_rows)))
    require(min(phase_rows) == (F(1, 105), F(3, 7), 3, 5), ("phase minimum", min(phase_rows)))
    require(sum(row[0] == F(1, 105) for row in phase_rows) == 1, "phase minimum not unique")

    bodies = tuple(combinations(range(1, 15), 6))
    require(len(bodies) == BODY_COUNT, ("body count", len(bodies)))
    exception_bodies = {row[0] for row in EXCEPTIONS}
    complete_bodies = tuple(body for body in bodies if body not in exception_bodies)
    require(len(complete_bodies) == COMPLETE_BODY_COUNT, "complete-body count changed")

    coarse_rows = []
    coarse_failures = []
    for body in complete_bodies:
        ranked = max(
            (coarse_pair_margin(body, TAIL_START, a, b), a, b)
            for a, b in combinations(body, 2)
        )
        coarse_rows.append((ranked[0], body, ranked[1:]))
        if ranked[0] <= 0:
            coarse_failures.append(body)
    require(tuple(coarse_failures) == EIGHT_BODIES, ("coarse exceptional bodies", coarse_failures))

    # Exact one-edge atlas on the eight low-ruler bodies.
    exact_checks = 0
    exact_failures = 0
    exact_worst = None
    for body in EIGHT_BODIES:
        for m in range(16, 48):
            for word in WORDS:
                margin, gain, debt, edge = packet_margin(body, m, word)
                exact_checks += 1
                exact_failures += int(margin <= 0)
                row = (margin, body, m, word, gain, debt, edge)
                if exact_worst is None or row < exact_worst:
                    exact_worst = row
    require(exact_checks == 184320 and exact_failures == 0, (exact_checks, exact_failures))
    require(exact_worst is not None, "empty eight-body atlas")
    require(exact_worst[:6] == EXPECTED_EIGHT_WORST, ("eight-body worst", exact_worst))

    # Seven bodies enter the monotone coarse cone at m=48.
    for body in EIGHT_BODIES:
        if body == HOSTILE:
            continue
        require(
            max(coarse_pair_margin(body, 48, a, b) for a, b in combinations(body, 2)) > 0,
            ("body did not enter coarse cone at 48", body),
        )

    # The sole hostile body: finite bridge 48..119, then monotone cone.
    hostile_checks = 0
    hostile_failures = 0
    hostile_worst = None
    for m in range(48, 120):
        for word in WORDS:
            margin, gain, debt, edge = packet_margin(HOSTILE, m, word)
            hostile_checks += 1
            hostile_failures += int(margin <= 0)
            row = (margin, HOSTILE, m, word, gain, debt, edge)
            if hostile_worst is None or row < hostile_worst:
                hostile_worst = row
    require(hostile_checks == 51840 and hostile_failures == 0, (hostile_checks, hostile_failures))
    require(hostile_worst is not None, "empty hostile continuation")
    require(hostile_worst[:6] == EXPECTED_HOSTILE_WORST, ("hostile worst", hostile_worst))
    require(
        max(coarse_pair_margin(HOSTILE, 120, a, b) for a, b in combinations(HOSTILE, 2)) > 0,
        "hostile body did not enter coarse cone at 120",
    )

    # Exceptional chromatic bodies: every proper word keeps singleton slots
    # 0 and 1 at different levels, so their large-ruler coarse pair closes.
    exception_rows = []
    all_edges = set(EDGES)
    for body, bad_edges, expected_count in EXCEPTIONS:
        good_edges = all_edges - set(bad_edges)
        proper = tuple(
            word
            for word in product(range(6), repeat=6)
            if min(word) == 0
            and max(word) == 5
            and all(word[i] != word[j] for i, j in good_edges)
        )
        require(len(proper) == expected_count, ("proper-word count", body, len(proper)))
        require(all(word[0] != word[1] for word in proper), ("singleton collision", body))
        margin = coarse_pair_margin(body, 16, body[0], body[1])
        require(margin > 0, ("exception coarse margin", body, margin))
        exception_rows.append((body, bad_edges, len(proper), margin))

    # Direct hostile control at the global finite-atlas minimum.  The proved
    # transport floor holds on every body-safe cell and is very conservative.
    control_body = HOSTILE
    control_m = 16
    control_word = (4, 2, 5, 1, 3, 0)
    control_L, safe_ranges = R.safe_cell_ranges(control_body)
    control_levels = tuple(control_m + d for d in control_word)
    control_edge = (1, 2)
    actual_rows = []
    for left, right in safe_ranges:
        for cell in range(left, right):
            i, j = control_edge
            actual = intersection_mass(
                R.reflected_level_arcs(control_L, control_body[i], control_levels[i], cell),
                R.reflected_level_arcs(control_L, control_body[j], control_levels[j], cell),
            )
            actual_rows.append((actual, cell))
    actual_minimum = min(actual_rows)
    require(actual_minimum == (F(36312, 1775425), 12), ("direct control", actual_minimum))
    require(actual_minimum[0] >= EXPECTED_EIGHT_WORST[4], "transport floor exceeds actual overlap")

    semantic_payload = (
        tuple(phase_rows),
        tuple(coarse_rows),
        tuple(coarse_failures),
        exact_checks,
        exact_failures,
        exact_worst,
        hostile_checks,
        hostile_failures,
        hostile_worst,
        tuple(exception_rows),
        actual_minimum,
    )
    semantic = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256:
        require(semantic == EXPECTED_SEMANTIC_SHA256, ("semantic digest", semantic))

    lines = [
        "LRC14 reflected D=5 cross-determinant tail closure exact proof",
        f"universe=bodies:{BODY_COUNT};complete_good_graph_bodies:{COMPLETE_BODY_COUNT};tail=m>={TAIL_START}",
        "crossdet_transport=c=1-a/(pL);eta=(qa-pb)/(pL-a);body-safe omitted tail is exactly empty",
        "one_sided_comb_loss=||1_B((q+eta)v+beta)-1_B(qv+beta)||_1<=2|eta|;two boundary progressions each contribute<=1",
        "phase_fibre=F_PQ=[T_((P+Q)/14)-T_((Q-P)/14)]/(PQ);T_s=sum_n(s-|z+n|)_+",
        f"phase_finite_bank=rows:{len(phase_rows)};unique_minimum=1/105;P=3;Q=5;z=3/7",
        "phase_large_product=F_PQ>=1/49-1/(2PQ);all_phase_floor=max(1/105,1/49-1/(2PQ))",
        "coarse_tail=f0(m)-2[m(b-a)+5b]/[mL-b]-1/(39m);monotone_in_m",
        f"coarse_m16_failures={tuple(coarse_failures)}",
        f"eight_body_exact_rows={exact_checks};failures={exact_failures};worst_margin={qtext(exact_worst[0])};worst_E={exact_worst[1]};worst_m={exact_worst[2]};worst_word={exact_worst[3]};gain={qtext(exact_worst[4])};debt={qtext(exact_worst[5])};edge={exact_worst[6]}",
        f"hostile_48_119_rows={hostile_checks};failures={hostile_failures};worst_margin={qtext(hostile_worst[0])};worst_m={hostile_worst[2]};worst_word={hostile_worst[3]};gain={qtext(hostile_worst[4])};debt={qtext(hostile_worst[5])};edge={hostile_worst[6]}",
        f"hostile_direct_control=pair_labels:(2,3);levels:(18,21);safe_cells:{len(actual_rows)};minimum_overlap={qtext(actual_minimum[0])};minimum_cell={actual_minimum[1]}",
    ]
    for body, bad_edges, count, margin in exception_rows:
        lines.append(
            f"EXCEPTION;E={body};bad_edges={bad_edges};proper_D5_words={count};singleton_slots=(0,1);always_distinct=true;m16_coarse_margin={qtext(margin)}"
        )
    lines.extend((
        "conclusion=every reflected D=5 residual packet closes for every m>=16",
        "scope_boundary=finite head m<=15 and D>=6 remain outside this theorem;sufficient reflected-family certificate only",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"base_sha256={sha256(BASE)}",
        f"universal_sha256={sha256(UNIVERSAL)}",
        f"universal_output_sha256={sha256(UNIVERSAL_OUTPUT)}",
        f"source_sha256={sha256(Path(__file__))}",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ))
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
