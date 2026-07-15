#!/usr/bin/env python3
"""Exact replay for THM-800, the Hamming-two common-scale theorem.

There are two logically different ingredients.

1.  A new *uniform* lifted-splice argument.  If

        A = (c*[12] - {c*r,c*s}) union {w_r,w_s},
        w_i == c*i (mod 13),  13 does not divide c,

    were exactly tight at ``1/13``, inspect the left-hand safe germ at every
    lifted missing-owner splice.  Endpoint phases eligible to cover that germ
    lie in the half-open arc ``(-1/13,1/13]``.  A deck of order ``D`` places at
    most ``ceil(2D/13)`` sheets in this arc.  Capacity and the exact order-two
    parity table force both deck orders ``c/gcd(c,w_i)`` to equal one.  Thus a
    tight packet would descend to a scale-one double lift.

2.  The sharp scale-one floor, first computed in the S52/HYP-4103 artifact.
    Every proper double lift

        ([12] - {r,s}) union {r+13*j,s+13*k},  j,k >= 1,

    has ``M >= 2/25``.  A periodic-danger measure bound reduces to the exact
    box ``w_b <= 258`` and ``w_a <= 24*w_b``.  This script replays all 600,756
    rows in that box without importing the old script (whose unrelated E3
    ladder section crashes).  Every row has either a missing-divisor witness
    or an explicitly checked rational ``2/25`` witness.  The bound is sharp at
    ``{4,6}->{17,19}``, whose exact maximin is ``2/25`` at ``6/25``.

Together these statements prove strict Hamming-two looseness at every AP
scale.  They do *not* prove an all-scale ``2/25`` floor: the deck descent uses
the oriented boundary germ specifically at equality ``M=1/13``.

Tournament Analysis / assumption challenge
------------------------------------------
The exact vertices are lifted splice sheets and the two replacement colours;
the incidence bit is membership of the oriented half-open phase arc.  Runners,
unoriented endpoints, and closed danger arcs all destroy the germ orientation.
For telemetry only, use the twelve missing-owner residue labels as vertices.
The pair observable says whether an order-two replacement of one label covers
the opposite parity at the other label.  A live pair is oriented provider to
recipient; silent pairs use increasing residue as the tie gauge.  Reversing
the live-pair switch flips 24 edges.  The tournament fingerprints are printed,
but the proof uses the labelled incidence relation and the fact that the two
opposite cross obligations cannot both be live.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import gcd
from typing import Sequence


DELTA = F(1, 13)
BETA = F(2, 25)
BASE = tuple(range(1, 13))


def fmt(x: F) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def signed_phase(x: F) -> F:
    """Representative in ``(-1/2,1/2]``."""

    z = F(x.numerator % x.denominator, x.denominator)
    return z - 1 if z > F(1, 2) else z


def norm_fraction(x: F) -> F:
    return abs(signed_phase(x))


def left_germ_eligible(x: F) -> bool:
    """A positive-speed tooth covers immediately to the left of this phase."""

    z = signed_phase(x)
    return -DELTA < z <= DELTA


def deck_capacity(order: int) -> int:
    """Maximum points of a uniform order-D grid in ``(-1/13,1/13]``."""

    return (2 * order + 12) // 13  # ceil(2D/13)


def deck_phase(c: int, owner: int, replacement: int, sheet: int) -> F:
    inverse = pow(owner, -1, 13)
    return F(replacement * (inverse + 13 * sheet), 13 * c)


def deck_and_parity_audit() -> None:
    print("ORIENTED_SPLICE_DECK_AUDIT")
    print("eligible_endpoint_arc=(-1/13,1/13]  orientation=left_core_safe_germ")

    # The uniform tail is algebraic: for D>=7, ceil(2D/13)<=D/3 because
    # 2D/13+1 <= D/3 once D>=39/7.  D=3..6 are the four explicit atoms.
    for order in range(3, 7):
        assert deck_capacity(order) * 3 <= order
    for order in range(7, 4097):
        assert deck_capacity(order) * 3 <= order
    assert deck_capacity(1) == 1
    assert deck_capacity(2) == 1
    print("capacity=ceil(2D/13); D=1 -> 1, D=2 -> 1, D>=3 -> capacity/D<=1/3")
    print("two-colour sheet cover => one D=1, or both D=2")

    order_two_cross = []
    complementary_rows = 0
    for owner in BASE:
        inverse = pow(owner, -1, 13)
        own_even = left_germ_eligible(F(1, 13))
        own_odd = left_germ_eligible(F(1, 13) + F(1, 2))
        assert own_even and not own_odd

        for provider in BASE:
            if provider == owner:
                continue
            ratio = provider * inverse % 13
            cross_even = left_germ_eligible(F(ratio, 13))
            cross_odd = left_germ_eligible(F(ratio, 13) + F(1, 2))
            assert not cross_even  # +1 is the own label; -1 is excluded.
            assert cross_odd == (ratio in {6, 7})
            if ratio == 12:
                complementary_rows += 1
                assert not cross_odd
            if cross_odd:
                order_two_cross.append((provider, owner, ratio))

    assert complementary_rows == 12
    assert len(order_two_cross) == 24
    for provider, owner, ratio in order_two_cross:
        reverse = owner * pow(provider, -1, 13) % 13
        assert ratio in {6, 7}
        assert reverse in {11, 2}
        assert reverse not in {6, 7}
    print("D=2 own colour covers even parity; cross covers odd parity iff provider/owner in {6,7}")
    print("reverse ratios are {6^-1,7^-1}={11,2}; simultaneous two-owner cover impossible")

    # An order-one replacement has constant phase +1/13 at its own owner.
    # At every distinct owner its ratio is never +1; the only boundary case is
    # -1/13, which the left-oriented half-open arc deliberately excludes.
    order_one_cross_hits = 0
    for label in BASE:
        for owner in BASE:
            ratio = label * pow(owner, -1, 13) % 13
            eligible = left_germ_eligible(F(ratio, 13))
            assert eligible == (label == owner)
            if label != owner and eligible:
                order_one_cross_hits += 1
    assert order_one_cross_hits == 0
    print("D=1 covers every sheet at its own owner and zero sheets at every distinct owner")
    print("therefore a mixed D=1/D>1 cover is impossible; exact tightness forces D_r=D_s=1")

    # Closed endpoint incidence alone misses the complementary-owner issue.
    # This explicit row is retained as a guardrail: wr covers all endpoints,
    # yet ws fails the oriented left-germ condition on one sheet.
    c, r, s, wr, ws = 2, 1, 12, 28, 11
    assert wr % 13 == (c * r) % 13
    assert ws % 13 == (c * s) % 13
    assert c // gcd(c, wr) == 1 and c // gcd(c, ws) == 2
    closed_masks = []
    oriented_masks = []
    for owner in (r, s):
        closed = []
        oriented = []
        for w in (wr, ws):
            phases = tuple(deck_phase(c, owner, w, ell) for ell in range(c))
            closed.append(tuple(norm_fraction(z) <= DELTA for z in phases))
            oriented.append(tuple(left_germ_eligible(z) for z in phases))
        closed_masks.append(tuple(closed))
        oriented_masks.append(tuple(oriented))
    assert all(
        a or b for owner_masks in closed_masks for a, b in zip(*owner_masks)
    )
    assert any(
        not (a or b) for owner_masks in oriented_masks for a, b in zip(*owner_masks)
    )
    print(
        "complementary_guardrail=(c=2,r=1,s=12,wr=28,ws=11) "
        f"closed_masks={closed_masks} oriented_masks={oriented_masks}"
    )
    print("interpretation=closed endpoint sheets give a false mixed-deck survivor; the half-open germ removes it")
    print()


def distance_mod(n: int, q: int) -> int:
    z = n % q
    return min(z, q - z)


def margin_num(speeds: Sequence[int], a: int, q: int) -> int:
    return min(distance_mod(a * v, q) for v in speeds)


def missing_divisor(speeds: Sequence[int]) -> int | None:
    for modulus in range(2, 13):
        if not any(v % modulus == 0 for v in speeds):
            return modulus
    return None


def scan_beta_certificate(speeds: Sequence[int], wa: int, wb: int) -> tuple[int, int, int] | None:
    qlist = (
        [13 * u for u in range(2, 21)]
        + [q for q in range(8, 41) if q % 13]
        + [25, 50, wa + wb, abs(wa - wb), wb + 1, wa + 1]
    )
    seen: set[int] = set()
    for q in qlist:
        if q < 2 or q in seen:
            continue
        seen.add(q)
        if any(v % q == 0 for v in speeds):
            continue
        for a in range(1, q // 2 + 1):
            if gcd(a, q) != 1:
                continue
            margin = margin_num(speeds, a, q)
            if 25 * margin >= 2 * q:
                return q, a, margin
    return None


def scale_one_floor_audit() -> None:
    print("SCALE_ONE_DOUBLE_LIFT_FLOOR")

    # Ten-runner core: M>=1/11.  Fatten a maximizing point down to beta.
    radius_ten = (F(1, 11) - BETA) / 12
    fee = radius_ten * (1 - 4 * BETA) / BETA
    wb_max = int(F(2, 1) / fee)
    assert radius_ten == F(1, 1100)
    assert fee == F(17, 2200)
    assert wb_max == 258
    assert F(2, 259) < fee

    # Once wb is fixed, the eleven-runner core has M>=1/12.  Its beta-safe
    # interval has length 1/(150 wb), while a wa tooth has length 4/(25 wa).
    # A single tooth can cover it only if wa<=24 wb.
    gap_eleven = F(1, 12) - BETA
    assert gap_eleven == F(1, 300)
    ratio_cap = 24
    print(
        f"beta={fmt(BETA)} ten_core_radius={fmt(radius_ten)} fee={fmt(fee)} "
        f"both_large_cutoff=259"
    )
    print(
        "periodic-danger bound: |D_w(beta) intersect I| <= 2*beta*|I| + 2*beta/w"
    )
    print("outer reduction: smaller replacement wb<=258 and larger replacement wa<=24*wb")

    stats = Counter(total=0, sieve=0, rational=0, failure=0)
    sieve_histogram: Counter[int] = Counter()
    q_histogram: Counter[int] = Counter()
    equality_certificates = 0
    minimum_certificate: F | None = None
    maximum_q = 0
    digest = sha256()

    for small_label in BASE:
        k = 1
        while True:
            wb = small_label + 13 * k
            if wb > wb_max:
                break
            for large_label in BASE:
                if large_label == small_label:
                    continue
                core = tuple(v for v in BASE if v not in (small_label, large_label))
                j_low = max(1, (wb - large_label + 12) // 13)
                j_high = (ratio_cap * wb - large_label) // 13
                for j in range(j_low, j_high + 1):
                    wa = large_label + 13 * j
                    if wa < wb:
                        continue
                    speeds = tuple(sorted((*core, wa, wb)))
                    assert len(speeds) == 12 and len(set(speeds)) == 12
                    assert wa % 13 == large_label and wb % 13 == small_label
                    stats["total"] += 1

                    modulus = missing_divisor(speeds)
                    if modulus is not None:
                        margin = F(margin_num(speeds, 1, modulus), modulus)
                        assert margin >= F(1, modulus) >= F(1, 12) > BETA
                        stats["sieve"] += 1
                        sieve_histogram[modulus] += 1
                        digest.update(
                            (",".join(map(str, speeds)) + f"|S|{modulus}\n").encode()
                        )
                        continue

                    certificate = scan_beta_certificate(speeds, wa, wb)
                    if certificate is None:
                        stats["failure"] += 1
                        digest.update((",".join(map(str, speeds)) + "|FAIL\n").encode())
                        continue

                    q, a, numerator = certificate
                    margin = F(numerator, q)
                    assert margin >= BETA
                    stats["rational"] += 1
                    q_histogram[q] += 1
                    maximum_q = max(maximum_q, q)
                    if margin == BETA:
                        equality_certificates += 1
                    if minimum_certificate is None or margin < minimum_certificate:
                        minimum_certificate = margin
                    digest.update(
                        (",".join(map(str, speeds)) + f"|R|{q}|{a}|{numerator}\n").encode()
                    )
            k += 1
        print(
            f"small_label={small_label:2d} cumulative_total={stats['total']:6d} "
            f"sieve={stats['sieve']:6d} rational={stats['rational']:6d} "
            f"fail={stats['failure']}"
        )

    assert stats == Counter(total=600756, sieve=393962, rational=206794, failure=0)
    assert minimum_certificate == BETA
    print(f"final_stats={dict(stats)}")
    print(f"sieve_modulus_histogram={dict(sorted(sieve_histogram.items()))}")
    print(f"rational_q_histogram={dict(sorted(q_histogram.items()))}")
    print(
        f"minimum_rational_certificate={fmt(minimum_certificate)} "
        f"equality_certificate_rows={equality_certificates} maximum_q={maximum_q}"
    )
    print(f"certificate_digest={digest.hexdigest()}")
    print("result=no scale-one proper double lift lies below 2/25")
    print()


def exact_maximin(speeds: Sequence[int]) -> tuple[F, tuple[F, ...]]:
    """Exact piecewise-linear maximum via all cusps and pair crossings."""

    denominators = {2 * v for v in speeds}
    for a, b in combinations(speeds, 2):
        denominators.add(a + b)
        denominators.add(abs(a - b))
    denominators.discard(0)

    best = F(0)
    witnesses: set[F] = set()
    for q in denominators:
        for a in range(1, q):
            value = F(margin_num(speeds, a, q), q)
            t = F(a, q)
            if value > best:
                best = value
                witnesses = {t}
            elif value == best:
                witnesses.add(t)
    return best, tuple(sorted(witnesses))


def sharpness_audit() -> None:
    speeds = tuple(sorted((set(BASE) - {4, 6}) | {17, 19}))
    value, witnesses = exact_maximin(speeds)
    assert value == BETA
    assert witnesses == (F(6, 25), F(19, 25))
    vector = tuple(distance_mod(6 * v, 25) for v in speeds)
    binders = tuple(v for v, z in zip(speeds, vector) if z == 2)
    assert binders == (8, 17)

    print("SHARPNESS")
    print(f"speeds={speeds}")
    print(
        f"exact_M={fmt(value)} maximizers={[fmt(t) for t in witnesses]} "
        f"binders_at_6/25={binders} vector={vector}"
    )
    print()


def build_residue_tournament(reverse_live: bool) -> tuple[list[list[bool]], int]:
    n = 12
    adjacency = [[False] * n for _ in range(n)]
    live_edges = 0
    for low in BASE:
        for high in range(low + 1, 13):
            high_covers_low = high * pow(low, -1, 13) % 13 in {6, 7}
            low_covers_high = low * pow(high, -1, 13) % 13 in {6, 7}
            assert not (high_covers_low and low_covers_high)
            if high_covers_low or low_covers_high:
                live_edges += 1
                provider = high if high_covers_low else low
                recipient = low if high_covers_low else high
                if reverse_live:
                    provider, recipient = recipient, provider
                adjacency[provider - 1][recipient - 1] = True
            else:
                adjacency[low - 1][high - 1] = True  # tie gauge
    return adjacency, live_edges


def tournament_fingerprint(adjacency: Sequence[Sequence[bool]]) -> tuple[dict[int, int], int, list[int], int, list[int]]:
    n = len(adjacency)
    scores = [sum(row) for row in adjacency]
    score_histogram = dict(sorted(Counter(scores).items()))

    triangles = 0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                if (
                    adjacency[i][j] and adjacency[j][k] and adjacency[k][i]
                ) or (
                    adjacency[i][k] and adjacency[k][j] and adjacency[j][i]
                ):
                    triangles += 1

    reach = [[i == j or adjacency[i][j] for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                reach[i][j] = reach[i][j] or (reach[i][k] and reach[k][j])
    unseen = set(range(n))
    scc_sizes: list[int] = []
    while unseen:
        root = min(unseen)
        component = {v for v in unseen if reach[root][v] and reach[v][root]}
        scc_sizes.append(len(component))
        unseen -= component
    scc_sizes.sort(reverse=True)

    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp[mask][last]
            if count == 0:
                continue
            for nxt in range(n):
                if not (mask >> nxt) & 1 and adjacency[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += count
    hamiltonian_paths = sum(dp[-1])

    # Redei insertion, with residue order as the deterministic tie path input.
    path: list[int] = []
    for vertex in range(n):
        position = 0
        while position < len(path) and adjacency[path[position]][vertex]:
            position += 1
        path.insert(position, vertex)
    assert all(adjacency[path[i]][path[i + 1]] for i in range(n - 1))
    return score_histogram, triangles, scc_sizes, hamiltonian_paths, [v + 1 for v in path]


def tournament_audit() -> None:
    forward, live = build_residue_tournament(False)
    reverse, reverse_live = build_residue_tournament(True)
    assert live == reverse_live == 24
    edge_flips = sum(
        forward[i][j] != reverse[i][j]
        for i in range(12)
        for j in range(i + 1, 12)
    )
    assert edge_flips == 24
    first = tournament_fingerprint(forward)
    second = tournament_fingerprint(reverse)

    print("TOURNAMENT_TELEMETRY")
    print("vertices=nonzero residue owners pair_observable=order-two cross-parity eligibility")
    print("orientation=provider_to_recipient tie_gauge=increasing_residue switch=reverse_live_edges")
    print(f"live_edges={live} edge_flips_under_switch={edge_flips}")
    print(
        f"forward_score_histogram={first[0]} directed_triangles={first[1]} "
        f"scc_sizes={first[2]} hamiltonian_paths={first[3]} tie_hamiltonian_path={first[4]}"
    )
    print(
        f"reverse_score_histogram={second[0]} directed_triangles={second[1]} "
        f"scc_sizes={second[2]} hamiltonian_paths={second[3]} tie_hamiltonian_path={second[4]}"
    )
    print("interpretation=tournament is a cyclic shadow; the proof needs labelled two-colour sheet incidence")
    print()


def main() -> None:
    print("LRC13 HAMMING-TWO COMMON-SCALE RIGIDITY — EXACT REPLAY")
    print("arithmetic=integer+fractions.Fraction floating_point=none")
    print("scope=exact tightness descent at 1/13 plus sharp scale-one floor 2/25")
    print()

    deck_and_parity_audit()
    scale_one_floor_audit()
    sharpness_audit()
    tournament_audit()

    print(
        "PASS: exact tightness forces common scale; the descended proper double lift has M>=2/25>1/13."
    )


if __name__ == "__main__":
    main()
