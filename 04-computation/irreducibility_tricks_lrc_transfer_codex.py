#!/usr/bin/env python3
"""Irreducibility-trick transfer atlas for the LRC(14) proof search.

This script continues HYP-2451/HYP-2452 and applies their polynomial
irreducibility grammar to the live LRC14 boundary in HYP-2470.

The guiding dictionary is:

    primitive polynomial              primitive speed row
    mod-p irreducibility test         Q27 / small-shell denominator gate
    Eisenstein/Newton valuation       prime-ideal carry/valuation channel
    hidden convolution lift           hidden obligation-cover lift
    Cohn/Singh factor capture         finite resource allocation budget

The concrete test case is the two HYP-2470 four-deletion addresses.  They are
Q27-feasible but open immediately at the missing plain shells 33 and 31.  The
useful pattern is almost Eisenstein-like: 12 of 13 speeds sit in the 7-ideal,
and the single primitive escape speed is a 13-clock token.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from functools import reduce
from itertools import combinations
from math import gcd
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
sys.path.append(str(HERE))

from lrc14_near_core_q27_setcover_codex import CORE_T, Q27, blocks, obligations
from lrc14_pisano_band_ladder_codex import (
    bprime_any,
    first_witness_in,
    min_witness_modulus,
    q_lattice,
)
from lrc14_lebesgue_wall_s676 import depth_sweep


PRIMES = (2, 3, 5, 7, 11, 13, 17, 23, 31, 37, 41)
MISSING_PLAIN_TO_41 = tuple(q for q in range(2, 42) if q not in set(Q27))


@dataclass(frozen=True)
class ExceptionPacket:
    label: str
    deleted: tuple[int, ...]
    added: tuple[int, ...]

    @property
    def retained(self) -> tuple[int, ...]:
        return tuple(v for v in CORE_T if v not in self.deleted)

    @property
    def row(self) -> tuple[int, ...]:
        return tuple(sorted(self.retained + self.added))


EXCEPTIONS = (
    ExceptionPacket("E1", (28, 42, 56, 84), (91, 322, 350, 504, 936)),
    ExceptionPacket("E2", (42, 56, 70, 84), (91, 119, 700, 1008, 1066)),
)


@dataclass(frozen=True)
class Trick:
    name: str
    polynomial_move: str
    lrc_transfer: str
    proof_use: str
    risk: str
    scores: tuple[int, int, int, int, int, int]


TRICKS = (
    Trick(
        "mod_p_convolution_blocker",
        "Find a prime p for which f mod p has no nontrivial split.",
        "Find a denominator family, such as Q27, whose safe-twist obligations cannot all be covered.",
        "Turns local survivor emptiness into a global no-counterexample certificate in a normalized class.",
        "Fails on ramified packets where residue data collapse.",
        (5, 5, 4, 5, 5, 3),
    ),
    Trick(
        "eisenstein_newton_valuation",
        "Use p-adic heights/Newton slopes when f mod p looks reducible.",
        "Use prime-ideal carry data: many speeds in one ideal, one primitive escape, next-shell witness.",
        "Explains the two HYP-2470 exceptions as ramified local packets rather than random misses.",
        "Needs a clean inequality replacing the scan over q=31,33.",
        (5, 4, 5, 5, 4, 5),
    ),
    Trick(
        "gauss_primitive_normal_form",
        "Strip content; prove irreducibility for the primitive part.",
        "Normalize rows by gcd, core retention, carry window, and Bprime/opening channels.",
        "Prevents wasting effort on nonprimitive 7-multiple shadows.",
        "Can hide large outside-window runners unless paired with a normalizer.",
        (5, 5, 3, 4, 5, 4),
    ),
    Trick(
        "hensel_recombination_ladder",
        "Factor locally, then ask whether local factors recombine globally.",
        "Lift Q27 residue/fiber packets through divisor shells and missing plain-shell layers.",
        "Suggests storing survivor ledgers instead of only first-witness denominators.",
        "May be too detailed unless the state space is compressed.",
        (4, 4, 4, 4, 4, 5),
    ),
    Trick(
        "cohn_perron_digit_dominance",
        "Large-base digit/value dominance can certify irreducibility.",
        "Large speed as a carry digit: either it opens Bprime or folds to a lower carry address.",
        "A route for outside-window normalization after the eight-core carry-window theorem.",
        "Dominance constants are the hard part.",
        (4, 5, 4, 5, 3, 5),
    ),
    Trick(
        "singh_factor_capture_budget",
        "A value f(m) with small Omega limits possible irreducible factors.",
        "Treat blocker obligations as tokens to allocate among added speeds; small budget forces failure.",
        "Turns HYP-2465/HYP-2470 MILP certificates into proof-shaped resource inequalities.",
        "The exact Singh theorem has hypotheses; the LRC use is an analogy until formalized.",
        (4, 4, 5, 5, 5, 4),
    ),
    Trick(
        "fixed_divisor_admissibility",
        "Bunyakovsky admissibility forbids one prime from dividing every value.",
        "No single residue/clock should be allowed to cover all twists without a compensating opening.",
        "Pushes the 13-clock residual toward a finite exceptional atlas.",
        "Admissibility is necessary, not sufficient.",
        (3, 4, 4, 4, 5, 3),
    ),
    Trick(
        "galois_cycle_type_sequence",
        "Factorization mod p gives Frobenius cycle-type data.",
        "First-witness denominator sequence gives cover-type data: which shells open which packets.",
        "Could classify below-eight-core packets by cycle-type-like signatures.",
        "Mostly descriptive until a Chebotarev-like replacement lemma is found.",
        (3, 3, 3, 4, 4, 3),
    ),
    Trick(
        "discriminant_ramification_map",
        "Bad primes divide the discriminant and mark ramified local behavior.",
        "Sparse unresolved deletion addresses mark ramification: add shell 31/33/41 or Bprime.",
        "Treats finite exceptions as data to resolve, not as theorem failure.",
        "Requires avoiding post-hoc overfitting.",
        (4, 4, 5, 5, 4, 4),
    ),
    Trick(
        "integral_convolution_lift_ilp",
        "Reducibility is existence of an integer coefficient-grid lift.",
        "Q-blocking is existence of an obligation-cover lift by added speeds.",
        "Makes exact infeasibility certificates the analogue of irreducibility certificates.",
        "MILP dual certificates must be extracted to become human proof.",
        (5, 5, 5, 5, 5, 4),
    ),
)


def prime_factorization(n: int) -> dict[int, int]:
    n = abs(n)
    out: dict[int, int] = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def factor_string(n: int) -> str:
    fac = prime_factorization(n)
    pieces = []
    for p, e in sorted(fac.items()):
        pieces.append(str(p) if e == 1 else f"{p}^{e}")
    return "*".join(pieces) if pieces else "1"


def row_gcd(row: tuple[int, ...]) -> int:
    return reduce(gcd, row, 0)


def v_p(n: int, p: int) -> int:
    n = abs(n)
    out = 0
    while n and n % p == 0:
        out += 1
        n //= p
    return out


def valuation_histogram(row: tuple[int, ...], p: int) -> dict[int, int]:
    return dict(sorted(Counter(v_p(v, p) for v in row).items()))


def safe_witnesses(row: tuple[int, ...], q: int) -> tuple[tuple[int, tuple[int, ...]], ...]:
    band = q // 14
    out = []
    for a in range(1, q):
        if gcd(a, q) != 1:
            continue
        distances = tuple(min((a * v) % q, q - (a * v) % q) for v in row)
        if all(d > band for d in distances):
            out.append((a, distances))
    return tuple(out)


def distance_profile(row: tuple[int, ...], q: int, a: int) -> tuple[tuple[int, int], ...]:
    return tuple((v, min((a * v) % q, q - (a * v) % q)) for v in row)


def cover_count_for_speed(v: int, obs: tuple[tuple[int, int], ...]) -> int:
    return sum(1 for q, a in obs if blocks(v, q, a))


def exception_diagnostics(packet: ExceptionPacket) -> dict[str, object]:
    row = packet.row
    q27 = first_witness_in(list(row), list(Q27))
    q41 = first_witness_in(list(row), q_lattice(41))
    plain_min_q = min_witness_modulus(list(row), 1200)
    sweep = depth_sweep(row)
    base_obs = obligations(packet.retained, Q27)
    missing_openings = []
    for q in MISSING_PLAIN_TO_41:
        witnesses = safe_witnesses(row, q)
        if witnesses:
            missing_openings.append((q, witnesses[0][0], len(witnesses)))
    non7 = tuple(v for v in packet.added if v % 7 != 0)
    added_cover = tuple((v, cover_count_for_speed(v, base_obs)) for v in packet.added)
    added_cover = tuple(sorted(added_cover, key=lambda vc: (-vc[1], vc[0])))
    return {
        "label": packet.label,
        "deleted": packet.deleted,
        "retained": packet.retained,
        "added": packet.added,
        "row": row,
        "gcd": row_gcd(row),
        "q27": q27,
        "q41": q41,
        "plain_min_q": plain_min_q,
        "bprime": bprime_any(list(row)),
        "p0": sweep.p0,
        "components": sweep.positive_safe_components,
        "safe_points": len(sweep.safe_points),
        "base_q27_obligations": len(base_obs),
        "added_cover": added_cover,
        "non7": non7,
        "ideal7_count": sum(1 for v in row if v % 7 == 0),
        "non7_all_13": all(v % 13 == 0 for v in non7),
        "valuation_hists": {p: valuation_histogram(row, p) for p in (7, 13, 31, 41)},
        "missing_openings": tuple(missing_openings),
    }


def print_trick_transfer_table() -> None:
    print("A. Irreducibility tricks as LRC proof carriers")
    for trick in TRICKS:
        print(f"  {trick.name}")
        print(f"    poly: {trick.polynomial_move}")
        print(f"    LRC : {trick.lrc_transfer}")
        print(f"    use : {trick.proof_use}")
        print(f"    risk: {trick.risk}")
    print()


def print_exception_section() -> None:
    print("B. HYP-2470 exceptions through the irreducibility lens")
    print(f"  Q27 size={len(Q27)}; missing plain shells through 41={MISSING_PLAIN_TO_41}")
    for packet in EXCEPTIONS:
        diag = exception_diagnostics(packet)
        print(f"  {packet.label}: deleted={diag['deleted']} added={diag['added']}")
        print(f"    row={diag['row']}")
        print(
            f"    gcd={diag['gcd']} Q27={diag['q27']} first_plain_q={diag['plain_min_q']} "
            f"Q41={diag['q41']} Bprime={diag['bprime']}"
        )
        print(
            f"    p0={diag['p0']} components={diag['components']} safe_points={diag['safe_points']}"
        )
        print(
            f"    base_Q27_obligations={diag['base_q27_obligations']} "
            f"added_Q27_cover_counts={diag['added_cover']}"
        )
        print(
            f"    7-ideal occupancy={diag['ideal7_count']}/13; "
            f"primitive_non7={diag['non7']}; non7_all_13_clock={diag['non7_all_13']}"
        )
        for v in packet.added:
            flags = []
            if v % 7 == 0:
                flags.append("7-ideal")
            if v % 13 == 0:
                flags.append("13-clock")
            if v % 7 != 0:
                flags.append("primitive-escape")
            print(f"      added {v:4d} = {factor_string(v):12s} {'/'.join(flags)}")
        print(f"    valuation histograms={diag['valuation_hists']}")
        print(f"    first missing-shell openings={diag['missing_openings'][:4]}")
        if diag["q41"] is not None:
            q, a = diag["q41"]
            print(f"    witness profile at q={q}, a={a}, band={q // 14}:")
            print(f"      {distance_profile(packet.row, q, a)}")
    print()


def tournament_fingerprint() -> dict[str, object]:
    names = [t.name for t in TRICKS]
    scores = {t.name: t.scores for t in TRICKS}
    out = {name: set() for name in names}
    edge_flips = 0
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            av = sum(x > y for x, y in zip(scores[a], scores[b]))
            bv = sum(y > x for x, y in zip(scores[a], scores[b]))
            if av >= bv:
                out[a].add(b)
            else:
                out[b].add(a)
                edge_flips += 1

    cycles = 0
    for a, b, c in combinations(names, 3):
        degs = {a: 0, b: 0, c: 0}
        for x, y in ((a, b), (a, c), (b, c)):
            if y in out[x]:
                degs[x] += 1
            else:
                degs[y] += 1
        if sorted(degs.values()) == [1, 1, 1]:
            cycles += 1

    index = {name: i for i, name in enumerate(names)}
    n = len(names)
    adj = [[False] * n for _ in range(n)]
    for a in names:
        for b in out[a]:
            adj[index[a]][index[b]] = True

    # Transitive closure for SCC sizes.
    reach = [row[:] for row in adj]
    for i in range(n):
        reach[i][i] = True
    for k in range(n):
        for i in range(n):
            if reach[i][k]:
                for j in range(n):
                    reach[i][j] = reach[i][j] or reach[k][j]
    seen = set()
    scc_sizes = []
    for i in range(n):
        if i in seen:
            continue
        comp = {j for j in range(n) if reach[i][j] and reach[j][i]}
        seen |= comp
        scc_sizes.append(len(comp))

    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            if not dp[mask][last]:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += dp[mask][last]
    hpaths = sum(dp[(1 << n) - 1])

    ranking = sorted(names, key=lambda name: (-len(out[name]), name))
    return {
        "score_hist": dict(sorted(Counter(len(out[name]) for name in names).items())),
        "directed_3cycles": cycles,
        "scc_sizes": sorted(scc_sizes, reverse=True),
        "hamiltonian_paths": hpaths,
        "edge_flips_vs_tie_path": edge_flips,
        "ranking": ranking,
        "out": out,
        "scores": scores,
    }


def print_tournament_analysis() -> None:
    fp = tournament_fingerprint()
    print("C. Tournament Analysis")
    print("  vertices: proof carriers/tricks, not runners or arcs")
    print(
        "  observable=(exactness, local_to_global, exception_power, "
        "LRC_leverage, computability, descent_power)"
    )
    print("  switch/gauge: primitive-normal-form first; scalar q-blocking is only the residue face")
    print("  tie Hamiltonian path: the listed trick order in this script")
    print(
        f"  fingerprint: score_hist={fp['score_hist']} "
        f"directed_3cycles={fp['directed_3cycles']} "
        f"scc_sizes={fp['scc_sizes']} "
        f"hamiltonian_paths={fp['hamiltonian_paths']} "
        f"edge_flips_vs_tie_path={fp['edge_flips_vs_tie_path']}/45"
    )
    for name in fp["ranking"]:
        print(f"    {name:34s} out={len(fp['out'][name])} scores={fp['scores'][name]}")
    print()


def print_proof_program() -> None:
    print("D. Proof-route hypotheses generated by the transfer")
    items = (
        (
            "Residue-face theorem",
            "HYP-2465+HYP-2470 already say: retaining >=8 core speeds in the carry window "
            "forces Q27 or q<=41.  The next proof artifact should extract reusable dual "
            "set-cover certificates, analogous to a clean mod-p irreducibility certificate.",
        ),
        (
            "Ramification/Eisenstein portal",
            "The only four-deletion Q27 exceptions have 12/13 speeds in the 7-ideal and "
            "a single non-7 primitive escape divisible by 13.  Prove this shape always "
            "opens at q in {31,33,41} or by Bprime, instead of relying on enumeration.",
        ),
        (
            "Factor-capture budget lemma",
            "Model each Q27 safe twist as a prime-token-like obligation.  With <=e+1 added "
            "speeds, the obligation allocation should be impossible for e<=4 except the "
            "ramified portal above; this is the LRC analogue of small Omega(f(m)).",
        ),
        (
            "Cohn/Perron outside-window normalizer",
            "A very large speed should behave like a dominant base digit: either it creates "
            "a Bprime interval, or its carry residue can be folded into the 1..1092 window.",
        ),
        (
            "Cycle-type signature atlas",
            "For below-eight-core rows, store the ordered first-witness shells, Q27 divisor "
            "fibers, 13-clock debt, and Bprime owner as a Frobenius-cycle-type analogue.",
        ),
    )
    for title, body in items:
        print(f"  {title}: {body}")
    print()

    print("E. Assumption challenge")
    print("  Alternate vertex sets considered: runners, gaps, fixed circle sections, section")
    print("  boundaries, wall-crossing events, residues, cover arcs, Fourier modes, matroid")
    print("  circuits, denominator obligations, added-speed cover masks, valuation primes,")
    print("  and proof tricks.  For this transfer the selected quotient is proof tricks plus")
    print("  denominator obligations.  It preserves exact local-to-global certificate logic")
    print("  and finite exception structure, but destroys continuous timing geometry and")
    print("  arbitrary outside-window interactions.  Challenged assumption: a tournament")
    print("  vertex need not be a runner; the proof-bearing vertex can be a surviving split")
    print("  obligation, exactly as polynomial irreducibility is about hidden factor grids.")


def main() -> None:
    print("=" * 78)
    print("Codex irreducibility tricks -> LRC14 proof-transfer atlas")
    print("=" * 78)
    print()
    print_trick_transfer_table()
    print_exception_section()
    print_tournament_analysis()
    print_proof_program()


if __name__ == "__main__":
    main()
