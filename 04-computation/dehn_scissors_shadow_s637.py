#!/usr/bin/env python3
"""
dehn_scissors_shadow_s637.py

Scout for the S637 Dehn/scissors side-channel addendum.

The prompt asks to compare rational/irrational with odd/even and
addition/multiplication, while also folding in cuboids versus simplices and the
repo's recent carrier-compression work.  Incoming HYP-2211 owns the two-shadow
field-descent packet; this addendum layers in Dehn invariants and
equidecomposability.  It repeats the rationality correction as a guardrail,
then treats parity, rationality, Dehn invariants, tournament H, LRC shells, and
unit-distance spine/bulk splits as examples of the same rule:

    a scalar quotient is useful only after naming the side channel it forgets.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations
from math import acos, pi
import os


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULT_PATH = os.path.join(
    ROOT, "05-knowledge", "results", "dehn_scissors_shadow_s637.out"
)


def hamiltonian_count(edge: dict[tuple[int, int], int], n: int) -> int:
    """Count directed Hamiltonian paths in a small tournament."""
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1

    for mask in range(1 << n):
        for last in range(n):
            ways = dp[mask][last]
            if not ways:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                a, b = sorted((last, nxt))
                bit = edge[(a, b)]
                last_beats = (bit == 1) if last == a else (bit == 0)
                if last_beats:
                    dp[mask | (1 << nxt)][nxt] += ways
    return sum(dp[-1])


def tournament_fingerprint(orderings: list[list[str]]) -> tuple[list[str], dict[str, int], int, int]:
    """Majority tournament from ranked proof lenses."""
    vertices = sorted({v for order in orderings for v in order})
    positions = [{v: i for i, v in enumerate(order)} for order in orderings]
    edge: dict[tuple[int, int], int] = {}
    scores = Counter({v: 0 for v in vertices})

    for i, a in enumerate(vertices):
        for j, b in enumerate(vertices):
            if i >= j:
                continue
            votes_a = sum(pos[a] < pos[b] for pos in positions)
            if votes_a > len(orderings) // 2:
                edge[(i, j)] = 1
                scores[a] += 1
            else:
                edge[(i, j)] = 0
                scores[b] += 1

    c3 = 0
    for i, j, k in combinations(range(len(vertices)), 3):
        def beats(x: int, y: int) -> bool:
            a, b = sorted((x, y))
            bit = edge[(a, b)]
            return (bit == 1) if x == a else (bit == 0)

        if (beats(i, j) and beats(j, k) and beats(k, i)) or (
            beats(i, k) and beats(k, j) and beats(j, i)
        ):
            c3 += 1

    return vertices, dict(scores), c3, hamiltonian_count(edge, len(vertices))


def shorten(n: float) -> str:
    return f"{n:.12f}"


def main() -> None:
    lines: list[str] = []
    add = lines.append

    add("S637 Dehn / Scissors Side-Channel Addendum")
    add("================================================")
    add("")
    add("1. Corrected e and pi statement")
    add("--------------------------------")
    add(
        "Let s=e+pi and p=e*pi.  Since e and pi are transcendental, they cannot both be algebraic: if s and p were algebraic, e and pi would be the two roots of x^2 - s*x + p, hence algebraic over the algebraic numbers, contradiction."
    )
    add("")
    add("Status table:")
    rows = [
        (
            "both s and p algebraic",
            "impossible",
            "quadratic trap x^2-s*x+p has e,pi as roots",
        ),
        (
            "s rational or merely algebraic",
            "would force p transcendental",
            "otherwise both symmetric coordinates are algebraic",
        ),
        (
            "p rational or merely algebraic",
            "would force s transcendental",
            "same symmetric-coordinate trap",
        ),
        (
            "both s and p transcendental/irrational",
            "open/compatible with known theorem",
            "the theorem proves at least one, not exactly one",
        ),
        (
            "exactly one rational",
            "not known",
            "possible only as a conditional branch, not a theorem",
        ),
    ]
    add("claim/status".ljust(40) + "reading".ljust(36) + "reason")
    for claim, status, reason in rows:
        add(claim.ljust(40) + status.ljust(36) + reason)
    add("")
    add(
        "Reading: addition and multiplication are the elementary symmetric coordinates of the unordered pair {e,pi}.  The theorem is not a parity statement; it is a retained-side-channel statement about algebraic dependence."
    )
    add("")

    add("2. Rationality is not parity")
    add("----------------------------")
    add("channel".ljust(28) + "closed quotient?".ljust(38) + "failure mode / useful replacement")
    channel_rows = [
        (
            "odd/even",
            "yes: Z -> Z/2",
            "parity is a genuine two-color ring quotient",
        ),
        (
            "rational/irrational",
            "no",
            "irrationals are not closed: sqrt(2)+(-sqrt(2)) and sqrt(2)*sqrt(2)",
        ),
        (
            "algebraic/transcendental",
            "better, but subtle",
            "algebraic numbers form a field; transcendence is not closed either",
        ),
        (
            "algebraic dependence",
            "best channel here",
            "(sum,product) can trap a two-element set over a base field",
        ),
        (
            "add/multiply",
            "bridged by symmetric polynomials",
            "for positive values logs also convert product to sum, with new side data",
        ),
    ]
    for name, closed, failure in channel_rows:
        add(name.ljust(28) + closed.ljust(38) + failure)
    add("")

    add("3. Cuboids, simplices, and Dehn side channels")
    add("----------------------------------------------")
    theta = acos(1 / 3)
    add(f"regular tetrahedron dihedral angle theta/pi ~= {shorten(theta / pi)}")
    add(
        "Niven check: if theta/pi were rational and cos(theta)=1/3, Niven's rational-angle theorem would force cos(theta) in {0,+/-1/2,+/-1}; contradiction.  Thus theta is nonzero in R/(pi Q)."
    )
    add("")
    add("object".ljust(40) + "volume channel".ljust(28) + "Dehn / retained channel")
    dehn_rows = [
        (
            "cuboid",
            "ordinary volume",
            "zero Dehn invariant: right angles vanish in R/(pi Q)",
        ),
        (
            "cube triangulated into simplices",
            "same total volume",
            "simplex Dehn pieces cancel in pairs/families",
        ),
        (
            "regular tetrahedron",
            "can match cuboid volume",
            "nonzero Dehn invariant from acos(1/3)",
        ),
        (
            "tetrahedra packed in a box",
            "volume inequality only",
            "packing is not scissors congruence unless boundary Dehn data cancels",
        ),
    ]
    for obj, volume, dehn in dehn_rows:
        add(obj.ljust(40) + volume.ljust(28) + dehn)
    add("")
    add(
        "Reading: volume is the scalar quotient; Dehn data is the side channel.  A simplex can sit inside a cuboid, and a cuboid can be triangulated into simplices, but equal volume alone does not give equidecomposability."
    )
    add("")

    add("4. Repo dictionary")
    add("------------------")
    add("problem".ljust(28) + "scalar quotient".ljust(28) + "side channel that must survive")
    repo_rows = [
        (
            "e+pi / e*pi",
            "two symmetric values",
            "algebraic-dependence trap over Qbar",
        ),
        (
            "cuboid vs simplex",
            "volume",
            "Dehn invariant in R/(pi Q)",
        ),
        (
            "tournament H",
            "Hamiltonian-path count",
            "strong-component multiset, beta/c3 packets, OCF fibers",
        ),
        (
            "unit-distance n=21",
            "57 edges",
            "20-edge unit spine plus 37-edge centered-hex bulk",
        ),
        (
            "LRC pair-sum sieve",
            "C=2n-1 shell modulus",
            "unit/nonunit/lift/carry/owner labels",
        ),
        (
            "SC tournaments",
            "fixed-count sequence",
            "Aut/Anti transporter and rooted-perspective flip",
        ),
        (
            "hard sequences",
            "next term",
            "fixed/merged/nonfixed/q-shadow/skinny/transporter packet",
        ),
    ]
    for problem, scalar, side in repo_rows:
        add(problem.ljust(28) + scalar.ljust(28) + side)
    add("")

    add("5. Tournament Analysis over proof lenses")
    add("-----------------------------------------")
    vertices = [
        "algebraic-dependence channel",
        "Dehn invariant",
        "symmetric-polynomial trap",
        "sequence-shadow packet",
        "carrier compression",
        "volume-only packing",
        "parity analogy",
        "raw exact-one rational claim",
        "log add-multiply bridge",
    ]
    certainty_order = [
        "algebraic-dependence channel",
        "symmetric-polynomial trap",
        "Dehn invariant",
        "carrier compression",
        "sequence-shadow packet",
        "log add-multiply bridge",
        "parity analogy",
        "volume-only packing",
        "raw exact-one rational claim",
    ]
    scissors_order = [
        "Dehn invariant",
        "volume-only packing",
        "carrier compression",
        "algebraic-dependence channel",
        "symmetric-polynomial trap",
        "sequence-shadow packet",
        "log add-multiply bridge",
        "parity analogy",
        "raw exact-one rational claim",
    ]
    repo_order = [
        "carrier compression",
        "sequence-shadow packet",
        "Dehn invariant",
        "algebraic-dependence channel",
        "symmetric-polynomial trap",
        "parity analogy",
        "log add-multiply bridge",
        "volume-only packing",
        "raw exact-one rational claim",
    ]
    got_vertices, scores, c3, ham = tournament_fingerprint(
        [certainty_order, scissors_order, repo_order]
    )
    score_hist = Counter(scores.values())
    add(
        "Vertices are proof lenses, not numbers.  Pairwise observable: which lens preserves more obstruction data before quotienting.  Gauges: algebraic certainty, scissors retention, and repo portability."
    )
    add(f"score_hist={dict(sorted(score_hist.items()))} directed_3cycles={c3} H={ham}")
    for v in sorted(got_vertices, key=lambda x: (-scores[x], x)):
        add(f"  score={scores[v]}: {v}")
    add("")

    add("6. Assumption Challenge")
    add("-----------------------")
    add(
        "Alternate vertices considered: numbers, operations, fields, proof obligations, Dehn pieces, simplex facets, cuboid boxes, unit-distance edges, LRC shells, Fourier/log modes, and tournament carriers.  Chosen vertices: proof lenses."
    )
    add(
        "Preserved predicate: whether a scalar equality/inequality remains valid after the forgotten side channel is reattached.  Destroyed data: exact embeddings, individual point coordinates, and representative-level isomorphism data."
    )
    add(
        "Challenged assumption: rational versus irrational behaves like odd versus even.  The better analogy is not a two-color quotient but a field-of-definition side channel."
    )
    add("")

    add("7. Working hypotheses")
    add("---------------------")
    hypotheses = [
        "For unit-distance proofs, deliberately construct small witnesses where edge count is unchanged but spine traceability, bulk shell, or Dehn-like packet changes; these are the geometry analogues of e+pi/e*pi side-channel traps.",
        "For LRC, C=2n-1 is the modulus identity only after the odd-shell side channel is named; unit/nonunit/lift labels play the role that Dehn invariants play for scissors congruence.",
        "For tournament H, the forbidden values 7 and 21 should be attacked as phantom scalar volumes: no strong-component scissors packet realizes them, while nearby values realize after adding the missing component-side data.",
        "For hard sequences such as A000568/SC, the useful recursion may live in the packet of shadows around the value, not in the value itself.",
    ]
    for i, hyp in enumerate(hypotheses, 1):
        add(f"{i}. {hyp}")

    print("\n".join(lines))


if __name__ == "__main__":
    main()
