"""S142: affine-depth grammar for the S140 unital block chain.

HYP-2941 kept the prompt's triangular/perfect-number idea as an affine-depth
carrier.  HYP-2942 supplied the calibrated q=3 unital path

    AP/GW --H12-- near/K33 --D9-- petal10.

This script composes those transfer components as words in

    a(x)=x/2, b(x)=x+1

and records the noncommutative depth signatures.  The main signal is that the
S140 linked order has component depths [3,4,1], whose affine suffix-depth sum is
14.  The same components in other orders do not have that LRC14 depth.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import permutations
from math import gcd
from typing import Iterable


def triangular(n: int) -> int:
    return n * (n + 1) // 2


def c27_depth(residue: int) -> int:
    """Depth from the gcd stratum in C=27.

    Units are depth 0, gcd-3 residues depth 1, gcd-9 residues depth 2.
    Zero is not used by the low-frontier marked transfers, but gets depth 3.
    """

    g = gcd(residue, 27)
    if g == 1:
        return 0
    if g == 3:
        return 1
    if g == 9:
        return 2
    return 3


@dataclass(frozen=True)
class Component:
    key: str
    label: str
    hole: int
    double: int
    block: tuple[str, ...]
    branch: str
    note: str

    @property
    def depth(self) -> int:
        # One base halving plus visible C27 depth of the two endpoints.
        return 1 + c27_depth(self.hole) + c27_depth(self.double)

    @property
    def word(self) -> str:
        return "b" + "a" * self.depth

    @property
    def strata(self) -> tuple[int, int]:
        return c27_depth(self.hole), c27_depth(self.double)


COMPONENTS = {
    "GW": Component(
        key="GW",
        label="GW 12->24",
        hole=12,
        double=3,
        block=("AP", "D3", "GW", "H12"),
        branch="tight-floor",
        note="AP/GW anchor block",
    ),
    "K33": Component(
        key="K33",
        label="near-miss 12->36",
        hole=12,
        double=9,
        block=("D10", "D9", "H12", "H9"),
        branch="K33-unit-excess",
        note="near/K33 block",
    ),
    "P10": Component(
        key="P10",
        label="petal 10->20",
        hole=10,
        double=7,
        block=("D7", "D9", "H10", "H13"),
        branch="C27-petal",
        note="unit petal10 block",
    ),
    "P13": Component(
        key="P13",
        label="petal 13->26",
        hole=13,
        double=1,
        block=("D1", "D11", "D12", "H13"),
        branch="C27-petal",
        note="unit petal13 block",
    ),
}


def affine_depth_profile(depths: Iterable[int]) -> tuple[str, list[int], Fraction, int]:
    """Return word, future-halving depths per b, beta, and depth sum."""

    ds = list(depths)
    word = "".join("b" + "a" * d for d in ds)
    suffixes: list[int] = []
    running = 0
    for depth in reversed(ds):
        running += depth
        suffixes.append(running)
    suffixes.reverse()
    beta = sum(Fraction(1, 2**s) for s in suffixes)
    return word, suffixes, beta, sum(suffixes)


def sequence_report(keys: tuple[str, ...]) -> dict[str, object]:
    comps = [COMPONENTS[k] for k in keys]
    depths = [c.depth for c in comps]
    word, suffixes, beta, depth_sum = affine_depth_profile(depths)
    union = set().union(*(set(c.block) for c in comps))
    intersections = []
    for i in range(len(comps)):
        for j in range(i + 1, len(comps)):
            intersections.append(tuple(sorted(set(comps[i].block) & set(comps[j].block))))
    return {
        "keys": keys,
        "depths": depths,
        "word": word,
        "suffixes": suffixes,
        "depth_sum": depth_sum,
        "beta": beta,
        "union_size": len(union),
        "intersections": tuple(intersections),
    }


@dataclass(frozen=True)
class Carrier:
    key: str
    name: str
    score: tuple[int, int, int, int, int, int]
    role: str


CARRIERS = [
    Carrier(
        "M",
        "exact M/Farey branch",
        (6, 6, 6, 5, 5, 6),
        "keeps p/q, q, and excess e=14p-q",
    ),
    Carrier(
        "unital",
        "C27 unital calibrated chain",
        (5, 6, 6, 5, 5, 6),
        "AP/GW--H12--near/K33--D9--petal10 packet",
    ),
    Carrier(
        "aff14",
        "affine depth-14 signature",
        (5, 5, 5, 5, 5, 6),
        "noncommutative depth sum 14 for linked S140 order",
    ),
    Carrier(
        "splice",
        "two-block splice strips",
        (5, 5, 5, 5, 4, 6),
        "P10+GW and P10+K33 triangular-strip depths",
    ),
    Carrier(
        "K33",
        "Kpq/K33 owner packet",
        (5, 4, 4, 4, 6, 5),
        "p>=3 incidence and forbidden-H state-lift address",
    ),
    Carrier(
        "tri",
        "triangular/perfect product lane",
        (4, 4, 4, 5, 3, 5),
        "T_{2x-1} and p*(14p-1) product/coimage side channel",
    ),
    Carrier(
        "raw",
        "raw scalar equality",
        (1, 1, 1, 2, 1, 1),
        "false equality if used without labels",
    ),
]


def tournament_edges(carriers: list[Carrier]) -> dict[tuple[str, str], str]:
    edges: dict[tuple[str, str], str] = {}
    for i, left in enumerate(carriers):
        for right in carriers[i + 1 :]:
            winner = left.key if left.score >= right.score else right.key
            edges[(left.key, right.key)] = winner
    return edges


def count_directed_triangles(carriers: list[Carrier], edges: dict[tuple[str, str], str]) -> int:
    def beats(a: str, b: str) -> bool:
        key = (a, b) if (a, b) in edges else (b, a)
        return edges[key] == a

    total = 0
    keys = [c.key for c in carriers]
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            for k in range(j + 1, len(keys)):
                a, b, c = keys[i], keys[j], keys[k]
                if (beats(a, b) and beats(b, c) and beats(c, a)) or (
                    beats(a, c) and beats(c, b) and beats(b, a)
                ):
                    total += 1
    return total


def score_hist(carriers: list[Carrier], edges: dict[tuple[str, str], str]) -> dict[int, int]:
    scores = {c.key: 0 for c in carriers}
    for winner in edges.values():
        scores[winner] += 1
    hist: dict[int, int] = {}
    for score in scores.values():
        hist[score] = hist.get(score, 0) + 1
    return dict(sorted(hist.items()))


def main() -> None:
    print("S142 LRC14 AFFINE-DEPTH UNITAL CHAIN")
    print("=" * 78)

    print("[0] Assumption challenge / quotient declaration")
    print("  considered vertices:")
    print("    runners, C27 residues, unital blocks, calibrated transfer pairs,")
    print("    affine words, suffix-depth profiles, triangular strips, K33 flags,")
    print("    and proof obligations.")
    print("  chosen quotient:")
    print("    calibrated S140 transfer components composed as affine a/b words.")
    print("  quotient preserves:")
    print("    C27 gcd depth, unital block intersections, branch order, and")
    print("    exact M/Farey labels through external metadata.")
    print("  quotient destroys:")
    print("    exact safe-time geometry and any canonical AP8 symmetry.")
    print("  challenged assumption:")
    print("    the triangular/perfect-number motif is only scalar numerology.")
    print("    Here it becomes a noncommutative depth-order certificate.")

    print("\n[1] Component depths from C27 strata")
    print("  depth = 1 + depth_gcd(hole) + depth_gcd(double),")
    print("  where unit,gcd3,gcd9 strata have depths 0,1,2.")
    for comp in COMPONENTS.values():
        print(
            "  {key:<3} {label:<18} strata={strata} depth={depth} "
            "block={block} note={note}".format(
                key=comp.key,
                label=comp.label,
                strata=comp.strata,
                depth=comp.depth,
                block=comp.block,
                note=comp.note,
            )
        )

    print("\n[2] S138 two-block splice affine strips")
    for keys in [("P10", "GW"), ("P10", "K33"), ("GW", "P10"), ("K33", "P10")]:
        rep = sequence_report(keys)
        print(
            "  {keys}: component_depths={depths} word={word} "
            "suffix_depths={suffixes} depth_sum={depth_sum} beta={beta} "
            "union_size={union_size} intersections={intersections}".format(**rep)
        )
    print("  readout:")
    print("    P10 then GW gives suffix depths [4,3] and sum 7.")
    print("    P10 then K33 gives [5,4] and sum 9, one depth level higher.")
    print("    Reversing the order changes the invariant.  The affine packet is")
    print("    noncommutative, as required by the S140 branch-local rule.")

    print("\n[3] S140 linked unital chain")
    linked = sequence_report(("GW", "K33", "P10"))
    print(
        "  linked order GW -> K33 -> P10: component_depths={depths} "
        "word={word} suffix_depths={suffixes} depth_sum={depth_sum} beta={beta} "
        "union_size={union_size} intersections={intersections}".format(**linked)
    )
    print("  all permutations of the same depth multiset [3,4,1]:")
    seen: set[tuple[int, ...]] = set()
    for depths in sorted(set(permutations([3, 4, 1]))):
        word, suffixes, beta, depth_sum = affine_depth_profile(depths)
        mark = "<-- S140 linked order" if depths == (3, 4, 1) else ""
        seen.add(depths)
        line = f"    depths={depths} suffix={suffixes} sum={depth_sum:2d} beta={beta}"
        if mark:
            line += f" {mark}"
        print(line)
    print("  readout:")
    print("    The S140 path order gives suffix-depth sum 14 exactly.")
    print("    The other five orders give 13,15,17,18,19.  Thus LRC14 appears")
    print("    as an order-sensitive affine-depth signature, not as a raw count.")
    print(f"    Triangular comparison: T_4={triangular(4)}, T_5={triangular(5)};")
    print("    depth 14 is the interior strip between nearby triangular totals,")
    print("    selected by the AP/GW--H12--K33--D9--petal chain.")

    print("\n[4] Tournament Analysis")
    print("  vertices: proof carriers, not runners.")
    print("  pairwise observable:")
    print("    (M/Farey retention, unital incidence retention, affine order retention,")
    print("     finite checkability, state-lift fit, anti-scalar guard).")
    print("  switch/gauge: lexicographically larger role score wins.")
    edges = tournament_edges(CARRIERS)
    for carrier in CARRIERS:
        print(f"    {carrier.name:<34} score={carrier.score} role={carrier.role}")
    print(
        "  fingerprint: score_hist={} c3={} hp=1".format(
            score_hist(CARRIERS, edges), count_directed_triangles(CARRIERS, edges)
        )
    )
    print("  Hamiltonian carrier order:")
    print("    " + " > ".join(c.name for c in CARRIERS))

    print("\n[5] Proof readout")
    print("  HYP-2941 said to attach affine-depth packets after exact labels.")
    print("  HYP-2942 supplies the calibrated labels.  Under the component-depth")
    print("  rule above, the linked K33 path has affine suffix-depth sum 14.")
    print("  Candidate lemma refinement:")
    print("    In the second-gap frontier, any residual that reaches the calibrated")
    print("    AP/GW--H12--near/K33--D9--petal10 path should carry the depth-14")
    print("    affine signature.  Unit-only splices stay in lower triangular strips")
    print("    and should discharge by C27 petal/two-swap rigidity; depth-14 linked")
    print("    packets should feed the HYP-2908 / THM-572 forbidden-H state lift.")
    print("  LRC14 is not proved here.  The contribution is a concrete packet")
    print("  grammar that turns the repeated triangular/perfect-number prompt into")
    print("  an order-sensitive invariant on the known low-frontier branches.")


if __name__ == "__main__":
    main()
