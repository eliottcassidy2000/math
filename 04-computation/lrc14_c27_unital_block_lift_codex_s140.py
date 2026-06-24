#!/usr/bin/env python3
"""S140: lift marked C27 transfers into q=3 unital 4-point blocks.

The q=3 Hermitian unital has parameters 2-(28,4,1).  The current LRC14
prompt asks whether HYP-2937's marked C=27 shell transfers can be lifted into
these 4-point blocks after AP/Goddyn-Wong labels are attached.

This script builds the Hermitian unital over GF(9), attaches the 28 labels

    AP, GW, H_1..H_13, D_1..D_13,

and calibrates one unital block to be

    {AP, GW, H12, D3},

so the Goddyn-Wong transfer H12 -> D3 becomes literally the AP/GW anchor block.
The remaining labels are noncanonical, by design: S107 already warned that
there is no natural S8-invariant AP8 pair-slot model.  What is canonical is the
pair-unique completion rule: every marked transfer pair lies in exactly one
4-point block.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations, permutations, product


Q = 3
FIELD_SIZE = 9


@dataclass(frozen=True, order=True)
class F9:
    """GF(9)=F3[w]/(w^2+1), stored as a+b*w."""

    a: int
    b: int = 0

    def __post_init__(self) -> None:
        object.__setattr__(self, "a", self.a % 3)
        object.__setattr__(self, "b", self.b % 3)

    def __add__(self, other: "F9") -> "F9":
        return F9(self.a + other.a, self.b + other.b)

    def __neg__(self) -> "F9":
        return F9(-self.a, -self.b)

    def __sub__(self, other: "F9") -> "F9":
        return self + (-other)

    def __mul__(self, other: "F9") -> "F9":
        # w^2 = -1 = 2 in F3.
        return F9(self.a * other.a + 2 * self.b * other.b, self.a * other.b + self.b * other.a)

    def __pow__(self, n: int) -> "F9":
        out = F9(1, 0)
        base = self
        while n:
            if n & 1:
                out *= base
            base *= base
            n >>= 1
        return out

    def inv(self) -> "F9":
        if self == F9(0, 0):
            raise ZeroDivisionError
        return self ** 7

    def __truediv__(self, other: "F9") -> "F9":
        return self * other.inv()

    def code(self) -> int:
        return self.a + 3 * self.b

    def short(self) -> str:
        if self.b == 0:
            return str(self.a)
        if self.a == 0:
            return "w" if self.b == 1 else "2w"
        return f"{self.a}+{self.b}w"


ZERO = F9(0, 0)
ONE = F9(1, 0)
ELEMENTS = tuple(F9(a, b) for a in range(3) for b in range(3))
NONZERO = tuple(x for x in ELEMENTS if x != ZERO)


Point = tuple[F9, F9, F9]


def point_key(point: Point) -> tuple[int, int, int]:
    return tuple(x.code() for x in point)


def canonical_projective(point: Point) -> Point:
    if all(x == ZERO for x in point):
        raise ValueError("zero projective point")
    return min((tuple(lam * x for x in point) for lam in NONZERO), key=point_key)


def hermitian_value(point: Point) -> F9:
    # For q=3, x^(q+1)=x^4 lies in F3.
    return point[0] ** 4 + point[1] ** 4 + point[2] ** 4


def dot(line: Point, point: Point) -> F9:
    return line[0] * point[0] + line[1] * point[1] + line[2] * point[2]


def build_unital() -> tuple[list[Point], list[frozenset[int]], dict[frozenset[int], int], Counter[int]]:
    projective = sorted(
        {canonical_projective(p) for p in product(ELEMENTS, repeat=3) if any(x != ZERO for x in p)},
        key=point_key,
    )
    points = [p for p in projective if hermitian_value(p) == ZERO]
    point_index = {p: i for i, p in enumerate(points)}

    blocks: set[frozenset[int]] = set()
    line_intersections: Counter[int] = Counter()
    for line in projective:
        hit = frozenset(point_index[p] for p in points if dot(line, p) == ZERO)
        line_intersections[len(hit)] += 1
        if len(hit) == 4:
            blocks.add(hit)

    block_list = sorted(blocks, key=lambda b: tuple(sorted(b)))
    pair_to_block: dict[frozenset[int], int] = {}
    for bi, block in enumerate(block_list):
        for pair in combinations(block, 2):
            key = frozenset(pair)
            if key in pair_to_block:
                raise AssertionError("pair occurs in two blocks")
            pair_to_block[key] = bi
    return points, block_list, pair_to_block, line_intersections


def all_labels() -> list[str]:
    labels = ["AP", "GW"]
    labels += [f"H{a}" for a in range(1, 14)]
    labels += [f"D{a}" for a in range(1, 14)]
    return labels


@dataclass(frozen=True)
class Labelling:
    point_label: dict[int, str]
    label_point: dict[str, int]
    anchor_block: int
    petal10_block: int
    near_block: int


def block_labels(blocks: list[frozenset[int]], labelling: Labelling, block_id: int) -> tuple[str, ...]:
    return tuple(sorted(labelling.point_label[p] for p in blocks[block_id]))


def choose_labelling(blocks: list[frozenset[int]], pair_to_block: dict[frozenset[int], int]) -> Labelling:
    """Pick a deterministic AP/GW-calibrated labelling with a useful splice chain."""

    anchor_id = 0
    anchor = set(blocks[anchor_id])
    all_points = set().union(*blocks)
    outside = sorted(all_points - anchor)

    for ap, gw, h12, d3 in permutations(sorted(anchor), 4):
        assigned = {ap: "AP", gw: "GW", h12: "H12", d3: "D3"}
        used = set(assigned)

        for h10, d7 in permutations(outside, 2):
            petal_block = pair_to_block[frozenset((h10, d7))]
            if blocks[petal_block] & anchor:
                continue
            for d9 in outside:
                if d9 in {h10, d7}:
                    continue
                near_block = pair_to_block[frozenset((h12, d9))]
                # The near block must be the 12-branch: it meets the AP/GW
                # anchor exactly at H12, and meets the unit petal block once.
                if blocks[near_block] & anchor != {h12}:
                    continue
                if len(blocks[near_block] & blocks[petal_block]) != 1:
                    continue

                for h13, d1 in permutations(outside, 2):
                    if len({h10, d7, d9, h13, d1}) < 5:
                        continue
                    petal13_block = pair_to_block[frozenset((h13, d1))]
                    if blocks[petal13_block] & anchor:
                        continue
                    if petal13_block in {petal_block, near_block}:
                        continue

                    assigned = {
                        ap: "AP",
                        gw: "GW",
                        h12: "H12",
                        d3: "D3",
                        h10: "H10",
                        d7: "D7",
                        d9: "D9",
                        h13: "H13",
                        d1: "D1",
                    }
                    remaining_labels = [
                        label
                        for label in all_labels()
                        if label not in set(assigned.values())
                    ]
                    remaining_points = [p for p in sorted(all_points) if p not in assigned]
                    point_label = dict(assigned)
                    for p, label in zip(remaining_points, remaining_labels):
                        point_label[p] = label
                    label_point = {label: point for point, label in point_label.items()}
                    return Labelling(point_label, label_point, anchor_id, petal_block, near_block)

    raise RuntimeError("no calibrated labelling found")


@dataclass(frozen=True)
class Transfer:
    name: str
    components: tuple[tuple[str, str], ...]
    branch: str
    M: str


TRANSFERS = (
    Transfer("GW 12->24", (("H12", "D3"),), "tight-floor", "1/14"),
    Transfer("near-miss 12->36", (("H12", "D9"),), "K33-unit-excess", "3/41"),
    Transfer("petal 10->20", (("H10", "D7"),), "C27-petal", "2/27"),
    Transfer("petal 13->26", (("H13", "D1"),), "C27-petal", "2/27"),
    Transfer("splice 10,12 -> 20,24", (("H10", "D7"), ("H12", "D3")), "C27 two-hole splice", "2/27"),
    Transfer("splice 10,12 -> 20,36", (("H10", "D7"), ("H12", "D9")), "C27/K33 two-hole splice", "2/27"),
)


def component_block(
    labelling: Labelling,
    pair_to_block: dict[frozenset[int], int],
    component: tuple[str, str],
) -> int:
    a, b = component
    return pair_to_block[frozenset((labelling.label_point[a], labelling.label_point[b]))]


def transfer_summary(
    blocks: list[frozenset[int]],
    pair_to_block: dict[frozenset[int], int],
    labelling: Labelling,
    transfer: Transfer,
) -> dict[str, object]:
    block_ids = tuple(component_block(labelling, pair_to_block, c) for c in transfer.components)
    union = set().union(*(blocks[b] for b in block_ids))
    intersections = []
    for a, b in combinations(block_ids, 2):
        intersections.append(tuple(sorted(labelling.point_label[p] for p in (blocks[a] & blocks[b]))))
    return {
        "name": transfer.name,
        "components": transfer.components,
        "branch": transfer.branch,
        "M": transfer.M,
        "blocks": block_ids,
        "block_labels": tuple(block_labels(blocks, labelling, b) for b in block_ids),
        "union_size": len(union),
        "component_intersections": tuple(intersections),
    }


def relation_score(summary: dict[str, object]) -> tuple[int, int, int, int, int]:
    labels = set()
    for block in summary["block_labels"]:
        labels.update(block)
    has_apgw = int({"AP", "GW"}.issubset(labels))
    has_unit = int(any(label in labels for label in ("H10", "D7", "H13", "D1")))
    has_nonunit = int(any(label in labels for label in ("H12", "D3", "D9")))
    splice_depth = len(summary["blocks"])
    k33 = int("K33" in summary["branch"])
    return (has_apgw, has_unit + has_nonunit, splice_depth, k33, -summary["union_size"])


def tournament_fingerprint(scores: list[tuple[int, ...]]) -> tuple[dict[int, int], int, int]:
    n = len(scores)
    wins = [0] * n
    edges: dict[tuple[int, int], int] = {}
    for i, j in combinations(range(n), 2):
        if scores[i] >= scores[j]:
            winner = i
        else:
            winner = j
        wins[winner] += 1
        edges[(i, j)] = winner

    def beats(a: int, b: int) -> bool:
        key = (a, b) if a < b else (b, a)
        return edges[key] == a

    c3 = 0
    for a, b, c in combinations(range(n), 3):
        if beats(a, b) and beats(b, c) and beats(c, a):
            c3 += 1
        if beats(a, c) and beats(c, b) and beats(b, a):
            c3 += 1
    hp = 1 if sorted(wins) == list(range(n)) and c3 == 0 else 0
    return dict(sorted(Counter(wins).items())), c3, hp


def main() -> None:
    points, blocks, pair_to_block, line_intersections = build_unital()
    labelling = choose_labelling(blocks, pair_to_block)

    print("S140 LRC14 C27 TRANSFERS INTO q=3 UNITAL BLOCKS")
    print("=" * 78)
    print("[0] q=3 Hermitian unital construction check")
    print(f"  GF(9) model: w^2=-1 over F3; Hermitian equation x^4+y^4+z^4=0.")
    print(f"  projective unital points={len(points)}")
    print(f"  non-tangent 4-point blocks={len(blocks)}")
    print(f"  line intersection histogram={dict(sorted(line_intersections.items()))}")
    block_size_hist = Counter(len(b) for b in blocks)
    point_rep = Counter()
    for b in blocks:
        for p in b:
            point_rep[p] += 1
    pair_rep = Counter(pair_to_block.values())
    print(f"  block_size_hist={dict(sorted(block_size_hist.items()))}")
    print(f"  point_replication_hist={dict(sorted(Counter(point_rep.values()).items()))}")
    print(f"  pair_to_block entries={len(pair_to_block)} expected={28*27//2}")
    print(f"  secant pair multiplicity per block={dict(sorted(Counter(pair_rep.values()).items()))}")

    print("\n[1] Assumption challenge / quotient declaration")
    print("  considered vertices:")
    print("    runners, C27 shell events, transfer pairs, unital points, unital blocks,")
    print("    AP/GW labels, two-swap splice components, S139 affine-depth packets,")
    print("    and HYP-2908 state-lift obligations.")
    print("  chosen quotient:")
    print("    q=3 unital points labelled by AP, GW, H_a, D_a; marked transfers")
    print("    are lifted by the unique 4-point block through the pair H_a,D_b.")
    print("  preserved predicate:")
    print("    pair-unique incidence, AP/GW anchor visibility, unit/nonunit shell")
    print("    label, and transfer-component intersections.")
    print("  destroyed predicate:")
    print("    canonical S8 symmetry and exact LRC time geometry; exact M/Farey")
    print("    branch remains attached outside the unital lift.")
    print("  challenged assumption:")
    print("    the q=3 unital should be a canonical AP8 pair-slot design.  S107 says")
    print("    no; here it is used as an AP/GW-calibrated pair-completion frame.")

    print("\n[2] AP/Goddyn-Wong anchor block")
    anchor_labels = block_labels(blocks, labelling, labelling.anchor_block)
    print(f"  anchor_block_id={labelling.anchor_block}")
    print(f"  anchor labels={anchor_labels}")
    print("  calibration: block(AP,GW) = block(H12,D3) = {AP,GW,H12,D3}.")
    print("  readout: GW's marked transfer is the AP/GW-labelled unital block.")

    print("\n[3] HYP-2937 / HYP-2940 transfer lifts")
    summaries = [transfer_summary(blocks, pair_to_block, labelling, t) for t in TRANSFERS]
    for s in summaries:
        print(f"  {s['name']}: M={s['M']} branch={s['branch']}")
        for comp, bid, labels in zip(s["components"], s["blocks"], s["block_labels"]):
            print(f"    component {comp[0]}->{comp[1]} lifts to block {bid}: {labels}")
        if len(s["blocks"]) > 1:
            print(f"    component block intersections={s['component_intersections']} union_size={s['union_size']}")

    print("\n[4] Splice readout")
    petal_block = block_labels(blocks, labelling, labelling.petal10_block)
    near_block = block_labels(blocks, labelling, labelling.near_block)
    print(f"  petal10 block={petal_block}")
    print(f"  near/K33 block={near_block}")
    print("  chosen calibration makes:")
    print("    GW block and petal10 block disjoint, so the 10+12->20+24 splice")
    print("    is a clean product of two primary defects.")
    print("    near/K33 block shares H12 with the AP/GW block and one latent point")
    print("    with the petal10 block, producing a 3-block chain")
    print("      AP/GW --H12-- near/K33 --latent-- petal10.")

    print("\n[5] Tournament Analysis on lifted transfer packets")
    print("  vertices: transfer packets, not runners.")
    print("  pairwise observable:")
    print("    AP/GW anchor visibility, unit/nonunit content, component depth,")
    print("    K33 flag, and smaller unital union as tie/pressure gauge.")
    names = [s["name"] for s in summaries]
    scores = [relation_score(s) for s in summaries]
    hist, c3, hp = tournament_fingerprint(scores)
    for name, score in sorted(zip(names, scores), key=lambda x: x[1], reverse=True):
        print(f"    {name:28s} score={score}")
    print(f"  fingerprint: score_hist={hist} c3={c3} hp={hp}")
    print("  Hamiltonian path:")
    print("    " + " > ".join(name for name, _ in sorted(zip(names, scores), key=lambda x: x[1], reverse=True)))

    print("\n[6] Proof readout")
    print("  The lift makes the q=3 unital useful in exactly the conservative sense:")
    print("    it supplies pair-unique 4-point completions for marked C27 transfer")
    print("    pairs after AP/GW labels are attached.")
    print("  It does not canonically identify the unital with AP8 pair slots.")
    print("  New local target:")
    print("    show that any low-gap non-AP/GW residual must lift to either the")
    print("    AP/GW anchor block, a unit-petal block, or the three-block")
    print("    AP/GW--near/K33--petal chain.  Then route the unit block to C27")
    print("    petal/two-swap discharge and the nonunit chain to HYP-2908.")
    print("  POKE contribution:")
    print("    treat q=3 unital blocks as a pair-completion forum for transfer")
    print("    packets, not as scalar evidence or a canonical AP8 design.")


if __name__ == "__main__":
    main()
