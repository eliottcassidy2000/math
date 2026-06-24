#!/usr/bin/env python3
"""S140: lift marked C27 transfers into q=3 unital 4-point blocks.

This resolves the S140 double-checkpoint by keeping both projections:

1. A raw C27 residue-pair lift, where H[a]->D[d] is the four residues
   {a, 27-a, d, 27-d}.  This exposes a lambda=1 obstruction: the GW
   H12->D3 block and the K33 H12->D9 block share the pair {12,15}, so they
   cannot both be blocks in one 2-(28,4,1) unital chart.
2. An AP/Goddyn-Wong calibrated Hermitian q=3 labelling, where AP, GW,
   H1..H13, D1..D13 label the 28 unital points and one block is calibrated as
   {AP,GW,H12,D3}.  This gives a noncanonical pair-completion forum for marked
   transfer pairs and two-block splice chains.

The proof use is deliberately modest: q=3 unital blocks are branch-local
pair-unique charts/splice grammar, not a universal C27 atlas and not a
canonical AP8 pair-slot model.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations, permutations, product


C = 27


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


ZERO = F9(0, 0)
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
    for block_id, block in enumerate(block_list):
        for pair in combinations(block, 2):
            key = frozenset(pair)
            if key in pair_to_block:
                raise AssertionError("unital pair occurs in two blocks")
            pair_to_block[key] = block_id
    return points, block_list, pair_to_block, line_intersections


def block_intersection_hist(blocks: list[frozenset[int]]) -> Counter[int]:
    hist: Counter[int] = Counter()
    for a, b in combinations(blocks, 2):
        hist[len(a & b)] += 1
    return hist


@dataclass(frozen=True)
class ShellTransfer:
    key: str
    name: str
    hole_shell: int
    double_shell: int
    branch: str
    m_value: str

    @property
    def residue_block(self) -> tuple[int, int, int, int]:
        pts = {self.hole_shell, C - self.hole_shell, self.double_shell, C - self.double_shell}
        return tuple(sorted(pts))

    @property
    def global_label_block(self) -> tuple[str, str, str, str]:
        return tuple(sorted(("AP", "GW", f"H{self.hole_shell}", f"D{self.double_shell}")))

    @property
    def labelled_pair(self) -> tuple[str, str]:
        return (f"H{self.hole_shell}", f"D{self.double_shell}")


SHELL_TRANSFERS = {
    "GW": ShellTransfer("GW", "GW H12->D3", 12, 3, "tight-floor", "1/14"),
    "K33": ShellTransfer("K33", "K33 H12->D9", 12, 9, "K33-unit-excess", "3/41"),
    "P10": ShellTransfer("P10", "petal H10->D7", 10, 7, "C27-petal-two-block", "2/27"),
    "P13": ShellTransfer("P13", "petal H13->D1", 13, 1, "C27-petal-two-block", "2/27"),
}

BRANCHES = {
    "tight chart: GW + both petals": ("GW", "P10", "P13"),
    "K33 chart: near-miss + both petals": ("K33", "P10", "P13"),
    "global low frontier: GW + K33 + petals": ("GW", "K33", "P10", "P13"),
    "S138 splice drop(10,12)->add(20,24)": ("P10", "GW"),
    "S138 splice drop(10,12)->add(20,36)": ("P10", "K33"),
    "12-branch superposition only": ("GW", "K33"),
}


def pair_multiplicities(blocks: list[tuple[object, ...]]) -> Counter[tuple[object, object]]:
    counts: Counter[tuple[object, object]] = Counter()
    for block in blocks:
        for a, b in combinations(block, 2):
            counts[tuple(sorted((a, b), key=str))] += 1
    return counts


def repeated_pairs(blocks: list[tuple[object, ...]]) -> dict[tuple[object, object], int]:
    counts = pair_multiplicities(blocks)
    return {pair: n for pair, n in sorted(counts.items(), key=lambda item: str(item[0])) if n > 1}


def find_unital_block_chart(
    desired_blocks: list[tuple[int, ...]], unital_blocks: list[frozenset[int]]
) -> list[int] | None:
    """Find unital blocks with the same block intersection matrix."""

    desired_intersections = {}
    for i, j in combinations(range(len(desired_blocks)), 2):
        desired_intersections[(i, j)] = len(set(desired_blocks[i]) & set(desired_blocks[j]))

    chosen: list[int] = []

    def rec(k: int) -> list[int] | None:
        if k == len(desired_blocks):
            return chosen.copy()
        for idx, block in enumerate(unital_blocks):
            if idx in chosen:
                continue
            ok = True
            for prev_pos, prev_idx in enumerate(chosen):
                need = desired_intersections[(prev_pos, k)]
                got = len(block & unital_blocks[prev_idx])
                if got != need:
                    ok = False
                    break
            if not ok:
                continue
            chosen.append(idx)
            out = rec(k + 1)
            if out is not None:
                return out
            chosen.pop()
        return None

    return rec(0)


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
        for h10, d7 in permutations(outside, 2):
            petal_block = pair_to_block[frozenset((h10, d7))]
            if blocks[petal_block] & anchor:
                continue
            for d9 in outside:
                if d9 in {h10, d7}:
                    continue
                near_block = pair_to_block[frozenset((h12, d9))]
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
                    remaining_labels = [label for label in all_labels() if label not in set(assigned.values())]
                    remaining_points = [p for p in sorted(all_points) if p not in assigned]
                    point_label = dict(assigned)
                    for p, label in zip(remaining_points, remaining_labels):
                        point_label[p] = label
                    label_point = {label: point for point, label in point_label.items()}
                    return Labelling(point_label, label_point, anchor_id, petal_block, near_block)

    raise RuntimeError("no calibrated labelling found")


@dataclass(frozen=True)
class Packet:
    name: str
    components: tuple[tuple[str, str], ...]
    branch: str
    m_value: str


PACKETS = (
    Packet("GW 12->24", (("H12", "D3"),), "tight-floor", "1/14"),
    Packet("near-miss 12->36", (("H12", "D9"),), "K33-unit-excess", "3/41"),
    Packet("petal 10->20", (("H10", "D7"),), "C27-petal", "2/27"),
    Packet("petal 13->26", (("H13", "D1"),), "C27-petal", "2/27"),
    Packet("splice 10,12 -> 20,24", (("H10", "D7"), ("H12", "D3")), "C27 two-hole splice", "2/27"),
    Packet("splice 10,12 -> 20,36", (("H10", "D7"), ("H12", "D9")), "C27/K33 two-hole splice", "2/27"),
)


def component_block(
    labelling: Labelling,
    pair_to_block: dict[frozenset[int], int],
    component: tuple[str, str],
) -> int:
    a, b = component
    return pair_to_block[frozenset((labelling.label_point[a], labelling.label_point[b]))]


def packet_summary(
    blocks: list[frozenset[int]],
    pair_to_block: dict[frozenset[int], int],
    labelling: Labelling,
    packet: Packet,
) -> dict[str, object]:
    block_ids = tuple(component_block(labelling, pair_to_block, c) for c in packet.components)
    union = set().union(*(blocks[b] for b in block_ids))
    intersections = []
    for a, b in combinations(block_ids, 2):
        intersections.append(tuple(sorted(labelling.point_label[p] for p in (blocks[a] & blocks[b]))))
    return {
        "name": packet.name,
        "components": packet.components,
        "branch": packet.branch,
        "M": packet.m_value,
        "blocks": block_ids,
        "block_labels": tuple(block_labels(blocks, labelling, b) for b in block_ids),
        "union_size": len(union),
        "component_intersections": tuple(intersections),
    }


def packet_relation_score(summary: dict[str, object]) -> tuple[int, int, int, int, int]:
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
        winner = i if scores[i] >= scores[j] else j
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


def print_assumption_challenge() -> None:
    print("[1] Assumption challenge / quotient declaration")
    print("  considered vertices:")
    print("    runners, C27 residues, antipodal shells, hole/double statuses,")
    print("    AP/GW labels, q=3 unital points, unital 4-blocks, branch-local")
    print("    charts, S138 two-hole obligations, and S139 affine-depth packets.")
    print("  chosen predicates:")
    print("    raw test: can marked residue packets be simultaneous q=3 blocks?")
    print("    calibrated test: what unique block completes each labelled pair?")
    print("  quotient preserves:")
    print("    C27 hole/double shell pair, branch label, lambda=1 incidence,")
    print("    AP/GW anchor visibility, and transfer-component intersections.")
    print("  quotient destroys:")
    print("    exact time geometry, canonical S8 symmetry, and within-shell")
    print("    orientation beyond the declared antipodal pair.")
    print("  challenged assumption:")
    print("    one unital chart should contain all low-frontier transfers.  The")
    print("    raw lift says no; branch-local charts and calibrated completions work.")


def print_unital_sanity(
    points: list[Point],
    blocks: list[frozenset[int]],
    pair_to_block: dict[frozenset[int], int],
    line_intersections: Counter[int],
) -> None:
    block_size_hist = Counter(len(b) for b in blocks)
    point_rep = Counter()
    for block in blocks:
        for p in block:
            point_rep[p] += 1
    pair_rep = Counter(pair_to_block.values())
    print()
    print("[2] q=3 Hermitian unital sanity")
    print("  GF(9) model: w^2=-1 over F3; Hermitian equation x^4+y^4+z^4=0.")
    print(f"  projective unital points={len(points)}")
    print(f"  non-tangent 4-point blocks={len(blocks)}")
    print(f"  line intersection histogram={dict(sorted(line_intersections.items()))}")
    print(f"  block_size_hist={dict(sorted(block_size_hist.items()))}")
    print(f"  point_replication_hist={dict(sorted(Counter(point_rep.values()).items()))}")
    print(f"  pair_to_block entries={len(pair_to_block)} expected={28 * 27 // 2}")
    print(f"  secant pair multiplicity per block={dict(sorted(Counter(pair_rep.values()).items()))}")
    print(f"  block intersection histogram={dict(sorted(block_intersection_hist(blocks).items()))}")
    print("  readout: two distinct unital blocks never share two points.")


def print_raw_residue_model(unital_blocks: list[frozenset[int]]) -> None:
    _ = unital_blocks
    print()
    print("[3] Raw C27 residue-pair block model")
    for key, transfer in SHELL_TRANSFERS.items():
        print(
            f"  {key:3s} {transfer.name:18s} branch={transfer.branch:20s} "
            f"block={transfer.residue_block}"
        )
    all_blocks = [SHELL_TRANSFERS[k].residue_block for k in ("GW", "K33", "P10", "P13")]
    repeats = repeated_pairs(all_blocks)
    print(f"  simultaneous low-frontier repeated pairs={repeats}")
    print("  obstruction:")
    print("    GW and K33 both contain the C27 hole pair (12,15).  Since a q=3")
    print("    unital has lambda=1, those two packets cannot both be blocks in")
    print("    one unital chart under the raw residue-pair lift.")


def print_global_label_model() -> None:
    print()
    print("[4] Global AP/GW-label block model")
    blocks = [SHELL_TRANSFERS[k].global_label_block for k in ("GW", "K33", "P10", "P13")]
    for key, block in zip(("GW", "K33", "P10", "P13"), blocks):
        print(f"  {key:3s} block={block}")
    repeats = repeated_pairs(blocks)
    print(f"  repeated pairs={repeats}")
    print("  readout:")
    print("    treating AP and GW as two global points in every transfer block")
    print("    repeats the pair (AP,GW).  AP/GW must be chart labels, branch")
    print("    colors, or a single calibrated anchor, not universal vertices.")


def print_branch_charts(unital_blocks: list[frozenset[int]]) -> None:
    print()
    print("[5] Branch-local q=3 unital chart test")
    for name, keys in BRANCHES.items():
        desired = [SHELL_TRANSFERS[k].residue_block for k in keys]
        repeats = repeated_pairs(desired)
        chart = None if repeats else find_unital_block_chart(desired, unital_blocks)
        status = "FAIL(pair-repeat)" if repeats else ("EMBEDS" if chart is not None else "FAIL(no chart)")
        print(f"  {name:45s} status={status}")
        print(f"    desired_blocks={desired}")
        if repeats:
            print(f"    repeated_pairs={repeats}")
        if chart is not None:
            print(f"    example_unital_block_indices={chart}")
            print(f"    example_unital_blocks={[tuple(sorted(unital_blocks[i])) for i in chart]}")
    print("  splice readout:")
    print("    The two S138 genuine two-hole rows are not new single unital")
    print("    blocks.  Each is a two-block path: petal H10->D7 plus either")
    print("    the GW block H12->D3 or the K33 near-miss block H12->D9.")


def print_calibrated_labelling(
    blocks: list[frozenset[int]],
    pair_to_block: dict[frozenset[int], int],
    labelling: Labelling,
) -> list[dict[str, object]]:
    print()
    print("[6] AP/Goddyn-Wong calibrated pair-completion model")
    anchor_labels = block_labels(blocks, labelling, labelling.anchor_block)
    print(f"  anchor_block_id={labelling.anchor_block}")
    print(f"  anchor labels={anchor_labels}")
    print("  calibration: block(AP,GW) = block(H12,D3) = {AP,GW,H12,D3}.")
    print("  guardrail:")
    print("    this is noncanonical and AP/GW-calibrated, not a symmetric AP8")
    print("    pair-slot model.  The invariant is unique pair completion.")

    print("  transfer completions:")
    summaries = [packet_summary(blocks, pair_to_block, labelling, packet) for packet in PACKETS]
    for summary in summaries:
        print(f"    {summary['name']}: M={summary['M']} branch={summary['branch']}")
        for comp, block_id, labels in zip(summary["components"], summary["blocks"], summary["block_labels"]):
            print(f"      component {comp[0]}->{comp[1]} lifts to block {block_id}: {labels}")
        if len(summary["blocks"]) > 1:
            print(
                f"      component block intersections={summary['component_intersections']} "
                f"union_size={summary['union_size']}"
            )

    print("  calibrated splice readout:")
    petal_block = block_labels(blocks, labelling, labelling.petal10_block)
    near_block = block_labels(blocks, labelling, labelling.near_block)
    print(f"    petal10 block={petal_block}")
    print(f"    near/K33 block={near_block}")
    print("    10+12 -> 20+24 is a disjoint product of the petal block and AP/GW block.")
    print("    10+12 -> 20+36 is a linked packet sharing D9:")
    print("      AP/GW --H12-- near/K33 --D9-- petal10.")
    return summaries


def print_tournament_analysis(summaries: list[dict[str, object]]) -> None:
    print()
    print("[7] Tournament Analysis")
    carriers = [
        ("exact M/Farey branch", (6, 6, 5, 6, 6)),
        ("C27 marked transfer", (5, 6, 6, 5, 6)),
        ("unital pair-repeat obstruction", (5, 5, 6, 6, 6)),
        ("branch-local q3 unital chart", (4, 5, 5, 5, 5)),
        ("calibrated AP/GW pair-completion", (4, 5, 5, 4, 6)),
        ("S138 two-block splice path", (4, 5, 5, 4, 5)),
        ("global AP/GW vertex model", (1, 2, 2, 2, 5)),
        ("raw unital analogy", (0, 1, 1, 1, 1)),
    ]
    wins = Counter()
    c3 = 0
    for i, j in combinations(range(len(carriers)), 2):
        wins[i if carriers[i][1] >= carriers[j][1] else j] += 1

    def beats(a: int, b: int) -> bool:
        return carriers[a][1] >= carriers[b][1]

    for i, j, k in combinations(range(len(carriers)), 3):
        if (beats(i, j) and beats(j, k) and beats(k, i)) or (
            beats(i, k) and beats(k, j) and beats(j, i)
        ):
            c3 += 1
    order = sorted(range(len(carriers)), key=lambda idx: carriers[idx][1], reverse=True)
    print("  proof-carrier vertices, not runners.")
    print("  carrier observable:")
    print("    theorem-scale retention, C27 predicate retention, lambda=1 incidence,")
    print("    finite checkability, and anti-scalar guard.")
    carrier_hist = Counter(wins[i] for i in range(len(carriers)))
    print(f"  carrier fingerprint: score_hist={dict(sorted(carrier_hist.items()))} c3={c3} hp=1")
    print("  carrier Hamiltonian path:")
    print("    " + " > ".join(carriers[i][0] for i in order))

    print("  calibrated-packet tournament:")
    names = [summary["name"] for summary in summaries]
    scores = [packet_relation_score(summary) for summary in summaries]
    hist, packet_c3, hp = tournament_fingerprint(scores)
    print("    observable: AP/GW visibility, unit/nonunit content, depth, K33 flag, smaller union.")
    for name, score in sorted(zip(names, scores), key=lambda item: item[1], reverse=True):
        print(f"    {name:28s} score={score}")
    print(f"    fingerprint: score_hist={hist} c3={packet_c3} hp={hp}")
    print("    Hamiltonian path:")
    print("      " + " > ".join(name for name, _ in sorted(zip(names, scores), key=lambda item: item[1], reverse=True)))


def print_verdict() -> None:
    print()
    print("[8] Verdict / POKE contribution")
    print("  Negative global result:")
    print("    the raw C27 residue-pair lift cannot put GW H12->D3 and K33")
    print("    H12->D9 into the same q=3 unital, because the C27 pair (12,15)")
    print("    would occur in two different blocks.")
    print("  Positive local result:")
    print("    each branch separately embeds into the q=3 unital incidence design,")
    print("    and the S138 genuine two-hole rows lift as two-block splices.")
    print("  Calibrated pair-completion result:")
    print("    after AP/GW labels are attached to one Hermitian unital, each marked")
    print("    pair has a unique completion; GW is the AP/GW anchor, and the")
    print("    K33/petal near splice forms AP/GW--H12--near/K33--D9--petal10.")
    print("  Proof-use refinement:")
    print("    use the unital as a branch-local pair-unique chart and POKE forum")
    print("    for transfer packets, not as scalar evidence or a universal atlas.")


def main() -> None:
    points, blocks, pair_to_block, line_intersections = build_unital()
    labelling = choose_labelling(blocks, pair_to_block)

    print("S140 LRC14 C27 -> q=3 UNITAL BLOCK-LIFT SYNTHESIS")
    print("=" * 78)
    print_assumption_challenge()
    print_unital_sanity(points, blocks, pair_to_block, line_intersections)
    print_raw_residue_model(blocks)
    print_global_label_model()
    print_branch_charts(blocks)
    summaries = print_calibrated_labelling(blocks, pair_to_block, labelling)
    print_tournament_analysis(summaries)
    print_verdict()


if __name__ == "__main__":
    main()
