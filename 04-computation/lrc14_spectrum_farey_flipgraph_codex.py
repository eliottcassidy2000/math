"""LRC14 spectrum walk: Farey node labels on tournament flip paths.

This scout tests the reframe:

    Sigma(S) = { iso(T(S,t)) : t in [0,1) }

where T(S,t) is the winding tournament on the speeds.  As t crosses a
breakpoint, one or more pairwise arcs flip, so the spectrum is a walk in the
tournament flip graph / metagraph G_n.  The LRC14 status is recorded by the
marked Farey node M(S):

    tight       -> marked at 1/14
    loose up    -> q(S) < 14, divisibility/coarser node
    loose down  -> Farey child p/q > 1/14, often det(1/14,p/q) = -1

The script keeps the class labels practical: exact graph isomorphism through
networkx with score/c3 prefilters.  It also records the raw walk statistics,
so the flip-graph statement can be checked independently of the class labels.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
from math import gcd

import networkx as nx


N = 14
THR = F(1, N)
AP = tuple(range(1, 14))


@dataclass(frozen=True)
class Row:
    label: str
    speeds: tuple[int, ...]
    note: str


def norm1(x: F) -> F:
    r = x % 1
    return min(r, 1 - r)


def candidate_taus(S: tuple[int, ...]) -> set[F]:
    """Exact LRC candidate times from runner walls and pair crossings."""
    S = tuple(sorted(set(S)))
    cands: set[F] = {F(1, 2)}
    for i, a in enumerate(S):
        k = 0
        while True:
            t = F(2 * k + 1, 2 * a)
            if t > F(1, 2):
                break
            cands.add(t)
            k += 1
        for b in S[i + 1 :]:
            for d in (a + b, b - a):
                if d <= 0:
                    continue
                k = 1
                while True:
                    t = F(k, d)
                    if t > F(1, 2):
                        break
                    cands.add(t)
                    k += 1
    return cands


def exact_M(S: tuple[int, ...]) -> tuple[F, tuple[F, ...]]:
    best = F(0)
    pts: list[F] = []
    for t in candidate_taus(S):
        val = min(norm1(s * t) for s in S)
        if val > best:
            best = val
            pts = [t]
        elif val == best:
            pts.append(t)
    return best, tuple(sorted(pts))


def q_threshold(S: tuple[int, ...]) -> int:
    q = 2
    while any(s % q == 0 for s in S):
        q += 1
    return q


def frac_part(x: F) -> F:
    return x - (x.numerator // x.denominator)


def bit_index(n: int, i: int, j: int) -> int:
    if i > j:
        i, j = j, i
    return i * n - i * (i + 1) // 2 + (j - i - 1)


def edge(mask: int, n: int, i: int, j: int) -> bool:
    if i == j:
        raise ValueError("loop")
    if i < j:
        return bool(mask & (1 << bit_index(n, i, j)))
    return not edge(mask, n, j, i)


def winding_mask(S: tuple[int, ...], t: F) -> int:
    """Bit 1 means lower-index/lower-speed vertex beats higher index."""
    speeds = tuple(sorted(S))
    n = len(speeds)
    mask = 0
    for i in range(n):
        for j in range(i + 1, n):
            rel = frac_part(F(speeds[i] - speeds[j]) * t)
            if rel == 0 or rel == F(1, 2):
                forward = True
            else:
                forward = rel < F(1, 2)
            if forward:
                mask |= 1 << bit_index(n, i, j)
    return mask


def tournament_graph(mask: int, n: int) -> nx.DiGraph:
    graph = nx.DiGraph()
    graph.add_nodes_from(range(n))
    for i in range(n):
        for j in range(i + 1, n):
            if edge(mask, n, i, j):
                graph.add_edge(i, j)
            else:
                graph.add_edge(j, i)
    return graph


def scores(mask: int, n: int) -> tuple[int, ...]:
    return tuple(sum(1 for j in range(n) if j != i and edge(mask, n, i, j)) for i in range(n))


def c3(mask: int, n: int) -> int:
    count = 0
    for a, b, c in combinations(range(n), 3):
        local = [sum(1 for y in (a, b, c) if x != y and edge(mask, n, x, y)) for x in (a, b, c)]
        if sorted(local) == [1, 1, 1]:
            count += 1
    return count


def fingerprint(mask: int, n: int) -> tuple[tuple[int, ...], int]:
    return tuple(sorted(scores(mask, n), reverse=True)), c3(mask, n)


class IsoClassifier:
    def __init__(self, n: int):
        self.n = n
        self.reps: list[tuple[int, tuple[tuple[int, ...], int], nx.DiGraph]] = []
        self.by_mask: dict[int, int] = {}

    def classify(self, mask: int) -> int:
        if mask in self.by_mask:
            return self.by_mask[mask]
        fp = fingerprint(mask, self.n)
        graph = tournament_graph(mask, self.n)
        for class_id, rep_fp, rep_graph in self.reps:
            if rep_fp == fp and nx.is_isomorphic(graph, rep_graph):
                self.by_mask[mask] = class_id
                return class_id
        class_id = len(self.reps)
        self.reps.append((class_id, fp, graph))
        self.by_mask[mask] = class_id
        return class_id


def orientation_breaks(S: tuple[int, ...]) -> dict[F, list[tuple[int, int, str]]]:
    """Breakpoints for the winding tournament on [0,1]."""
    speeds = tuple(sorted(S))
    events: dict[F, list[tuple[int, int, str]]] = defaultdict(list)
    for i, j in combinations(range(len(speeds)), 2):
        d = abs(speeds[j] - speeds[i])
        for k in range(d + 1):
            events[F(k, d)].append((i, j, "zero"))
        for k in range(d):
            events[F(2 * k + 1, 2 * d)].append((i, j, "anti"))
    return events


def spectrum_walk(S: tuple[int, ...], extra_marks: tuple[F, ...]) -> dict[str, object]:
    n = len(S)
    events = orientation_breaks(S)
    breaks = sorted(events)
    intervals: list[tuple[F, F, F, int]] = []
    for left, right in zip(breaks, breaks[1:]):
        if left == right:
            continue
        mid = (left + right) / 2
        intervals.append((left, right, right - left, winding_mask(S, mid)))

    unique_masks = {mask for _, _, _, mask in intervals}
    unique_masks.update(winding_mask(S, t) for t in extra_marks)

    classifier = IsoClassifier(n)
    for mask in sorted(unique_masks):
        classifier.classify(mask)

    class_measure: Counter[int] = Counter()
    for _, _, length, mask in intervals:
        class_measure[classifier.classify(mask)] += length

    flip_hist: Counter[int] = Counter()
    event_hist: Counter[int] = Counter()
    internal_breaks = breaks[1:-1]
    left_masks = [item[3] for item in intervals[:-1]]
    right_masks = [item[3] for item in intervals[1:]]
    for t, lm, rm in zip(internal_breaks, left_masks, right_masks):
        changed = (lm ^ rm).bit_count()
        pair_count = len({(i, j) for i, j, _ in events[t]})
        flip_hist[changed] += 1
        event_hist[pair_count] += 1

    return {
        "breaks": breaks,
        "intervals": intervals,
        "classes": classifier,
        "class_measure": class_measure,
        "flip_hist": flip_hist,
        "event_hist": event_hist,
    }


def farey_status(M: F, qdiv: int) -> str:
    if M == THR:
        return "tight: marked apex 1/14"
    if qdiv < N:
        return f"loose-up: q(S)={qdiv} coarser divisibility node"
    excess = N * M.numerator - M.denominator
    det = M.denominator - N * M.numerator
    if excess == 1:
        return f"loose-down: unit Farey child, det={det}"
    return f"loose-down: nonunit excess e={excess}, det={det}"


def packet_status(value: F) -> str:
    a, b = value.numerator, value.denominator
    graph = f"K_{{{a},{b}}}"
    state = "nonplanar/K33-minor" if min(a, b) >= 3 else "bipartite-planar"
    return f"{graph}, sum={a + b}, product={a * b}, {state}"


def child_ladder() -> None:
    print("[Farey child ladder above 1/14]")
    print("  p/q = p/(14p-1); determinant q-14p = -1")
    for p in range(1, 7):
        node = F(p, N * p - 1)
        label = "right parent" if p == 1 else "child"
        print(f"  p={p}: {str(node):<6} {label:<12} {packet_status(node)}")
    print("  First child with a K_{3,3} carrier is 3/41, because 2/27 is K_{2,27}.")
    print()


def rows() -> list[Row]:
    def repl(old: int, new: int) -> tuple[int, ...]:
        return tuple(sorted((set(AP) - {old}) | {new}))

    return [
        Row("AP", AP, "known tight"),
        Row("GW 12->24", repl(12, 24), "known tight, n*2 petal"),
        Row("11->24", repl(11, 24), "coarser q=11 test"),
        Row("12->26", repl(12, 26), "coarser q=12 residue-liar"),
        Row("13->26", repl(13, 26), "q=14 first child test"),
        Row("12->36", repl(12, 36), "near-miss 3/41"),
        Row("12->96", repl(12, 96), "large Farey descendant"),
        Row("8->16", repl(8, 16), "small n*2 petal"),
        Row("10->20", repl(10, 20), "small n*2 petal"),
    ]


def print_row(row: Row) -> None:
    M, pts = exact_M(row.speeds)
    qdiv = q_threshold(row.speeds)
    marks = tuple(sorted(set((THR, M) + pts)))
    walk = spectrum_walk(row.speeds, marks)
    classes: IsoClassifier = walk["classes"]  # type: ignore[assignment]
    class_measure: Counter[int] = walk["class_measure"]  # type: ignore[assignment]
    apex_id = classes.classify(winding_mask(row.speeds, THR))
    opt_id = classes.classify(winding_mask(row.speeds, M))
    opt_measure = class_measure.get(opt_id, F(0))
    top = sorted(class_measure.items(), key=lambda kv: (-kv[1], kv[0]))[:3]
    top_s = ", ".join(f"C{cid}:{measure}" for cid, measure in top)
    flip_hist: Counter[int] = walk["flip_hist"]  # type: ignore[assignment]
    event_hist: Counter[int] = walk["event_hist"]  # type: ignore[assignment]
    intervals = walk["intervals"]  # type: ignore[assignment]

    print(f"[{row.label}] {row.note}")
    print(f"  S={row.speeds}")
    print(f"  q(S)={qdiv}; M={M}; opt-times={pts}; {farey_status(M, qdiv)}")
    print(f"  marked packet: {packet_status(M)}")
    print(
        f"  Sigma intervals={len(intervals)} classes={len(classes.reps)} "
        f"apex-class=C{apex_id} opt-class=C{opt_id} opt-class-measure={opt_measure}"
    )
    print(f"  top class measures: {top_s}")
    print(f"  flip histogram across breakpoints: {dict(sorted(flip_hist.items()))}")
    print(f"  event multiplicity histogram: {dict(sorted(event_hist.items()))}")
    print()


def main() -> None:
    print("LRC14 SPECTRUM AS FAREY-LABELLED TOURNAMENT FLIP WALK")
    print("=" * 78)
    child_ladder()
    for row in rows():
        print_row(row)
    print("[Readout]")
    print("  Class counts here use full directed tournament isomorphism over [0,1).")
    print("  Folded, complement-merged, or endpoint-inclusive conventions may count differently.")
    print("  Tight rows keep the marked node at 1/14.")
    print("  Divisibility failures migrate upward to q(S)<14.")
    print("  Farey-neighbor failures migrate downward to p/q with q-14p=-1.")
    print("  The spectrum walk confirms the metagraph hook: phase breakpoints are arc-flip events,")
    print("  generically single flips and at symmetric phases multi-flips.")


if __name__ == "__main__":
    main()
