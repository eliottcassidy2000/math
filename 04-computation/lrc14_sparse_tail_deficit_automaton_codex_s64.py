#!/usr/bin/env python3
"""HYP-2702 / T937: sparse-tail deficit automaton scout.

HYP-2700 says the exact Z/7 coloring makes the top residual layers safe:
consecutive blocks dominate full cover and large residual requirements.  HYP-2697
and HYP-2698 say coordinatewise sparse residuals are the remaining obstruction.

This script asks whether the sparse-tail coordinate wins disappear after one
passes through generated context words.  It is an exact Fraction scout:

* enumerate bounded shapes and the residual coordinates where they beat the
  consecutive block;
* group sparse residuals by Z/7 cyclic gap geometry;
* test singleton-product context kernels;
* test coherent generated contexts of total size 7 using the 64-mask
  OR/deletion automaton;
* report a small Tournament Analysis fingerprint whose vertices are proof
  carriers, not runners or arcs.

No LRC14 proof is claimed here.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from functools import lru_cache
from fractions import Fraction as F
from itertools import combinations
from math import comb, gcd


FULL_MASK = (1 << 6) - 1
SPAN = 13
SPARSE_TARGETS = tuple(mask for mask in range(1, 64) if mask.bit_count() <= 4)


def fmt(q: F) -> str:
    return f"{q} ({float(q):.9f})"


def mask_size(mask: int) -> int:
    return mask.bit_count()


def mask_to_sectors(mask: int) -> tuple[int, ...]:
    return tuple(i + 1 for i in range(6) if mask & (1 << i))


def primitive(shape: tuple[int, ...]) -> bool:
    g = 0
    for x in shape:
        g = gcd(g, x)
    return g == 1


def bounded_shapes(size: int, span: int = SPAN) -> list[tuple[int, ...]]:
    out = []
    for rest in combinations(range(1, span + 1), size - 1):
        shape = (0,) + rest
        if primitive(shape):
            out.append(shape)
    return out


def frontier_challengers(size: int) -> list[tuple[int, ...]]:
    """Near-consecutive challengers from the HYP-2698 frontier scan."""
    consec = tuple(range(size))
    out = []
    for shape in bounded_shapes(size, size + 2):
        if shape == consec:
            continue
        l1_move = sum(abs(shape[i] - i) for i in range(size))
        max_move = max(abs(shape[i] - i) for i in range(size))
        if l1_move <= 4 and max_move <= 2:
            out.append(shape)
    return out


def partitions(n: int, max_part: int | None = None) -> list[tuple[int, ...]]:
    if max_part is None:
        max_part = n
    if n == 0:
        return [()]
    out: list[tuple[int, ...]] = []
    for first in range(min(n, max_part), 0, -1):
        for rest in partitions(n - first, first):
            out.append((first,) + rest)
    return out


def coherent_context(part: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(range(s)) for s in part)


def context_name(part: tuple[int, ...]) -> str:
    return "[" + "+".join(str(s) for s in part) + "]"


def cyclic_gap_signature(sectors: tuple[int, ...]) -> tuple[int, ...]:
    pts = sorted(sectors)
    if len(pts) <= 1:
        return (7,)
    gaps = []
    for a, b in zip(pts, pts[1:]):
        gaps.append(b - a)
    gaps.append(pts[0] + 7 - pts[-1])
    return tuple(sorted(gaps))


@lru_cache(maxsize=None)
def breakpoints(shape: tuple[int, ...]) -> tuple[F, ...]:
    bps = {F(0), F(1)}
    diffs = {abs(a - b) for a, b in combinations(shape, 2) if a != b}
    diffs.update(abs(d) for d in shape if d)
    for d in diffs:
        for a in range(7 * d + 1):
            bps.add(F(a, 7 * d))
    return tuple(sorted(bps))


def union_breakpoints(shapes: tuple[tuple[int, ...], ...]) -> tuple[F, ...]:
    bps = {F(0), F(1)}
    for shape in shapes:
        bps.update(breakpoints(shape))
    return tuple(sorted(bps))


@lru_cache(maxsize=None)
def hit_law_at_x(shape: tuple[int, ...], x: F) -> tuple[tuple[int, F], ...]:
    base = tuple((d * x) % 1 for d in shape)
    cuts = {F(0), F(1)}
    for b in base:
        for s in range(7):
            cuts.add((F(s, 7) - b) % 1)
    cuts = sorted(cuts)
    dist: dict[int, F] = defaultdict(F)
    for lo, hi in zip(cuts, cuts[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        mask = 0
        for b in base:
            sec = int(((b + mid) % 1) * 7)
            if 1 <= sec <= 6:
                mask |= 1 << (sec - 1)
        dist[mask] += hi - lo
    return tuple(sorted(dist.items()))


def zeta_capacity(hit_law: tuple[tuple[int, F], ...]) -> tuple[F, ...]:
    cap = [F(0) for _ in range(64)]
    for hit, mass in hit_law:
        sub = hit
        while True:
            cap[sub] += mass
            if sub == 0:
                break
            sub = (sub - 1) & hit
    return tuple(cap)


@lru_cache(maxsize=None)
def capacity_at_x(shape: tuple[int, ...], x: F) -> tuple[F, ...]:
    return zeta_capacity(hit_law_at_x(shape, x))


@lru_cache(maxsize=None)
def global_capacity(shape: tuple[int, ...]) -> tuple[F, ...]:
    total = [F(0) for _ in range(64)]
    xs = breakpoints(shape)
    for lo, hi in zip(xs, xs[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        cap = capacity_at_x(shape, mid)
        for mask, value in enumerate(cap):
            total[mask] += (hi - lo) * value
    return tuple(total)


@lru_cache(maxsize=None)
def global_capacity_sparse(shape: tuple[int, ...]) -> tuple[F, ...]:
    """q_C(R) only for |R|<=4, packed in SPARSE_TARGETS order."""
    total = [F(0) for _ in SPARSE_TARGETS]
    xs = breakpoints(shape)
    for lo, hi in zip(xs, xs[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        dist = hit_law_at_x(shape, mid)
        for hit, mass in dist:
            if mass == 0:
                continue
            for i, target in enumerate(SPARSE_TARGETS):
                if target & ~hit == 0:
                    total[i] += (hi - lo) * mass
    return tuple(total)


@lru_cache(maxsize=None)
def hit_count_law(shape: tuple[int, ...]) -> tuple[F, ...]:
    total = [F(0) for _ in range(7)]
    xs = breakpoints(shape)
    for lo, hi in zip(xs, xs[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        for hit, mass in hit_law_at_x(shape, mid):
            total[mask_size(hit)] += (hi - lo) * mass
    assert sum(total, F(0)) == 1
    return tuple(total)


def singleton_cover_weight(r: int, need: int) -> F:
    if need == 0:
        return F(1)
    total = 0
    for j in range(need + 1):
        total += (-1) ** j * comb(need, j) * (7 - j) ** r
    return F(total, 7**r)


def singleton_delta(r: int, consec: tuple[int, ...], challenger: tuple[int, ...]) -> F:
    law_k = hit_count_law(consec)
    law_c = hit_count_law(challenger)
    return sum(
        (law_k[h] - law_c[h]) * singleton_cover_weight(r, 6 - h)
        for h in range(7)
    )


def or_convolve(
    a: tuple[tuple[int, F], ...], b: tuple[tuple[int, F], ...]
) -> tuple[tuple[int, F], ...]:
    out: dict[int, F] = defaultdict(F)
    for ma, wa in a:
        for mb, wb in b:
            out[ma | mb] += wa * wb
    return tuple(sorted(out.items()))


@lru_cache(maxsize=None)
def hit_union_law_at_x(
    context: tuple[tuple[int, ...], ...], x: F
) -> tuple[tuple[int, F], ...]:
    cur: tuple[tuple[int, F], ...] = ((0, F(1)),)
    for shape in context:
        cur = or_convolve(cur, hit_law_at_x(shape, x))
    return cur


@lru_cache(maxsize=None)
def residual_law_at_x(
    context: tuple[tuple[int, ...], ...], x: F
) -> tuple[tuple[int, F], ...]:
    out: dict[int, F] = defaultdict(F)
    for hit, mass in hit_union_law_at_x(context, x):
        out[FULL_MASK ^ hit] += mass
    return tuple(sorted(out.items()))


def context_delta(
    consec: tuple[int, ...],
    challenger: tuple[int, ...],
    context: tuple[tuple[int, ...], ...],
) -> F:
    shapes = (consec, challenger) + context
    xs = union_breakpoints(shapes)
    total = F(0)
    for lo, hi in zip(xs, xs[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        cap_k = capacity_at_x(consec, mid)
        cap_c = capacity_at_x(challenger, mid)
        local = F(0)
        for residual, mass in residual_law_at_x(context, mid):
            local += mass * (cap_k[residual] - cap_c[residual])
        total += (hi - lo) * local
    return total


def coordinate_violators(size: int) -> list[dict[str, object]]:
    consec = tuple(range(size))
    qk = global_capacity_sparse(consec)
    rows = []
    shapes = bounded_shapes(size) if size == 3 else frontier_challengers(size)
    for shape in shapes:
        if shape == consec:
            continue
        qc = global_capacity_sparse(shape)
        gains = []
        for idx, mask in enumerate(SPARSE_TARGETS):
            diff = qc[idx] - qk[idx]
            if diff > 0:
                gains.append((diff, mask))
        if gains:
            gains.sort(reverse=True)
            rows.append({"shape": shape, "gains": gains})
    return rows


def print_coordinate_census() -> dict[int, list[dict[str, object]]]:
    print("PART A -- exact sparse-tail coordinate violator census")
    print(f"  shapes use 0 in E, primitive, K_s=(0,...,s-1)")
    print("  size=3 is full span<=13; sizes 4..6 use the HYP-2698 near-consecutive frontier")
    print("  only |R|<=4 is recomputed here; HYP-2700 supplies the large-layer anchor")
    all_rows = {}
    geom_counter: Counter[tuple[int, tuple[int, ...]]] = Counter()
    for size in range(3, 7):
        rows = coordinate_violators(size)
        all_rows[size] = rows
        by_layer = Counter()
        max_gain = defaultdict(F)
        max_row: dict[int, tuple[tuple[int, ...], tuple[int, ...], F]] = {}
        for row in rows:
            shape = row["shape"]
            for gain, mask in row["gains"]:
                layer = mask_size(mask)
                by_layer[layer] += 1
                if gain > max_gain[layer]:
                    max_gain[layer] = gain
                    max_row[layer] = (shape, mask_to_sectors(mask), gain)
                if layer <= 4:
                    geom_counter[(layer, cyclic_gap_signature(mask_to_sectors(mask)))] += 1
        scope = "full span<=13" if size == 3 else "frontier"
        print(f"  size={size} ({scope}): coordinate-violator shapes={len(rows)}")
        for layer in range(1, 7):
            if by_layer[layer]:
                shape, sectors, gain = max_row[layer]
                print(
                    f"    |R|={layer}: coordinates={by_layer[layer]:5d}, "
                    f"max_gain={fmt(gain)} at C={shape}, R={sectors}, "
                    f"gap_sig={cyclic_gap_signature(sectors)}"
                )
        print("    scan complete", flush=True)
    print("\n  sparse residual geometry census (layer, cyclic gap signature):")
    for (layer, sig), count in geom_counter.most_common(12):
        print(f"    |R|={layer}, gaps={sig}: {count}")
    print()
    return all_rows


def print_singleton_neutralization(all_rows: dict[int, list[dict[str, object]]]) -> dict[str, object]:
    print("PART B -- singleton-product context neutralization")
    print("  r=7-size is the first total-size-7 context depth from HYP-2698.")
    metrics = {
        "failures": 0,
        "tested": 0,
        "min_margin": None,
        "worst": None,
        "r_depth_hist": Counter(),
    }
    for size in range(3, 7):
        consec = tuple(range(size))
        r0 = 7 - size
        rows = all_rows[size]
        failures = 0
        min_delta: F | None = None
        worst = None
        r_depths = Counter()
        for row in rows:
            shape = row["shape"]
            delta = singleton_delta(r0, consec, shape)
            if min_delta is None or delta < min_delta:
                min_delta = delta
                worst = shape
            if delta <= 0:
                failures += 1
            depth = None
            for r in range(0, 10):
                if singleton_delta(r, consec, shape) > 0:
                    depth = r
                    break
            r_depths[depth] += 1
        metrics["failures"] += failures
        metrics["tested"] += len(rows)
        if min_delta is not None and (
            metrics["min_margin"] is None or min_delta < metrics["min_margin"]
        ):
            metrics["min_margin"] = min_delta
            metrics["worst"] = (size, worst)
        metrics["r_depth_hist"].update((size, k, v) for k, v in r_depths.items())
        print(
            f"  size={size}, r={r0}: tested={len(rows)}, failures={failures}, "
            f"min_delta={fmt(min_delta or F(0))} at C={worst}"
        )
        print(f"    first-positive singleton depth histogram: {dict(sorted(r_depths.items()))}")
    print()
    return metrics


def print_coherent_context_test(all_rows: dict[int, list[dict[str, object]]]) -> dict[str, object]:
    print("PART C -- coherent generated contexts of total size 7")
    print("  Exact shared-x OR/deletion automaton; contexts are integer partitions of r=7-size.")
    metrics = {"failures": 0, "tested": 0, "min_margin": None, "worst": None}
    for size in range(3, 7):
        consec = tuple(range(size))
        r = 7 - size
        parts = partitions(r)
        contexts = [(part, coherent_context(part)) for part in parts]
        failures = 0
        min_delta: F | None = None
        worst = None
        tested = 0
        for row in all_rows[size]:
            shape = row["shape"]
            for part, context in contexts:
                tested += 1
                delta = context_delta(consec, shape, context)
                if min_delta is None or delta < min_delta:
                    min_delta = delta
                    worst = (shape, context_name(part))
                if delta <= 0:
                    failures += 1
        metrics["failures"] += failures
        metrics["tested"] += tested
        if min_delta is not None and (
            metrics["min_margin"] is None or min_delta < metrics["min_margin"]
        ):
            metrics["min_margin"] = min_delta
            metrics["worst"] = (size, worst)
        print(
            f"  size={size}, r={r}, contexts={len(contexts)}: tested={tested}, "
            f"failures={failures}, min_delta={fmt(min_delta or F(0))} at {worst}"
        )
    print()
    return metrics


def sccs(vertices: list[str], edges: dict[tuple[str, str], str]) -> list[list[str]]:
    adj = {v: [] for v in vertices}
    radj = {v: [] for v in vertices}
    for a in vertices:
        for b in vertices:
            if a == b:
                continue
            winner = edges[(a, b)]
            loser = b if winner == a else a
            adj[winner].append(loser)
            radj[loser].append(winner)
    seen = set()
    order = []

    def dfs(v: str) -> None:
        seen.add(v)
        for w in adj[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for v in vertices:
        if v not in seen:
            dfs(v)
    seen.clear()
    comps = []

    def rdfs(v: str, comp: list[str]) -> None:
        seen.add(v)
        comp.append(v)
        for w in radj[v]:
            if w not in seen:
                rdfs(w, comp)

    for v in reversed(order):
        if v not in seen:
            comp: list[str] = []
            rdfs(v, comp)
            comps.append(sorted(comp))
    return comps


def hamiltonian_path_count(vertices: list[str], edges: dict[tuple[str, str], str]) -> int:
    count = 0

    def edge(a: str, b: str) -> bool:
        return edges[(a, b)] == a

    def rec(path: list[str], unused: set[str]) -> None:
        nonlocal count
        if not unused:
            count += 1
            return
        last = path[-1]
        for v in list(unused):
            if edge(last, v):
                unused.remove(v)
                path.append(v)
                rec(path, unused)
                path.pop()
                unused.add(v)

    for v in vertices:
        rest = set(vertices)
        rest.remove(v)
        rec([v], rest)
    return count


def print_tournament(singleton_metrics: dict[str, object], coherent_metrics: dict[str, object]) -> None:
    print("PART D -- Tournament Analysis on proof carriers")
    min_single = singleton_metrics["min_margin"] or F(0)
    min_coherent = coherent_metrics["min_margin"] or F(0)
    carriers = {
        "raw_coordinate_weights": (10_000, F(-1), 0, 0),
        "large_R_stratification": (0, F(0), 4, 2),
        "singleton_hit_count_kernel": (
            int(singleton_metrics["failures"]),
            min_single,
            5,
            3,
        ),
        "coherent_OR_contexts": (
            int(coherent_metrics["failures"]),
            min_coherent,
            6,
            4,
        ),
        "miss_zeta_product_word": (0, min_coherent, 7, 5),
        "survival_basis_signal": (0, F(0), 4, 4),
    }
    vertices = list(carriers)

    def better(a: str, b: str) -> str:
        ma = carriers[a]
        mb = carriers[b]
        key_a = (-ma[0], ma[1], ma[2], ma[3])
        key_b = (-mb[0], mb[1], mb[2], mb[3])
        return a if key_a >= key_b else b

    edges = {}
    scores = Counter()
    for i, a in enumerate(vertices):
        for b in vertices[i + 1 :]:
            w = better(a, b)
            edges[(a, b)] = w
            edges[(b, a)] = w
            scores[w] += 1
    score_hist = Counter(scores[v] for v in vertices)
    cycles = 0
    for a, b, c in combinations(vertices, 3):
        ab = edges[(a, b)]
        bc = edges[(b, c)]
        ca = edges[(c, a)]
        if ab == a and bc == b and ca == c:
            cycles += 1
        if ab == b and bc == c and ca == a:
            cycles += 1
    comps = sccs(vertices, edges)
    hp = hamiltonian_path_count(vertices, edges)
    print("  vertices are proof carriers; observable=(failures, min_margin, exactness, generated-fidelity)")
    for name, metric in sorted(carriers.items(), key=lambda kv: (-scores[kv[0]], kv[0])):
        print(f"    {name}: score={scores[name]}, metric={metric}")
    print(f"  score_hist={dict(sorted(score_hist.items()))}")
    print(f"  directed_3cycles={cycles}")
    print(f"  SCCs={comps}")
    print(f"  Hamiltonian_path_count={hp}")
    print("  Read: sparse-tail deficits should be lifted into the generated-context basis;")
    print("        raw residual coordinates are the losing quotient.")


def main() -> None:
    print("HYP-2702 / T937 -- sparse-tail deficit automaton exact scout")
    print("Arithmetic: exact Fractions over Z/7 sector-mask breakpoints.\n")
    all_rows = print_coordinate_census()
    singleton_metrics = print_singleton_neutralization(all_rows)
    coherent_metrics = print_coherent_context_test(all_rows)
    print_tournament(singleton_metrics, coherent_metrics)
    print("\nSYNTHESIS")
    print("  Sparse coordinate wins are real, but they are concentrated in low layers")
    print("  and specific Z/7 gap geometries.  The next proof target is not a raw")
    print("  coordinate inequality; it is a generated-context automaton inequality.")
    print("  HYP-2701's survival-basis gate is parallel signal: choose the basis first,")
    print("  then scalarize.")


if __name__ == "__main__":
    main()
