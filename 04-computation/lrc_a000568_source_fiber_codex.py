#!/usr/bin/env python3
"""Expose the A000568 source fiber hidden in LRC clocks.

The earlier S509/S512/S535 computations showed that raw unlabelled
tournament classes are the wrong final invariant for LRC: good and bad
states mix inside the same A000568 bucket.  This follow-up isolates the
clean positive statement.

For a time state, keep the moving-runner half-turn tournament, but attach the
stationary observer by a threshold edge:

    observer -> runner  iff  dist(runner, observer) >= 1/N.

Then the LRC witness predicate is exactly "the observer is a source."  At such
a state, deleting the observer leaves an ordinary unlabelled tournament class
on the moving runners.  Conversely, adding a new source to any unlabelled
tournament gives a rooted observer-source class.  Thus the LRC target slice is
a source-cone copy of the A000568(m) base, where m=N-1 is the number of moving
runners.  The arithmetic clock only visits a restricted subimage of that copy.

Tournament Analysis declaration:
    primary vertices: quotient layers, not necessarily runners.
    pairwise observable: profile tuple
        (source purity, A000568 proximity, n14 transfer, non-tautology,
         compression, low label cost).
    switch/gauge: lexicographic profile dominance.
    tie Hamiltonian path: layer declaration order.

Alternate LRC vertex sets explicitly considered:
    runners, gaps, fixed circle sections, section boundaries, wall-crossing
    events, residues, cover arcs, Fourier modes, matroid circuits, and proof
    obligations.  The selected quotient preserves the observer-source
    predicate and the deleted moving-runner class; it destroys exact gap
    lengths, labels, owner intervals, and multiplicative cover loads.  The
    challenged assumption is that vertices must be runners or that raw
    unmarked A000568 classes should already decide loneliness.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations
from math import gcd, log2


ONE = Fraction(1, 1)
HALF = Fraction(1, 2)

Adj = tuple[tuple[int, ...], ...]
ClassKey = tuple[tuple[int, ...], tuple[int, ...]]

A000568 = {
    1: 1,
    2: 1,
    3: 2,
    4: 4,
    5: 12,
    6: 56,
    7: 456,
    8: 6880,
    9: 191536,
    10: 9733056,
}


@dataclass(frozen=True)
class Family:
    name: str
    moving: tuple[int, ...]


@dataclass(frozen=True)
class State:
    speeds: tuple[int, ...]
    t: Fraction
    good: bool
    pos: tuple[Fraction, ...]
    runner_adj: Adj
    source_adj: Adj


@dataclass(frozen=True)
class Layer:
    name: str
    source_purity: int
    a000568_proximity: int
    n14_transfer: int
    non_tautology: int
    compression: int
    low_label_cost: int

    @property
    def profile(self) -> tuple[int, ...]:
        return (
            self.source_purity,
            self.a000568_proximity,
            self.n14_transfer,
            self.non_tautology,
            self.compression,
            self.low_label_cost,
        )


def mod1(x: Fraction) -> Fraction:
    return x % ONE


def dist0(x: Fraction) -> Fraction:
    x = mod1(x)
    return min(x, ONE - x)


def fmt_frac(x: Fraction) -> str:
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


def positions(speeds: tuple[int, ...], t: Fraction) -> tuple[Fraction, ...]:
    return tuple(mod1(Fraction(v) * t) for v in speeds)


def primitive(moving: tuple[int, ...]) -> bool:
    g = 0
    for v in moving:
        g = gcd(g, v)
    return g == 1


def tournament_from_winner(n: int, winner) -> Adj:
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        w = winner(i, j)
        if w == i:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return tuple(tuple(row) for row in adj)


def phase_adj(pos: tuple[Fraction, ...], vertices: tuple[int, ...]) -> Adj:
    def winner(i: int, j: int) -> int:
        a = vertices[i]
        b = vertices[j]
        delta = mod1(pos[b] - pos[a])
        if delta == 0 or delta == HALF:
            return i
        return i if delta < HALF else j

    return tournament_from_winner(len(vertices), winner)


def source_lift_adj(pos: tuple[Fraction, ...]) -> Adj:
    n = len(pos)
    threshold = Fraction(1, n)
    runner = phase_adj(pos, tuple(range(1, n)))
    adj = [[0] * n for _ in range(n)]
    for j in range(1, n):
        if dist0(pos[j] - pos[0]) >= threshold:
            adj[0][j] = 1
        else:
            adj[j][0] = 1
    for a in range(1, n):
        for b in range(1, n):
            if a != b:
                adj[a][b] = runner[a - 1][b - 1]
    return tuple(tuple(row) for row in adj)


@lru_cache(maxsize=None)
def perms(n: int) -> tuple[tuple[int, ...], ...]:
    return tuple(permutations(range(n)))


@lru_cache(maxsize=None)
def canonical_colored(adj: Adj, colors: tuple[int, ...] | None = None) -> ClassKey:
    n = len(adj)
    if colors is None:
        colors = (0,) * n
    best: ClassKey | None = None
    for p in perms(n):
        pc = tuple(colors[p[i]] for i in range(n))
        bits = tuple(adj[p[i]][p[j]] for i in range(n) for j in range(i + 1, n))
        key = (pc, bits)
        if best is None or key < best:
            best = key
    assert best is not None
    return best


def delete_root(adj: Adj) -> Adj:
    return tuple(tuple(adj[i][j] for j in range(1, len(adj))) for i in range(1, len(adj)))


def source_cone(adj: Adj) -> Adj:
    m = len(adj)
    out = [[0] * (m + 1) for _ in range(m + 1)]
    for j in range(1, m + 1):
        out[0][j] = 1
    for i in range(m):
        for j in range(m):
            if i != j:
                out[i + 1][j + 1] = adj[i][j]
    return tuple(tuple(row) for row in out)


def all_tournaments(n: int) -> list[Adj]:
    pairs = tuple(combinations(range(n), 2))
    out: list[Adj] = []
    for mask in range(1 << len(pairs)):
        adj = [[0] * n for _ in range(n)]
        for bit, (i, j) in enumerate(pairs):
            if (mask >> bit) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        out.append(tuple(tuple(row) for row in adj))
    return out


def class_representatives(n: int) -> dict[ClassKey, Adj]:
    reps: dict[ClassKey, Adj] = {}
    for adj in all_tournaments(n):
        reps.setdefault(canonical_colored(adj), adj)
    return reps


def source_cone_audit(max_m: int = 6) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for m in range(1, max_m + 1):
        reps = class_representatives(m)
        cones: dict[ClassKey, ClassKey] = {}
        deletion_failures = 0
        for base_key, adj in reps.items():
            cone = source_cone(adj)
            cone_key = canonical_colored(cone, (1,) + (0,) * m)
            cones[cone_key] = base_key
            deleted = canonical_colored(delete_root(cone))
            if deleted != base_key:
                deletion_failures += 1
        rows.append(
            {
                "m": m,
                "base_classes": len(reps),
                "A000568": A000568.get(m),
                "source_cone_rooted_classes": len(cones),
                "collisions": len(reps) - len(cones),
                "deletion_failures": deletion_failures,
            }
        )
    return rows


def event_times(speeds: tuple[int, ...]) -> list[Fraction]:
    n = len(speeds)
    threshold = Fraction(1, n)
    out: set[Fraction] = {Fraction(0), Fraction(1)}

    for v in speeds[1:]:
        for k in range(v):
            out.add(mod1((Fraction(k) + threshold) / v))
            out.add(mod1((Fraction(k) - threshold) / v))

    for i, j in combinations(range(1, n), 2):
        d = abs(speeds[i] - speeds[j])
        if d == 0:
            continue
        for k in range(2 * d + 1):
            out.add(Fraction(k, 2 * d))

    return sorted(x for x in out if 0 <= x <= 1)


def sampled_states(family: Family) -> list[State]:
    assert primitive(family.moving)
    speeds = (0,) + family.moving
    walls = event_times(speeds)
    times: set[Fraction] = set(walls)
    for left, right in zip(walls, walls[1:]):
        if left < right:
            times.add((left + right) / 2)

    out: list[State] = []
    threshold = Fraction(1, len(speeds))
    for t in sorted(times):
        pos = positions(speeds, t)
        good = all(dist0(pos[i] - pos[0]) >= threshold for i in range(1, len(pos)))
        runner = phase_adj(pos, tuple(range(1, len(pos))))
        out.append(
            State(
                speeds=speeds,
                t=t,
                good=good,
                pos=pos,
                runner_adj=runner,
                source_adj=source_lift_adj(pos),
            )
        )
    return out


def mixed_count(class_to_goodness: dict[ClassKey, Counter[bool]]) -> int:
    return sum(1 for c in class_to_goodness.values() if c[True] and c[False])


def hamiltonian_paths(adj: Adj) -> int:
    n = len(adj)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp[mask][last]
            if not count:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += count
    return sum(dp[full])


def directed_triangles(adj: Adj) -> int:
    total = 0
    for a, b, c in combinations(range(len(adj)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            total += 1
    return total


def scc_sizes(adj: Adj) -> tuple[int, ...]:
    n = len(adj)
    reach = [[bool(adj[i][j]) for j in range(n)] for i in range(n)]
    for i in range(n):
        reach[i][i] = True
    for k in range(n):
        for i in range(n):
            if reach[i][k]:
                for j in range(n):
                    reach[i][j] = reach[i][j] or reach[k][j]
    seen: set[int] = set()
    sizes = []
    for i in range(n):
        if i in seen:
            continue
        comp = {j for j in range(n) if reach[i][j] and reach[j][i]}
        seen |= comp
        sizes.append(len(comp))
    return tuple(sorted(sizes, reverse=True))


def score_hist(adj: Adj) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(Counter(sum(row) for row in adj).items()))


def analyze_states(label: str, states: list[State]) -> dict[str, object]:
    n = len(states[0].speeds)
    m = n - 1
    runner_classes: dict[ClassKey, Counter[bool]] = defaultdict(Counter)
    rooted_classes: dict[ClassKey, Counter[bool]] = defaultdict(Counter)
    unrooted_source_classes: dict[ClassKey, Counter[bool]] = defaultdict(Counter)
    good_runner_classes: set[ClassKey] = set()
    good_rooted_classes: set[ClassKey] = set()
    source_cone_mismatches = 0
    deletion_mismatches = 0
    mixed_example: tuple[Fraction, Fraction] | None = None
    mixed_witness: dict[ClassKey, dict[bool, Fraction]] = defaultdict(dict)

    for state in states:
        runner_key = canonical_colored(state.runner_adj)
        rooted_key = canonical_colored(state.source_adj, (1,) + (0,) * m)
        unrooted_key = canonical_colored(state.source_adj)
        runner_classes[runner_key][state.good] += 1
        rooted_classes[rooted_key][state.good] += 1
        unrooted_source_classes[unrooted_key][state.good] += 1
        mixed_witness[runner_key].setdefault(state.good, state.t)

        if state.good:
            good_runner_classes.add(runner_key)
            good_rooted_classes.add(rooted_key)
            cone_key = canonical_colored(source_cone(state.runner_adj), (1,) + (0,) * m)
            if cone_key != rooted_key:
                source_cone_mismatches += 1
            deleted = canonical_colored(delete_root(state.source_adj))
            if deleted != runner_key:
                deletion_mismatches += 1

    for data in mixed_witness.values():
        if True in data and False in data:
            mixed_example = (data[True], data[False])
            break

    source_cone_image = {
        canonical_colored(source_cone(rep), (1,) + (0,) * m)
        for rep_key, rep in class_representatives_for_keys(good_runner_classes, m).items()
    }

    return {
        "label": label,
        "N": n,
        "moving": m,
        "states": len(states),
        "good_states": sum(1 for s in states if s.good),
        "runner_classes": len(runner_classes),
        "runner_good_classes": len(good_runner_classes),
        "runner_mixed": mixed_count(runner_classes),
        "rooted_source_classes": len(rooted_classes),
        "rooted_good_classes": len(good_rooted_classes),
        "rooted_mixed": mixed_count(rooted_classes),
        "unrooted_source_mixed": mixed_count(unrooted_source_classes),
        "source_cone_image": len(source_cone_image),
        "source_cone_exact": good_rooted_classes == source_cone_image,
        "source_cone_mismatches": source_cone_mismatches,
        "deletion_mismatches": deletion_mismatches,
        "a000568_m": A000568.get(m),
        "source_codim_bits": (
            log2(A000568[m] / len(good_runner_classes))
            if m in A000568 and good_runner_classes
            else None
        ),
        "mixed_example": mixed_example,
    }


def class_representatives_for_keys(keys: set[ClassKey], n: int) -> dict[ClassKey, Adj]:
    reps: dict[ClassKey, Adj] = {}
    if not keys:
        return reps
    for adj in all_tournaments(n):
        key = canonical_colored(adj)
        if key in keys and key not in reps:
            reps[key] = adj
            if len(reps) == len(keys):
                break
    return reps


def primitive_speed_sets(total_n: int, max_speed: int) -> list[tuple[int, ...]]:
    out: list[tuple[int, ...]] = []
    for moving in combinations(range(1, max_speed + 1), total_n - 1):
        if primitive(moving):
            out.append(moving)
    return out


def batch_audit(total_n: int, max_speed: int) -> dict[str, object]:
    rows = primitive_speed_sets(total_n, max_speed)
    all_states: list[State] = []
    sets_with_good = 0
    for moving in rows:
        states = sampled_states(Family(str(moving), moving))
        all_states.extend(states)
        if any(s.good for s in states):
            sets_with_good += 1
    summary = analyze_states(f"batch N={total_n} max={max_speed}", all_states)
    summary["speed_sets"] = len(rows)
    summary["sets_with_good"] = sets_with_good
    return summary


def unit_witness(row: tuple[int, ...], qmax: int = 80) -> tuple[int, int] | None:
    n = len(row) + 1
    threshold = Fraction(1, n)
    for q in range(2, qmax + 1):
        for a in range(1, q):
            if gcd(a, q) != 1:
                continue
            if all(dist0(Fraction(a * v, q)) >= threshold for v in row):
                return q, a
    return None


def row_fingerprint(row: tuple[int, ...], t: Fraction) -> dict[str, object]:
    speeds = (0,) + row
    pos = positions(speeds, t)
    threshold = Fraction(1, len(speeds))
    good = all(dist0(pos[i] - pos[0]) >= threshold for i in range(1, len(pos)))
    runner = phase_adj(pos, tuple(range(1, len(pos))))
    source = source_lift_adj(pos)
    dists = [dist0(pos[i] - pos[0]) for i in range(1, len(pos))]
    return {
        "good": good,
        "root_outdegree": sum(source[0]),
        "min_dist": min(dists),
        "tight": sum(1 for d in dists if d == threshold),
        "score_hist": score_hist(runner),
        "c3": directed_triangles(runner),
        "scc": scc_sizes(runner),
        "H": hamiltonian_paths(runner),
    }


def meta_tournament(layers: list[Layer]) -> Adj:
    n = len(layers)
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        a = layers[i].profile
        b = layers[j].profile
        if a > b or (a == b and i < j):
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return tuple(tuple(row) for row in adj)


def print_source_cone_audit() -> None:
    print("A. SOURCE-CONE COPY OF A000568")
    print("=" * 88)
    print("m  A000568  enumerated  rooted_source_classes  collisions  deletion_failures")
    for row in source_cone_audit():
        print(
            "{m:>1} {A000568:>8} {base_classes:>11} {source_cone_rooted_classes:>22} "
            "{collisions:>10} {deletion_failures:>18}".format(**row)
        )
    print()
    print(
        "Conclusion: adding a distinguished source is faithful, and deleting that "
        "source returns the base A000568 class."
    )
    print()


def print_state_summary(row: dict[str, object]) -> None:
    codim = row["source_codim_bits"]
    codim_text = "-" if codim is None else f"{float(codim):.3f}"
    example = row["mixed_example"]
    if example is None:
        example_text = "-"
    else:
        example_text = f"good_t={fmt_frac(example[0])}, bad_t={fmt_frac(example[1])}"
    print(
        "{:<22s} N={:>2} states={:>6} good={:>5} "
        "runner_cls={:>4} good_runner={:>4} mixed_runner={:>3} "
        "rooted_good={:>4} rooted_mixed={:>2} cone_exact={} "
        "codim_bits={} mixed_example={}".format(
            str(row["label"])[:22],
            row["N"],
            row["states"],
            row["good_states"],
            row["runner_classes"],
            row["runner_good_classes"],
            row["runner_mixed"],
            row["rooted_good_classes"],
            row["rooted_mixed"],
            row["source_cone_exact"],
            codim_text,
            example_text,
        )
    )


def print_lrc_audits() -> None:
    print("B. EXACT LRC CLOCK AUDITS")
    print("=" * 88)
    families = [
        Family("N4 consecutive", (1, 2, 3)),
        Family("N4 sparse", (1, 3, 4)),
        Family("N5 prime-ish", (2, 3, 5, 7)),
        Family("N6 mixed", (1, 3, 4, 7, 9)),
        Family("N7 mixed", (1, 4, 6, 9, 10, 15)),
    ]
    for fam in families:
        print_state_summary(analyze_states(fam.name, sampled_states(fam)))
    print()
    print("Batch scans over primitive speed sets:")
    for total_n, max_speed in ((4, 8), (5, 7), (6, 7), (7, 7)):
        print_state_summary(batch_audit(total_n, max_speed))
    print()
    print(
        "Read: runner classes alone still mix good and bad states, but the "
        "observer-source lift has zero rooted mixed fibers in these audits. "
        "At every good state, rooted class = source_cone(deleted runner class)."
    )
    print()


def print_lrc14_snapshots() -> None:
    print("C. LRC14 HARD-ROW SNAPSHOTS")
    print("=" * 88)
    rows = [
        ("AP13", tuple(range(1, 14))),
        ("one-stranger-611", tuple(sorted((7, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77, 84, 611)))),
        ("HYP2470-E1-q33", (7, 14, 21, 35, 49, 63, 70, 77, 91, 322, 350, 504, 936)),
        ("HYP2470-E2-q31", (7, 14, 21, 28, 35, 49, 63, 77, 91, 119, 700, 1008, 1066)),
        ("THM497-band2-refuter", (39, 56, 108, 128, 164, 168, 200, 266, 323, 341, 518, 594, 598)),
    ]
    for label, row in rows:
        witness = unit_witness(row, 80)
        if witness is None:
            print(f"{label:<22s} no plain unit witness through q=80")
            continue
        q, a = witness
        fp = row_fingerprint(row, Fraction(a, q))
        print(
            "{:<22s} first_witness={:>2}/{:<2} source={} root_outdeg={:>2} "
            "min_dist={} tight={:>2} H={} c3={} scc={} score_hist={}".format(
                label,
                a,
                q,
                fp["good"],
                fp["root_outdegree"],
                fmt_frac(fp["min_dist"]),
                fp["tight"],
                fp["H"],
                fp["c3"],
                fp["scc"],
                fp["score_hist"],
            )
        )
    print()
    print(
        "For N=14, full canonical A000568(13) enumeration is deliberately not "
        "attempted here.  The witness state still lands in the same formal "
        "source-cone slice: the observer has outdegree 13 and deleting it "
        "leaves a 13-vertex runner tournament fingerprint."
    )
    print()


def print_meta_tournament() -> None:
    print("D. META-TOURNAMENT OVER QUOTIENT LAYERS")
    print("=" * 88)
    layers = [
        Layer("raw_runner_A000568", 1, 5, 2, 5, 5, 5),
        Layer("observer_marked_phase", 2, 4, 2, 4, 4, 4),
        Layer("observer_source_lift", 5, 4, 4, 3, 3, 3),
        Layer("source_deleted_fiber", 5, 5, 4, 4, 5, 5),
        Layer("gap_threshold_stack", 5, 3, 4, 4, 3, 2),
        Layer("band_cover_unit_fiber", 4, 2, 5, 5, 2, 3),
        Layer("blocking_height_dominance", 3, 2, 5, 4, 4, 4),
    ]
    adj = meta_tournament(layers)
    scores = [sum(row) for row in adj]
    print("vertices: quotient layers")
    print("observable: (purity, A000568 proximity, n14 transfer, non-tautology, compression, low label cost)")
    print("switch: lexicographic dominance; tie path is declaration order")
    print(f"score_hist={tuple(sorted(Counter(scores).items()))}")
    print(f"directed_3cycles={directed_triangles(adj)} scc={scc_sizes(adj)} H={hamiltonian_paths(adj)}")
    for layer, score in sorted(zip(layers, scores), key=lambda x: (-x[1], x[0].name)):
        print(f"  score={score} {layer.name:<28s} profile={layer.profile}")
    print()
    print(
        "Interpretation: source_deleted_fiber wins because it is both close to "
        "A000568 and exactly source-pure.  The band-cover and blocking-height "
        "layers are the LRC14 transfer side, not replacements for the source slice."
    )
    print()


def main() -> None:
    print("LRC A000568 source-fiber audit -- codex-2026-06-13")
    print("=" * 88)
    print(
        "Goal: make precise where unlabelled tournament isomorphism classes are "
        "hiding in LRC.  The answer is not raw phase class; it is the deleted "
        "class below an observer-source cone."
    )
    print()
    print_source_cone_audit()
    print_lrc_audits()
    print_lrc14_snapshots()
    print_meta_tournament()
    print("SYNTHESIS")
    print("=" * 88)
    print("1. A000568 appears exactly as the base of the observer-source cone.")
    print("2. LRC asks whether the arithmetic clock reaches that source-cone slice.")
    print("3. Raw moving-runner classes mix good and bad states, so projection is lossy.")
    print("4. THM-497 says scalar cover counts cannot prove LRC14; the missing data is")
    print("   alignment/fiber structure.  The source-deleted A000568 fiber is one such")
    print("   structural coordinate, while Q27/Q31 and blocking height are the LRC14")
    print("   multiplicative-cover coordinates.")
    print("5. New proof target: show that any long LRC14 blocked-band walk either enters")
    print("   a source-cone deleted class, or its failure to do so forces the balanced")
    print("   cover congruences feeding the Q31/7-ideal/13-clock portal.")


if __name__ == "__main__":
    main()
