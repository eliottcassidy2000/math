#!/usr/bin/env python3
"""
LRC(14) covering residual classifier.

This is a proof-hunt script, not a counterexample search.  It starts from the
current gap-side reductions:

  * THM-523: a counterexample must be q-covering, i.e. contain a multiple of
    every q=2..14.
  * THM-524: the exact gap M(S) is attained at a binding-pair crossing.
  * THM-526: if S=A union {14m} and the widest level-1/14 safe arc of A has
    width W(A)>1/(98m), then S is certified loose by arc width alone.

The goal is to classify where the remaining proof pressure sits.  In particular
we keep separate:

  * q-covering obligations (THM-523),
  * unit-residue coverage modulo 14 (a stronger, different condition),
  * arc-width certificates (THM-526),
  * binding-crossing slack (THM-524), and
  * private-obligation pressure of parked runners.

Tournament Analysis:
  Vertices are sampled q-covering rows.  The pairwise observable is "which row is
  harder" under the exact M-slack gauge, and the comparison gauge is whether the
  private/blocking-height proxy predicts the same orientation.  Ties are broken
  by the fixed Hamiltonian path of lexicographic row order.  This tests the
  phrase "dominance grows with blocking height" in a concrete LRC14 quotient:
  not runner dominance, but obligation-pressure dominance.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache, reduce
from math import gcd
import itertools
import random


N = 14
LEVEL = F(1, N)
QS = tuple(range(2, N + 1))
AP13 = tuple(range(1, 14))
RNG = random.Random(20260617)


def lcm(a: int, b: int) -> int:
    return a * b // gcd(a, b)


def gcd_all(values: tuple[int, ...]) -> int:
    return reduce(gcd, values, 0)


def norm1(x: F) -> F:
    r = x - int(x)
    if r < 0:
        r += 1
    return r if r <= F(1, 2) else 1 - r


def g_value(S: tuple[int, ...], tau: F) -> F:
    return min(norm1(v * tau) for v in S)


@lru_cache(maxsize=None)
def candidate_taus(S: tuple[int, ...]) -> frozenset[F]:
    S = tuple(sorted(set(S)))
    out: set[F] = {F(1, 2)}
    for v in S:
        k = 0
        while True:
            t = F(2 * k + 1, 2 * v)
            if t > F(1, 2):
                break
            out.add(t)
            k += 1
    for a, b in itertools.combinations(S, 2):
        for d in (a + b, abs(b - a)):
            if d <= 0:
                continue
            k = 1
            while True:
                t = F(k, d)
                if t > F(1, 2):
                    break
                out.add(t)
                k += 1
    return frozenset(out)


@lru_cache(maxsize=None)
def exact_M(S: tuple[int, ...]) -> tuple[F, list[F]]:
    best = F(0)
    ats: list[F] = []
    for t in candidate_taus(S):
        val = g_value(S, t)
        if val > best:
            best = val
            ats = [t]
        elif val == best:
            ats.append(t)
    return best, sorted(ats)


def binding_records(S: tuple[int, ...], tau: F) -> list[dict[str, object]]:
    val = g_value(S, tau)
    binders = [v for v in S if norm1(v * tau) == val]
    records: list[dict[str, object]] = []
    for a, b in itertools.combinations(binders, 2):
        for kind, d in (("sum", a + b), ("diff", abs(b - a))):
            if d <= 0:
                continue
            if norm1(d * tau) == 0:
                j = val * d
                records.append(
                    {
                        "pair": (a, b),
                        "kind": kind,
                        "D": d,
                        "j": j,
                        "slack": 14 * j - d,
                    }
                )
    return records


def add_wrapped_interval(out: list[tuple[F, F]], lo: F, hi: F) -> None:
    while lo < 0:
        lo += 1
        hi += 1
    while lo >= 1:
        lo -= 1
        hi -= 1
    if hi <= 1:
        out.append((lo, hi))
    else:
        out.append((lo, F(1)))
        out.append((F(0), hi - 1))


def merged(intervals: list[tuple[F, F]]) -> list[tuple[F, F]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    res: list[tuple[F, F]] = []
    lo, hi = intervals[0]
    for a, b in intervals[1:]:
        if a <= hi:
            if b > hi:
                hi = b
        else:
            res.append((lo, hi))
            lo, hi = a, b
    res.append((lo, hi))
    return res


def safe_arcs(S: tuple[int, ...], level: F = LEVEL) -> list[tuple[F, F, F]]:
    danger: list[tuple[F, F]] = []
    for v in set(S):
        hw = level / v
        for k in range(v):
            add_wrapped_interval(danger, F(k, v) - hw, F(k, v) + hw)
    dz = merged(danger)
    if not dz:
        return [(F(0), F(1), F(1))]
    arcs: list[tuple[F, F, F]] = []
    for i, (_, hi) in enumerate(dz):
        next_lo = dz[(i + 1) % len(dz)][0]
        if i == len(dz) - 1:
            next_lo += 1
        width = next_lo - hi
        if width > 0:
            arcs.append((hi % 1, next_lo % 1, width))
    return arcs


@lru_cache(maxsize=None)
def widest_safe_arc(S: tuple[int, ...]) -> F:
    arcs = safe_arcs(S)
    return max((w for _, _, w in arcs), default=F(0))


def is_q_covering(S: tuple[int, ...]) -> bool:
    return all(any(v % q == 0 for v in S) for q in QS)


def unit_residue_coverage(S: tuple[int, ...]) -> tuple[int, tuple[int, ...]]:
    units = tuple(r for r in range(1, 14) if gcd(r, 14) == 1)
    residues = {v % 14 for v in S if v % 14 != 0}
    missing = tuple(r for r in units if r not in residues)
    return len(units) - len(missing), missing


def cover_counts(S: tuple[int, ...]) -> dict[int, int]:
    return {q: sum(1 for v in S if v % q == 0) for q in QS}


def parked_profiles(S: tuple[int, ...]) -> list[dict[str, object]]:
    counts = cover_counts(S)
    out = []
    for w in S:
        if w % 14 != 0:
            continue
        private = tuple(q for q in QS if w % q == 0 and counts[q] == 1)
        covered = tuple(q for q in QS if w % q == 0)
        core = tuple(v for v in S if v != w)
        W = widest_safe_arc(core)
        m = w // 14
        margin = 98 * m * W - 1
        out.append(
            {
                "w": w,
                "m": m,
                "private": private,
                "covered": covered,
                "W": W,
                "arc_margin": margin,
            }
        )
    return sorted(out, key=lambda p: (p["arc_margin"], len(p["private"])), reverse=True)


def fragility_score(S: tuple[int, ...]) -> F:
    counts = cover_counts(S)
    return sum(F(1, counts[q]) for q in QS)


def principal_family(max_k: int = 24) -> dict[tuple[int, ...], str]:
    rows: dict[tuple[int, ...], str] = {}
    for q in range(2, 14):
        L = lcm(q, 14)
        core = tuple(v for v in AP13 if v != q)
        for k in range(1, max_k + 1):
            w = L * k
            S = tuple(sorted(set(core + (w,))))
            if len(S) == 13 and gcd_all(S) == 1 and is_q_covering(S):
                rows[S] = f"principal_drop_q{q}_k{k}"
    return rows


def two_drop_core_family(max_x: int = 60, max_m: int = 3, cap: int = 35) -> dict[tuple[int, ...], str]:
    rows: dict[tuple[int, ...], str] = {}
    candidates: list[tuple[int, tuple[int, ...], str]] = []
    for a, b in itertools.combinations(AP13, 2):
        core0 = tuple(v for v in AP13 if v not in (a, b))
        for x in range(14, max_x + 1):
            if x in core0:
                continue
            core = tuple(sorted(set(core0 + (x,))))
            if len(core) != 12:
                continue
            for m in range(1, max_m + 1):
                w = 14 * m
                S = tuple(sorted(set(core + (w,))))
                if len(S) == 13 and gcd_all(S) == 1 and is_q_covering(S):
                    label = f"two_drop_{a}_{b}_x{x}_m{m}"
                    candidates.append((max(S), S, label))
    for _, S, label in sorted(candidates)[:cap]:
        rows[S] = label
    return rows


def seeded_random_covering(count: int = 20, bound: int = 70) -> dict[tuple[int, ...], str]:
    rows: dict[tuple[int, ...], str] = {}
    attempts = 0
    while len(rows) < count and attempts < 20000:
        attempts += 1
        vals: set[int] = set()
        for q in QS:
            vals.add(q * RNG.randint(1, bound // q))
        vals.add(1)
        while len(vals) < 13:
            vals.add(RNG.randint(1, bound))
        if len(vals) > 13:
            vals = set(RNG.sample(sorted(vals), 13))
        S = tuple(sorted(vals))
        if len(S) == 13 and gcd_all(S) == 1 and is_q_covering(S):
            rows[S] = "seeded_random_covering"
    return rows
    

def classify_rows(rows: dict[tuple[int, ...], str]) -> list[dict[str, object]]:
    records = []
    for idx, (S, source) in enumerate(sorted(rows.items(), key=lambda kv: (max(kv[0]), kv[0])), 1):
        M, taus = exact_M(S)
        profiles = parked_profiles(S)
        best_arc = profiles[0]["arc_margin"] if profiles else None
        best_private = max((len(p["private"]) for p in profiles), default=0)
        bind = []
        for tau in taus:
            bind.extend(binding_records(S, tau))
        sum_bind = [b for b in bind if b["kind"] == "sum"]
        D_min = min((b["D"] for b in sum_bind), default=None)
        unit_count, unit_missing = unit_residue_coverage(S)
        records.append(
            {
                "idx": idx,
                "S": S,
                "source": source,
                "M": M,
                "taus": taus,
                "M_slack": M - LEVEL,
                "arc_margin": best_arc,
                "arc_cert": bool(best_arc is not None and best_arc > 0),
                "private14": best_private,
                "fragility": fragility_score(S),
                "min_cover_count": min(cover_counts(S).values()),
                "unit_count": unit_count,
                "unit_missing": unit_missing,
                "bindings": bind,
                "min_sum_D": D_min,
                "parked": profiles,
            }
        )
    return records


def route_summary(records: list[dict[str, object]]) -> None:
    print("=" * 78)
    print("LRC(14) covering residual classifier")
    print("=" * 78)
    print(f"sampled q-covering primitive 13-sets: {len(records)}")
    print(f"arc-width certified by THM-526: {sum(r['arc_cert'] for r in records)}")
    print(f"arc-width residuals: {sum(not r['arc_cert'] for r in records)}")
    print(f"zero LRC breaks found: {sum(r['M'] < LEVEL for r in records)}")
    print()

    print("Hardest rows by exact gap M(S) (smallest first):")
    for r in sorted(records, key=lambda x: (x["M"], x["S"]))[:16]:
        S = r["S"]
        source = r["source"]
        print(
            f"  M={r['M']} ({float(r['M']):.6f}), slack={r['M_slack']}, "
            f"arc_cert={r['arc_cert']}, arc_margin={r['arc_margin']}, "
            f"priv14={r['private14']}, unit_missing={r['unit_missing']}, "
            f"source={source}, S={S}"
        )
        sums = [b for b in r["bindings"] if b["kind"] == "sum"]
        if sums:
            b = sorted(sums, key=lambda x: (x["D"], x["pair"]))[0]
            print(f"    first sum-binding pair={b['pair']} D={b['D']} j={b['j']} 14j-D={b['slack']}")
    print()

    by_source = defaultdict(list)
    for r in records:
        by_source[str(r["source"]).split("_k")[0]].append(r)
    print("By source prefix:")
    for src, group in sorted(by_source.items()):
        best = min(group, key=lambda r: r["M"])
        print(
            f"  {src:28s} count={len(group):4d} best_M={str(best['M']):>8} "
            f"arc_cert={sum(r['arc_cert'] for r in group):4d}/{len(group):4d}"
        )


def tournament_fingerprint(records: list[dict[str, object]], limit: int = 96) -> None:
    """Compare exact hardness order with blocking-height proxy order."""
    hard = sorted(records, key=lambda r: (r["M"], r["arc_margin"] if r["arc_margin"] is not None else F(-10), r["S"]))[:limit]
    n = len(hard)
    if n < 3:
        return

    def exact_cmp(a: dict[str, object], b: dict[str, object]) -> int:
        ka = (a["M"], a["arc_margin"] if a["arc_margin"] is not None else F(-10), a["S"])
        kb = (b["M"], b["arc_margin"] if b["arc_margin"] is not None else F(-10), b["S"])
        return -1 if ka < kb else 1

    def pressure_key(r: dict[str, object]) -> tuple[object, ...]:
        # Larger fragility/private pressure is predicted to be harder; lower arc
        # margin and lower unit coverage are also pressure signals.  Lex order is
        # the fixed Hamiltonian-path tie-break.
        arc = r["arc_margin"] if r["arc_margin"] is not None else F(-10)
        return (-r["fragility"], -r["private14"], arc, r["unit_count"], r["S"])

    agreements = 0
    flips = 0
    scores = Counter()
    adj = [[False] * n for _ in range(n)]
    for i, j in itertools.combinations(range(n), 2):
        a = hard[i]
        b = hard[j]
        exact = exact_cmp(a, b)
        proxy = -1 if pressure_key(a) < pressure_key(b) else 1
        if exact == proxy:
            agreements += 1
        else:
            flips += 1
        # Tournament orientation: exact harder row points to easier row.
        if exact < 0:
            adj[i][j] = True
            scores[i] += 1
        else:
            adj[j][i] = True
            scores[j] += 1

    # Directed 3-cycles under the exact hardness tournament should be zero unless
    # ties escaped the fixed path.  We still report it as a fingerprint.
    cycles = 0
    for i, j, k in itertools.combinations(range(n), 3):
        if adj[i][j] and adj[j][k] and adj[k][i]:
            cycles += 1
        if adj[i][k] and adj[k][j] and adj[j][i]:
            cycles += 1
    score_hist = Counter(scores[i] for i in range(n))
    total_pairs = n * (n - 1) // 2

    print()
    print("Tournament Analysis: exact hardness vs blocking-height proxy")
    print(f"  vertices: {n} hardest sampled q-covering rows")
    print("  pairwise observable: exact M-slack hardness; tie path = lexicographic S")
    print("  comparison gauge: fragility=sum_q 1/count_q, private parked obligations, arc margin")
    print(f"  pressure-order agreements: {agreements}/{total_pairs} = {agreements/total_pairs:.3f}")
    print(f"  edge flips vs pressure proxy: {flips}/{total_pairs} = {flips/total_pairs:.3f}")
    print(f"  exact tournament score histogram: {dict(sorted(score_hist.items()))}")
    print(f"  directed 3-cycles in exact tournament: {cycles}")
    print("  readout: exact hardness is nearly scalar/transitive by construction; the useful")
    print("          signal is the flip rate, i.e. how much blocking pressure fails to predict M.")


def principal_tower_table(records: list[dict[str, object]]) -> None:
    print()
    print("Principal single-drop towers: current q-covering floor")
    print("  These rows are AP13 with q removed and lcm(q,14)k inserted.")
    print("  The table separates q-covering from the stronger unit-residue coverage condition.")
    print(f"  {'q':>2} {'k':>2} {'w':>4} {'M':>10} {'slack':>10} {'unit_missing':>16} {'arc?':>5}")
    best_for_q: dict[int, dict[str, object]] = {}
    for r in records:
        source = str(r["source"])
        if not source.startswith("principal_drop_q"):
            continue
        q = int(source.split("_q")[1].split("_")[0])
        if q not in best_for_q or r["M"] < best_for_q[q]["M"]:
            best_for_q[q] = r
    for q in sorted(best_for_q):
        r = best_for_q[q]
        missing = r["unit_missing"]
        w = max(set(r["S"]) - (set(AP13) - {q}))
        k = w // lcm(q, 14)
        print(f"  {q:2d} {k:2d} {w:4d} {str(r['M']):>10} {str(r['M_slack']):>10} {str(missing):>16} {str(r['arc_cert']):>5}")


def assumption_challenge(records: list[dict[str, object]]) -> None:
    print()
    print("Assumption challenge ledger")
    print("  Runner vertices preserve exact speed identity but destroy endpoint/component adjacency.")
    print("  q-obligation vertices preserve THM-523 covering necessity but destroy timing and binders.")
    print("  Parked-runner private obligations preserve the small-m resonance source but miss")
    print("    multi-parked rows where deletion leaves a covering core.")
    print("  Safe arcs/endpoints preserve constructive witnesses (THM-526) but W(A) is non-uniform.")
    print("  Binding-pair vertices preserve the exact LRC predicate (THM-524) but the crossing")
    print("    index j>=D/14 is tautological unless j is derived from non-M arithmetic.")
    residuals = [r for r in records if not r["arc_cert"]]
    no_private = [r for r in residuals if r["private14"] == 0]
    print(f"  In this sample, arc residuals with no parked private obligation: {len(no_private)}/{len(residuals)}")
    if no_private:
        print("  First no-private residuals (suggesting the next quotient must recurse beyond one parked runner):")
        for r in sorted(no_private, key=lambda x: (x["M"], x["S"]))[:5]:
            print(f"    M={r['M']} source={r['source']} S={r['S']}")


def main() -> None:
    rows: dict[tuple[int, ...], str] = {}
    for source in (principal_family(max_k=4), two_drop_core_family(), seeded_random_covering()):
        rows.update(source)

    records = classify_rows(rows)
    route_summary(records)
    principal_tower_table(records)
    tournament_fingerprint(records)
    assumption_challenge(records)

    print()
    print("Hypothesis generated:")
    print("  HYP-2579 (proposed): after THM-523/524/526, the remaining proof is a")
    print("  private-obligation recursion.  Arc-width discharges large parked runners;")
    print("  principal towers collapse to k=1; the real residual consists of rows where")
    print("  deleting any parked runner either leaves another q-covering core or gives")
    print("  a small-m binding crossing whose private q-obligations determine the flank.")
    print("  The immediate correction is that q-covering and unit-residue coverage are")
    print("  different floors: {1..12,182} is a q-covering principal row with M=14/183,")
    print("  closer to 1/14 than the older {1..11,13,84} cover-all-unit champion.")


if __name__ == "__main__":
    main()
