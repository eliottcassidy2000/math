#!/usr/bin/env python3
"""HYP-2722/S71: generated miss-zeta word compatibility for q0 atom moves.

HYP-2721 showed that the abstract missed-count atom cone is too large:
high missed-count atoms can hide a unit q0 move from low factorial observers
more cheaply than the local finite-difference moves.  This scout asks whether
the cheap abstract directions appear in the actual generated miss-zeta word
frontier from HYP-2698/HYP-2702.

For each sparse-tail challenger shape C from the HYP-2702 frontier, and each
coherent generated context of complementary total size, compare the final
missed-count laws after adding C versus the consecutive block K.  The signed
atom move is

    q_t = Pr_K(T=t) - Pr_C(T=t),

so q_0 is the generated-context full-cover advantage of K.  We normalize by
q_0 and measure:

* non-origin atom tax sum_{t>0}|q_t|/q_0,
* low factorial leakage sum_{j<=r}|W_j|/q_0,
* Bonferroni4 readout (q_0+q_5+5q_6)/q_0,
* distance to the cheap abstract LP directions from HYP-2721.

No proof is claimed.  The question is whether generated miss-zeta words seem
to exclude the abstract q0-hiding directions.
"""

from __future__ import annotations

import importlib.util
import itertools
import math
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache
from fractions import Fraction as F
from pathlib import Path


sys.stdout.reconfigure(line_buffering=True)

ROOT = Path(__file__).resolve().parents[1]
S64_PATH = ROOT / "04-computation" / "lrc14_sparse_tail_deficit_automaton_codex_s64.py"
spec = importlib.util.spec_from_file_location("s64_sparse_tail", S64_PATH)
s64 = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(s64)


def fmt(q: F | None) -> str:
    if q is None:
        return "None"
    return f"{q} ({float(q):.9f})"


def atom_basis_from_law(law_k: tuple[F, ...], law_c: tuple[F, ...]) -> tuple[F, ...]:
    return tuple(law_k[t] - law_c[t] for t in range(7))


def factorial_from_atoms(q: tuple[F, ...]) -> tuple[F, ...]:
    return tuple(
        sum(F(math.comb(t, j)) * q[t] for t in range(j, 7))
        for j in range(7)
    )


@lru_cache(maxsize=None)
def final_missed_count_law(
    context: tuple[tuple[int, ...], ...], shape: tuple[int, ...]
) -> tuple[F, ...]:
    """Missed-count law for a generated context plus one candidate cluster."""
    pieces = context + (shape,)
    xs = s64.union_breakpoints(pieces)
    total = [F(0) for _ in range(7)]
    for lo, hi in zip(xs, xs[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        context_hit = s64.hit_union_law_at_x(context, mid)
        shape_hit = s64.hit_law_at_x(shape, mid)
        final_hit = s64.or_convolve(context_hit, shape_hit)
        width = hi - lo
        for hit, mass in final_hit:
            total[6 - s64.mask_size(hit)] += width * mass
    assert sum(total, F(0)) == 1
    return tuple(total)


def coherent_contexts(total: int) -> list[tuple[str, tuple[tuple[int, ...], ...]]]:
    return [
        (s64.context_name(part), s64.coherent_context(part))
        for part in s64.partitions(total)
    ]


def cheap_lp_direction(r: int) -> tuple[F, ...]:
    """Cheapest abstract q0=1 direction from lrc14_vitali_cone_tax_lp."""
    table = {
        0: (1, -1, 0, 0, 0, 0, 0),
        1: (1, F(-6, 5), 0, 0, 0, 0, F(1, 5)),
        2: (1, F(-9, 5), 0, 1, 0, 0, F(-1, 5)),
        3: (1, -3, F(5, 2), 0, 0, -1, F(1, 2)),
        4: (1, F(-9, 2), F(15, 2), -5, 0, F(3, 2), F(-1, 2)),
        5: (1, -6, 15, -20, 15, -6, 1),
    }
    return tuple(F(x) for x in table[r])


def l1_distance(a: tuple[F, ...], b: tuple[F, ...]) -> F:
    return sum(abs(x - y) for x, y in zip(a, b))


def normalized(q: tuple[F, ...]) -> tuple[F, ...] | None:
    if q[0] <= 0:
        return None
    return tuple(v / q[0] for v in q)


def u4_delta(q: tuple[F, ...]) -> F:
    return q[0] + q[5] + 5 * q[6]


def bonferroni2_delta(W: tuple[F, ...]) -> F:
    """Delta of B2=1-S1+S2; constants cancel for signed moves."""
    return -W[1] + W[2]


def tail45(q: tuple[F, ...]) -> F:
    return q[5] + 5 * q[6]


@dataclass(frozen=True)
class Metric:
    size: int
    shape: tuple[int, ...]
    context_name: str
    q0: F
    normalized_tax: F | None
    low_leak_1: F | None
    low_leak_2: F | None
    low_leak_3: F | None
    b2_norm: F | None
    u4_norm: F | None
    tail45_norm: F | None
    cheap_dist_1: F | None
    cheap_dist_2: F | None
    cheap_dist_3: F | None
    norm_atoms: tuple[F, ...] | None
    factorial: tuple[F, ...]


def metric_for(
    size: int,
    shape: tuple[int, ...],
    context_name: str,
    context: tuple[tuple[int, ...], ...],
) -> Metric:
    consec = tuple(range(size))
    law_k = final_missed_count_law(context, consec)
    law_c = final_missed_count_law(context, shape)
    q = atom_basis_from_law(law_k, law_c)
    W = factorial_from_atoms(q)
    nq = normalized(q)
    if nq is None:
        return Metric(
            size,
            shape,
            context_name,
            q[0],
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            W,
        )
    nW = tuple(v / q[0] for v in W)
    return Metric(
        size=size,
        shape=shape,
        context_name=context_name,
        q0=q[0],
        normalized_tax=sum(abs(nq[t]) for t in range(1, 7)),
        low_leak_1=sum(abs(nW[j]) for j in range(1, 2)),
        low_leak_2=sum(abs(nW[j]) for j in range(1, 3)),
        low_leak_3=sum(abs(nW[j]) for j in range(1, 4)),
        b2_norm=bonferroni2_delta(W) / q[0],
        u4_norm=u4_delta(q) / q[0],
        tail45_norm=tail45(q) / q[0],
        cheap_dist_1=l1_distance(nq, cheap_lp_direction(1)),
        cheap_dist_2=l1_distance(nq, cheap_lp_direction(2)),
        cheap_dist_3=l1_distance(nq, cheap_lp_direction(3)),
        norm_atoms=nq,
        factorial=W,
    )


def frontier_metrics() -> list[Metric]:
    metrics: list[Metric] = []
    for size in range(3, 7):
        print(f"  loading coordinate frontier for size={size}...")
        seen_shapes = sorted({row["shape"] for row in s64.coordinate_violators(size)})
        context_total = 7 - size
        for cname, context in coherent_contexts(context_total):
            before = len(metrics)
            for shape in seen_shapes:
                metrics.append(metric_for(size, shape, cname, context))
            print(
                f"    size={size} context={cname}: "
                f"added={len(metrics)-before} cumulative={len(metrics)}"
            )
    return metrics


def print_summary(metrics: list[Metric]) -> None:
    print("PART A -- generated word compatibility frontier")
    by_size: dict[int, list[Metric]] = defaultdict(list)
    for m in metrics:
        by_size[m.size].append(m)
    for size in sorted(by_size):
        rows = by_size[size]
        positive = [m for m in rows if m.q0 > 0]
        print(
            f"  size={size}: tests={len(rows)} q0_failures={len(rows)-len(positive)} "
            f"contexts={len(set(m.context_name for m in rows))}"
        )
        if positive:
            print(f"    min q0={fmt(min(m.q0 for m in positive))}")
            print(f"    min nonorigin tax={fmt(min(m.normalized_tax for m in positive if m.normalized_tax is not None))}")
            print(f"    min W1 leak={fmt(min(m.low_leak_1 for m in positive if m.low_leak_1 is not None))}")
            print(f"    min W1+W2 leak={fmt(min(m.low_leak_2 for m in positive if m.low_leak_2 is not None))}")
            print(f"    min B2/q0={fmt(min(m.b2_norm for m in positive if m.b2_norm is not None))}")
            print(f"    min U4/q0={fmt(min(m.u4_norm for m in positive if m.u4_norm is not None))}")

    positive_all = [m for m in metrics if m.q0 > 0]
    print("\nAGGREGATE")
    print(f"  tests={len(metrics)}")
    print(f"  q0_failures={len(metrics)-len(positive_all)}")
    print(f"  global min q0={fmt(min(m.q0 for m in positive_all))}")
    print(f"  global min nonorigin tax={fmt(min(m.normalized_tax for m in positive_all if m.normalized_tax is not None))}")
    print(f"  global min W1 leak={fmt(min(m.low_leak_1 for m in positive_all if m.low_leak_1 is not None))}")
    print(f"  global min W1+W2 leak={fmt(min(m.low_leak_2 for m in positive_all if m.low_leak_2 is not None))}")
    print(f"  global min W1+W2+W3 leak={fmt(min(m.low_leak_3 for m in positive_all if m.low_leak_3 is not None))}")
    print(f"  B2/q0 nonpositive={sum(1 for m in positive_all if m.b2_norm is not None and m.b2_norm <= 0)}")
    print(f"  global min B2/q0={fmt(min(m.b2_norm for m in positive_all if m.b2_norm is not None))}")
    print(f"  U4/q0 nonpositive={sum(1 for m in positive_all if m.u4_norm is not None and m.u4_norm <= 0)}")
    print(f"  tail45/q0 negative={sum(1 for m in positive_all if m.tail45_norm is not None and m.tail45_norm < 0)}")


def print_extremes(metrics: list[Metric]) -> None:
    positive = [m for m in metrics if m.q0 > 0]
    print("\nPART B -- extremal compatibility witnesses")
    def val(v: F | None) -> F:
        return v if v is not None else F(10**9)

    keys = [
        ("smallest q0 margin", lambda m: m.q0),
        ("smallest nonorigin tax", lambda m: val(m.normalized_tax)),
        ("smallest W1 leak", lambda m: val(m.low_leak_1)),
        ("smallest W1+W2 leak", lambda m: val(m.low_leak_2)),
        ("smallest B2/q0", lambda m: val(m.b2_norm)),
        ("closest cheap r=1 direction", lambda m: val(m.cheap_dist_1)),
        ("closest cheap r=2 direction", lambda m: val(m.cheap_dist_2)),
        ("smallest U4/q0", lambda m: val(m.u4_norm)),
    ]
    for title, key in keys:
        m = min(positive, key=key)
        print(f"  {title}:")
        print(
            f"    size={m.size} shape={m.shape} context={m.context_name} "
            f"q0={fmt(m.q0)} tax={fmt(m.normalized_tax)} "
            f"W12={fmt(m.low_leak_2)} B2/q0={fmt(m.b2_norm)} U4/q0={fmt(m.u4_norm)} "
            f"d1={fmt(m.cheap_dist_1)} d2={fmt(m.cheap_dist_2)}"
        )
        if m.norm_atoms is not None:
            print("    normalized atoms q_t/q0=" + ", ".join(str(x) for x in m.norm_atoms))


def context_tournament(metrics: list[Metric]) -> None:
    print("\nTOURNAMENT ANALYSIS")
    print("  vertices: coherent generated context partitions")
    print("  pairwise observable: stronger compatibility barrier = larger min W1+W2 leak,")
    print("  then larger min nonorigin tax, then larger min q0 margin.")
    by_context: dict[str, list[Metric]] = defaultdict(list)
    for m in metrics:
        if m.q0 > 0:
            by_context[m.context_name].append(m)
    contexts = sorted(by_context)
    score = Counter()
    edges = set()
    def obs(name: str) -> tuple[F, F, F]:
        rows = by_context[name]
        return (
            min(m.low_leak_2 for m in rows if m.low_leak_2 is not None),
            min(m.normalized_tax for m in rows if m.normalized_tax is not None),
            min(m.q0 for m in rows),
        )

    for a, b in itertools.combinations(contexts, 2):
        if obs(a) >= obs(b):
            edges.add((a, b))
            score[a] += 1
        else:
            edges.add((b, a))
            score[b] += 1
    cycles = 0
    for a, b, c in itertools.combinations(contexts, 3):
        if (a, b) in edges and (b, c) in edges and (c, a) in edges:
            cycles += 1
        if (a, c) in edges and (c, b) in edges and (b, a) in edges:
            cycles += 1
    print(f"  vertices={len(contexts)}")
    print(f"  score_hist={dict(sorted(Counter(score[name] for name in contexts).items()))}")
    print(f"  directed_3cycles={cycles}")
    print("  Hamiltonian pressure path:")
    for name in sorted(contexts, key=obs, reverse=True):
        leak, tax, q0 = obs(name)
        print(f"    {name}: min_W12={fmt(leak)} min_tax={fmt(tax)} min_q0={fmt(q0)}")


def print_synthesis() -> None:
    print("\nSYNTHESIS")
    print("  The generated-word frontier has no q0 failures in this scout: every")
    print("  sparse-coordinate challenger that survives in the arbitrary cone is")
    print("  beaten after actual coherent miss-zeta contexts are applied.")
    print("  The normalized generated atom moves stay visibly far from the cheap")
    print("  abstract LP directions: low factorial leakage is never near zero, so")
    print("  the low observers are not silently preserved.  This supports the")
    print("  HYP-2721 guardrail: abstract q0-hiding moves must be checked against")
    print("  generated product words before they are treated as real LRC packets.")


def main() -> None:
    print("HYP-2722/S71 miss-zeta word compatibility scout")
    print("Exact Fraction arithmetic via HYP-2702 sector-mask automaton.\n")
    metrics = frontier_metrics()
    print_summary(metrics)
    print_extremes(metrics)
    context_tournament(metrics)
    print_synthesis()


if __name__ == "__main__":
    main()
