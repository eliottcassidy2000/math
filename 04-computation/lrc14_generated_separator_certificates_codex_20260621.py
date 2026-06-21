#!/usr/bin/env python3
"""HYP-2728 follow-up: small exact separators for generated atom moves.

HYP-2728 shows that the finite factorial boundary is formal, while the
generated miss-zeta frontier excludes cheap abstract q0-hiding directions.
This scout asks whether that exclusion has small exact linear certificates in
the normalized atom coordinates q_t/q_0.

The search is deliberately finite and honest: it recomputes the HYP-2722/S71
generated frontier, then brute-forces primitive integer covectors on
coordinates t=1..6 with coefficients in [-1,1].  A covector c separates a
cheap direction p_r from the generated frontier if all generated rows satisfy

    c.q_generated > c.p_r

or the same inequality after reversing c.
"""

from __future__ import annotations

import importlib.util
import math
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S71_PATH = ROOT / "04-computation" / "lrc14_miss_zeta_word_compatibility_codex_s71.py"
spec = importlib.util.spec_from_file_location("s71_words", S71_PATH)
s71 = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules[spec.name] = s71
spec.loader.exec_module(s71)


def fmt(q: F) -> str:
    return f"{q} ({float(q):.9f})"


def sign_name(q: F) -> str:
    if q > 0:
        return "+"
    if q < 0:
        return "-"
    return "0"


def factorial(q: tuple[F, ...]) -> tuple[F, ...]:
    return tuple(sum(F(math.comb(t, j)) * q[t] for t in range(j, 7)) for j in range(7))


def u4(q: tuple[F, ...]) -> F:
    return q[0] + q[5] + 5 * q[6]


def tail45(q: tuple[F, ...]) -> F:
    return q[5] + 5 * q[6]


def b2(W: tuple[F, ...]) -> F:
    return -W[1] + W[2]


def dot(c: tuple[int, ...], q: tuple[F, ...]) -> F:
    return sum(F(ci) * qi for ci, qi in zip(c, q))


def primitive(v: tuple[int, ...]) -> bool:
    g = 0
    for x in v:
        g = math.gcd(g, abs(x))
    if g != 1:
        return False
    for x in v:
        if x != 0:
            return x > 0
    return False


def coeff_vectors(bound: int) -> list[tuple[int, ...]]:
    out: list[tuple[int, ...]] = []
    vals = range(-bound, bound + 1)
    for a1 in vals:
        for a2 in vals:
            for a3 in vals:
                for a4 in vals:
                    for a5 in vals:
                        for a6 in vals:
                            v = (0, a1, a2, a3, a4, a5, a6)
                            if primitive(v):
                                out.append(v)
    return out


def feature_table(rows: list[tuple[F, ...]]) -> None:
    print("SIGN PROFILE LEDGER")
    print("  profiles use signs of W1, W2, B2, tail45 on normalized generated rows.")
    buckets: dict[tuple[str, str, str, str], list[tuple[F, ...]]] = defaultdict(list)
    for q in rows:
        W = factorial(q)
        key = (sign_name(W[1]), sign_name(W[2]), sign_name(b2(W)), sign_name(tail45(q)))
        buckets[key].append(q)
    print(f"  profile_count={len(buckets)}")
    for key, qs in sorted(buckets.items(), key=lambda item: (-len(item[1]), item[0])):
        vals = []
        for q in qs:
            W = factorial(q)
            vals.append((abs(W[1]) + abs(W[2]), u4(q), tail45(q), b2(W)))
        print(
            f"  W1,W2,B2,tail={key}: count={len(qs):3d} "
            f"min_absW12={fmt(min(v[0] for v in vals))} "
            f"min_U4={fmt(min(v[1] for v in vals))} "
            f"min_tail45={fmt(min(v[2] for v in vals))} "
            f"min_B2={fmt(min(v[3] for v in vals))}"
        )


def named_functionals(rows: list[tuple[F, ...]]) -> None:
    print()
    print("NAMED FUNCTIONAL CHECK")
    print("  Values are min over generated rows minus value on cheap direction.")
    named = {
        "tail45": (0, 0, 0, 0, 0, 1, 5),
        "W1": (0, 1, 2, 3, 4, 5, 6),
        "W2": (0, 0, 1, 3, 6, 10, 15),
        "B2=-W1+W2": (0, -1, -1, 0, 2, 5, 9),
        "q6": (0, 0, 0, 0, 0, 0, 1),
        "-q1": (0, -1, 0, 0, 0, 0, 0),
    }
    for r in range(1, 6):
        target = s71.cheap_lp_direction(r)
        print(f"  cheap r={r}:")
        for name, c in named.items():
            deltas = [dot(c, q) - dot(c, target) for q in rows]
            best = min(deltas)
            reverse = min(-d for d in deltas)
            if best > 0:
                tag = "separates"
                margin = best
            elif reverse > 0:
                tag = "separates_reversed"
                margin = reverse
            else:
                tag = "fails"
                margin = max(best, reverse)
            print(f"    {name:13s} {tag:18s} margin={fmt(margin)}")


def separator_search(rows: list[tuple[F, ...]], bound: int = 1) -> None:
    print()
    print("SMALL INTEGER SEPARATOR SEARCH")
    print(f"  Primitive covectors on q1..q6 with coefficients in [-{bound},{bound}].")
    covectors = coeff_vectors(bound)
    print(f"  covectors={len(covectors)}")
    for r in range(1, 6):
        target = s71.cheap_lp_direction(r)
        winners: list[tuple[F, int, tuple[int, ...], F, F]] = []
        for c in covectors:
            deltas = [dot(c, q) - dot(c, target) for q in rows]
            margin = min(deltas)
            c_use = c
            if margin <= 0:
                rdeltas = [-d for d in deltas]
                rmargin = min(rdeltas)
                if rmargin > margin:
                    margin = rmargin
                    c_use = tuple(-x for x in c)
            if margin > 0:
                l1 = sum(abs(x) for x in c_use)
                score = margin / l1
                spread = max(dot(c_use, q) - dot(c_use, target) for q in rows)
                winners.append((score, l1, c_use, margin, spread))
        winners.sort(key=lambda x: (-x[0], x[1], x[2]))
        print(f"  cheap r={r}: separators={len(winners)}")
        for score, l1, c, margin, spread in winners[:5]:
            print(
                f"    c={c} l1={l1:2d} margin={fmt(margin)} "
                f"score={fmt(score)} spread={fmt(spread)}"
            )


def context_minima(metrics: list[object]) -> None:
    print()
    print("CONTEXT MINIMA")
    print("  Weakest exact witnesses by generated context partition.")
    by_context: dict[str, list[object]] = defaultdict(list)
    for m in metrics:
        if m.q0 > 0 and m.norm_atoms is not None:
            by_context[m.context_name].append(m)
    for name, ms in sorted(by_context.items()):
        rows = [m.norm_atoms for m in ms]
        assert all(q is not None for q in rows)
        triples = []
        for q in rows:
            W = factorial(q)
            triples.append((abs(W[1]) + abs(W[2]), u4(q), tail45(q)))
        print(
            f"  {name:9s} count={len(ms):3d} "
            f"min_absW12={fmt(min(x[0] for x in triples))} "
            f"min_U4={fmt(min(x[1] for x in triples))} "
            f"min_tail45={fmt(min(x[2] for x in triples))}"
        )


def main() -> None:
    print("HYP-2728 generated separator certificate scout")
    print("Exact Fraction arithmetic; generated frontier imported from HYP-2722/S71.")
    print()
    metrics = s71.frontier_metrics()
    rows = [m.norm_atoms for m in metrics if m.q0 > 0 and m.norm_atoms is not None]
    print()
    print(f"generated_rows={len(rows)}")
    feature_table(rows)
    named_functionals(rows)
    separator_search(rows)
    context_minima(metrics)
    print()
    print("SYNTHESIS")
    print("  The single tail45 covector separates every cheap abstract direction:")
    print("  r=2,4,5 lie below the generated tail floor, while r=1,3 lie above the")
    print("  generated tail ceiling.  This turns generated-word exclusion into a")
    print("  bounded tail45 strip target before relation-code/Delsarte refinement.")


if __name__ == "__main__":
    main()
