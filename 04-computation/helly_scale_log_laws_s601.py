#!/usr/bin/env python3
"""S601: Helly scale as the logarithm of certificate search entropy.

S599 extracted minimal empty subfamilies from the bounded two-block determinant
automaton.  S600 proposed that this certificate size can be a source of log,
loglog, or logloglog factors.  This script makes that currency explicit.

For M live component languages and certificate rank H, the search budget is

    B_H(M) = sum_{h<=H} binom(M, h)
    Lambda_H(M) = log B_H(M).

If H is bounded, Lambda_H(M) ~= H log M.  Thus the same local Helly theorem
has ordinary-log, loglog, or logloglog cost depending only on whether the live
component count M is polynomial, logarithmic, or doubly-logarithmic in the
ambient parameter.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import comb, exp, log

import lrc_twoblock_helly_s599 as s599


EMPTY_CLASSES = {
    "singleton_empty": 1,
    "pair_empty": 2,
    "triple_empty": 3,
    "quad_empty": 4,
}


def log_budget(m: int, h: int) -> float:
    if h <= 0 or m <= 0:
        return 0.0
    h = min(h, m)
    return log(sum(comb(m, r) for r in range(1, h + 1)))


def fmt(x: float) -> str:
    return f"{x:8.3f}"


@dataclass
class HellyStats:
    rows: int = 0
    preempt: int = 0
    high: int = 0
    live: int = 0
    rank_counts: Counter[int] | None = None
    component_sum: int = 0
    component_max: int = 0
    lambda_sum: float = 0.0
    lambda_max: float = 0.0
    wlog_sum: float = 0.0
    examples: dict[int, dict[str, object]] | None = None

    def __post_init__(self) -> None:
        if self.rank_counts is None:
            self.rank_counts = Counter()
        if self.examples is None:
            self.examples = {}

    @property
    def certified(self) -> int:
        return sum(self.rank_counts.values())

    @property
    def nonpreempt(self) -> int:
        return self.certified + self.high + self.live

    def add(self, info: dict[str, object], n: int) -> None:
        cls = str(info["class"])
        self.rows += 1
        h = EMPTY_CLASSES.get(cls)
        if h is None:
            if cls == "high_order_empty":
                self.high += 1
            elif cls == "bounded_live":
                self.live += 1
            else:
                self.preempt += 1
            return

        comps = info["components"]
        m = len(comps)
        lam = log_budget(m, h)
        self.rank_counts[h] += 1
        self.component_sum += m
        self.component_max = max(self.component_max, m)
        self.lambda_sum += lam
        self.lambda_max = max(self.lambda_max, lam)
        self.wlog_sum += log(max(int(info["w_bound"]), 1))
        self.examples.setdefault(h, info)

    def avg_m(self) -> float:
        return self.component_sum / self.certified if self.certified else 0.0

    def avg_lambda(self) -> float:
        return self.lambda_sum / self.certified if self.certified else 0.0

    def avg_wlog(self) -> float:
        return self.wlog_sum / self.certified if self.certified else 0.0


def collect(n: int, regime: str, trials: int, seed: int) -> HellyStats:
    rng = s599.random.Random(seed + 1009 * n + (17 if regime == "full_stack" else 0))
    out = HellyStats()
    for _ in range(trials):
        row = s599.row_from_sample(n, rng)
        if row is None:
            continue
        out.add(s599.classify_row(row, n, regime), n)
    return out


def print_budget_atlas() -> None:
    print("CERTIFICATE ENTROPY ATLAS")
    print("  Lambda_H(M) = log sum_{h<=H} binom(M,h)")
    print("  Fixed H gives Lambda_H(M) ~ H log M.")
    print()
    print(" M       H=1      H=2      H=3      H=ceil(log M)")
    for m in (4, 8, 16, 32, 64, 128):
        hlog = max(1, int(round(log(m))))
        print(
            f"{m:<5d} {fmt(log_budget(m,1))} {fmt(log_budget(m,2))} "
            f"{fmt(log_budget(m,3))} {fmt(log_budget(m,hlog))}"
        )
    print()
    print("AMBIENT SCALE TRANSLATION")
    print(" N            log N     loglog N   logloglog N   Lambda_2(log N)   Lambda_2(loglog N)")
    for n in (10**6, 10**12, 10**60, 10**600):
        l1 = log(n)
        l2 = log(l1)
        l3 = log(l2)
        m1 = max(2, int(round(l1)))
        m2 = max(2, int(round(l2)))
        print(
            f"1e{len(str(n))-1:<9d} {fmt(l1)} {fmt(l2)} {fmt(l3)} "
            f"{fmt(log_budget(m1,2))} {fmt(log_budget(m2,2))}"
        )
    print()


def print_sample_audit() -> None:
    print("S599/S601 TWO-BLOCK HELLY-SCALE AUDIT")
    print("  Trials are deterministic samples of Cprime rows with one multiple of n.")
    print("  lambda = Lambda_h(M) for the extracted minimal empty h-subfamily.")
    print()
    print("regime       n  rows preempt  h1   h2   h3   h4 high live  avgM maxM  avgLam  maxLam avglogW")
    combined = Counter()
    for regime in ("BC_only", "full_stack"):
        for n in (6, 8, 10, 12, 14):
            st = collect(n, regime, trials=1800, seed=601)
            combined.update(
                {
                    "h1": st.rank_counts[1],
                    "h2": st.rank_counts[2],
                    "h3": st.rank_counts[3],
                    "h4": st.rank_counts[4],
                    "high": st.high,
                    "live": st.live,
                    "preempt": st.preempt,
                }
            )
            print(
                f"{regime:<11} {n:2d} {st.rows:5d} {st.preempt:7d} "
                f"{st.rank_counts[1]:3d} {st.rank_counts[2]:4d} {st.rank_counts[3]:4d} "
                f"{st.rank_counts[4]:4d} {st.high:4d} {st.live:4d} "
                f"{st.avg_m():5.2f} {st.component_max:4d} {st.avg_lambda():7.3f} "
                f"{st.lambda_max:7.3f} {st.avg_wlog():7.3f}"
            )
        print()

    print("AGGREGATE HELLY RANKS")
    total_empty = sum(combined[f"h{i}"] for i in range(1, 5))
    total_seen = total_empty + combined["high"] + combined["live"]
    print(f"  certified_empty={total_empty} high_order_empty={combined['high']} bounded_live={combined['live']}")
    print(f"  nonpreempted_seen={total_seen} preempted={combined['preempt']}")
    if total_empty:
        for h in range(1, 5):
            print(f"  h={h}: {combined[f'h{h}']} ({combined[f'h{h}']/total_empty:.3f} of extracted certificates)")
    print()


def print_laws() -> None:
    print("HELLY-LOG LAWS")
    print("1. Certificate-entropy law:")
    print("   Lambda_H(M)=log sum_{h<=H} C(M,h).  For fixed H,")
    print("   Lambda_H(M)=H log M + O_H(1).")
    print("2. Ambient translation:")
    print("   If M(N)~N^a, fixed-H Helly costs ordinary log N.")
    print("   If M(N)~log N, fixed-H Helly costs loglog N.")
    print("   If M(N)~loglog N, fixed-H Helly costs logloglog N.")
    print("3. Tax/dividend split:")
    print("   A Helly tax is a union bound over B_H(M) possible witnesses.")
    print("   A Helly dividend is an early empty subfamily that avoids the global CRT modulus.")
    print("4. S599/S601 reading:")
    print("   The sampled two-block determinant residual is almost all h=1, with a few")
    print("   h=2 rows at n=6/8.  Thus the current n=14 determinant branch has not")
    print("   started paying a large Helly logarithm; it is earning bounded local")
    print("   certificates before the global CRT scale appears.")
    print()
    print("HARMONIC HELLY LADDER")
    print("   If rank-h certificates remove alpha/h of the remaining residual,")
    print("   R_H <= R_1 exp(-alpha sum_{h<=H} 1/h) <= R_1 H^{-alpha}.")
    print("   If H itself is log N, this is the Helly-scale source of loglog N.")
    print()


def tournament_fingerprint() -> None:
    vertices = [
        ("bounded_helly_certificate", 6, 1),
        ("component_count_log", 5, 2),
        ("certificate_entropy_log", 4, 3),
        ("bounded_w_window_log", 3, 4),
        ("prime_power_CRT_modulus", 2, 5),
        ("ambient_denominator_height", 1, 6),
    ]

    def beats(a: tuple[str, int, int], b: tuple[str, int, int]) -> bool:
        return (a[1], -a[2], a[0]) > (b[1], -b[2], b[0])

    scores = Counter({v[0]: 0 for v in vertices})
    cycles = 0
    for i, a in enumerate(vertices):
        for b in vertices[i + 1 :]:
            scores[(a if beats(a, b) else b)[0]] += 1
    for i, a in enumerate(vertices):
        for j, b in enumerate(vertices[i + 1 :], i + 1):
            for c in vertices[j + 1 :]:
                ab, bc, ca = beats(a, b), beats(b, c), beats(c, a)
                if (ab and bc and ca) or ((not ab) and (not bc) and (not ca)):
                    cycles += 1
    path = " > ".join(v[0] for v in sorted(vertices, key=lambda x: (x[1], -x[2], x[0]), reverse=True))
    print("TOURNAMENT ANALYSIS")
    print("  vertices: proof currencies, not runners or arcs")
    print("  observable: (certificate locality, compression depth, search entropy)")
    print("  switch/gauge: earlier bounded certificates beat larger global log currencies")
    print(f"  score_hist={dict(sorted(Counter(scores.values()).items()))}")
    print(f"  directed_3_cycles={cycles}; SCCs={(1,) * len(vertices)}; Hamiltonian_paths=1")
    print(f"  tie_HP={path}")
    print()


def assumption_challenge() -> None:
    print("ASSUMPTION CHALLENGE")
    print("  Candidate vertices considered: runners, safe arcs, component rows, owner")
    print("  pairs, w residues, prime powers, determinant rows, and proof obligations.")
    print("  Chosen vertices: Helly certificate budgets B_H(M).")
    print("  Preserved predicate: whether a bounded common w can survive the chosen")
    print("  component-language subfamily.")
    print("  Destroyed information: exact phase geometry, row labels outside the")
    print("  bounded w window, and the full CRT modulus.")
    print("  Challenged assumption: the logarithm need not come from denominator height;")
    print("  it can come from the number of possible local certificates.")


def main() -> None:
    print("S601 Helly-scale logarithm atlas")
    print()
    print_budget_atlas()
    print_sample_audit()
    print_laws()
    tournament_fingerprint()
    assumption_challenge()


if __name__ == "__main__":
    main()
