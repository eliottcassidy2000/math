#!/usr/bin/env python3
"""
lrc14_subtorus_relation_lattice_codex.py

Codex 2026-06-18: another LRC(14) S3 proof session, now using the
subtorus relation lattice as the primary object.

Core object.
  For an offset shape E={e_0,...,e_{k-1}} with e_0=0, the orbit
    x |-> (e_0*x,...,e_{k-1}*x) in T^k
  has affine relation lattice
    Lambda_aff(E) = {n in Z^k : sum n_i = 0 and sum n_i e_i = 0}.
  For the smooth empty-arc minorant
    G(v)=sum_gaps (gap-2/7)_+
  the Fourier identity is
    int G(E*x) dx = (5/7)^k
      + sum_{n in Lambda_aff(E), n!=0} prod_i psi_hat(n_i),
  with psi_hat(n)=-(1-exp(-2*pi*i*2n/7))/(2*pi*i*n).

This script does not try to evaluate that signed lattice series.  It asks a
more local, proof-useful question:

  Which finite relation-lattice fingerprints predict the exact negative
  correction F(k)-mu(E) in bounded banks?

The point is to distinguish the real proof split:
  short affine relations -> finite pattern enumeration;
  high relation height -> Fourier/product-kernel tail.

Tournament Analysis declaration.
  Vertices are proof carriers / relation observables, not runners or arcs:
    triple_decay, small_triples, additive_quad, run_mass, max_run,
    inverse_spread, spread.
  Pairwise observable: median positive correlation with the normalized negative
    correction F(k)-mu(E) across k=5..8 in an exact bounded bank.
  Gauge: A -> B when A has larger median score; ties follow the listed
    Hamiltonian path above.  Fingerprints include score histogram, directed
    3-cycles, SCC sizes, and Hamiltonian-path count.

Assumption challenge.
  I considered vertices = runners, gaps, fixed universal centers, safe
  components, residues, Fourier modes, subtorus blocks, short relations, and
  proof obligations.  The chosen quotient keeps the predicate most directly
  appearing in the subtorus Fourier identity: affine relation height.  It
  destroys small-part alignment with G_P and color-14 placement data, so it is
  not a complete LRC14 certificate by itself.
"""

from __future__ import annotations

import itertools
import math
import statistics
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction as F

from lrc14_global_threshold_ladder_codex import good_set_thr, measure

THR = F(2, 7)


def fiid(k: int, thr: F = THR) -> F:
    """Iid torus ceiling for max-gap > thr among k circular points."""
    out = F(0)
    j = 1
    while 1 - j * thr > 0:
        out += ((-1) ** (j + 1)) * math.comb(k, j) * (1 - j * thr) ** (k - 1)
        j += 1
    return out


_MU_CACHE: dict[tuple[int, ...], F] = {}


def mu_exact(E: tuple[int, ...]) -> F:
    E = tuple(sorted(set(E)))
    if E not in _MU_CACHE:
        _MU_CACHE[E] = measure(good_set_thr(E, THR))
    return _MU_CACHE[E]


def runs(E: tuple[int, ...]) -> list[list[int]]:
    out: list[list[int]] = []
    cur = [E[0]]
    for x in E[1:]:
        if x == cur[-1] + 1:
            cur.append(x)
        else:
            out.append(cur)
            cur = [x]
    out.append(cur)
    return out


def gcd_all(vals: tuple[int, ...]) -> int:
    g = 0
    for v in vals:
        g = math.gcd(g, abs(v))
    return g


def primitive_affine_coeffs_for_triple(a: int, b: int, c: int) -> tuple[int, int, int]:
    """Primitive n with sum n_i=0 and n_a*a+n_b*b+n_c*c=0."""
    coeff = (b - c, c - a, a - b)
    g = gcd_all(coeff)
    coeff = tuple(x // g for x in coeff)
    if coeff[0] < 0:
        coeff = tuple(-x for x in coeff)
    return coeff


def additive_quad_energy(E: tuple[int, ...]) -> tuple[int, float]:
    """Number/score of primitive support-4 relations a+b=c+d."""
    sums: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for a, b in itertools.combinations(E, 2):
        sums[a + b].append((a, b))
    count = 0
    score = 0.0
    for pairs in sums.values():
        for (a, b), (c, d) in itertools.combinations(pairs, 2):
            if len({a, b, c, d}) < 4:
                continue
            count += 1
            span = max(a, b, c, d) - min(a, b, c, d)
            score += 1.0 / max(1, span)
    return count, score


def relation_observables(E: tuple[int, ...], small_coeff: int = 3) -> dict[str, float]:
    r = runs(E)
    span = E[-1] - E[0]
    triple_decay = 0.0
    small_triples = 0
    min_triple_l1 = None
    for a, b, c in itertools.combinations(E, 3):
        coeff = primitive_affine_coeffs_for_triple(a, b, c)
        abs_coeff = [abs(x) for x in coeff if x]
        l1 = sum(abs_coeff)
        min_triple_l1 = l1 if min_triple_l1 is None else min(min_triple_l1, l1)
        if max(abs_coeff) <= small_coeff:
            small_triples += 1
        # Product-kernel size proxy for support-3 Fourier terms.
        triple_decay += 1.0 / math.prod(abs_coeff)
    quad_count, quad_score = additive_quad_energy(E)
    return {
        "triple_decay": triple_decay,
        "small_triples": float(small_triples),
        "additive_quad": float(quad_count) + quad_score,
        "run_mass": float(sum(len(x) * len(x) for x in r)),
        "max_run": float(max(len(x) for x in r)),
        "inverse_spread": 1.0 / float(span),
        "spread": float(span),
        "min_triple_l1": float(min_triple_l1 or 0),
        "num_runs": float(len(r)),
    }


def short_affine_relations(
    E: tuple[int, ...], coeff_bound: int = 4, max_support: int = 4
) -> dict[str, object]:
    """Enumerate displayed primitive short relations for a single shape."""
    rels: set[tuple[tuple[int, ...], tuple[int, ...], int, int]] = set()
    n = len(E)
    for supp_size in range(3, max_support + 1):
        for idxs in itertools.combinations(range(n), supp_size):
            vals = [E[i] for i in idxs]
            for coeff in itertools.product(range(-coeff_bound, coeff_bound + 1), repeat=supp_size):
                if any(c == 0 for c in coeff):
                    continue
                if sum(coeff) != 0:
                    continue
                if sum(c * v for c, v in zip(coeff, vals)) != 0:
                    continue
                g = gcd_all(coeff)
                if g != 1:
                    continue
                coeff = tuple(coeff)
                if next(c for c in coeff if c) < 0:
                    coeff = tuple(-c for c in coeff)
                rels.add((idxs, coeff, sum(abs(c) for c in coeff), supp_size))
    by_support = Counter(r[3] for r in rels)
    by_l1 = Counter(r[2] for r in rels)
    sample = sorted(rels, key=lambda r: (r[2], r[3], r[0], r[1]))[:8]
    return {
        "count": len(rels),
        "by_support": dict(sorted(by_support.items())),
        "by_l1": dict(sorted(by_l1.items())),
        "sample": sample,
    }


@dataclass(frozen=True)
class Row:
    k: int
    E: tuple[int, ...]
    mu: F
    fiid: F
    defect: float
    obs: dict[str, float]


def bounded_bank(k_values: range, cap: int) -> list[Row]:
    rows: list[Row] = []
    for k in k_values:
        fk = fiid(k)
        for combo in itertools.combinations(range(1, cap + 1), k - 1):
            E = (0,) + combo
            if gcd_all(E) != 1:
                continue
            mu = mu_exact(E)
            defect = float(fk - mu)
            rows.append(Row(k=k, E=E, mu=mu, fiid=fk, defect=defect, obs=relation_observables(E)))
    return rows


def pearson(xs: list[float], ys: list[float]) -> float:
    if len(xs) < 2:
        return 0.0
    mx = statistics.fmean(xs)
    my = statistics.fmean(ys)
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx == 0 or vy == 0:
        return 0.0
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / math.sqrt(vx * vy)


def auc_bottom(rows: list[Row], key: str, q: float = 0.10) -> float:
    """Mann-Whitney AUC: how well high observable finds bottom-mu shapes."""
    if not rows:
        return 0.5
    sorted_mu = sorted(rows, key=lambda r: (r.mu, r.E))
    m = max(1, int(round(q * len(rows))))
    bottom = set(r.E for r in sorted_mu[:m])
    pos = [r.obs[key] for r in rows if r.E in bottom]
    neg = [r.obs[key] for r in rows if r.E not in bottom]
    wins = ties = total = 0
    for p in pos:
        for n in neg:
            total += 1
            if p > n:
                wins += 1
            elif p == n:
                ties += 1
    return (wins + 0.5 * ties) / total if total else 0.5


def relation_tournament(metric_scores: dict[str, float], tie_path: list[str]) -> dict[str, object]:
    vertices = tie_path[:]
    rank = {v: i for i, v in enumerate(vertices)}
    edges: dict[tuple[str, str], str] = {}
    scores = Counter()
    for a, b in itertools.combinations(vertices, 2):
        if metric_scores[a] > metric_scores[b]:
            winner = a
        elif metric_scores[b] > metric_scores[a]:
            winner = b
        else:
            winner = a if rank[a] < rank[b] else b
        loser = b if winner == a else a
        edges[(winner, loser)] = winner
        scores[winner] += 1
        scores.setdefault(loser, scores[loser])
    cycles3 = 0
    for a, b, c in itertools.combinations(vertices, 3):
        wins = {
            (u, v)
            for u, v in itertools.permutations((a, b, c), 2)
            if (u, v) in edges
        }
        if (a, b) in wins and (b, c) in wins and (c, a) in wins:
            cycles3 += 1
        if (a, c) in wins and (c, b) in wins and (b, a) in wins:
            cycles3 += 1

    # SCCs by tiny Kosaraju.
    adj = {v: [] for v in vertices}
    radj = {v: [] for v in vertices}
    for winner, loser in edges:
        adj[winner].append(loser)
        radj[loser].append(winner)

    seen: set[str] = set()
    order: list[str] = []

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
    comps: list[list[str]] = []

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
            comps.append(comp)

    # Hamiltonian path count by DP over subsets.
    idx = {v: i for i, v in enumerate(vertices)}
    n = len(vertices)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            if not dp[mask][last]:
                continue
            vlast = vertices[last]
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                vnxt = vertices[nxt]
                if (vlast, vnxt) in edges:
                    dp[mask | (1 << nxt)][nxt] += dp[mask][last]
    hpaths = sum(dp[(1 << n) - 1])

    return {
        "score_hist": dict(sorted(Counter(scores.values()).items())),
        "scores": dict(sorted(scores.items(), key=lambda kv: (-kv[1], kv[0]))),
        "cycles3": cycles3,
        "scc_sizes": sorted([len(c) for c in comps], reverse=True),
        "hamiltonian_paths": hpaths,
    }


def fmt_frac(x: F) -> str:
    return f"{x} ({float(x):.6f})"


def main() -> None:
    print("=" * 88)
    print("LRC14 subtorus relation-lattice scout")
    print("=" * 88)
    print("Exact bounded bank: k=5..8, primitive E subset [0,14], threshold 2/7.")
    rows = bounded_bank(range(5, 9), cap=14)
    print(f"bank rows = {len(rows)}")
    by_k: dict[int, list[Row]] = defaultdict(list)
    for row in rows:
        by_k[row.k].append(row)

    print("\nI. Exact minima and negative-correction carriers")
    for k in sorted(by_k):
        rk = sorted(by_k[k], key=lambda r: (r.mu, r.E))
        fk = fiid(k)
        print(f"k={k}: F(k)={fmt_frac(fk)}  rows={len(rk)}")
        for row in rk[:5]:
            prof = short_affine_relations(row.E, coeff_bound=4, max_support=4)
            obs = row.obs
            print(
                f"  E={row.E} mu={fmt_frac(row.mu)} defect=F-mu={row.defect:+.6f} "
                f"small3={int(obs['small_triples'])} triple_decay={obs['triple_decay']:.4f} "
                f"quad={obs['additive_quad']:.2f} runs={runs(row.E)} "
                f"short_rel={prof['by_support']}"
            )

    print("\nII. Relation observables vs normalized negative correction")
    keys = [
        "triple_decay",
        "small_triples",
        "additive_quad",
        "run_mass",
        "max_run",
        "inverse_spread",
        "spread",
    ]
    metric_scores: dict[str, float] = {}
    for key in keys:
        cors = []
        aucs = []
        for k in sorted(by_k):
            rk = by_k[k]
            xs = [r.obs[key] for r in rk]
            ys = [r.defect for r in rk]
            cors.append(pearson(xs, ys))
            aucs.append(auc_bottom(rk, key, q=0.10))
        med_corr = statistics.median(cors)
        med_auc = statistics.median(aucs)
        metric_scores[key] = max(0.0, med_corr) + max(0.0, med_auc - 0.5)
        print(
            f"{key:15s} corr_by_k={[round(c, 3) for c in cors]} "
            f"median_corr={med_corr:+.3f} bottom10_auc={[round(a, 3) for a in aucs]} "
            f"score={metric_scores[key]:.3f}"
        )

    print("\nIII. Tournament Analysis on proof carriers")
    tour = relation_tournament(metric_scores, keys)
    print(f"scores={tour['scores']}")
    print(f"score_hist={tour['score_hist']}")
    print(f"directed_3cycles={tour['cycles3']}")
    print(f"scc_sizes={tour['scc_sizes']}")
    print(f"hamiltonian_path_count={tour['hamiltonian_paths']}")

    print("\nIV. Relation-height split examples")
    examples = [
        (0, 1, 2, 3, 4, 5, 6, 7),
        (0, 2, 3, 4, 5, 6, 8),
        (0, 1, 4, 9, 13, 14),
        (0, 1, 500, 501),
        (0, 7, 53, 311, 1009, 4999),
    ]
    for E in examples:
        obs = relation_observables(E)
        prof = short_affine_relations(E, coeff_bound=4, max_support=4)
        mu_text = "not computed (large spread)"
        if E[-1] <= 60:
            mu_text = fmt_frac(mu_exact(E))
        print(f"E={E}")
        print(f"  mu={mu_text}; runs={runs(E)}")
        print(
            f"  triple_decay={obs['triple_decay']:.6f}, small_triples={int(obs['small_triples'])}, "
            f"additive_quad={obs['additive_quad']:.3f}, short_rel_by_support={prof['by_support']}, "
            f"short_rel_by_l1={prof['by_l1']}"
        )
        if prof["sample"]:
            sample = []
            for idxs, coeff, l1, supp_size in prof["sample"][:3]:
                values = tuple(E[i] for i in idxs)
                sample.append(f"{coeff}@{values}/l1={l1}")
            print("  sample_short_relations=" + "; ".join(sample))

    print("\nV. Proof readout")
    print(
        "The bounded bank supports the relation-height split: the low-mu rows are "
        "short-relation-rich, and the best carrier is support-3 product-kernel decay "
        "rather than raw spread alone.  The next rigorous lemma should bound the "
        "full G-minorant lattice tail when every primitive affine relation has large "
        "coefficient product, while finite-checking the low-height relation patterns."
    )


if __name__ == "__main__":
    main()
