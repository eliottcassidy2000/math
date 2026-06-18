#!/usr/bin/env python3
"""
lrc14_finite_endpoint_feasibility_codex.py

Codex 2026-06-18: feasibility scout for the finite endpoint left by the
colored-resonance LRC(14) route.

Context.
  HYP-2595 says that a uniform Sigma floor plus the conjectural colored
  resonance discrepancy bound would certify all V >= 711.  This leaves a
  finite endpoint V < 711.  The question is whether that endpoint is feasible
  to compute directly in this session.

Findings encoded here.
  * Naive finite enumeration of S3 shapes below V=711 is astronomically large.
  * The tempting q=14V CRT finite check is not sufficient at low V: there are
    primitive q-covering rows with no q=14V witness, but exact M(S)>1/14.
  * Existing k=2 finite-core work remains the correct low-end model: prove a
    structural finite core first, then exact-check that core.
  * For k>=3 the immediate proof target should be a new structural endpoint
    reduction: resonance bound + Sigma floor for large V, and arc-width /
    endpoint-protection / exact-M cores for small V.

Tournament Analysis declaration.
  Vertex set: proof/end-point strategies, not runners:
    naive_shape_enumeration, q14V_milp, k2_finite_core, arc_width_dropmax,
    colored_resonance_largeV, endpoint_exactM_core.
  Pairwise observable: (proved scope, endpoint cost, remaining proof gap), with
    the lexicographic gauge larger scope, then smaller cost, then smaller gap.
  Tie Hamiltonian path: the order above.
  Fingerprints: score histogram, directed 3-cycles, SCCs, Hamiltonian paths.

Assumption challenge.
  I considered runners, large-cluster offsets E, q=14V residues, forbidden
  residue classes, exact THM-524 crossings, and proof strategies as vertices.
  The q=14V residue quotient preserves the constructive CRT witness predicate
  but destroys exact LRC witnesses at other denominators; V=15 already exposes
  that loss.  The challenged assumption is that the colored large-V modulus can
  also serve as the whole low-V finite endpoint.
"""

from __future__ import annotations

import itertools
import math
import statistics
import time
from collections import Counter, deque
from dataclasses import dataclass
from fractions import Fraction as F

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_array


C = F(1, 14)
QS = range(2, 15)


def nrm(x: F) -> F:
    r = x - int(x)
    if r < 0:
        r += 1
    return r if r <= F(1, 2) else 1 - r


def candidate_times(S: tuple[int, ...]) -> set[F]:
    S = tuple(sorted(set(S)))
    out: set[F] = {F(1, 2)}
    for v in S:
        k = 0
        while F(2 * k + 1, 2 * v) <= F(1, 2):
            out.add(F(2 * k + 1, 2 * v))
            k += 1
    for i, a in enumerate(S):
        for b in S[i + 1 :]:
            for d in (a + b, b - a):
                if d <= 0:
                    continue
                k = 1
                while F(k, d) <= F(1, 2):
                    out.add(F(k, d))
                    k += 1
    return out


def mexact(S: tuple[int, ...]) -> tuple[F, F]:
    best = F(0)
    arg = F(0)
    for t in candidate_times(S):
        m = min(nrm(v * t) for v in S)
        if m > best:
            best = m
            arg = t
    return best, arg


def mfloat(S: tuple[int, ...]) -> float:
    S = tuple(sorted(set(S)))
    cs: set[float] = {0.5}
    for v in S:
        k = 0
        while 2 * k + 1 <= v:
            cs.add((2 * k + 1) / (2.0 * v))
            k += 1
    for i, a in enumerate(S):
        for b in S[i + 1 :]:
            for d in (a + b, b - a):
                if d <= 0:
                    continue
                k = 1
                while 2 * k <= d:
                    cs.add(k / float(d))
                    k += 1
    best = 0.0
    for t in cs:
        m = 1.0
        for v in S:
            r = (v * t) % 1.0
            r = min(r, 1.0 - r)
            if r < m:
                m = r
            if m < best:
                break
        if m > best:
            best = m
    return best


def covering(S: tuple[int, ...]) -> bool:
    return all(any(v % q == 0 for v in S) for q in QS)


def primitive(S: tuple[int, ...]) -> bool:
    return math.gcd(*S) == 1


def case_s3(S: tuple[int, ...]) -> bool:
    S = tuple(sorted(S))
    k = sum(1 for v in S if v > 13)
    return k >= 2 and S[-1] >= 13 * S[0]


def endpoint_shape_counts(Vmax: int = 710) -> list[tuple[int, int]]:
    """
    Count raw S3 coordinates S=P union {V-e:e in E}, with 0 in E.
    This deliberately ignores covering, primitivity, and exact-M pruning.
    """
    rows = []
    for k in range(1, 14):
        total = 0
        for V in range(14, Vmax + 1):
            if V - 14 >= k - 1:
                total += math.comb(13, 13 - k) * math.comb(V - 14, k - 1)
        rows.append((k, total))
    return rows


def primes_upto(n: int) -> list[int]:
    out = []
    for x in range(2, n + 1):
        if all(x % p for p in out if p * p <= x):
            out.append(x)
    return out


@dataclass
class MilpResult:
    V: int
    status: int
    message: str
    seconds: float
    S: tuple[int, ...] | None
    exact_M: F | None = None
    exact_tau: F | None = None


def q14V_cover_milp(V: int, seconds: float = 5.0) -> MilpResult:
    """
    Search for a primitive q-covering 13-set S containing V whose forbidden
    residue classes cover Z/(14V).  Feasible means q=14V is NOT a universal
    finite endpoint for this V.  Infeasible is a strong q=14V certificate.
    """
    t0 = time.time()
    q = 14 * V
    rows: list[int] = []
    cols: list[int] = []
    data: list[int] = []
    lb: list[float] = []
    ub: list[float] = []
    row = 0

    # Cover every residue by at least one chosen speed's danger preimage.
    for a in range(q):
        for v in range(1, V + 1):
            r = (a * v) % q
            if 14 * min(r, q - r) < q:
                rows.append(row)
                cols.append(v - 1)
                data.append(1)
        lb.append(1.0)
        ub.append(np.inf)
        row += 1

    # Exactly 13 speeds.
    for v in range(1, V + 1):
        rows.append(row)
        cols.append(v - 1)
        data.append(1)
    lb.append(13.0)
    ub.append(13.0)
    row += 1

    # The endpoint/max speed is V.
    rows.append(row)
    cols.append(V - 1)
    data.append(1)
    lb.append(1.0)
    ub.append(1.0)
    row += 1

    # THM-523 q-covering obligations.
    for d in QS:
        for v in range(d, V + 1, d):
            rows.append(row)
            cols.append(v - 1)
            data.append(1)
        lb.append(1.0)
        ub.append(np.inf)
        row += 1

    # Primitive: for each prime p, at least one chosen speed is not divisible
    # by p.  This is equivalent to gcd(S)=1 for integer speeds <= V.
    for p in primes_upto(V):
        for v in range(1, V + 1):
            if v % p != 0:
                rows.append(row)
                cols.append(v - 1)
                data.append(1)
        lb.append(1.0)
        ub.append(np.inf)
        row += 1

    A = coo_array((data, (rows, cols)), shape=(row, V)).tocsr()
    constraints = LinearConstraint(A, np.array(lb), np.array(ub))
    res = milp(
        c=np.zeros(V),
        integrality=np.ones(V),
        bounds=Bounds(0, 1),
        constraints=constraints,
        options={"time_limit": seconds, "mip_rel_gap": 0},
    )
    S = None
    exact_M = None
    exact_tau = None
    if res.x is not None:
        S = tuple(i + 1 for i, x in enumerate(res.x) if x > 0.5)
        if len(S) == 13 and max(S) <= 62:
            exact_M, exact_tau = mexact(S)
    return MilpResult(V, int(res.status), str(res.message), time.time() - t0, S, exact_M, exact_tau)


def exact_k3_mini_core(B: int = 42) -> dict:
    """
    Small exact k=3 core scan.  This is not the endpoint; it is a calibrated
    model of what a structural finite core would need to produce.
    """
    t0 = time.time()
    total = 0
    exact_checked = 0
    below = 0
    worst: tuple[F, tuple[int, ...], F] | None = None
    threshold = 1 / 14
    large = tuple(range(14, B + 1))

    for P in itertools.combinations(range(1, 14), 10):
        miss = {q for q in QS if not any(p % q == 0 for p in P)}
        cover_bits = {v: {q for q in miss if v % q == 0} for v in large}
        for L in itertools.combinations(large, 3):
            if miss and set().union(*(cover_bits[v] for v in L)) != miss:
                continue
            S = tuple(sorted(P + L))
            if not primitive(S) or not case_s3(S):
                continue
            total += 1
            mf = mfloat(S)
            if mf < threshold + 0.012:
                exact_checked += 1
                m, tau = mexact(S)
                if worst is None or m < worst[0]:
                    worst = (m, S, tau)
                if m < C:
                    below += 1

    return {
        "B": B,
        "total": total,
        "exact_checked": exact_checked,
        "below": below,
        "worst": worst,
        "seconds": time.time() - t0,
    }


def count_directed_3cycles(adj: list[list[bool]]) -> int:
    total = 0
    for a, b, c in itertools.combinations(range(len(adj)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (adj[a][c] and adj[c][b] and adj[b][a]):
            total += 1
    return total


def scc_sizes(adj: list[list[bool]]) -> list[int]:
    n = len(adj)
    radj = [[adj[j][i] for j in range(n)] for i in range(n)]
    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for w, ok in enumerate(adj[v]):
            if ok and w not in seen:
                dfs(w)
        order.append(v)

    for v in range(n):
        if v not in seen:
            dfs(v)

    seen.clear()
    out = []
    for start in reversed(order):
        if start in seen:
            continue
        stack = [start]
        seen.add(start)
        size = 0
        while stack:
            v = stack.pop()
            size += 1
            for w, ok in enumerate(radj[v]):
                if ok and w not in seen:
                    seen.add(w)
                    stack.append(w)
        out.append(size)
    return sorted(out, reverse=True)


def hamiltonian_paths_count(adj: list[list[bool]]) -> int:
    n = len(adj)
    dp = [Counter() for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v, count in list(dp[mask].items()):
            for w in range(n):
                if count and not (mask & (1 << w)) and adj[v][w]:
                    dp[mask | (1 << w)][w] += count
    return sum(dp[-1].values())


def strategy_tournament() -> dict:
    names = [
        "naive_shape_enumeration",
        "q14V_milp",
        "k2_finite_core",
        "arc_width_dropmax",
        "colored_resonance_largeV",
        "endpoint_exactM_core",
    ]
    # Manual but explicit score triples from this session and prior canon:
    # larger scope is better; lower cost/gap is better after sign conversion.
    # scope: proven or credibly targeted endpoint coverage.
    # cost: inverse computational tractability.
    # gap: amount still conjectural.
    metrics = {
        "naive_shape_enumeration": (1, -6, -6),
        "q14V_milp": (3, -3, -4),
        "k2_finite_core": (6, -1, 0),
        "arc_width_dropmax": (5, -2, -2),
        "colored_resonance_largeV": (5, -1, -3),
        "endpoint_exactM_core": (4, -2, -2),
    }
    n = len(names)
    adj = [[False] * n for _ in range(n)]
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if i == j:
                continue
            if metrics[a] > metrics[b] or (metrics[a] == metrics[b] and i < j):
                adj[i][j] = True
    scores = {names[i]: sum(adj[i]) for i in range(n)}
    return {
        "names": names,
        "scores": scores,
        "score_hist": dict(sorted(Counter(scores.values()).items())),
        "cycles3": count_directed_3cycles(adj),
        "sccs": scc_sizes(adj),
        "hamiltonian_paths": hamiltonian_paths_count(adj),
        "leader": max(scores, key=scores.get),
    }


def main() -> None:
    print("=" * 78)
    print("LRC14 finite endpoint feasibility scout")
    print("=" * 78)
    counts = endpoint_shape_counts(710)
    total = sum(v for _, v in counts)
    print("\n1. Raw endpoint shape count below V=711")
    print(f"  total raw S3 coordinate shapes: {total}")
    for k, value in counts:
        note = "  (k=2 already structurally finite-core proved)" if k == 2 else ""
        print(f"  k={k:2d}: {value}{note}")
    print("  Conclusion: a literal all-shapes exact-M endpoint is not feasible by brute force.")

    print("\n2. q=14V covering-system MILP scout")
    print("  status 0=feasible obstruction, 1=time limit, 2=infeasible certificate")
    probes = [14, 15, 16, 20, 28, 42, 56, 70, 84, 98, 120]
    milp_rows = [q14V_cover_milp(V, seconds=5.0) for V in probes]
    for row in milp_rows:
        if row.S is None:
            print(f"  V={row.V:3d}: status={row.status}, seconds={row.seconds:5.2f}, S=None")
        else:
            mtxt = ""
            if row.exact_M is not None:
                mtxt = f", exact M={row.exact_M} at tau={row.exact_tau}"
            print(
                f"  V={row.V:3d}: status={row.status}, seconds={row.seconds:5.2f}, "
                f"S={row.S}{mtxt}"
            )
    print(
        "  Readout: q=14V is a strong large-V colored tool, but V=15 gives a "
        "primitive q-covering row with no q=14V witness while exact M=2/23. "
        "So the low endpoint cannot be closed by q=14V alone."
    )

    print("\n3. Small exact k=3 core calibration")
    mini = exact_k3_mini_core(B=42)
    print(
        f"  k=3, all large speeds <= {mini['B']}: scanned {mini['total']} "
        f"covering primitive S3 rows; exact-checked {mini['exact_checked']} near-threshold rows."
    )
    if mini["worst"] is not None:
        m, S, tau = mini["worst"]
        print(f"  worst exact checked M={m} at tau={tau}, S={S}")
    print(f"  # exact-checked rows below 1/14: {mini['below']}  ({mini['seconds']:.1f}s)")
    print(
        "  This is useful calibration, not the endpoint.  It says exact-M cores "
        "are workable after structural pruning; it does not license the V<711 brute force."
    )

    print("\n4. Proof-strategy tournament")
    tour = strategy_tournament()
    print(f"  scores: {tour['scores']}")
    print(
        f"  score_hist={tour['score_hist']}, cycles3={tour['cycles3']}, "
        f"sccs={tour['sccs']}, Hamiltonian_paths={tour['hamiltonian_paths']}, "
        f"leader={tour['leader']}"
    )
    print(
        "  The tournament ranks k2_finite_core as the proven model and "
        "colored_resonance_largeV as the best large-V lever; the missing bridge "
        "is a structural finite-core reduction for k>=3."
    )

    print("\n5. Next proof directions")
    print("  A. Promote HYP-2595: prove Delta <= 8*(k+cGP)+1 by bounding color-compatible V-resonances.")
    print("  B. Promote HYP-2593: prove the uniform Sigma floor, likely by bounded-spread plus Weyl tails.")
    print("  C. For V<711, avoid all-shape enumeration: classify q=14V obstructions, then exact-M only the residual cores.")
    print("  D. Generalize the k=2 finite-core pattern: drop one large speed, prove W(A) floor except a bounded hard core.")


if __name__ == "__main__":
    main()
