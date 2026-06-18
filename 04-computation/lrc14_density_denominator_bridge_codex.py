#!/usr/bin/env python3
"""
lrc14_density_denominator_bridge_codex.py

Codex 2026-06-18: bridge the current LRC(14) S3 endgame

    rho*(P,E) = meas{x in G_P : maxgap({e*x}) > 2/7} >= c0 > 0

with the gap-side/binding-denominator language

    M(S) = j/D > 1/14, so D = 14*j - r with remainder debt r = 14*j - D.

The point of this scout is not to re-prove THM-527/528.  It tests the
translation layer:

  * exact rho*(P,E) for several known hard shapes;
  * discretized CRT density rho_q at q=14*Vmax for actual covering sets S;
  * exact M(S), binding pair denominators, and the small-r profile;
  * a Tournament Analysis quotient on blocking dominance.

Tournament Analysis declaration.
  Vertex set chosen here: runners v in S.
  Pairwise observable: with all other runners fixed, how many of the remaining
    residues mod q=14*Vmax does v newly forbid, compared with w?
  Gauge: orient v -> w if v has larger conditional blocking count.
  Tie Hamiltonian path: increasing numerical runner order.
  Fingerprints: score histogram, directed 3-cycles, SCC sizes, Hamiltonian paths.

Assumption challenge.
  I considered runner, residue, four-window, safe-component, denominator-remainder,
  and proof-obligation vertices.  Runner-vertices preserve the exact CRT covering
  obstruction and expose blocking dominance, but they destroy the geometry of the
  Good_E interval components.  So this is a deliberately partial quotient: useful
  if small-r denial is a blocker-competition phenomenon, incomplete if the real
  proof lives in interval-component or denominator-event vertices.
"""

from __future__ import annotations

import itertools
import random
from collections import Counter, defaultdict, deque
from fractions import Fraction as F
from functools import reduce
from math import gcd

H = F(1, 14)
THR = F(2, 7)
RNG = random.Random(140618)


# ---------------------------------------------------------------------------
# Exact rho* machinery (adapted from the THM-527/Angle-E exact engine)

def danger_arcs(u: int, h: F = H):
    out = []
    for j in range(u):
        c = F(j, u)
        a = (c - h / u) % 1
        b = (c + h / u) % 1
        if a < b:
            out.append((a, b))
        else:
            out.append((a, F(1)))
            out.append((F(0), b))
    return out


def merge(arcs):
    arcs = sorted(arcs)
    out = []
    for a, b in arcs:
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def safe_set(A, h: F = H):
    if not A:
        return [(F(0), F(1))]
    dz = merge([iv for u in A for iv in danger_arcs(u, h)])
    safe = []
    prev = F(0)
    for a, b in dz:
        if a > prev:
            safe.append((prev, a))
        prev = max(prev, b)
    if prev < 1:
        safe.append((prev, F(1)))
    return safe


def intersect(A, B):
    out = []
    i = j = 0
    A = sorted(A)
    B = sorted(B)
    while i < len(A) and j < len(B):
        a0, a1 = A[i]
        b0, b1 = B[j]
        lo = max(a0, b0)
        hi = min(a1, b1)
        if lo < hi:
            out.append((lo, hi))
        if a1 < b1:
            i += 1
        else:
            j += 1
    return out


def meas(arcs):
    return sum(b - a for a, b in arcs)


def good_set_exact(E):
    es = sorted(set(E))
    bps = {F(0), F(1)}
    for e in es:
        for n in range(1, e):
            bps.add(F(n, e))
    for i, ei in enumerate(es):
        for ej in es[i + 1:]:
            d = ej - ei
            for n in range(1, d):
                bps.add(F(n, d))
    bp = sorted(bps)
    good = []
    for lo, hi in zip(bp, bp[1:]):
        mid = (lo + hi) / 2
        valmap = {}
        for e in es:
            v = (e * mid) % 1
            if v not in valmap:
                valmap[v] = (e, e * mid - v)
        keys = sorted(valmap)
        slot_e = [valmap[v][0] for v in keys]
        slot_c = [valmap[v][1] for v in keys]
        m = len(keys)
        if m == 1:
            good.append((lo, hi))
            continue
        sub = []
        for k in range(m):
            if k < m - 1:
                alpha = slot_e[k + 1] - slot_e[k]
                beta = -(slot_c[k + 1] - slot_c[k])
            else:
                alpha = slot_e[0] - slot_e[m - 1]
                beta = F(1) - (slot_c[0] - slot_c[m - 1])
            if alpha == 0:
                if beta > THR:
                    sub.append((lo, hi))
            else:
                xb = (THR - beta) / alpha
                if alpha > 0:
                    a = max(lo, xb)
                    if a < hi:
                        sub.append((a, hi))
                else:
                    b = min(hi, xb)
                    if lo < b:
                        sub.append((lo, b))
        good.extend(merge(sub))
    return merge(good)


def rho_star_exact(P, E):
    good = good_set_exact(E)
    gp = safe_set(P)
    both = intersect(good, gp)
    return meas(both), meas(gp), meas(good), both


# ---------------------------------------------------------------------------
# Actual finite covering sets and exact M/binding-denominator profile

def circ_norm(x: F):
    r = x - int(x)
    if r < 0:
        r += 1
    return r if r <= F(1, 2) else 1 - r


def gcd_all(vals):
    return reduce(gcd, vals, 0)


def is_covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))


def classify(S):
    return "S3" if max(S) >= 13 * min(S) and sum(v > 13 for v in S) >= 2 else "not-S3"


def build_S(P, E, Vmax):
    L = [Vmax - e for e in E]
    S = sorted(set(P) | set(L))
    ok = (
        len(S) == 13
        and all(v > 13 for v in L)
        and gcd_all(S) == 1
        and is_covering(S)
        and classify(S) == "S3"
    )
    return S, ok


def candidate_taus(S):
    out = {F(1, 2)}
    S = sorted(S)
    for v in S:
        k = 0
        while True:
            t = F(2 * k + 1, 2 * v)
            if t > F(1, 2):
                break
            out.add(t)
            k += 1
    for a, b in itertools.combinations(S, 2):
        for D in (a + b, abs(b - a)):
            if D <= 0:
                continue
            k = 1
            while True:
                t = F(k, D)
                if t > F(1, 2):
                    break
                out.add(t)
                k += 1
    return out


def exact_M(S):
    best = F(0)
    ats = []
    for t in candidate_taus(S):
        val = min(circ_norm(v * t) for v in S)
        if val > best:
            best = val
            ats = [t]
        elif val == best:
            ats.append(t)
    return best, sorted(ats)


def binding_records(S, tau):
    M = min(circ_norm(v * tau) for v in S)
    binders = [v for v in S if circ_norm(v * tau) == M]
    rows = []
    for a, b in itertools.combinations(binders, 2):
        for kind, D in (("sum", a + b), ("diff", abs(b - a))):
            if D <= 0:
                continue
            if (D * tau).denominator != 1:
                continue
            j = M * D
            if j.denominator == 1:
                jj = int(j)
                rows.append({
                    "pair": (a, b),
                    "kind": kind,
                    "D": D,
                    "j": jj,
                    "r": 14 * jj - D,
                })
    return M, binders, rows


def is_safe_residue(r, q):
    d = min(r % q, (-r) % q)
    return 14 * d >= q


def good_residues(S, q):
    return [a for a in range(q) if all(is_safe_residue(v * a, q) for v in S)]


def circ_frac(x: F):
    r = x - int(x)
    if r < 0:
        r += 1
    return r


def maxgap_frac(points):
    pts = sorted(set(circ_frac(p) for p in points))
    if len(pts) <= 1:
        return F(1)
    gaps = [pts[i + 1] - pts[i] for i in range(len(pts) - 1)]
    gaps.append(pts[0] + 1 - pts[-1])
    return max(gaps)


def rho_star_grid_residues(P, E, Vmax, q):
    """Residues sampling the limiting rho* carrier on the q-grid.

    These are not automatically actual CRT witnesses at the same residue:
    maxgap(E*x)>2/7 says a suitable fast phase exists, while a/q fixes the
    finite fast phase Vmax*x.  This is exactly the integer-vs-real phase
    alignment that THM-527 leaves as part of the tail/rate problem.
    """
    out = []
    for a in range(q):
        x = F(a, q)
        if circ_norm(Vmax * x) < H:
            continue
        if any(circ_norm(p * x) < H for p in P):
            continue
        if maxgap_frac([e * x for e in E]) > THR:
            out.append(a)
    return out


# ---------------------------------------------------------------------------
# Blocking tournament: conditional marginal blocker dominance at q=14*Vmax.

def bad_set(v, q):
    return {a for a in range(q) if not is_safe_residue(v * a, q)}


def scc_sizes(edge):
    n = len(edge)
    seen = [False] * n
    order = []

    def dfs(v):
        seen[v] = True
        for w in range(n):
            if edge[v][w] and not seen[w]:
                dfs(w)
        order.append(v)

    for v in range(n):
        if not seen[v]:
            dfs(v)
    redge = [[edge[j][i] for j in range(n)] for i in range(n)]
    seen = [False] * n
    sizes = []

    def rdfs(v):
        seen[v] = True
        total = 1
        for w in range(n):
            if redge[v][w] and not seen[w]:
                total += rdfs(w)
        return total

    for v in reversed(order):
        if not seen[v]:
            sizes.append(rdfs(v))
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(edge):
    n = len(edge)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp[mask][last]
            if not val:
                continue
            rem = ((1 << n) - 1) ^ mask
            while rem:
                bit = rem & -rem
                nxt = bit.bit_length() - 1
                if edge[last][nxt]:
                    dp[mask | bit][nxt] += val
                rem ^= bit
    return sum(dp[-1])


def blocking_tournament(S, q):
    S = list(S)
    n = len(S)
    universe = set(range(q))
    bad = {v: bad_set(v, q) for v in S}
    edge = [[False] * n for _ in range(n)]
    margins = {}
    for i, v in enumerate(S):
        for j, w in enumerate(S):
            if i >= j:
                continue
            other_bad = set()
            for u in S:
                if u != v and u != w:
                    other_bad |= bad[u]
            context = universe - other_bad
            bv = len(bad[v] & context)
            bw = len(bad[w] & context)
            if bv > bw or (bv == bw and v < w):
                edge[i][j] = True
                margins[(v, w)] = bv - bw
            else:
                edge[j][i] = True
                margins[(w, v)] = bw - bv
    scores = [sum(edge[i][j] for j in range(n)) for i in range(n)]
    cycles = 0
    for a, b, c in itertools.combinations(range(n), 3):
        degs = [
            int(edge[a][b]) + int(edge[a][c]),
            int(edge[b][a]) + int(edge[b][c]),
            int(edge[c][a]) + int(edge[c][b]),
        ]
        if sorted(degs) == [1, 1, 1]:
            cycles += 1
    leaders = sorted(zip(scores, S), reverse=True)[:5]
    return {
        "score_hist": dict(sorted(Counter(scores).items())),
        "cycles3": cycles,
        "scc_sizes": scc_sizes(edge),
        "ham_paths": hamiltonian_path_count(edge),
        "leaders": leaders,
    }


# ---------------------------------------------------------------------------
# Scenario reports

def fstr(x: F):
    return f"{x} = {float(x):.7f}"


def report_scenario(name, P, E, Vs, do_M_until=160):
    print("\n" + "=" * 88)
    print(f"SCENARIO: {name}")
    print("=" * 88)
    rho, gp, mu, intervals = rho_star_exact(P, E)
    print(f"  P={P}")
    print(f"  E={E}  k={len(E)} spread={max(E)}")
    print(f"  exact mu(E)       = {fstr(mu)}")
    print(f"  exact meas(G_P)   = {fstr(gp)}")
    print(f"  exact rho*(P,E)   = {fstr(rho)}  intervals={len(intervals)}")
    if intervals:
        print("  first witness intervals:")
        for a, b in intervals[:4]:
            print(f"    ({a}, {b}) width={b-a}")
    for V in Vs:
        S, ok = build_S(P, E, V)
        print(f"\n  Vmax={V}: ok={ok} S={S}")
        if not ok:
            continue
        qmod = 14 * V
        A_global = good_residues(S, qmod)
        rho_global = F(len(A_global), qmod)
        A_grid = rho_star_grid_residues(P, E, V, qmod)
        rho_grid = F(len(A_grid), qmod)
        print(f"    global CRT witness density = {len(A_global)}/{qmod} = {float(rho_global):.7f}")
        print(f"    rho* grid-opportunity density = {len(A_grid)}/{qmod} = {float(rho_grid):.7f}"
              f"  (grid/rho*={float(rho_grid / rho):.3f})")
        tour = blocking_tournament(S, qmod)
        print(f"    blocking tournament: score_hist={tour['score_hist']}"
              f" cycles3={tour['cycles3']} scc={tour['scc_sizes']}"
              f" Hpaths={tour['ham_paths']}")
        print(f"      leaders(score,runner)={tour['leaders']}")
        if V <= do_M_until:
            M, ats = exact_M(S)
            print(f"    exact M(S) = {fstr(M)} at {len(ats)} tau(s); first tau={ats[0]}")
            M0, binders, rows = binding_records(S, ats[0])
            rows = sorted(rows, key=lambda r: (r["r"], r["D"], r["pair"]))
            print(f"      binders={binders}")
            if rows:
                best = rows[0]
                print("      binding denominator records (smallest remainder first):")
                for r in rows[:5]:
                    print(f"        {r['kind']:>4s} pair={r['pair']} D={r['D']} "
                          f"j={r['j']}  D=14*j-{r['r']}")
                print(f"      smallest M-binding remainder debt r={best['r']}")
            else:
                print("      no pair-denominator record found at first optimum tau")


def gen_random_covering_sets(target=80, vmax_hi=130):
    out = []
    seen = set()
    tries = 0
    while len(out) < target and tries < 500_000:
        tries += 1
        k = RNG.randint(3, 10)
        psize = 13 - k
        if psize < 1:
            continue
        P = RNG.sample(range(1, 14), psize)
        V = RNG.randint(max(25, k + 20), vmax_hi)
        spread = RNG.randint(k - 1, min(3 * k + 8, V - 14))
        E = sorted(set([0, spread] + RNG.sample(range(1, spread), k - 2)))
        if len(E) != k:
            continue
        S, ok = build_S(P, E, V)
        t = tuple(S)
        if ok and t not in seen:
            seen.add(t)
            out.append((S, P, E, V))
    return out


def random_remainder_census():
    print("\n" + "=" * 88)
    print("RANDOM CENSUS: M-binding small-remainder debt vs CRT density")
    print("=" * 88)
    rows = []
    for S, P, E, V in gen_random_covering_sets():
        M, ats = exact_M(S)
        _, binders, recs = binding_records(S, ats[0])
        if not recs:
            continue
        rec = min(recs, key=lambda r: (r["r"], r["D"]))
        qmod = 14 * V
        rho_global = F(len(good_residues(S, qmod)), qmod)
        rho_grid = F(len(rho_star_grid_residues(P, E, V, qmod)), qmod)
        rows.append((rec["r"], M, rec["D"], rec["j"], rho_global, rho_grid, S, P, E, V, binders))
    rows.sort(key=lambda t: (t[0], t[1]))
    print(f"  exact rows with a pair-binding record: {len(rows)}")
    print(f"  LRC breaks below 1/14: {sum(1 for _, M, *_ in rows if M < H)}")
    print(f"  remainder histogram r<=30: {dict(sorted(Counter(r for r, *_ in rows if r <= 30).items()))}")
    if rows:
        print(f"  min global CRT density in census = {min(r[4] for r in rows)}"
              f" = {float(min(r[4] for r in rows)):.7f}")
        print(f"  min rho* grid-opportunity density in census = {min(r[5] for r in rows)}"
              f" = {float(min(r[5] for r in rows)):.7f}")
        print("  12 smallest-remainder rows:")
        for r, M, D, j, rho_global, rho_grid, S, P, E, V, binders in rows[:12]:
            print(f"    r={r:3d}  M={M}={float(M):.7f}  D={D}=14*{j}-{r}"
                  f"  rho_global={rho_global}={float(rho_global):.5f}"
                  f"  rho_grid={rho_grid}={float(rho_grid):.5f} V={V}"
                  f"  binders={binders} P={P} E={E}")
    print("\n  Interpretation:")
    print("    Small r is the finite gap-side avatar of a near-1/14 optimum.  In this")
    print("    bounded S3 census, the same rows still have positive CRT witness density;")
    print("    the obstruction is therefore not a single pair denominator by itself, but")
    print("    a covering-system alignment that would have to erase every rho_q witness.")


def main():
    print("=" * 88)
    print("LRC(14) density-floor / binding-denominator bridge scout")
    print("=" * 88)
    print("All exact decisions use Fraction arithmetic except the displayed decimals.")

    scenarios = [
        (
            "Angle-E k=5 exact rho* minimizer (scaled consecutive resonance)",
            [1, 2, 3, 4, 7, 9, 12, 13],
            [0, 2, 4, 6, 8],
            [28, 42, 70, 140, 280],
        ),
        (
            "Angle-E k=6 exact rho* minimizer (raw consecutive)",
            [1, 2, 3, 7, 9, 11, 13],
            [0, 1, 2, 3, 4, 5],
            [42, 70, 112, 140, 280],
        ),
        (
            "Angle-B broad-scan 1/90 witness shape",
            [1, 2, 9],
            [0, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            [56, 84, 140, 280],
        ),
    ]
    for args in scenarios:
        report_scenario(*args)
    random_remainder_census()

    print("\n" + "=" * 88)
    print("TAKEAWAY")
    print("=" * 88)
    print("  The density formulation gives an interval-residue reservoir: rho*(P,E)>0")
    print("  turns into linearly many grid opportunities; the actual CRT witnesses are")
    print("  the aligned subset where the finite fast phase also lands correctly.")
    print("  The small-r denominator language is equivalent only after all other runners")
    print("  clear at the same crossing.  The useful next lemma should therefore bound")
    print("  the ability of the conditional blockers to cover the rho* reservoir, rather")
    print("  than trying to prove a private pair-denominator inequality in isolation.")


if __name__ == "__main__":
    main()
