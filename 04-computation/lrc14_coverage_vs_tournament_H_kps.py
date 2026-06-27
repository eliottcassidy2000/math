"""
lrc14_coverage_vs_tournament_H_kps.py  (kind-pasteur-2026-06-27-S31ag)

NEW PROJECT-NATIVE ANGLE: is the LRC coverage extremality "consec/AP maximizes
P(N=0)" the SAME as the project's PROVED "regular maximizes H = #Ham paths"?

The AP cluster's winding tournament is R_k (regular rotational, the H-maximizer).
If coverage-max <=> H-max, the LRC extremality follows from THM-027 (regular max H).

For each k-cluster E={0,e_1,...,e_{k-1}}, compute:
  - P(N=0) = meas{x: all 7 sectors hit}  (the coverage / cover-bound LHS)
  - H = #Hamiltonian paths of the winding tournament T(E): i->j iff (e_i-e_j) mod 14 in {1..6}
    (the project's apex winding tournament, HYP-2922; magnitude-blind, uses e_i mod 14)
Check: does consec maximize BOTH? Do P(N=0) and H correlate? Does the magnitude
(same residues, different scale) change P(N=0) while H is fixed (the known blindness)?
"""
import sys, itertools, random
from fractions import Fraction as F

def sector_of(p): return int((p % 1) * 7)

def P_N0(E):
    """meas{x: all 7 sectors hit} = q_0, exact."""
    E = sorted(set(E)); b = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for m in range(0, 7 * e + 1): b.add(F(m, 7 * e))
    b = sorted(b); q0 = F(0)
    for i in range(len(b) - 1):
        x0, x1 = b[i], b[i + 1]
        if x1 <= x0: continue
        cov = set(sector_of(e * ((x0 + x1) / 2)) for e in E)
        if len(cov) == 7: q0 += x1 - x0
    return q0

def winding_tournament(E):
    """adjacency: adj[i][j]=True if i->j, i.e. (e_i-e_j) mod 14 in {1..6}. Returns None if a tie (mod14 in {0,7})."""
    k = len(E)
    adj = [[False]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i == j: continue
            d = (E[i] - E[j]) % 14
            if d in (0, 7):
                return None  # tie / antipodal -> not a tournament
            adj[i][j] = d in (1,2,3,4,5,6)
    return adj

def ham_paths(adj):
    """count Hamiltonian paths (directed) via Held-Karp DP. k <= ~12."""
    k = len(adj)
    if k == 0: return 1
    from functools import lru_cache
    full = (1<<k) - 1
    # dp[mask][last] = # paths covering mask ending at last
    dp = [[0]*k for _ in range(1<<k)]
    for v in range(k): dp[1<<v][v] = 1
    for mask in range(1<<k):
        for last in range(k):
            c = dp[mask][last]
            if c == 0: continue
            for nxt in range(k):
                if mask & (1<<nxt): continue
                if adj[last][nxt]:
                    dp[mask | (1<<nxt)][nxt] += c
    return sum(dp[full][v] for v in range(k))

def H_of_cluster(E):
    adj = winding_tournament(E)
    if adj is None: return None
    return ham_paths(adj)

if __name__ == "__main__":
    sys.stdout.reconfigure(line_buffering=True)
    random.seed(27)
    for k in (8, 9, 10):
        consec = tuple(range(k))
        pN0_c = P_N0(consec); H_c = H_of_cluster(consec)
        print(f"\n=== k={k} ===  consec: P(N=0)={float(pN0_c):.5f}  H={H_c}")
        # sample clusters (with 0), compute P(N=0) and H, check extremality + correlation
        data = [("consec", consec, pN0_c, H_c)]
        for _ in range(150):
            cfg = tuple(sorted([0] + random.sample(range(1, 28), k-1)))
            if len(set(cfg)) != k: continue
            H = H_of_cluster(cfg)
            if H is None: continue
            p = P_N0(cfg)
            data.append(("rand", cfg, p, H))
        # extremality: does consec maximize P(N=0)? maximize H?
        max_p = max(d[2] for d in data); max_H = max(d[3] for d in data)
        argmax_p = [d for d in data if d[2]==max_p][0]
        argmax_H = [d for d in data if d[3]==max_H][0]
        print(f"  max P(N=0) = {float(max_p):.5f} at {argmax_p[0]} (consec? {argmax_p[1]==consec})")
        print(f"  max H      = {max_H} at {argmax_H[0]} (consec? {argmax_H[1]==consec})")
        # correlation (only tournament-valid clusters; consec itself may be degenerate)
        import statistics
        valid = [d for d in data if d[3] is not None]
        print(f"  tournament-valid clusters: {len(valid)}/{len(data)} "
              f"(consec valid? {H_c is not None} -- antipodal-mod-14 tie at apex 7)")
        if len(valid) < 3:
            print("  too few valid for correlation"); continue
        ps = [float(d[2]) for d in valid]; Hs = [float(d[3]) for d in valid]
        # Pearson
        mp = statistics.mean(ps); mH = statistics.mean(Hs)
        cov = sum((p-mp)*(h-mH) for p,h in zip(ps,Hs))
        sp = (sum((p-mp)**2 for p in ps))**0.5; sH = (sum((h-mH)**2 for h in Hs))**0.5
        pear = cov/(sp*sH) if sp*sH>0 else 0
        print(f"  Pearson corr(P(N=0), H) = {pear:+.3f}  over {len(data)} clusters")
        # magnitude test: same residues mod 14, different scale -> same H, different P(N=0)?
        base = consec
        scaled = tuple(e + 14*(i%3) for i,e in enumerate(base))  # shift some by mult of 14 (same residue)
        if len(set(scaled))==k:
            print(f"  magnitude test: scaled {scaled} H={H_of_cluster(scaled)} (vs {H_c}) "
                  f"P(N=0)={float(P_N0(scaled)):.5f} (vs {float(pN0_c):.5f})")
