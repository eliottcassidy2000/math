#!/usr/bin/env python3
"""
blue_as_tournament_s20cr.py — The blue subgraph as an emerging tournament
kind-pasteur-2026-03-23-S20cr

KEY INSIGHT (from human):
  The blue subgraph of G_n/Z_2 approaches the complete graph as n grows
  (E_blue/E_total -> 1). A complete undirected graph on V vertices is
  the underlying graph of a tournament on V vertices.

  QUESTION: How does the blue subgraph approach completeness, and
  how does this line up to tournaments of different sizes?

ANALYSES:
  1. Blue density (E_blue / C(V_merged, 2)) — how complete is the blue graph?
  2. H-oriented meta-tournament: orient blue edges by H-gradient to get a
     directed graph. Is it a tournament? What are its properties?
  3. Score sequence of the meta-tournament vs actual tournaments
  4. V_merged values vs tournament count sequence A000568
  5. The missing edges: how many blue pairs are NOT connected?
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations, combinations
from collections import defaultdict, Counter, deque
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  THE BLUE SUBGRAPH AS AN EMERGING TOURNAMENT")
print("  kind-pasteur-2026-03-23-S20cr")
print("=" * 80)

# ============================================================================
# TOURNAMENT HELPERS
# ============================================================================
def tadj(n, bits):
    a = [[0]*n for _ in range(n)]; idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): a[i][j] = 1
            else: a[j][i] = 1
            idx += 1
    return a

def Hdp(a, n):
    dp = [0]*((1<<n)*n)
    for v in range(n): dp[(1<<v)*n+v] = 1
    for S in range(1, 1<<n):
        if bin(S).count('1') >= n: continue
        for v in range(n):
            if not(S&(1<<v)): continue
            val = dp[S*n+v]
            if val == 0: continue
            for u in range(n):
                if S&(1<<u): continue
                if a[v][u]: dp[(S|(1<<u))*n+u] += val
    return sum(dp[((1<<n)-1)*n+v] for v in range(n))

def canon(a, n):
    sc = [sum(a[i][j] for j in range(n)) for i in range(n)]
    sg = defaultdict(list)
    for v in range(n): sg[sc[v]].append(v)
    gs = [sg[s] for s in sorted(set(sc))]
    if all(len(g)==1 for g in gs):
        p = [g[0] for g in gs]
        return tuple(a[p[i]][p[j]] for i in range(n) for j in range(n))
    best = None
    def gp(gs):
        if not gs: yield []; return
        for p in permutations(gs[0]):
            for r in gp(gs[1:]): yield list(p)+r
    for p in gp(gs):
        f = tuple(a[p[i]][p[j]] for i in range(n) for j in range(n))
        if best is None or f < best: best = f
    return best

def comp(a, n):
    return [[1-a[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]

def ssq(a, n):
    return tuple(sorted(sum(a[i][j] for j in range(n)) for i in range(n)))

def c3c(a, n):
    c = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if a[i][j] and a[j][k] and a[k][i]: c += 1
                if a[i][k] and a[k][j] and a[j][i]: c += 1
    return c

def build(n):
    m = comb(n,2); total = 1<<m; t0 = time.time()
    iso = defaultdict(list)
    for b in range(total):
        a = tadj(n, b); iso[canon(a, n)].append(b)
    cl = []; c2c = {}
    for idx, (cn, mem) in enumerate(sorted(iso.items())):
        a = tadj(n, mem[0])
        cl.append({'cid':idx,'canon':cn,'mem':mem,'adj':a,
                   'H':Hdp(a,n),'sc':canon(comp(a,n),n)==cn,
                   'score':ssq(a,n),'c3':c3c(a,n),
                   'aut':factorial(n)//len(mem),'size':len(mem)})
        c2c[cn] = idx
    for d in cl:
        cc = canon(comp(d['adj'],n),n); d['comp_cid'] = c2c.get(cc,-1)
    edges = set(); ecol = {}
    for d in cl:
        ci = d['cid']
        for b in d['mem']:
            for ai in range(m):
                f = b^(1<<ai); fa = tadj(n, f); fc = canon(fa, n)
                nb = c2c.get(fc)
                if nb is not None and nb != ci:
                    e = (min(ci,nb),max(ci,nb)); edges.add(e)
                    if e not in ecol:
                        ecol[e] = 'blue' if d['sc']==cl[nb]['sc'] else 'black'
    mi = {}; mid = 0
    for d in cl:
        ci = d['cid']
        if ci in mi: continue
        cp = d['comp_cid']; mi[ci] = mid
        if cp != ci and cp >= 0: mi[cp] = mid
        mid += 1
    Vm = mid
    me = set(); mec = {}
    for (a,b) in edges:
        ma, mb = mi[a], mi[b]
        if ma != mb:
            e = (min(ma,mb),max(ma,mb)); me.add(e)
            c = ecol.get((min(a,b),max(a,b)),'blue')
            if e not in mec: mec[e] = c
            elif mec[e] != c: mec[e] = 'mixed'
    mc = []; seen = set()
    for d in cl:
        mv = mi[d['cid']]
        if mv not in seen:
            seen.add(mv)
            mc.append({'mid':mv,'H':d['H'],'sc':d['sc'],'score':d['score'],
                       'c3':d['c3'],'aut':d['aut'],'size':d['size']})
    mc.sort(key=lambda x: x['mid'])
    # Build adjacency
    Ab = np.zeros((Vm,Vm),dtype=int); Ak = np.zeros((Vm,Vm),dtype=int)
    for e, c in mec.items():
        a,b = e
        if c in ('blue','mixed'): Ab[a][b]=1; Ab[b][a]=1
        if c in ('black','mixed'): Ak[a][b]=1; Ak[b][a]=1
    print(f"  Built n={n}: V_merged={Vm}, E_total={len(me)}, "
          f"E_blue={sum(1 for c in mec.values() if c=='blue')}, "
          f"E_black={sum(1 for c in mec.values() if c=='black')} [{time.time()-t0:.1f}s]")
    return Vm, mc, me, mec, Ab, Ak


# ============================================================================
# ANALYSIS
# ============================================================================

all_data = {}

for n in [3, 4, 5, 6]:
    print(f"\n{'#'*80}")
    print(f"  n = {n}")
    print(f"{'#'*80}")

    Vm, mc, me, mec, Ab, Ak = build(n)
    h_map = {d['mid']: d['H'] for d in mc}
    sc_map = {d['mid']: d['sc'] for d in mc}

    E_blue = sum(1 for c in mec.values() if c == 'blue')
    E_black = sum(1 for c in mec.values() if c == 'black')
    E_total = len(me)
    max_edges = comb(Vm, 2)  # complete graph

    # ================================================================
    # 1. APPROACH TO COMPLETENESS
    # ================================================================
    print(f"\n  --- APPROACH TO COMPLETENESS ---")
    print(f"    V_merged = {Vm}")
    print(f"    E_blue = {E_blue}")
    print(f"    C(V,2) = {max_edges}")
    blue_density = E_blue / max_edges if max_edges > 0 else 0
    total_density = E_total / max_edges if max_edges > 0 else 0
    blue_frac = E_blue / E_total if E_total > 0 else 0
    print(f"    Blue density E_blue/C(V,2) = {blue_density:.4f}")
    print(f"    Total density E_total/C(V,2) = {total_density:.4f}")
    print(f"    Blue fraction E_blue/E_total = {blue_frac:.4f}")
    print(f"    Missing blue edges = {max_edges - E_blue}")
    print(f"    Avg blue degree = {2*E_blue/Vm:.2f} (vs complete: {Vm-1})")

    # ================================================================
    # 2. H-ORIENTED META-TOURNAMENT
    # ================================================================
    print(f"\n  --- H-ORIENTED META-TOURNAMENT ---")

    # Orient blue edges by H: higher H "beats" lower H
    # This creates a directed graph. Count:
    # - How many pairs have DISTINCT H (can be oriented)
    # - How many pairs have SAME H (level edges, can't orient)
    # - Is it a tournament (all pairs oriented)?

    oriented = 0; level = 0; missing = 0
    meta_adj = [[0]*Vm for _ in range(Vm)]  # directed: meta_adj[i][j]=1 means i beats j (H_i > H_j)

    for i in range(Vm):
        for j in range(i+1, Vm):
            if Ab[i][j]:  # blue edge exists
                hi, hj = h_map[i], h_map[j]
                if hi > hj:
                    meta_adj[i][j] = 1; oriented += 1
                elif hi < hj:
                    meta_adj[j][i] = 1; oriented += 1
                else:
                    level += 1
            else:
                missing += 1

    print(f"    Oriented edges: {oriented}")
    print(f"    Level edges (same H, blue): {level}")
    print(f"    Missing edges (no blue): {missing}")
    print(f"    Total pairs: {max_edges} = {oriented} + {level} + {missing}")
    print(f"    Tournament completeness: {oriented}/{max_edges} = {oriented/max_edges:.4f}" if max_edges > 0 else "")

    # Is it transitive? (DAG = transitive tournament)
    # Check: are there any 3-cycles?
    cycles_3 = 0
    for i in range(Vm):
        for j in range(Vm):
            if meta_adj[i][j] == 0: continue
            for k in range(Vm):
                if meta_adj[j][k] and meta_adj[k][i]:
                    cycles_3 += 1
    cycles_3 //= 3  # each cycle counted 3 times
    print(f"    Directed 3-cycles: {cycles_3}")
    print(f"    {'TRANSITIVE (DAG)' if cycles_3 == 0 else 'NOT transitive'}")

    # Score sequence of meta-tournament
    scores = sorted([sum(meta_adj[i][j] for j in range(Vm)) for i in range(Vm)])
    print(f"    Meta-tournament scores: min={scores[0]}, max={scores[-1]}, "
          f"median={scores[Vm//2]}")
    if Vm <= 20:
        print(f"    Full score sequence: {scores}")

    # Is it regular? (constant out-degree)
    score_set = set(scores)
    print(f"    Distinct scores: {len(score_set)}")
    if len(score_set) == 1:
        print(f"    REGULAR meta-tournament!")
    else:
        print(f"    Score variance: {np.var(scores):.2f}")

    # ================================================================
    # 3. COMPARISON WITH ACTUAL TOURNAMENTS ON V_merged VERTICES
    # ================================================================
    print(f"\n  --- COMPARISON WITH TOURNAMENTS ON {Vm} VERTICES ---")

    # A000568 values (non-isomorphic tournaments)
    a000568 = {1:1, 2:1, 3:2, 4:4, 5:12, 6:56, 7:456, 8:6880, 9:191536, 10:9733056}

    if Vm in a000568:
        print(f"    A000568({Vm}) = {a000568[Vm]} non-iso tournaments on {Vm} vertices")
    else:
        print(f"    A000568({Vm}) not precomputed")

    # The meta-tournament (H-oriented blue graph) is one specific "partial tournament"
    # on V_merged vertices. What fraction of edges does it have?
    t_edges = comb(Vm, 2)  # complete tournament has this many edges
    meta_edges = oriented  # how many we actually have
    print(f"    Complete tournament: {t_edges} directed edges")
    print(f"    Meta-tournament: {meta_edges} directed edges ({meta_edges/t_edges:.1%} complete)")

    # ================================================================
    # 4. THE MISSING EDGES: WHO IS NOT BLUE-CONNECTED?
    # ================================================================
    print(f"\n  --- MISSING BLUE EDGES (non-adjacent in blue) ---")

    non_adj_blue = []
    for i in range(Vm):
        for j in range(i+1, Vm):
            if Ab[i][j] == 0:
                non_adj_blue.append((i, j))

    print(f"    Missing blue pairs: {len(non_adj_blue)}")

    if non_adj_blue and Vm <= 40:
        # Classify missing pairs by type
        sc_sc_miss = 0; ns_ns_miss = 0; sc_ns_miss = 0
        for i, j in non_adj_blue:
            if sc_map[i] and sc_map[j]: sc_sc_miss += 1
            elif not sc_map[i] and not sc_map[j]: ns_ns_miss += 1
            else: sc_ns_miss += 1

        print(f"    Missing SC-SC: {sc_sc_miss}")
        print(f"    Missing NS-NS: {ns_ns_miss}")
        print(f"    Missing SC-NS: {sc_ns_miss} (these are impossible — would be black)")

        # What's the H-distance between missing pairs?
        if non_adj_blue:
            h_diffs = [abs(h_map[i] - h_map[j]) for i, j in non_adj_blue if sc_map[i] == sc_map[j]]
            if h_diffs:
                print(f"    Same-type missing pairs: {len(h_diffs)}")
                print(f"    H-difference: min={min(h_diffs)}, max={max(h_diffs)}, avg={np.mean(h_diffs):.1f}")

        # Compare with actual blue pairs
        adj_blue = [(i,j) for i in range(Vm) for j in range(i+1,Vm) if Ab[i][j]]
        if adj_blue:
            h_diffs_adj = [abs(h_map[i] - h_map[j]) for i, j in adj_blue]
            print(f"    Blue-connected pairs: {len(adj_blue)}")
            print(f"    H-difference: min={min(h_diffs_adj)}, max={max(h_diffs_adj)}, avg={np.mean(h_diffs_adj):.1f}")

    # ================================================================
    # 5. BLUE SUBGRAPH AS EVOLVING GRAPH SEQUENCE
    # ================================================================
    print(f"\n  --- BLUE GRAPH PROPERTIES ---")

    # Compute blue graph properties that are tournament-relevant
    # Independence number (tournament analog: independent set in underlying graph)
    # Clique number (tournament analog: transitive sub-tournament in underlying graph)

    # Blue degree sequence
    blue_deg = [int(np.sum(Ab[i])) for i in range(Vm)]
    print(f"    Blue degree sequence: min={min(blue_deg)}, max={max(blue_deg)}, avg={np.mean(blue_deg):.1f}")

    # Vertex types by blue degree
    sc_deg = [int(np.sum(Ab[i])) for i in range(Vm) if sc_map[i]]
    ns_deg = [int(np.sum(Ab[i])) for i in range(Vm) if not sc_map[i]]

    if sc_deg:
        print(f"    SC blue degree: avg={np.mean(sc_deg):.1f}, range [{min(sc_deg)},{max(sc_deg)}]")
    if ns_deg:
        print(f"    NS blue degree: avg={np.mean(ns_deg):.1f}, range [{min(ns_deg)},{max(ns_deg)}]")

    all_data[n] = {
        'V': Vm, 'E_blue': E_blue, 'E_total': E_total,
        'density_blue': blue_density, 'density_total': total_density,
        'blue_frac': blue_frac, 'oriented': oriented,
        'level': level, 'missing': missing, 'cycles_3': cycles_3,
        'scores': scores
    }


# ============================================================================
# CROSS-n SYNTHESIS
# ============================================================================

print(f"\n\n{'='*80}")
print("  CROSS-n SYNTHESIS: APPROACH TO TOURNAMENT COMPLETENESS")
print(f"{'='*80}")

# Known values for n=7,8 from opus data
known_extra = {
    7: {'V': 272, 'E_blue': 1573, 'E_total': 2123, 'SC': 88, 'NS': 184},
    8: {'V': 3528, 'E_blue': 43656, 'E_total': 45550, 'SC': 176, 'NS': 3352}
}

print(f"\n  {'n':>3} {'V_m':>6} {'E_blue':>8} {'C(V,2)':>10} {'blu_dens':>10} {'tot_dens':>10} "
      f"{'blu_frac':>10} {'avg_bdeg':>10} {'V-1':>6}")
print(f"  {'-'*80}")

for n in [3, 4, 5, 6]:
    d = all_data[n]
    V = d['V']
    max_e = comb(V, 2)
    avg_bdeg = 2*d['E_blue']/V if V > 0 else 0
    print(f"  {n:3d} {V:6d} {d['E_blue']:8d} {max_e:10d} {d['density_blue']:10.4f} "
          f"{d['density_total']:10.4f} {d['blue_frac']:10.4f} {avg_bdeg:10.2f} {V-1:6d}")

for n in [7, 8]:
    d = known_extra[n]
    V = d['V']; Eb = d['E_blue']; Et = d['E_total']
    max_e = comb(V, 2)
    avg_bdeg = 2*Eb/V
    bd = Eb/max_e; td = Et/max_e; bf = Eb/Et
    print(f"  {n:3d} {V:6d} {Eb:8d} {max_e:10d} {bd:10.6f} "
          f"{td:10.6f} {bf:10.4f} {avg_bdeg:10.2f} {V-1:6d}")

# ============================================================================
# KEY INSIGHT: THE TWO SCALES
# ============================================================================

print(f"\n\n{'='*80}")
print("  KEY INSIGHT: TWO RATES OF APPROACH")
print(f"{'='*80}")
print("""
  There are TWO senses in which the blue subgraph approaches a "tournament":

  (A) FRACTION: E_blue/E_total -> 1
      Among edges that EXIST in G_n/Z_2, almost all are blue.
      Sequence: 1.000, 0.333, 0.619, 0.685, 0.741, 0.958
      This converges to 1 because SC classes become negligible (NS/V -> 1)
      and NS nodes connect almost exclusively by blue edges.

  (B) DENSITY: E_blue/C(V,2) -> 0
      The blue graph is SPARSE relative to the complete graph.
      Sequence: 1.000, 0.333, 0.289, 0.175, 0.043, 0.007
      This decreases because V_merged grows much faster than E_blue.

  So the blue graph simultaneously becomes:
  - MORE DOMINANT (as a fraction of all meta-edges)
  - MORE SPARSE (relative to completeness)

  The AVG BLUE DEGREE tells the real story:
  Sequence: 1.00, 0.67, 2.60, 5.76, 11.57, 24.74
  This is INCREASING, roughly doubling each step.
  The ratio avg_degree/(V-1) = blue_density -> 0, but slowly.
""")

# ============================================================================
# THE META-TOURNAMENT INTERPRETATION
# ============================================================================

print(f"\n{'='*80}")
print("  THE META-TOURNAMENT: H-ORIENTED BLUE GRAPH")
print(f"{'='*80}")
print("""
  Since G_n/Z_2 is a DAG under H (0 downhill edges), the H-oriented
  blue subgraph is a TRANSITIVE partial tournament.

  A transitive tournament on V vertices has exactly ONE Hamiltonian path.
  The meta-tournament has:
    - V_merged vertices (iso classes modulo complement)
    - Edges oriented by H: higher H "wins"
    - Missing edges where blue doesn't connect

  This partial tournament is SPARSE but carries key information:
  - Its score sequence is the H-rank distribution
  - Its connected components are the blue components (SC + NS separate)
  - Its reachability (transitive closure) defines the H-dominance order

  THE DEEP ANALOGY:
  At each n, the meta-graph G_n/Z_2 is itself a "tournament landscape"
  describing how n-vertex tournaments relate to each other via arc flips.
  As n grows, this landscape becomes increasingly dominated by blue
  (within-type) connections, making it more tournament-like.

  V_merged values: 2, 3, 10, 34, 272, 3528, ...
  These are the "effective sizes" of the tournament landscape.
  The sequence grows as ~A000568(n)/2 = 2^{C(n,2)-1}/n!
""")

# ============================================================================
# COMPARISON: V_merged vs TOURNAMENT COUNTS
# ============================================================================

print(f"\n{'='*80}")
print("  V_MERGED vs A000568 (TOURNAMENT COUNTS)")
print(f"{'='*80}")

a000568 = {2:1, 3:2, 4:4, 5:12, 6:56, 7:456, 8:6880, 9:191536}
sc_count = {3:2, 4:2, 5:8, 6:12, 7:88, 8:176}  # from project data

print(f"\n  {'n':>3} {'A000568(n)':>12} {'SC(n)':>8} {'V_merged':>10} {'A000568/2':>12} {'ratio':>8}")
print(f"  {'-'*60}")
for n in range(3, 9):
    a = a000568.get(n, '?')
    s = sc_count.get(n, '?')
    if isinstance(a, int) and isinstance(s, int):
        vm = (a + s) // 2
        a2 = a / 2
        ratio = vm / a2 if a2 > 0 else '?'
        print(f"  {n:3d} {a:12d} {s:8d} {vm:10d} {a2:12.1f} {ratio:8.4f}")
    else:
        print(f"  {n:3d} {'?':>12} {'?':>8}")

print("""
  V_merged = (A000568 + SC) / 2 converges to A000568/2 as SC/A000568 -> 0.
  SC fraction: 1.00, 0.50, 0.67, 0.21, 0.19, 0.026
  SC becomes negligible for large n, so V_merged ~ A000568(n)/2.

  The "tournament of tournaments" has half the classes (complement identification).
""")

# ============================================================================
# TOURNAMENTS THAT MATCH V_MERGED
# ============================================================================

print(f"\n{'='*80}")
print("  WHICH TOURNAMENT SIZE HAS ~V_MERGED VERTICES?")
print(f"{'='*80}")

print(f"\n  For each n, what k has A000568(k) closest to V_merged(n)?")
for n in range(3, 9):
    a = a000568.get(n, None)
    s = sc_count.get(n, None)
    if a is None or s is None: continue
    vm = (a + s) // 2
    # Find k such that A000568(k) is closest to vm
    best_k = None; best_diff = float('inf')
    for k, ak in a000568.items():
        diff = abs(ak - vm)
        if diff < best_diff:
            best_diff = diff; best_k = k
    print(f"  n={n}: V_merged={vm}, closest A000568({best_k})={a000568[best_k]} (diff={best_diff})")


print(f"\n  DONE.")
print("=" * 80)
