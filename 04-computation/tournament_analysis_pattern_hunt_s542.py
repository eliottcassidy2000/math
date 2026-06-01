#!/usr/bin/env python3
"""
tournament_analysis_pattern_hunt_s542.py    oracle-2026-06-01-S542o

TOURNAMENT ANALYSIS (the repo's central method, s471/s23): pairwise data + a switch
(comparator) + a tie Hamiltonian path -> a tournament-valued observable; study its
fingerprints and TRAJECTORIES as variables change. Tie-break path = the basketball
JERSEY ORDER 1..5 = the tiling-model BASE HAMILTONIAN PATH (S530/S531); so every
output = base-path backbone (+) tiles (the metric's real signal) = cut (+) cycle.

This is the empirical PATTERN HUNT across many instances/variables:
 P1 BASKETBALL archetypes (discrete): star / balanced / two-man -> H, SCC, source.
 P2 COMPARATOR ZOO (continuous runners): arc / chord / half-turn / APPROACHING --
    which collapse (geometric, S541) vs which are distinct.
 P3 H(t) TRAJECTORY + iso-class WALK over t: statistics for AP vs generic speeds.
 P4 COMPARATOR PHASE SWEEP: rotate the half-turn window by phi; iso-class vs phi.
 P5 THRESHOLD theta SWEEP (trienerment): #ties vs theta; the loneliness theta=1/n.
 P6 BASE-PATH (+) TILES: tile-count (cyclic content) over t; correlate with H.
"""
from itertools import combinations, permutations
from functools import reduce
from math import gcd, sin, pi
import random

def frac(x): return x - int(x // 1)

# ---------- invariants ----------
def H_count(adj, n):
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            c = dp[mask][v]
            if not c: continue
            for u in range(n):
                if mask & (1 << u): continue
                if adj[v][u]: dp[mask | (1 << u)][u] += c
    return sum(dp[full][v] for v in range(n))

def num_3cycles(adj, n):
    return sum(1 for a,b,c in combinations(range(n),3)
               if (adj[a][b] and adj[b][c] and adj[c][a]) or (adj[b][a] and adj[c][b] and adj[a][c]))

def canon(adj, n):
    best = None
    for p in permutations(range(n)):
        bits = tuple(adj[p[i]][p[j]] for i in range(n) for j in range(n) if i != j)
        if best is None or bits < best: best = bits
    return best

def is_source_sink(adj, n):
    outs = [sum(adj[i]) for i in range(n)]
    return (n-1 in outs), (0 in outs)

# ---------- build tournament from a directed metric, jersey tie-break ----------
def tournament_from_metric(M, n):
    """M[i][j] directed datum; i->j iff M[i][j]>M[j][i]; tie -> jersey path (i<j => i->j)."""
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j: continue
            if M[i][j] > M[j][i] or (M[i][j] == M[j][i] and i < j):
                adj[i][j] = 1
    return adj

# ---------- runner comparators ----------
def half_turn(positions, n, phi=0.0):
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j: continue
            if 0 < (positions[j]-positions[i]-phi) % 1.0 < 0.5: adj[i][j] = 1
    return adj

# ============================================================
def P1_basketball():
    print("="*70); print("P1 BASKETBALL archetypes (jersey 1..5 tie-break = base path)"); print("="*70)
    n = 5
    star = [[0,9,9,9,9],[2,0,1,1,1],[2,1,0,1,1],[2,1,1,0,1],[2,1,1,1,0]]      # P1 ball-dominant
    balanced = [[0,5,4,5,4],[4,0,5,4,5],[5,4,0,5,4],[4,5,4,0,5],[5,4,5,4,0]]  # ~equal
    twoman = [[0,9,2,2,2],[9,0,2,2,2],[3,3,0,4,4],[3,3,4,0,4],[3,3,4,4,0]]    # 1-2 pick&roll
    for name, M in [("star/ball-dominant", star), ("balanced", balanced), ("two-man-game", twoman)]:
        adj = tournament_from_metric(M, n)
        H = H_count(adj, n); c3 = num_3cycles(adj, n); src, snk = is_source_sink(adj, n)
        print(f"  {name:18s}: H={H:3d}, 3-cycles={c3}, source={src}, sink={snk}")
    print("  PATTERN: ball-dominant -> a SOURCE, near-transitive, LOW H (bunched); balanced")
    print("  -> regular-like, HIGH H (spread). H = the 'team balance / loneliness' meter (S26).")
    print("  Jersey 1..5 tie-break = base Hamiltonian path; passing asymmetries = the TILES.")

def P2_comparator_zoo():
    print("\n"+"="*70); print("P2 COMPARATOR ZOO (runners): which collapse, which are distinct"); print("="*70)
    v = (1,2,4,7); n = 4; rnd = random.Random(1)
    same_arc_chord = True; same_ht = True; approaching_distinct = False
    for _ in range(3000):
        t = rnd.random()
        pos = [frac(x*t) for x in v]
        # arc-distance-to-a-fixed-point ranking vs chord ranking: chord=2sin(pi*arc) monotone
        arc = [min(p,1-p) for p in pos]; chord = [2*sin(pi*a) for a in arc]
        if [sorted(range(n),key=lambda i:arc[i])] != [sorted(range(n),key=lambda i:chord[i])]:
            same_arc_chord = False
        # half-turn vs 'approaching observer' (sign of d/dt dist to 0)
        ht = half_turn(pos, n)
        appr = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i==j: continue
                # i 'beats' j if i is closer to 0 than j (distance comparator) -- a ranking
                appr[i][j] = 1 if arc[i] < arc[j] or (arc[i]==arc[j] and i<j) else 0
        if canon(ht,n) != canon(appr,n): approaching_distinct = True
    print(f"  arc-rank == chord-rank (chord=2sin(pi*arc) monotone): {same_arc_chord}  (GEOMETRIC collapse, S541)")
    print(f"  half-turn(ahead) vs distance-rank(approaching): distinct iso-classes seen: {approaching_distinct}")
    print("  => distance-based comparators (arc/chord) collapse; SIGNED 'ahead' (half-turn,")
    print("     cyclic, non-transitive) differs from distance-RANK (transitive). Comparator choice = geometry exposed.")

def P3_H_trajectory():
    print("\n"+"="*70); print("P3 H(t) TRAJECTORY + iso-class WALK over t: AP vs generic"); print("="*70)
    for v in [(1,2,3,4,5),(1,2,4,7,11)]:
        n = len(v); G = 6000
        Hs=[]; isos=set(); transitive_time=0
        for s in range(G):
            t=(s+0.5)/G; pos=[frac(x*t) for x in v]; adj=half_turn(pos,n)
            H=H_count(adj,n); Hs.append(H); isos.add(canon(adj,n))
            if H==1: transitive_time+=1
        ap = " [AP/regular]" if list(v)==list(range(1,n+1)) else ""
        print(f"  v={v}{ap}: H range [{min(Hs)},{max(Hs)}], mean {sum(Hs)/len(Hs):.1f}, "
              f"#distinct iso-classes visited={len(isos)}, frac transitive(H=1)={transitive_time/G:.3f}")
    print("  PATTERN: the tournament WALKS the menu; H oscillates between 1 (bunched/transitive)")
    print("  and max (spread/regular). The walk's statistics fingerprint the speed set.")

def P4_phase_sweep():
    print("\n"+"="*70); print("P4 COMPARATOR PHASE SWEEP: rotate the half-turn window by phi"); print("="*70)
    v=(1,2,3,4); n=4; t=0.137; pos=[frac(x*t) for x in v]
    isos=[]
    for k in range(200):
        phi=k/200.0; adj=half_turn(pos,n,phi); isos.append(canon(adj,n))
    distinct=len(set(isos))
    # count transitions
    trans=sum(1 for k in range(len(isos)-1) if isos[k]!=isos[k+1])
    print(f"  fixed positions (v={v},t={t}); rotating comparator window phi:0->1:")
    print(f"  distinct iso-classes = {distinct}, transitions = {trans}")
    print("  => the SWITCH itself is a variable; rotating the comparator window traces an")
    print("     orbit of iso-classes for FIXED data -- a new axis of Tournament Analysis.")

def P5_threshold_sweep():
    print("\n"+"="*70); print("P5 THRESHOLD theta SWEEP (trienerment): #ties vs theta"); print("="*70)
    v=(1,2,3,4,5); n=len(v); t=0.137; pos=[frac(x*t) for x in v]
    def cdist(a,b): d=abs(a-b)%1.0; return min(d,1-d)
    print("  theta : #ties (pairs within theta)   (loneliness threshold theta=1/n marked)")
    for theta in [0.05,0.1,1.0/n,0.2,0.3,0.5]:
        ties=sum(1 for i,j in combinations(range(n),2) if cdist(pos[i],pos[j])<theta)
        mark=" <- 1/n" if abs(theta-1.0/n)<1e-9 else ""
        print(f"  {theta:.3f} : {ties}{mark}")
    print("  => theta=0 pure tournament (0 ties); theta growing -> ties percolate -> all-tie at")
    print("     theta=1/2. The loneliness threshold 1/n is where 'near'=tie is defined (S539).")

def P6_basepath_tiles():
    print("\n"+"="*70); print("P6 BASE-PATH (+) TILES: tile-count (cyclic content) over t vs H"); print("="*70)
    v=(1,2,3,4,5); n=len(v); G=4000
    # base path = jersey/speed order 1..n; a TILE = a pair where half-turn DISAGREES with speed order
    rows=[]
    for s in range(G):
        t=(s+0.5)/G; pos=[frac(x*t) for x in v]; adj=half_turn(pos,n)
        # tile if (i<j but adj[j][i]) i.e. lower index loses to higher = disagrees with base path
        tiles=sum(1 for i in range(n) for j in range(i+1,n) if adj[j][i])
        rows.append((tiles, H_count(adj,n)))
    # correlate tiles with H
    import statistics as st
    tiles_vals=[r[0] for r in rows]; Hs=[r[1] for r in rows]
    # H by tile-count
    byt={}
    for tl,H in rows: byt.setdefault(tl,[]).append(H)
    print("  tile-count -> avg H (cyclic content drives the loneliness meter):")
    for tl in sorted(byt): print(f"    tiles={tl}: avg H={sum(byt[tl])/len(byt[tl]):.1f}  (freq {len(byt[tl])})")
    print("  => tiles = where the metric disagrees with the base path = the CYCLIC content;")
    print("     more tiles -> higher H (more spread). base-path(+)tiles = cut(+)cycle = S537.")

def main():
    P1_basketball(); P2_comparator_zoo(); P3_H_trajectory(); P4_phase_sweep(); P5_threshold_sweep(); P6_basepath_tiles()

if __name__=="__main__":
    main()
