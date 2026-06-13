#!/usr/bin/env python3
"""
lrc_repeated_addition_hyperoperation_s548.py    oracle-2026-06-01-S548o

MULTIPLICATION = REPEATED ADDITION, taken wildly + recursively.

THE SEED. A runner at speed k, time t, sits at k*t = t+t+...+t (k times). So:
   r_0 = 0 (observer),   r_k = r_{k-1} + t.
The runner system at AP speeds {1..N} IS the REPEATED-ADDITION ORBIT {k t} of t -- a
single rotation orbit. 'Multiplication k*t' is literally 't added k times.'

THE HYPEROPERATION LADDER (the recursive synthesis):
   level 0 succ (+1) -> level 1 ADDITION (repeated succ) = the walk step t
   -> level 2 MULTIPLICATION (repeated addition) = the runner position k*t; AP speeds
      = consecutive multiples = repeated addition
   -> level 3 EXPONENTIATION (repeated multiplication) = GEOMETRIC speeds r^k (lacunary)
LRC lives at level 2 (positions k*t); the cascade PRODUCT of clearances (S545) is
level 3 (repeated x) over the runners (level 2); threshold 1/n & doubling x2 are level-2.

TWO PAYOFFS (computed):
 (1) THREE-GAP RIGIDITY is the signature of repeated addition. AP speeds {1..N} = a
     single rotation orbit {kt} -> at most 3 distinct gaps (Steinhaus), recursively
     from the continued fraction of t (= repeated SUBTRACTION, the Euclidean inverse
     of repeated addition). GENERAL speeds are NOT a single repeated-addition orbit ->
     they break three-gap (S538). So the additive (AP) structure is exactly what gives
     the rigid recursive gap structure.
 (2) THE HYPEROPERATION LEVEL of the speed set controls LRC spread/entropy (S543):
     level-2 (AP, repeated +) -> even orbit, HIGH H-entropy (the tight regular polygon);
     level-3 (geometric, repeated x) -> lacunary, LOW H-entropy. The recursive ladder
     of speed construction = the LRC difficulty spectrum.
"""
from itertools import combinations
from functools import reduce
from math import gcd, log2
import random, statistics as st

def frac(x): return x - int(x // 1)

# ---------- (A) repeated-addition orbit ----------
def part_A():
    print("="*70); print("(A) RUNNERS = the repeated-addition orbit: r_k = r_{k-1} + t = k*t"); print("="*70)
    t = 0.31831  # ~1/pi
    r = 0.0; seq = []
    for k in range(8):
        seq.append(round(r % 1.0, 4)); r += t
    print(f"  t={t}: r_k (= t added k times) = {seq}")
    print(f"  check k*t mod 1: {[round((k*t)%1.0,4) for k in range(8)]}  (identical)")
    print("  => 'multiplication k*t' IS 't repeated-added k times'. AP speeds {1..N} = the orbit.")
    print()

# ---------- (B) three-gap rigidity = repeated-addition signature ----------
def distinct_gaps(points):
    pts = sorted(p % 1.0 for p in points)
    gaps = [round(pts[(i+1) % len(pts)] - pts[i], 6) % 1.0 for i in range(len(pts)-1)]
    gaps.append(round(1 - pts[-1] + pts[0], 6))
    return sorted(set(round(g, 5) for g in gaps))

def part_B():
    print("="*70); print("(B) THREE-GAP rigidity: AP {1..N} (repeated-add orbit) <=3 gaps; general breaks it"); print("="*70)
    rnd = random.Random(7)
    for N in (6, 9, 12):
        t = rnd.random()
        ap_pts = [0.0] + [frac(k*t) for k in range(1, N)]        # AP speeds 1..N-1 = single orbit
        gen = tuple(sorted(rnd.sample(range(1, 8*N), N-1)))
        gen_pts = [0.0] + [frac(v*t) for v in gen]
        ng_ap = len(distinct_gaps(ap_pts)); ng_gen = len(distinct_gaps(gen_pts))
        print(f"  N={N}, t={t:.4f}: AP {{1..{N-1}}} orbit -> {ng_ap} distinct gaps (<=3: {ng_ap<=3}); "
              f"general speeds {gen} -> {ng_gen} distinct gaps")
    print("  => the SINGLE repeated-addition orbit (AP) obeys three-gap (<=3, recursive CF");
    print("     structure); general speeds are NOT one additive orbit -> >3 gaps (S538). The")
    print("     additive (x=repeated+) structure is EXACTLY the three-gap rigidity.")
    print()

# ---------- (C) hyperoperation level of speeds vs LRC entropy ----------
def H_count(adj, n):
    full=(1<<n)-1; dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        for v in range(n):
            c=dp[mask][v]
            if not c: continue
            for u in range(n):
                if mask&(1<<u): continue
                if adj[v][u]: dp[mask|(1<<u)][u]+=c
    return sum(dp[full][v] for v in range(n))
def half_turn(pos,n):
    adj=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i!=j and 0<(pos[j]-pos[i])%1.0<0.5: adj[i][j]=1
    return adj
def mean_SH(v, G=3000):
    n=len(v); s=0.0; c=0
    for k in range(G):
        t=(k+0.5)/G; H=H_count(half_turn([frac(x*t) for x in v],n),n)
        if H>0: s+=log2(H); c+=1
    return s/c if c else 0.0

def part_C():
    print("="*70); print("(C) HYPEROPERATION LEVEL of speeds vs LRC mean H-entropy (S543 spectrum)"); print("="*70)
    n=6
    fams={
        "level-2 ADD: AP 1,2,3,4,5,6 (repeated +)": tuple(range(1,7)),
        "level-2 ADD: Fibonacci 1,2,3,5,8,13 (s+s)": (1,2,3,5,8,13),
        "level-3 EXP: geometric 1,2,4,8,16,32 (x2)": tuple(2**k for k in range(6)),
        "level-3 EXP: geometric 1,3,9,27,81,243 (x3)": tuple(3**k for k in range(6)),
    }
    for name,v in fams.items():
        v=tuple(sorted(set(v)))
        if reduce(gcd,v)!=1:
            g=reduce(gcd,v); v=tuple(x//g for x in v)
        print(f"  {name:46s} v={v}: mean H-entropy={mean_SH(v):.3f}")
    print("  => ADDITIVE (AP / repeated +) speeds -> HIGH mean H-entropy (even orbit, the tight")
    print("     regular polygon); EXPONENTIAL (geometric / repeated x) speeds -> LOW (lacunary).")
    print("     The hyperoperation LEVEL of the speed recursion = the LRC difficulty/entropy spectrum.")
    print()

def main():
    part_A(); part_B(); part_C()
    print("="*70)
    print("RECURSIVE SYNTHESIS: x = repeated + builds the hyperoperation ladder")
    print("succ -> + -> x -> ^. The runner position k*t is level-2 (t repeated-added k times);")
    print("AP speeds = the pure repeated-addition orbit = three-gap-rigid (CF recursion =")
    print("repeated SUBTRACTION, the Euclidean inverse) = HIGH entropy = the tight regular")
    print("polygon. The cascade PRODUCT of clearances (S545) is level-3 (repeated x) over the")
    print("runners. The hyperoperation LEVEL at which the speeds are generated sets the LRC")
    print("regime: level-2 additive = even/rigid/tight; level-3 multiplicative = lacunary/loose.")
    print("LRC = whether the level-2 repeated-addition orbit's apex gap clears 1/n.")
    print("="*70)

if __name__=="__main__":
    main()
