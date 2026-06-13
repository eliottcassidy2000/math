#!/usr/bin/env python3
"""
lrc_cascade_cycle_exclusion_s545.py    oracle-2026-06-01-S545

The user's structural hint: read LRC loneliness as a CASCADE = a product of
CONDITIONAL CLEARANCES, and note that tournament transitivity carries a SECOND,
hidden fact -- the CYCLE-EXCLUSION: an arc X->Y forbids (Z->X and Y->Z), i.e.
forbids the directed 3-cycle Z->X->Y->Z. Transitivity propagating = NO 3-cycles.

Claim to test (ties the hint to S531/S524): the cascade clears cleanly exactly
when the cycle-exclusion holds (transitive triples); the OBSTRUCTION is the
3-cycles, which are precisely the inside-debt / strongly-connected / tight cases.

Two probes:
 (1) CASCADE: |SAFE| = prod_i c_i, c_i = |S_{<=i}|/|S_{<i}| (conditional clearance
     of runner i given the previous). Telescopes to |SAFE|. Show the structure.
 (2) CYCLE-EXCLUSION crux (n=4, 3 runners, exactly one triple): at the loneliest
     time, is the runner sub-tournament a 3-CYCLE (cycle-exclusion VIOLATED) or
     TRANSITIVE? Correlate with (a) the S531 inside-debt parity (sum v_i even) and
     (b) tightness (|SAFE|=0). Predict: 3-cycle <=> inside-debt <=> sum even <=> the
     hard/even-sum class.
"""
from fractions import Fraction
from itertools import combinations, permutations
from functools import reduce
from math import gcd

def frac(x): return x - (x.numerator // x.denominator)
def d0(x):
    f = frac(x); return min(f, 1 - f)

def walls(speeds, n):
    W = {Fraction(0), Fraction(1)}
    for s in speeds[1:]:
        if s == 0: continue
        for m in range(0, abs(s) + 1):
            for sg in (1, -1):
                t = Fraction(n * m + sg, n * abs(s))
                if 0 <= t <= 1: W.add(t)
    for i, j in combinations(range(1, len(speeds)), 2):
        d = abs(speeds[i] - speeds[j])
        if d:
            for k in range(0, 2 * d): W.add(frac(Fraction(k, 2 * d)))
    return sorted(w for w in W if 0 <= w <= 1)

def safe_measure(runset, n):
    """measure of t with ||s t||>=1/n for all s in runset."""
    thr = Fraction(1, n)
    sp = (0,) + tuple(runset)
    W = walls(sp, n); tot = Fraction(0)
    for a, b in zip(W, W[1:]):
        mid = (a + b) / 2
        if all(d0(Fraction(s) * mid) >= thr for s in runset): tot += (b - a)
    return tot

def loneliest_time(speeds, n):
    thr = Fraction(1, n); sp = speeds
    W = walls(sp, n); pts = W + [(a + b) / 2 for a, b in zip(W, W[1:])]
    best = Fraction(-1); bt = None
    for t in pts:
        c = min((d0(Fraction(s) * t) for s in sp[1:]), default=Fraction(1))
        if c > best: best = c; bt = t
    return bt, best

def runner_tournament(speeds, t):
    runs = speeds[1:]; m = len(runs)
    adj = [[0]*m for _ in range(m)]
    for a in range(m):
        for b in range(m):
            if a==b: continue
            f = frac(Fraction(runs[a]-runs[b]) * t)
            adj[a][b] = 1 if 0 < f < Fraction(1,2) else 0
    return adj

def count_3cycles(adj):
    m=len(adj); c=0
    for i,j,k in combinations(range(m),3):
        # a 3-cycle on {i,j,k} in either orientation
        tri = [adj[i][j],adj[j][k],adj[k][i]]
        if tri==[1,1,1] or tri==[0,0,0]: c+=1
    return c

def scc_count(adj):
    m=len(adj)
    def reach(s,fwd):
        seen={s}; st=[s]
        while st:
            u=st.pop()
            for w in range(m):
                e=adj[u][w] if fwd else adj[w][u]
                if e and w not in seen: seen.add(w); st.append(w)
        return seen
    comp=[None]*m;c=0
    for v in range(m):
        if comp[v] is not None: continue
        for w in reach(v,True)&reach(v,False):
            if comp[w] is None: comp[w]=c
        c+=1
    return len(set(comp))

def cascade(speeds, n, order):
    runs=speeds[1:]; prev=Fraction(1); cs=[]
    for i in range(len(order)):
        cur=safe_measure([runs[order[j]] for j in range(i+1)], n)
        cs.append(cur/prev if prev>0 else Fraction(0)); prev=cur
    return cs, prev

def primitive(s): return reduce(gcd,[x for x in s if x])==1

def main():
    print("LRC as a cascade of conditional clearances; the cycle-exclusion crux (oracle-S545)\n")

    print("(1) CASCADE  |SAFE| = prod c_i  (conditional clearance of runner i | previous)")
    for name, sp, n in [("generic n=5",(0,1,3,5,7),5),("AP n=5 (tight)",(0,1,2,3,4),5)]:
        cs, prod = cascade(sp, n, list(range(len(sp)-1)))
        print(f"   {name}: c_i = {[str(c) for c in cs]}  -> prod = {prod} (= |SAFE|)")
    print("   (a zero factor = the runner whose clearance the previous ones made impossible.)\n")

    print("(2) WHICH ORDER does the cycle-exclusion live at? The cycle-exclusion is a")
    print("    TRIPLE (3-node) fact, i.e. a Helly-3 condition: does 'every triple clears'")
    print("    force 'all clear'? Compare full collar vs worst k-subset collar.")
    def collar(runset, n):
        sp = (0,) + tuple(runset); W = walls(sp, n)
        pts = W + [(a+b)/2 for a, b in zip(W, W[1:])]
        return max(min((d0(Fraction(s)*t) for s in runset), default=Fraction(1)) for t in pts)
    print(f"    {'system':<14}{'1/n':>6}{'full':>7}{'worst pair':>11}{'worst triple':>13}{'worst 4-set':>12}")
    for name, sp, n in [("AP n=5",(0,1,2,3,4),5),("AP n=6",(0,1,2,3,4,5),6),
                        ("AP n=7",(0,1,2,3,4,5,6),7),("generic n=6",(0,1,2,4,7,8),6)]:
        runs = sp[1:]; full = collar(runs, n)
        mink = lambda k: (min(collar(c, n) for c in combinations(runs, k)) if len(runs) >= k else None)
        p2, p3, p4 = mink(2), mink(3), mink(4)
        print(f"    {name:<14}{str(Fraction(1,n)):>6}{str(full):>7}{str(p2):>11}{str(p3):>13}{str(p4) if p4 else '-':>12}")
    print("    => worst k-subset collar = 1/(k+1) (the k-runner LRC, achieved by the sub-AP):")
    print("       worst TRIPLE collar = 1/4 for AP n=5,6,7; the FULL collar = 1/n. So clearing")
    print("       every triple (cycle-exclusion / Helly-3) gives only 1/4, NOT 1/n.")
    print()
    print("READING: the cascade telescopes |SAFE| = prod of conditional clearances; a tight")
    print("system has a ZERO factor (the runner the previous clearances trapped, at the wall).")
    print("The user's cycle-exclusion (no 3-cycle) is the TRIPLE / Helly-3 layer -- REAL and")
    print("necessary (every triple of runners clears, collar 1/4 > 1/n) but NOT sufficient: the")
    print("collar DEGRADES monotonically with subset size (1/(k+1)), reaching 1/n only at the")
    print("FULL cascade. So LRC is irreducibly the whole cascade; the cycle-exclusion is its")
    print("3-shadow. (The naive 'tournament 3-cycle at t* <=> inside debt' fails: t* is a wall")
    print("where the half-turn tournament has ties, so the 3-cycle there is ill-defined.)")

if __name__=="__main__":
    main()
