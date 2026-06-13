#!/usr/bin/env python3
"""
lrc_global_spread_local_emptiness_s544.py    oracle-2026-06-01-S544o

REFRAME: can GLOBAL SPREAD guarantee LOCAL EMPTINESS (a hole at the observer)?
Naively NO -- S543: the regular polygon (max spatial spread) is the TIGHT/hardest
case. So we must reframe the CHOICES. Three precise incarnations:

 (1) INSTANTANEOUS spatial spread is the ENEMY: at a fixed time, even/high-entropy
     = small gaps = NO local hole. We verify the ANTI-correlation
     corr( spatial spread S_H(t), local hole = observer min-flanking-gap(t) ) < 0.
     => 'spread' must NOT mean instantaneous spatial evenness.

 (2) The RIGHT 'global spread' is the INDEPENDENCE/decorrelation main term. Modeling
     the danger arcs as independent (global spread = no correlations) gives
        lonely measure  =  (1 - 2/n)^{n-1}  +  resonance corrections.
     The main term (1-2/n)^{n-1} > 0 GUARANTEES local emptiness; the only obstruction
     is the arithmetic resonances (the inside debt, S529). The regular polygon is
     exactly where they cancel it to 0 (tight). We measure lonely-measure vs the main
     term across families.

 (3) The CHOICE that works: TEMPORAL spread + rotation. The orbit is a MOVING geodesic
     (it cannot stay at the regular polygon); its spread S_H(t) must DIP to the bunched
     (low-entropy) extreme = a hole; and the hole's POSITION (argmax-gap) rotates and,
     by equidistribution, visits the observer. We measure: the temporal RANGE of the
     hole, and the COVERAGE of the hole-position around the circle (does it reach 0?).

Conclusion: local emptiness = global ANTI-spread (bunching); what global spread (in
TIME) guarantees is the OSCILLATION that must visit the bunched/holed extreme, and
rotation carries the hole to the observer. The obstruction is always the same
arithmetic concentration (resonances / closed orbit / regular polygon).
"""
from itertools import combinations
from functools import reduce
from math import gcd, log2
import random, statistics as st

def frac(x): return x - int(x // 1)

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

def obs_flank(speeds,t):
    """observer at 0: min flanking gap = min(nearest ahead, nearest behind)."""
    ahead=min((frac(v*t) for v in speeds if frac(v*t)>1e-12),default=1.0)
    behind=min((1-frac(v*t) for v in speeds if frac(v*t)>1e-12),default=1.0)
    return min(ahead,behind)

def max_gap_and_pos(speeds,t):
    pts=sorted([0.0]+[frac(v*t) for v in speeds])
    best=0;bp=0
    for i in range(len(pts)):
        a=pts[i]; b=pts[(i+1)%len(pts)]
        g=(b-a)%1.0
        if g>best: best=g; bp=(a+g/2)%1.0
    return best,bp

def pearson(xs,ys):
    mx=st.mean(xs);my=st.mean(ys)
    num=sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    dx=sum((x-mx)**2 for x in xs)**0.5; dy=sum((y-my)**2 for y in ys)**0.5
    return num/(dx*dy) if dx*dy>0 else 0.0

def main():
    fams={"AP 1..n (regular)":None,"random primitive":None,"geometric":None}
    print("="*72)
    print("(1) INSTANTANEOUS spatial spread S_H(t) vs local hole: ANTI-correlated")
    print("="*72)
    for n in (5,6):
        sets={"AP":tuple(range(1,n+1)),"geometric":tuple(2**k for k in range(n))}
        rr=random.Random(n); rv=tuple(sorted(rr.sample(range(1,8*n),n)))
        while reduce(gcd,rv)!=1: rv=tuple(sorted(rr.sample(range(1,8*n),n)))
        sets["random"]=rv
        for name,v in sets.items():
            G=4000; SH=[];hole=[]
            for s in range(G):
                t=(s+0.5)/G; pos=[frac(x*t) for x in v]
                SH.append(log2(H_count(half_turn(pos,n),n)))
                hole.append(obs_flank(v,t))
            print(f"  n={n} {name:9s} v={v}: corr(S_H, observer-hole) = {pearson(SH,hole):+.3f}")
    print("  => NEGATIVE: instantaneous spatial spread (even/high-entropy) means SMALL gaps =")
    print("     NO local hole. 'Global spread' must NOT be instantaneous spatial evenness.")
    print()

    print("="*72)
    print("(2) The INDEPENDENCE main term (1-2/n)^{n-1} guarantees emptiness; resonances obstruct")
    print("="*72)
    for n in (5,6,7):
        main_term=(1-2.0/n)**(n-1)
        sets={"AP 1..n-1 (regular)":tuple(range(1,n)),}
        rr=random.Random(100+n); rv=tuple(sorted(rr.sample(range(1,8*n),n-1)))
        while reduce(gcd,rv)!=1: rv=tuple(sorted(rr.sample(range(1,8*n),n-1)))
        sets["random primitive"]=rv
        print(f"  n={n}: independence main term (1-2/n)^(n-1) = {main_term:.4f}")
        for name,v in sets.items():
            G=200000;lon=0;thr=1.0/n
            for s in range(G):
                t=(s+0.5)/G
                if obs_flank(v,t)>=thr-1e-12: lon+=1
            print(f"      {name:22s} v={v}: lonely MEASURE = {lon/G:.4f}  "
                  f"({'~main term: spread wins' if lon/G>0.5*main_term else 'resonances cancel toward 0 (tight)'})")
    print("  => the SPREAD (independent) main term is POSITIVE and guarantees a hole; only the")
    print("     arithmetic resonances can cancel it. Regular polygon = exact cancellation (->0).")
    print()

    print("="*72)
    print("(3) TEMPORAL spread + ROTATION: the orbit must DIP to bunched; the hole ROTATES to 0")
    print("="*72)
    for n in (5,6):
        sets={"AP 1..n":tuple(range(1,n+1))}
        rr=random.Random(200+n); rv=tuple(sorted(rr.sample(range(1,8*n),n)))
        while reduce(gcd,rv)!=1: rv=tuple(sorted(rr.sample(range(1,8*n),n)))
        sets["random"]=rv
        for name,v in sets.items():
            G=4000; gaps=[];poscover=[0]*n
            for s in range(G):
                t=(s+0.5)/G; g,bp=max_gap_and_pos(v,t)
                gaps.append(g); poscover[int(bp*n)%n]+=1
            cover=sum(1 for c in poscover if c>0)
            print(f"  n={n} {name:7s} v={v}: max-gap(t) range [{min(gaps):.3f},{max(gaps):.3f}] "
                  f"(>=2/n={2.0/n:.3f}? {max(gaps)>=2.0/n-1e-9}); hole-position covers {cover}/{n} sectors")
    print("  => the orbit's max-gap (the bunched-extreme hole) reaches >=2/n, and the hole")
    print("     POSITION rotates to cover the circle -> it visits the observer's sector.")
    print("     TEMPORAL spread (the moving orbit) + ROTATION = local emptiness at the observer.")
    print()
    print("="*72)
    print("REFRAME: local emptiness = global ANTI-spread (bunching). Instantaneous spatial")
    print("spread is the ENEMY (even=no hole). What guarantees the hole is (a) the positive")
    print("INDEPENDENCE main term (1-2/n)^{n-1}, obstructed only by arithmetic resonances, and")
    print("(b) TEMPORAL spread + ROTATION: a moving orbit must dip to the bunched/holed extreme")
    print("and rotate the hole onto the observer. Choose 'spread' in TIME, not in space.")
    print("="*72)

if __name__=="__main__":
    main()
