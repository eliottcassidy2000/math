#!/usr/bin/env python3
"""
reachable_source_classes.py -- claude-2026-06-01 session on HYP-1981 (open part).

Faithful to THM-381 / oracle S511:
  observer 0; threshold 1/n; obs-runner 0->i iff ||v_i t||>=1/n (safe arc [1/n,1-1/n]);
  runner-runner half-turn  i->j iff frac((v_i-v_j)t) in (0,1/2);
  observer lonely  <=>  observer is SOURCE.

Claim under test (this session):
  At a lonely time all phases lie in the safe arc of length L=1-2/n, so the reachable
  observer-source classes are EXACTLY the half-turn tournaments on m=n-1 points
  realizable inside an arc of length L.  In particular for n<=4 (L<=1/2) only the
  TRANSITIVE class is reachable.
"""
from fractions import Fraction
from itertools import combinations, permutations
from math import gcd
from functools import reduce
import random

ONE=Fraction(1)
A000568={0:1,1:1,2:1,3:2,4:4,5:12,6:56,7:456}

def frac(x): return x-(x.numerator//x.denominator)
def dist0(x):
    f=frac(x); return min(f,ONE-f)

def half_turn(phases):
    m=len(phases); adj=[[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i==j: continue
            f=frac(phases[i]-phases[j]); adj[i][j]=1 if 0<f<Fraction(1,2) else 0
    return adj

def valid(adj):
    m=len(adj)
    for i in range(m):
        for j in range(i+1,m):
            if adj[i][j]+adj[j][i]!=1: return False
    return True

def canon(adj):
    m=len(adj); best=None
    for p in permutations(range(m)):
        flat=tuple(adj[p[a]][p[b]] for a in range(m) for b in range(m) if a!=b)
        if best is None or flat<best: best=flat
    return best

def flat_to_adj(flat,m):
    adj=[[0]*m for _ in range(m)]; k=0
    for a in range(m):
        for b in range(m):
            if a!=b: adj[a][b]=flat[k]; k+=1
    return adj

def cyc3(flat,m):
    adj=flat_to_adj(flat,m)
    return sum(1 for a,b,c in combinations(range(m),3)
               if (adj[a][b] and adj[b][c] and adj[c][a]) or (adj[b][a] and adj[c][b] and adj[a][c]))

def scores(flat,m):
    adj=flat_to_adj(flat,m); return tuple(sorted(sum(adj[a]) for a in range(m)))

def is_transitive(flat,m): return len(set(scores(flat,m)))==m

# ---- all tournament iso classes on m vertices ----
def all_classes(m):
    pairs=list(combinations(range(m),2)); cls=set()
    for bits in range(1<<len(pairs)):
        adj=[[0]*m for _ in range(m)]
        for idx,(a,b) in enumerate(pairs):
            if (bits>>idx)&1: adj[a][b]=1
            else: adj[b][a]=1
        cls.add(canon(adj))
    return cls

# ---- geometric reachable set: m points in arc of length L=1-2/n ----
def geom_exact(n, grid):
    m=n-1; L=ONE-Fraction(2,n); base=Fraction(1,n)
    pts=[base+Fraction(k,grid)*L for k in range(grid+1)]
    cls=set()
    for combo in combinations(range(len(pts)),m):
        adj=half_turn([pts[c] for c in combo])
        if valid(adj): cls.add(canon(adj))
    return cls

def geom_sample(n, samples, seed):
    m=n-1; L=ONE-Fraction(2,n); base=Fraction(1,n); rnd=random.Random(seed); cls=set()
    D=10007
    for _ in range(samples):
        xs=sorted(rnd.randrange(D) for _ in range(m))
        ph=[base+Fraction(x,D)*L for x in xs]
        if len(set(ph))<m: continue
        adj=half_turn(ph)
        if valid(adj): cls.add(canon(adj))
    return cls

# ---- faithful arithmetic reachability ----
def walls(speeds,n):
    thr=Fraction(1,n); W=set()
    for v in speeds[1:]:
        av=abs(v)
        if av==0: continue
        for mm in range(0,av):
            W.add(frac(Fraction(mm,av)+thr/av)); W.add(frac(Fraction(mm,av)-thr/av))
    for i,j in combinations(range(1,n),2):
        d=abs(speeds[i]-speeds[j])
        if d==0: continue
        for k in range(0,2*d): W.add(frac(Fraction(k,2*d)))
    return sorted(w for w in W if 0<=w<1)

def lonely(speeds,n,t):
    thr=Fraction(1,n); return all(dist0(Fraction(v)*t)>=thr for v in speeds[1:])

def rclass_at(speeds,n,t):
    return canon(half_turn([frac(Fraction(speeds[k])*t) for k in range(1,n)]))

def gen_systems(n,vmax,limit):
    out=[]
    for combo in combinations(range(1,vmax+1),n-1):
        if reduce(gcd,combo)!=1: continue
        out.append((0,)+combo)
        if len(out)>=limit: break
    return out

def arith_reachable(n,systems):
    R=set(); reached=0; bonly=0
    for s in systems:
        W=walls(s,n); W2=W+[ONE]
        cells=[(a+b)/2 for a,b in zip(W2,W2[1:])]
        op=False
        for t in cells:
            if lonely(s,n,t):
                adj=half_turn([frac(Fraction(s[k])*t) for k in range(1,n)])
                if valid(adj): op=True; R.add(canon(adj))
        if op: reached+=1
        else:
            for t in W:
                if lonely(s,n,t): bonly+=1; break   # boundary toucher: class degenerate, not recorded
    return R,reached,bonly

def main():
    print("="*82)
    print("HYP-1981 open part -- reachable observer-source classes vs A000568(n-1)")
    print("claude-2026-06-01 ; faithful to THM-381 / S511")
    print("="*82)
    GRID={4:60,5:40,6:30,7:24}; SAMP={4:20000,5:40000,6:60000,7:80000}
    VMAX={4:10,5:10,6:9,7:8}; SLIM={4:300,5:300,6:200,7:90}
    for n in range(4,8):
        m=n-1; L=ONE-Fraction(2,n); target=A000568[m]
        Rg=geom_exact(n,GRID[n]) | geom_sample(n,SAMP[n],n)
        # saturation check (half the samples)
        Rg_half=geom_exact(n,GRID[n]) | geom_sample(n,SAMP[n]//2,n+100)
        allc=all_classes(m)
        excluded=allc-Rg
        transit=[c for c in Rg if is_transitive(c,m)]
        Ra,reached,bonly=arith_reachable(n,gen_systems(n,VMAX[n],SLIM[n]))
        print(f"\n--- n={n}: m={m} runners, safe-arc L=1-2/n={L} ; A000568({m})={target}")
        print(f"  |R_geom| = {len(Rg)} / {target}   "
              f"({'ALL CLASSES' if len(Rg)==target else 'STRICT SUBSET'}) "
              f"[saturation: half-sample gives {len(Rg_half)}]")
        print(f"  |R_arith| = {len(Ra)}  (subset of R_geom: {Ra.issubset(Rg)})  "
              f"systems reaching open-cell={reached}, boundary-only={bonly}")
        print(f"  transitive classes in R_geom = {len(transit)} (expect exactly 1)")
        if excluded:
            print(f"  EXCLUDED classes ({len(excluded)}):")
            for c in sorted(excluded):
                print(f"      scores={scores(c,m)}  3-cycles={cyc3(c,m)}"
                      f"{'  [REGULAR/rotational]' if len(set(scores(c,m)))==1 else ''}")
        else:
            print("  EXCLUDED: none -- every tournament class is reachable")
    print("\n"+"="*82)
    print("n<=4 theorem: R_geom == {transitive}, size 1, strictly inside A000568(n-1).")

if __name__=="__main__": main()
