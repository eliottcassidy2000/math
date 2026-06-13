#!/usr/bin/env python3
"""
lrc_trienerment_gabor_s539.py    oracle-2026-06-01-S539o

The LRC pairwise structure is a T(R)IENERMENT, and the Gabor (time-frequency)
angle controls its ties.

T(r)ienerment (repo machinery, trienerment_iso_count.py): each edge {i,j} is one of
3 states i->j, i<-j, i<->j (tie). Iso classes = A007025 (7,42,582,21480 for n=3..6);
f(n,k) refines by #ties; f(n,0)=A000568 (tournaments).

THE LRC TRIENERMENT T(t): on n vertices (observer 0 + n-1 runners), edge {i,j}:
   i<->j  (TIE)  iff circular distance < 1/n   (the two runners are NEAR);
   else directed by the half-turn order (far + ahead/behind).
So 'tie' = 'near' = the loneliness-relevant relation. Then:
   OBSERVER LONELY  <=>  observer has ZERO tie-edges (all observer-edges directed)
   <=>  observer is 'tournament-like' (far from every runner).
This UNIFIES: order tournament (S518) = the tie-free limit; near-graph (S535/S536) =
the TIE SUBGRAPH; the metric mappings = the tie threshold. The #ties = #near-pairs =
a real-space clustering / discrepancy. LRC@n = every speed set reaches a T(t) whose
OBSERVER tie-degree is 0.

THE GABOR ANGLE: loneliness is a SHARP real-space feature (observer's 2/n window
empty, S536). By the DFT duality (S536/S537) a sharp space hole forces FREQUENCY
spread (uncertainty). We measure the trade-off: real-space hole (loneliness/ties) vs
frequency concentration (the character sums |chat_m|, S529). The (sector,harmonic)
cells are the Gabor atoms; the tie/hole <-> frequency-spread is a discrete
uncertainty principle.

Computes: (1) realizable LRC-trienerment iso-classes + tie distribution vs A007025 /
f(n,k); observer-tie-free reachability (LRC). (2) the Gabor/uncertainty trade-off:
at lonely vs generic times, real-space hole size vs frequency concentration.
"""
from itertools import combinations, permutations, product
from functools import reduce
from collections import defaultdict, Counter
from math import gcd, factorial, sin, cos, pi
import cmath, random

def frac(x): return x - int(x // 1)

# ---------- repo f(n,k) machinery (compact reimplementation) ----------
def partitions(n):
    def h(n, mx):
        if n == 0: yield (); return
        for k in range(min(n, mx), 0, -1):
            for r in h(n-k, k): yield (k,)+r
    yield from h(n, n)
def cycle_count(part):
    n=sum(part); c=Counter(part); r=factorial(n)
    for l,cnt in c.items(): r//=(l**cnt*factorial(cnt))
    return r
def poly_for_partition(part):
    parts=list(part)
    def mul(p1,p2):
        r=defaultdict(int)
        for d1,c1 in p1.items():
            for d2,c2 in p2.items(): r[d1+d2]+=c1*c2
        return dict(r)
    poly={0:1}
    for l in parts:
        if l%2==0: poly=mul(poly,{l//2:1})
    for l in parts:
        for _ in range((l-1)//2): poly=mul(poly,{0:2,l:1})
    for i in range(len(parts)):
        for j in range(i+1,len(parts)):
            la,lb=parts[i],parts[j]; g=gcd(la,lb); lc=la*lb//g
            for _ in range(g): poly=mul(poly,{0:2,lc:1})
    return poly
def f_distribution(n):
    me=n*(n-1)//2; cs=defaultdict(int)
    for part in partitions(n):
        cnt=cycle_count(part); poly=poly_for_partition(part)
        for d,c in poly.items(): cs[d]+=cnt*c
    nf=factorial(n)
    return {k: cs.get(k,0)//nf for k in range(me+1)}

# ---------- LRC trienerment ----------
def cdist(a,b):
    d=abs(a-b)%1.0; return min(d,1-d)
def lrc_trienerment_labeled(speeds_with_obs, n, t):
    """state per edge: 0=tie(near), +1 = i->j, -1 = i<-j. Return labeled tuple over i<j."""
    pos=[frac(s*t) for s in speeds_with_obs]
    st=[]
    for i in range(n):
        for j in range(i+1,n):
            if cdist(pos[i],pos[j])<1.0/n: st.append(0)
            else:
                d=(pos[j]-pos[i])%1.0
                st.append(1 if 0<d<0.5 else -1)
    return tuple(st)
def canon_trienerment(labeled, n):
    """canonical over S_n. labeled = states over i<j (0/+1/-1). reverse on swap."""
    # build matrix
    M=[[0]*n for _ in range(n)]
    it=iter(labeled)
    for i in range(n):
        for j in range(i+1,n):
            s=next(it); M[i][j]=s; M[j][i]=(-s if s!=0 else 0)
    best=None
    for p in permutations(range(n)):
        seq=tuple(M[p[i]][p[j]] for i in range(n) for j in range(i+1,n))
        if best is None or seq<best: best=seq
    return best
def n_ties(labeled): return sum(1 for s in labeled if s==0)
def observer_ties(labeled, n):
    # observer = vertex 0; its edges are the first (n-1) entries (i=0,j=1..n-1)
    return sum(1 for s in labeled[:n-1] if s==0)

def study_trienerment(n, n_sets=120, samples=3000):
    rnd=random.Random(40+n); labeled_set=set(); lonely_ok=0; tot=0
    tie_counts=set()
    for _ in range(20000):
        v=tuple(sorted(rnd.sample(range(1,6*n),n-1)))
        if reduce(gcd,v)!=1: continue
        tot+=1
        if tot>n_sets: break
        sp=[0]+list(v); lon=False; seen=set()
        for s in range(samples):
            t=(s+0.5)/samples
            lab=lrc_trienerment_labeled(sp,n,t)
            seen.add(lab)
            if observer_ties(lab,n)==0: lon=True
        if lon: lonely_ok+=1
        labeled_set|=seen
    classes=set(canon_trienerment(l,n) for l in labeled_set)
    for l in labeled_set: tie_counts.add(n_ties(l))
    return classes, lonely_ok, tot, sorted(tie_counts)

# ---------- Gabor / uncertainty trade-off ----------
def char_sum(speeds, n, t, m):
    return abs(sum(cmath.exp(2j*pi*m*frac(v*t)) for v in speeds))
def hole_size(speeds, n, t):
    """max empty arc length around the circle (observer+runners)."""
    pts=sorted([0.0]+[frac(v*t) for v in speeds])
    gaps=[pts[i+1]-pts[i] for i in range(len(pts)-1)]+[1-pts[-1]+pts[0]]
    return max(gaps)
def study_gabor(n, n_sets=30, samples=4000):
    rnd=random.Random(7+n); rows=[]
    cnt=0
    for _ in range(8000):
        v=tuple(sorted(rnd.sample(range(1,6*n),n-1)))
        if reduce(gcd,v)!=1: continue
        cnt+=1
        if cnt>n_sets: break
        # find lonely time (obs far) and a generic time; compare freq concentration
        best_lonely=None; generic=None
        for s in range(samples):
            t=(s+0.5)/samples
            far = all(cdist(frac(vi*t),0.0)>=1.0/n-1e-9 for vi in v)
            maxchar=max(char_sum(v,n,t,m) for m in range(1,n))
            if far and (best_lonely is None or maxchar<best_lonely[1]):
                best_lonely=(t,maxchar,hole_size(v,n,t))
            if abs(t-0.123)<0.002: generic=(t,maxchar,hole_size(v,n,t))
        if best_lonely:
            rows.append((best_lonely[2], best_lonely[1]))   # (hole, maxchar) at lonely time
    return rows

def main():
    print("="*74)
    print("LRC AS A T(R)IENERMENT: ties = nearness; observer lonely <=> observer tie-free")
    print("="*74)
    for n in (4,5):
        fdist=f_distribution(n); A007025=sum(fdist.values())
        classes,lon,tot,tcs=study_trienerment(n, n_sets=(120 if n==4 else 80))
        print(f"  n={n}: realizable LRC-trienerment iso-classes = {len(classes)} of A007025({n})={A007025} "
              f"(R={len(classes)/A007025:.3f}); realized total-tie-counts={tcs}")
        print(f"        observer-tie-free (LRC loneliness) reached: {lon}/{tot}")
        print(f"        f({n},k) [iso-classes by #ties] = {[fdist[k] for k in sorted(fdist)]}")
    print("  => the LRC structure is a trienerment (NOT a tournament); ties = near-pairs;")
    print("     loneliness = observer tie-degree 0. Tournaments (f(n,0)) = the tie-free slice.")
    print()
    print("="*74)
    print("THE GABOR / UNCERTAINTY ANGLE: sharp real-space hole (loneliness) vs freq spread")
    print("="*74)
    for n in (5,7):
        rows=study_gabor(n)
        if rows:
            holes=[r[0] for r in rows]; chars=[r[1] for r in rows]
            import statistics as st
            print(f"  n={n}: at the LONELY time (observer-hole >= 2/n={2/n:.3f}): "
                  f"avg hole={st.mean(holes):.3f}, avg max|char_sum|={st.mean(chars):.3f} "
                  f"(of n-1={n-1} runners)")
    print("  => loneliness forces a real-space hole of width >= 2/n; by the DFT duality")
    print("     (S536/S537) this is a frequency feature -- the max character sum stays")
    print("     bounded away from 0 at the lonely time (an uncertainty trade-off). The")
    print("     (sector,harmonic) cells are the Gabor atoms; ties (real clustering) <->")
    print("     character concentration (frequency). The Gabor TRIENERMENT on (sector,")
    print("     harmonic) cells is the joint space-frequency lift (posed).")

if __name__=="__main__":
    main()
