#!/usr/bin/env python3
"""
lrc14_route_D_compression_metagraph_opus_s63.py   (opus-2026-06-20-S63)

ROUTE D: the COMPRESSION METAGRAPH with U4-gradient (HYP-2693/2694 obligation #1).

GOAL.  Prove that consec_k = {0,1,...,k-1} MAXIMIZES the empty-sector moment
functional U4(E) = Ly(E,k) = sum_r y_r S_r(E) over ALL bounded offset shapes E
(0 in E, |E|=k).  This would close the LRC(14) compression/extremality step.

We mirror the project's G_n H-gradient metagraph:
   - VERTICES = bounded k-shapes E (scale-reduced so E ~ dE identified).
   - EDGES    = elementary "compression moves" between shapes.
   - ORIENT each edge by sign(Delta U4).  Want: a DAG with UNIQUE SINK = consec,
     i.e. every move that reduces "spread" raises U4 monotonically.

HYP-2603 WARNING (a real dead end, do NOT re-walk):  the naive one-step move
"shift one offset toward its neighbor by 1" can DECREASE U4.  So a naive
nearest-neighbor compression is NOT monotone.  This script tests SEVERAL
candidate operators and reports, exactly and over a thorough bank, which (if
any) is monotone toward consec.  We are brutally honest about leverage.

The moment functional (THM-534 / THM-556 duals, k-dependent):
   k=8:  Ly = sum y_r S_r,  y = (1, -1, 1, -9/10, 3/5, 0, 0)
   k=9,10: y = (1, -13/18, 4/9, -1/6, 0, 0, 0)
   k>=11: g(t)=(t-3)(t-4)/12  =>  y = (1, -7/12, 1/6, ...)  [U4 = 1-S1+S2-S3+S4 form]
S_r(E) = E[C(N,r)] = sum_{|A|=r subset {1..6}} J(A,E),  J = avoid-measure (exact PL).

Engine (avoid_measure / S_r / Ly / cap) is reused from
lrc14_hyp2607_moment_extremal_macmini_0619s1.py (the canonical U4 engine).

python3 stdlib only; exact Fractions.
"""
import sys, itertools, functools
from math import comb
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

# ---------- canonical U4 engine (reused) ----------
@functools.lru_cache(maxsize=None)
def avoid_measure(Etuple, A):
    """meas{x in [0,1): orbit {frac(e x): e in E} avoids every sector in A subset {1..6}}."""
    E=sorted(set(Etuple))
    forb=[(F(j,7),F(j+1,7)) for j in A]
    bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for (lo,hi) in forb:
            for t in (lo,hi):
                for m in range(e): bps.add((t+m)/e)
    bps=sorted(z for z in bps if 0<=z<=1); tot=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; ok=True
        for e in E:
            p=(e*xm)%1
            for (lo,hi) in forb:
                if lo<=p<hi: ok=False; break
            if not ok: break
        if ok: tot+=x1-x0
    return tot

def S_r(E, r):
    Et=tuple(sorted(set(E)))
    if r==0: return F(1)
    return sum(avoid_measure(Et,A) for A in itertools.combinations(range(1,7),r))

def yvec(k):
    if k==8: g=lambda t:F((t-1)*(t-2)*(t-4)*(t-5),40)
    elif k in (9,10): g=lambda t:F(-(t-2)*(t-3)*(t-6),36)
    else: g=lambda t:F((t-3)*(t-4),12)
    y=[F(0)]*7
    for r in range(7):
        s=F(0)
        for i in range(r+1): s+=F((-1)**(r-i)*comb(r,i))*g(i)
        y[r]=s
    return y

@functools.lru_cache(maxsize=None)
def Ly(Etuple,k):
    y=yvec(k); E=list(Etuple)
    return sum(y[r]*S_r(E,r) for r in range(7) if y[r]!=0)

# ---------- cap_k (the threshold to clear) ----------
H=F(1,14)
def danger(u):
    iv=[]
    for j in range(u):
        c=F(j,u); a=(c-H/u)%1; b=(c+H/u)%1
        if a<b: iv.append((a,b))
        else: iv.append((a,F(1))); iv.append((F(0),b))
    return iv
def mgm(iv):
    iv=sorted(iv); o=[]
    for a,b in iv:
        if o and a<=o[-1][1]: o[-1]=(o[-1][0],max(o[-1][1],b))
        else: o.append((a,b))
    return o
def measGP(P):
    dz=mgm([iv for u in P for iv in danger(u)]); s=F(0); prev=F(0)
    for a,b in dz:
        if a>prev: s+=a-prev
        prev=max(prev,b)
    if prev<1: s+=1-prev
    return s
@functools.lru_cache(None)
def cap(k): return min(measGP(P) for P in itertools.combinations(range(1,14),13-k))

# ---------- shape normalization (scale + translation invariance) ----------
def norm(E):
    """canonical rep: translate so min=0, divide by gcd (dilation invariant E ~ dE)."""
    from math import gcd
    from functools import reduce
    E=sorted(set(E)); E=[e-E[0] for e in E]
    g=reduce(gcd,E,0)
    if g>1: E=[e//g for e in E]
    return tuple(E)

def span(E): return max(E)-min(E)

# =====================================================================
# COMPRESSION OPERATORS.  Each takes a shape E and yields candidate
# "more compressed" shapes E' (smaller spread, same size, 0 in E').
# We then ORIENT by sign(Ly(E')-Ly(E)) and look for monotonicity.
# =====================================================================

def op_neighbor_pull(E):
    """O1 (NAIVE, HYP-2603 baseline): move ONE interior offset 1 step toward
       the centroid, if the target slot is free.  Records ALL such moves."""
    E=sorted(set(E)); cen=F(sum(E),len(E)); out=[]
    for i,e in enumerate(E):
        step = 1 if e<cen else (-1 if e>cen else 0)
        if step==0: continue
        ne=e+step
        if ne in E: continue
        cand=norm([x for x in E if x!=e]+[ne])
        if len(cand)==len(E): out.append(('pull',e,ne,cand))
    return out

def op_gap_close(E):
    """O2: find a GAP (consecutive offsets differing by >1); shift the entire
       upper block down by 1 to close the gap by one unit.  This is a 'block
       compression' that preserves internal structure of each side."""
    E=sorted(set(E)); out=[]
    for i in range(len(E)-1):
        if E[i+1]-E[i]>1:
            upper=set(E[i+1:])
            cand=norm([x for x in E if x<=E[i]]+[x-1 for x in E if x>E[i]])
            if len(cand)==len(E): out.append(('gapclose',E[i],E[i+1],cand))
    return out

def op_endpoint_pull(E):
    """O3: pull the EXTREME offset (max, or the shifted min after norm) inward
       to the nearest free slot."""
    E=sorted(set(E)); out=[]
    # pull max down
    top=E[-1]
    for nv in range(top-1, E[-2], -1):
        if nv not in E:
            cand=norm(E[:-1]+[nv]);
            if len(cand)==len(E): out.append(('endpull-max',top,nv,cand)); break
    return out

def op_steiner(E):
    """O4 (Steiner-symmetrization analogue): reflect the shape to make it more
       'balanced' AND simultaneously contract spread by 1 by closing the largest
       interior gap, choosing the side (upper/lower) whose block move keeps it
       most balanced.  A two-sided gap close."""
    E=sorted(set(E)); out=[]
    # find the largest interior gap
    gaps=[(E[i+1]-E[i],i) for i in range(len(E)-1)]
    gaps.sort(reverse=True)
    if not gaps or gaps[0][0]<=1: return out
    _,i=gaps[0]
    # close it: move upper block down by 1
    cand=norm([x for x in E if x<=E[i]]+[x-1 for x in E if x>E[i]])
    if len(cand)==len(E): out.append(('steiner',i,E[i],cand))
    return out

def op_difference_compress(E):
    """O5 (the DIFFERENCE-SET / relation-lattice operator).  The tournament T(x)
       lives on the DIFFERENCES e_i-e_j.  Compressing the difference multiset
       toward {1,2,...} is the natural 'relation-lattice densification'.
       MOVE: pick the offset whose removal+reinsertion at a free slot MOST
       reduces the spread of the difference multiset (sum of |diffs|), among
       moves that do not increase span."""
    E=sorted(set(E)); out=[]
    def diffmass(S):
        S=sorted(S); return sum(S[j]-S[i] for i in range(len(S)) for j in range(i+1,len(S)))
    base=diffmass(E)
    for i,e in enumerate(E):
        rest=[x for x in E if x!=e]
        lo,hi=min(rest),max(rest)
        for nv in range(lo, hi+1):
            if nv in rest: continue
            cand_full=sorted(rest+[nv])
            if span(cand_full)>span(E): continue
            if diffmass(cand_full)<base:
                cand=norm(cand_full)
                if len(cand)==len(E): out.append(('diffcomp',e,nv,cand))
    return out

OPS={'O1_neighbor_pull':op_neighbor_pull,
     'O2_gap_close':op_gap_close,
     'O3_endpoint_pull':op_endpoint_pull,
     'O4_steiner':op_steiner,
     'O5_difference_compress':op_difference_compress}

# =====================================================================
# TEST 1: is each operator MONOTONE (every move raises Ly)?
#   For a thorough bank of bounded shapes, count moves that DECREASE Ly.
#   A monotone operator => 0 decreasing moves AND every non-consec shape has
#   at least one outgoing move (so consec is the unique sink/global max).
# =====================================================================
def all_shapes(k, max_span):
    """all normalized bounded shapes of size k, 0=min, span<=max_span."""
    seen=set()
    for rest in itertools.combinations(range(1,max_span+1),k-1):
        E=norm([0]+list(rest))
        if span(E)<=max_span and E not in seen:
            seen.add(E); yield E

def test_operator(k, max_span, opname, opf):
    consec=tuple(range(k))
    lyc=Ly(consec,k)
    dec=0          # moves that DECREASE Ly  (violations of monotonicity)
    eq=0           # neutral moves
    inc=0          # good moves
    stuck=[]       # non-consec shapes with NO outgoing move (would-be extra sinks)
    nshapes=0
    for E in all_shapes(k,max_span):
        nshapes+=1
        lyE=Ly(E,k)
        moves=opf(E)
        # only count moves that strictly reduce spread OR are the canonical compress
        outgoing=0
        for tag,a,b,cand in moves:
            if cand==E: continue
            outgoing+=1
            d=Ly(cand,k)-lyE
            if d<0: dec+=1
            elif d>0: inc+=1
            else: eq+=1
        if E!=consec and outgoing==0:
            stuck.append(E)
    return dict(nshapes=nshapes, inc=inc, dec=dec, eq=eq, stuck=stuck, consec_ly=lyc)

if __name__=="__main__":
    print("="*72)
    print("ROUTE D: compression metagraph, U4-gradient.  consec = unique sink?")
    print("="*72)
    for k in (8,9):
        ms = 11 if k==8 else 12
        print(f"\n### k={k}, shapes with span<={ms}  (cap_{k}={float(cap(k)):.5f}) ###")
        consec=tuple(range(k))
        print(f"    Ly(consec)={float(Ly(consec,k)):.6f}  (the claimed global MAX)")
        for opname,opf in OPS.items():
            r=test_operator(k,ms,opname,opf)
            verdict=("MONOTONE-DAG: all moves raise Ly, consec unique sink"
                     if r['dec']==0 and not r['stuck'] else
                     f"NOT monotone (dec={r['dec']}, stuck={len(r['stuck'])})")
            print(f"  {opname:24s} shapes={r['nshapes']:4d} inc={r['inc']:4d} "
                  f"dec={r['dec']:4d} eq={r['eq']:4d} stuck={len(r['stuck']):2d}  -> {verdict}")
            if r['stuck'] and len(r['stuck'])<=4:
                print(f"        stuck shapes (no outgoing move): {r['stuck']}")
    print("\nDONE TEST 1.")
