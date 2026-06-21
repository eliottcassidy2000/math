#!/usr/bin/env python3
"""
lrc14_route_D_ascent_and_order_opus_s63.py   (opus-2026-06-20-S63)

ROUTE D, part 2.  TEST 1 showed NO local compression operator is edge-monotone
(every operator has many Ly-DECREASING edges).  This part asks the weaker /
sharper questions that can still yield a proof:

  (A) STEEPEST-ASCENT REACHABILITY.  From every bounded shape, does greedy
      steepest Ly-ascent (over the union of ALL operator moves) terminate at
      consec?  If consec is the unique greedy attractor, the metagraph still has
      consec as its unique "peak", just not via a fixed move rule.

  (B) GLOBAL Ly-DAG.  Orient an edge E->E' iff Ly(E')>Ly(E) over ALL elementary
      moves (any operator).  Is consec the unique SINK (no outgoing edge)?  Are
      there OTHER sinks (local maxima that are not consec)?  A unique sink whose
      value is the global max is what we want.

  (C) CONVEX-ORDER / MAJORIZATION TEST (the HYP-2607 angle, since the
      extremality is JOINT not per-moment).  Does the empty-sector count
      distribution N_consec convex-dominate N_E for every bounded E?  i.e. is
      E[phi(N_consec)] >= E[phi(N_E)] for every convex phi?  Equivalently:
      do the partial sums of the SORTED-DESCENDING tail probabilities of N
      majorize?  We test the precise statement:  for every t,
            sum_{s>=t} p_s(consec)  vs  sum_{s>=t} p_s(E)
      and the convex-order statement via all convex test functions on {0..6}.
      The DUAL g (THM-534) is convex on the relevant range -- if consec is the
      convex-order MAX, then Ly(E)=E[g(N_E)]<=E[g(N_consec)]=Ly(consec).

  (D) MOMENT-VECTOR DOMINANCE.  The cleanest possible certificate: is there a
      SINGLE linear functional / ordering of the moment vector (S_1,...,S_4)
      under which consec is extremal AND which forces Ly?  TEST 1 already
      showed per-moment fails; here we test the 2-D / 3-D facet structure.

python3 stdlib only; exact Fractions.
"""
import sys, itertools, functools
from math import comb, gcd
from functools import reduce
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

@functools.lru_cache(maxsize=None)
def avoid_measure(Etuple, A):
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

@functools.lru_cache(maxsize=None)
def pdist(Etuple):
    """exact p_0..p_6: measure of x with exactly N empty sectors among {1..6}."""
    E=sorted(set(Etuple))
    bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for j in range(7):
            for t in (F(j,7),F(j+1,7)):
                for m in range(e): bps.add((t+m)/e)
    bps=sorted(z for z in bps if 0<=z<=1)
    P=[F(0)]*7
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        occ=[False]*7
        for e in E:
            occ[int(((e*xm)%1)*7)]=True
        N=sum(1 for j in range(1,7) if not occ[j])
        P[N]+=x1-x0
    return tuple(P)

def yvec(k):
    if k==8: g=lambda t:F((t-1)*(t-2)*(t-4)*(t-5),40)
    elif k in (9,10): g=lambda t:F(-(t-2)*(t-3)*(t-6),36)
    else: g=lambda t:F((t-3)*(t-4),12)
    return [g(t) for t in range(7)]  # g evaluated on N=0..6 directly

@functools.lru_cache(maxsize=None)
def Ly(Etuple,k):
    g=yvec(k); P=pdist(Etuple)
    return sum(g[t]*P[t] for t in range(7))

def norm(E):
    E=sorted(set(E)); E=[e-E[0] for e in E]
    g=reduce(gcd,E,0)
    if g>1: E=[e//g for e in E]
    return tuple(E)
def span(E): return max(E)-min(E)

# elementary moves = single offset relocation to ANY free slot keeping span<=current
def all_compress_moves(E):
    E=sorted(set(E)); out=[]
    lo,hi=0,span(E)
    for i,e in enumerate(E):
        rest=[x for x in E if x!=e]
        for nv in range(min(rest),max(rest)+1):
            if nv in rest: continue
            cand=norm(rest+[nv])
            if len(cand)!=len(E): continue
            if span(cand)>span(E): continue
            out.append(cand)
    return list(set(out))

def all_shapes(k, max_span):
    seen=set()
    for rest in itertools.combinations(range(1,max_span+1),k-1):
        E=norm([0]+list(rest))
        if span(E)<=max_span and E not in seen:
            seen.add(E); yield E

# ---------------- (A) steepest-ascent reachability ----------------
def steepest_ascent(E,k,max_steps=60):
    cur=norm(E)
    for _ in range(max_steps):
        ly=Ly(cur,k); best=None; bv=ly
        for nb in all_compress_moves(cur):
            v=Ly(nb,k)
            if v>bv: bv=v; best=nb
        if best is None: return cur  # local max
        cur=best
    return cur

# ---------------- (C) convex order ----------------
def convex_test_funcs():
    """basis of convex functions on {0..6}: hinges max(0, t-c) for c=0..5, plus
       constant and identity (affine part) -- convex order = nonneg on all hinges
       AND equal means; but our distributions need not share mean, so we test the
       full convex cone:  phi convex  <=>  second differences >=0.
       Generators of the convex cone on a discrete grid {0..6} are the hinge
       functions h_c(t)=max(0,t-c), c=-1..5  (and +/- affine).  Convex DOMINANCE
       p >=_cx q  means E_p[h_c] >= E_q[h_c] for all c AND E_p[t]=E_q[t]
       (same mean) OR (increasing convex order) E_p[h_c]>=E_q[h_c] all c."""
    funcs=[]
    for c in range(-1,6):
        funcs.append(('hinge>=%d'%(c+1), lambda t,c=c: F(max(0,t-c))))
    return funcs

if __name__=="__main__":
    print("="*72)
    print("ROUTE D part 2: ascent reachability, global DAG sinks, convex order")
    print("="*72)
    for k in (8,9):
        ms = 11 if k==8 else 12
        consec=tuple(range(k)); lyc=Ly(consec,k)
        shapes=list(all_shapes(k,ms))
        print(f"\n### k={k}, span<={ms}, {len(shapes)} shapes.  Ly(consec)={float(lyc):.6f} ###")

        # (A) steepest ascent
        landed=[steepest_ascent(E,k) for E in shapes]
        not_consec=[(E,L) for E,L in zip(shapes,landed) if L!=consec]
        print(f"(A) steepest-ascent: {len(shapes)-len(not_consec)}/{len(shapes)} land at consec")
        if not_consec:
            # report distinct non-consec attractors and whether they beat consec
            attr={}
            for E,L in not_consec: attr.setdefault(L,0); attr[L]+=1
            for L,c in sorted(attr.items(), key=lambda kv:-kv[1])[:6]:
                print(f"    attractor {L}  Ly={float(Ly(L,k)):.6f}  basin={c}  "
                      f"{'>=consec!' if Ly(L,k)>=lyc else '<consec'}")

        # (B) global DAG: count sinks (no Ly-increasing move)
        sinks=[]
        for E in shapes:
            ly=Ly(E,k)
            if all(Ly(nb,k)<=ly for nb in all_compress_moves(E)):
                sinks.append(E)
        consec_is_top = all(Ly(E,k)<=lyc for E in shapes)
        print(f"(B) global Ly-DAG: #local-max sinks={len(sinks)}; consec global top={consec_is_top}")
        nonconsec_sinks=[s for s in sinks if s!=consec]
        for s in nonconsec_sinks[:6]:
            print(f"    extra sink {s} Ly={float(Ly(s,k)):.6f} (vs consec {float(lyc):.6f})")

        # (C) convex order on N distribution
        funcs=convex_test_funcs()
        Pc=pdist(consec)
        cx_dom=0; cx_fail=0; mean_eq=0; fails=[]
        for E in shapes:
            if E==consec: continue
            Pe=pdist(E)
            ok=True
            for nm,phi in funcs:
                Ec=sum(Pc[t]*phi(t) for t in range(7))
                Ee=sum(Pe[t]*phi(t) for t in range(7))
                if Ec<Ee: ok=False; break
            if ok: cx_dom+=1
            else: cx_fail+=1; fails.append(E)
        print(f"(C) increasing-convex order N_consec >=_icx N_E: holds {cx_dom}, "
              f"FAILS {cx_fail}/{len(shapes)-1}")
        if fails:
            E=fails[0]
            print(f"    first failure: E={E}  meanN(consec)={float(sum(t*Pc[t] for t in range(7))):.4f} "
                  f"meanN(E)={float(sum(t*pdist(E)[t] for t in range(7))):.4f}")

        # Is the DUAL g convex on 0..6?  (if so, convex-order would give Ly)
        g=yvec(k)
        d2=[g[t+2]-2*g[t+1]+g[t] for t in range(5)]
        print(f"    dual g(N) second-differences on 0..6: {[str(x) for x in d2]}  "
              f"({'CONVEX' if all(x>=0 for x in d2) else 'NOT convex'})")
    print("\nDONE part 2.")
