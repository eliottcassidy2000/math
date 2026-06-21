#!/usr/bin/env python3
"""
lrc_periodmax_exact_full_macmini_0621s7.py  (mac-mini-2026-06-21-S7)
THREAD 1: COMPLETE THM-563 general bounded-base finite check.
period-max(B) <= 15*margin(B) for ALL primitive bounded bases B subset [0,14], 0 in B, |B|=k-1,
k=9..13.  (k=8 already done.)

period-max(B) = max over one period [0,P), P=7*lcm(B), of  f(w) = w*Delta_w
             =  sum over one-miss sector-j arcs (j,a,b) of  [Sc_j(w*b) - Sc_j(w*a)]   (THM-563 id).
Sc_j = centered sawtooth ( |S_j|<=3/49 , subtract mean ).

ALGORITHM (fast-then-exact, all candidates verified in EXACT Fractions):
  1. f(w) is periodic in w with period P=7*lcm(B). f depends on w only through w mod P.
  2. NUMPY float scan over all w in [0,P): compute f(w) for every residue (vectorized sawtooth).
     This is exact-shaped (piecewise linear in frac(w*t)); floats locate the argmax(es).
  3. Take ALL residues w0 within a tiny epsilon of the float max -> candidate set.
  4. Recompute f(w0) in EXACT Fractions for each candidate; period-max = exact max.
     (f(w) only takes finitely many rational values; the float max identifies the right residue
      class, and we then certify with exact arithmetic. We also exact-check the float runner-ups
      to guard against float ties.)
  5. PASS iff period-max < 15*margin. Report any FAIL exactly.

Optimization to bound work: the float scan needs O(P*#arcs). With #arcs ~ V/2 ~ 80 and P up to
2.5e6 this is ~2e8 float ops per base -> too slow for 12000 bases naively. So we FIRST screen with
a float estimate of period-max via a COARSE bound, then only do the full numpy scan for bases whose
coarse estimate is within a safety factor of 15*margin. The coarse screen: float period-max estimate
from a numpy scan at REDUCED period P_small = 7*lcm(B with elements capped) ... NO -- instead we use
the honest two-stage: (a) cheap a-priori sumR bound; if sumR < 15*margin PASS trivially (rare here);
(b) else full numpy float scan over [0,P) to get the float max, then exact-verify the candidates.
We run k by k and report progress; for each base we record (k,B,plat,margin,period_max_exact,ratio,PASS).
"""
import sys, itertools, time
from fractions import Fraction as F
from math import gcd
import numpy as np
sys.stdout.reconfigure(line_buffering=True)
def lcm(a,b): return a*b//gcd(a,b)
def sector_of(p): return int((p%1)*7)

def setup(E):
    """plat, arcs=[(j,a,b)], V, sumR; canonical breakpoint method."""
    E=sorted(set(E)); b=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): b.add(F(m,7*e))
    b=sorted(b); arcs=[];cur={};p0=F(0);p1=F(0)
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        secs=set(sector_of(e*((x0+x1)/2)) for e in E); miss=[j for j in range(1,7) if j not in secs]
        if len(secs)==7: p0+=x1-x0
        mj=miss[0] if len(miss)==1 else None
        if mj is not None: p1+=x1-x0
        for j in range(1,7):
            active=(mj==j)
            if active and j not in cur: cur[j]=x0
            if (not active) and j in cur: arcs.append((j,cur.pop(j),x0))
    for j,a in cur.items(): arcs.append((j,a,F(1)))
    eps=set()
    for (j,a,bb) in arcs: eps.add(a); eps.add(bb)
    return p0+p1/7, arcs, len(eps), F(6,49)*len(arcs)

# exact centered sawtooth
def Sj_raw(t,j):
    t=t-int(t)
    if t<0:t+=1
    a=F(j,7);b=F(j+1,7)
    if t<=a: return -t/7
    elif t<=b: return -a/7+F(6,7)*(t-a)
    else: return -a/7+F(6,7)*(b-a)-F(1,7)*(t-b)
def meanSj(j):
    a=F(j,7);b=F(j+1,7);pts=[F(0),a,b,F(1)];v=[Sj_raw(x,j) for x in pts];I=F(0)
    for i in range(3): I+=(v[i]+v[i+1])/2*(pts[i+1]-pts[i])
    return I
MEAN={j:meanSj(j) for j in range(1,7)}
MEANf={j:float(MEAN[j]) for j in range(1,7)}
def Sc(t,j): return Sj_raw(t,j)-MEAN[j]
def f_exact(arcs,w):
    tot=F(0)
    for (j,a,bb) in arcs: tot+=Sc(w*bb,j)-Sc(w*a,j)
    return tot

# float centered sawtooth, vectorized over an array of frac-values u in [0,1)
def Sc_vec(u,j):  # u = frac(w*t) array, returns Sc_j(u) float
    a=j/7.0; b=(j+1)/7.0; out=np.empty_like(u)
    m1=u<=a; m2=(~m1)&(u<=b); m3=(~m1)&(~m2)
    out[m1]=-u[m1]/7.0
    out[m2]=-a/7.0+(6.0/7.0)*(u[m2]-a)
    out[m3]=-a/7.0+(6.0/7.0)*(b-a)-(1.0/7.0)*(u[m3]-b)
    return out-MEANf[j]

def term_period_array(t,j,P):
    """Sc_j(frac(w*t)) for w=0..P-1, computed on the q-residue cycle then tiled to length P.
    t=p/q in lowest terms with q | P (since q | 7e and P=7*lcm(B))."""
    q=t.denominator
    p=t.numerator % q
    r=np.arange(q,dtype=np.int64)
    u=((r*p) % q).astype(np.float64)/q   # frac(w*t) for w=0..q-1
    base=Sc_vec(u,j)                      # length-q period
    reps=P//q
    return np.tile(base,reps)

def period_max_float_argmax(arcs,P):
    """Vectorized float f(w) for all w in [0,P) using per-endpoint period tiling.
    Return (float_max, candidate residues within tol)."""
    tot=np.zeros(P,dtype=np.float64)
    for (j,a,bb) in arcs:
        tot+=term_period_array(bb,j,P)
        if a!=0:
            tot-=term_period_array(a,j,P)
        # if a==0: Sc_j(0) is a constant (w-independent) -> shifts all entries equally, drop it
        else:
            tot-=float(Sc(F(0),j))
    mx=tot.max()
    tol=1e-9
    cand=np.where(tot>=mx-tol)[0]
    return mx, cand.tolist()

caps={9:F(1979,4004),10:F(55,91),11:F(66,91),12:F(6,7),13:F(1)}
ABS=F(6,49)

def process_k(k, Pcap):
    cap=caps[k]; t0=time.time()
    n=0; passed=0; fails=[]; n_trivial=0; n_scan=0; n_skip_big=0
    worst_ex=F(-10); worst_ex_B=None        # worst EXACT period-max ratio (scanned bases)
    worst_tr=F(0); worst_tr_B=None          # worst a-priori sumR ratio (trivial bases)
    skipped_bases=[]
    for combo in itertools.combinations(range(1,15),k-2):
        B=(0,)+combo; n+=1
        plat,arcs,V,sumR=setup(B)
        margin=cap-plat
        if margin<=0:
            fails.append(('OVERCAP',B,float(plat),float(margin)))
            continue
        thr=15*margin
        # cheap a-priori: period-max <= sumR (THM-563 range bound). If sumR<thr, PASS (no scan).
        if sumR<thr:
            n_trivial+=1; passed+=1
            r=sumR/margin
            if r>worst_tr: worst_tr=r; worst_tr_B=B
            continue
        L=1
        for e in B:
            if e>0: L=lcm(L,e)
        P=7*L
        if P>Pcap:
            n_skip_big+=1; skipped_bases.append((B,plat,margin,P,sumR)); continue
        n_scan+=1
        mxf,cands=period_max_float_argmax(arcs,P)
        pm=F(-10)
        for w0 in cands:
            v=f_exact(arcs,int(w0))
            if v>pm: pm=v
        r=pm/margin
        if r>worst_ex: worst_ex=r; worst_ex_B=B
        if pm>=thr:
            fails.append(('PERIODMAX',B,float(plat),float(margin),float(pm),str(pm),float(r)))
        else:
            passed+=1
    dt=time.time()-t0
    print(f"\n=== k={k} cap={cap}={float(cap):.5f}  {n} bases  ({dt:.1f}s) ===")
    print(f"  trivial(sumR<15m, no scan): {n_trivial}   scanned(exact pm): {n_scan}   skipped(P>{Pcap}): {n_skip_big}")
    print(f"  PASS: {passed}   FAIL: {len([x for x in fails])}")
    print(f"  worst EXACT period-max/margin (scanned bases) = {float(worst_ex):.4f} at B={worst_ex_B} (need <15)")
    print(f"  worst a-priori sumR/margin (trivial bases)    = {float(worst_tr):.4f} at B={worst_tr_B} (these PASS since sumR<15m => pm<=sumR<15m)")
    if fails:
        print("  FAILS:")
        for x in fails[:40]: print("    ",x)
    return n_skip_big, fails, skipped_bases

if __name__=="__main__":
    import sys
    ks=[int(x) for x in sys.argv[1:2]] if len(sys.argv)>1 else [9,10,11,12,13]
    Pcap=int(sys.argv[2]) if len(sys.argv)>2 else 200000
    allskip=0; allfail=0
    for k in ks:
        sk,fl,skb=process_k(k,Pcap)
        allskip+=sk; allfail+=len(fl)
    print(f"\nTOTAL skipped (P>{Pcap}, deferred to stage 2): {allskip}   TOTAL fails: {allfail}")
