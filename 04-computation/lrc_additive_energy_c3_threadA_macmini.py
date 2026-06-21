#!/usr/bin/env python3
"""
lrc_additive_energy_c3_threadA_macmini.py
THREAD A: additive energy / Schur triples  <->  tournament c3 (THM-559 Ising).

GROUNDED OBJECT (LRC measS7 / multi-dim Weyl):
  Offset set E = {e_1,...,e_k} subset Z_{>=0} (the LRC "offsets" / freqs).
  Color of offset e at point x in [0,1):  c(e,x) = floor(7 * frac(e*x))  in Z/7.
  measS7(E)  = meas{ x in [0,1) : the k colors c(e,x) hit ALL 7 values }  (= COVER).
  iid_k      = 7! S(k,7)/7^k = surjective-occupancy prob for k iid-uniform colors
             = 1 - sum_{j} (-1)^j C(7,j) (1-j/7)^k    (decorrelated / single-particle).
  corr(E)    = measS7(E) - iid_k    (the carrier / Weyl error; SIGNED).

FOURIER / OFFSET-RELATION-LATTICE picture (derived in docstring below):
  The all-7-hit measure is, by inclusion-exclusion over which sectors are MISSED,
     measS7(E) = sum_{T subset Z/7} (-1)^|T| * meas{x : no color lands in T}.
  meas{x: all colors avoid sector-subset T} = integral_0^1 prod_e g_T(e x) dx,
  where g_T(t)=1[floor(7 frac t) not in T] has Fourier series sum_n a_n(T) e(n t).
  The integral picks out ONLY frequency combinations summing to 0:
     = sum_{ (n_1,...,n_k): sum n_e e_e = 0 } prod_e a_{n_e}(T).
  The n=0 term (all n_e=0) = prod a_0(T) = ((7-|T|)/7)^k -> gives exactly iid_k after IE.
  Therefore:
     corr(E) = sum_{ 0 != (n_e) in Lambda(E) } [ IE-weighted Fourier product ],
  where Lambda(E) = { n in Z^k : sum_e n_e e_e = 0 } is the OFFSET RELATION LATTICE
  (the LRC analog of the GF(2) cycle space of K_n).

  SHORTEST nonzero relations supported on >=3 offsets are ADDITIVE TRIANGLES:
  Schur triples a+b=c  (relation +1.a +1.b -1.c = 0). Support-2 relations need
  e_i = e_j (none, offsets distinct) or n_i e_i + n_j e_j=0 (rational ratios).

CLAIM TO TEST:
  |corr(E)| is controlled by the offset set's ADDITIVE ENERGY
     E_+(E) = #{(a,b,c,d) in E^4 : a+b = c+d}
  (and/or Schur-triple count ST(E)=#{(a,b,c): a+b=c, all in E}),
  the LRC analog of c3 (THM-559: c3 is a 2-body line-graph Ising energy; regular = max).
  consec (max additive structure) <-> regular (max c3): is consec the max |corr| = "max-cover deviation"?

This script:
  (i)  correlates |corr(E)| with E_+(E) and ST(E) over many E (k=8,9), exact Fractions;
  (ii) checks whether consec maximizes additive energy AND |corr| (parallel to regular=max c3);
  (iii) reports whether additive energy SPLITS structured vs dissociated, i.e. whether
        low cross-block E_+ -> small |corr| (the lever for the multi-block carrier error).
"""
import sys, itertools, random
from fractions import Fraction as F
from math import comb
sys.stdout.reconfigure(line_buffering=True)

# ---------- exact measS7 (cover = all 7 sectors hit) ----------
def measS7(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; secs=set()
        for e in E:
            if e==0: continue
            secs.add(int(((e*xm)%1)*7))
        # offset 0 (if present) always sits in sector 0; include it:
        if 0 in E: secs.add(0)
        if len(secs)==7: total+=x1-x0
    return total

def iid_k(k):
    # iid floor of COVER (all 7 sectors hit) = surjective-occupancy prob = 7! S(k,7)/7^k
    # = sum_{j=0}^{7} (-1)^j C(7,j) (1-j/7)^k  (inclusion-exclusion over MISSED sectors).
    return sum((-1)**j*comb(7,j)*F(7-j,7)**k for j in range(8))

def corr(E):
    return measS7(E)-iid_k(len(set(E)))

# ---------- additive energy and Schur triples ----------
def additive_energy(E):
    """E_+(E) = #{(a,b,c,d) in E^4 : a+b=c+d}.  Via r(s)=#{(a,b): a+b=s}, E_+=sum r(s)^2."""
    from collections import Counter
    r=Counter()
    for a in E:
        for b in E:
            r[a+b]+=1
    return sum(v*v for v in r.values())

def schur_triples(E):
    """ST(E) = #ordered (a,b,c) in E^3 with a+b=c (a Schur/additive triangle)."""
    Eset=set(E); cnt=0
    for a in E:
        for b in E:
            if a+b in Eset: cnt+=1
    return cnt

def schur_triples_distinct(E):
    """unordered {a,b,c} distinct with a+b=c (true additive triangles, support-3 relations)."""
    Eset=set(E); cnt=0
    for a,b in itertools.combinations(sorted(E),2):
        if a!=b and a+b in Eset and (a+b)!=a and (a+b)!=b:
            cnt+=1
    return cnt

# ---------- (i) correlation over a population of E ----------
def pearson(xs, ys):
    n=len(xs); mx=sum(xs)/n; my=sum(ys)/n
    sxx=sum((x-mx)**2 for x in xs); syy=sum((y-my)**2 for y in ys)
    sxy=sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    if sxx==0 or syy==0: return float('nan')
    return sxy/(sxx*syy)**0.5

def spearman(xs, ys):
    def ranks(v):
        order=sorted(range(len(v)), key=lambda i:v[i])
        rk=[0]*len(v)
        i=0
        while i<len(v):
            j=i
            while j+1<len(v) and v[order[j+1]]==v[order[i]]: j+=1
            avg=(i+j)/2.0
            for t in range(i,j+1): rk[order[t]]=avg
            i=j+1
        return rk
    return pearson(ranks(xs), ranks(ys))

def gen_population(k, n_random=120, spread_max=40, seed=1):
    """A varied population of k-offset sets E (always containing 0, distinct positive offsets)."""
    rng=random.Random(seed)
    pop=[]
    # canonical structured shapes
    pop.append(("consec", list(range(k))))
    pop.append(("AP_step2", [2*i for i in range(k)]))
    pop.append(("AP_step3", [3*i for i in range(k)]))
    pop.append(("geom2", [0]+[2**i for i in range(k-1)]))   # dissociated-ish
    # Sidon-ish (small additive energy): Mian-Chowla
    mc=[0]; s=set()
    cand=1
    while len(mc)<k:
        ok=True
        for x in mc:
            if (cand+x) in s: ok=False; break
        if ok:
            for x in mc: s.add(cand+x)
            s.add(2*cand); mc.append(cand)
        cand+=1
    pop.append(("MianChowla", mc))
    # two-block (multi-block carrier!) shapes: small base + far-shifted copy
    half=k//2
    for M in (50, 200, 1000):
        base=list(range(half))
        blk=("twoblock_M%d"%M, base+[M+i for i in range(k-half)])
        pop.append(blk)
    # random dissociated-ish sets
    for _ in range(n_random):
        S={0}
        while len(S)<k:
            S.add(rng.randint(1,spread_max))
        pop.append(("rand", sorted(S)))
    # dedup by tuple
    seen=set(); out=[]
    for name,E in pop:
        t=tuple(sorted(set(E)))
        if len(t)!=k: continue
        if t in seen: continue
        seen.add(t); out.append((name,list(t)))
    return out

def run_correlation(k):
    print("="*92)
    print(f"(i) CORRELATION  |corr(E)|  vs  additive energy E_+  and Schur-triples  (k={k})")
    print("="*92)
    pop=gen_population(k)
    rows=[]
    for name,E in pop:
        c=corr(E); ac=abs(float(c))
        ep=additive_energy(E); st=schur_triples(E); std=schur_triples_distinct(E)
        rows.append((name,E,ac,ep,st,std,float(c)))
    # correlations over the whole population
    acs=[r[2] for r in rows]; eps=[float(r[3]) for r in rows]
    sts=[float(r[4]) for r in rows]; stds=[float(r[5]) for r in rows]
    print(f"  population size = {len(rows)}")
    print(f"  Pearson  r(|corr|, E_+)        = {pearson(acs,eps):+.4f}")
    print(f"  Spearman rho(|corr|, E_+)      = {spearman(acs,eps):+.4f}")
    print(f"  Pearson  r(|corr|, Schur)      = {pearson(acs,sts):+.4f}")
    print(f"  Spearman rho(|corr|, Schur)    = {spearman(acs,sts):+.4f}")
    print(f"  Spearman rho(|corr|, SchurDst) = {spearman(acs,stds):+.4f}")
    # show the extremes
    rows_by_corr=sorted(rows, key=lambda r:-r[2])
    print("\n  Top 6 by |corr| (largest carrier error):")
    for r in rows_by_corr[:6]:
        print(f"    {r[0]:<14} |corr|={r[2]:.4f}  E_+={r[3]:>5}  Schur={r[4]:>4}  E={r[1]}")
    print("  Bottom 6 by |corr| (smallest carrier error):")
    for r in rows_by_corr[-6:]:
        print(f"    {r[0]:<14} |corr|={r[2]:.4f}  E_+={r[3]:>5}  Schur={r[4]:>4}  E={r[1]}")
    return rows

# ---------- (ii) is consec the max additive energy AND max |corr|? ----------
def run_extremality(k):
    print("="*92)
    print(f"(ii) EXTREMALITY: does consec maximize E_+ AND |corr|?  (parallel: regular=max c3)  k={k}")
    print("="*92)
    consec=list(range(k))
    print(f"  consec E={consec}")
    print(f"    E_+(consec)   = {additive_energy(consec)}")
    print(f"    Schur(consec) = {schur_triples(consec)}")
    print(f"    |corr|(consec)= {abs(float(corr(consec))):.6f}   measS7={float(measS7(consec)):.6f}")
    # exhaustive over all k-subsets of {0..k+W} for small window W, find max E_+ and max |corr|
    W = {8:4, 9:3}.get(k, 2)
    universe=list(range(k+W))
    best_ep=(-1,None); best_corr=(-1.0,None)
    cnt=0
    for combo in itertools.combinations(universe, k):
        if combo[0]!=0:  # normalize: contain 0 (translation-invariant for E_+; corr depends on actual offsets so keep 0)
            continue
        E=list(combo); cnt+=1
        ep=additive_energy(E)
        if ep>best_ep[0]: best_ep=(ep,E)
        ac=abs(float(corr(E)))
        if ac>best_corr[0]: best_corr=(ac,E)
    print(f"  searched {cnt} subsets of [0,{k+W-1}] containing 0:")
    print(f"    MAX E_+   = {best_ep[0]} at E={best_ep[1]}  (consec? {best_ep[1]==consec})")
    print(f"    MAX |corr|= {best_corr[0]:.6f} at E={best_corr[1]}  (consec? {best_corr[1]==consec})")
    return best_ep, best_corr

# ---------- (iii) cross-block additive energy: the multi-block lever ----------
def cross_block_energy(base, blk):
    """E_+ restricted to quadruples (a,b,c,d) NOT all in one block: the CROSS-BLOCK additive energy."""
    E=base+blk
    full=additive_energy(E)
    within=additive_energy(base)+additive_energy(blk)
    # within over-counts? additive_energy(base)+additive_energy(blk) counts quadruples all in base
    # plus all in blk. cross = full - those-all-in-one-block. But a+b=c+d with a,b in base,
    # c,d in base is "all in base". Mixed quadruples are cross. So:
    return full - within, full, within

def run_multiblock(k):
    print("="*92)
    print(f"(iii) MULTI-BLOCK: cross-block additive energy vs |corr| as blocks separate  (k={k})")
    print("="*92)
    half=k//2
    base=list(range(half))
    print(f"  base block = {base} (size {half}); far block = M + range({k-half})")
    print(f"  {'M':>6}{'crossE_+':>10}{'|corr|':>10}{'measS7':>10}{'iid':>10}")
    iv=float(iid_k(k))
    last=None
    for M in [half, half+1, half+3, 10, 30, 100, 300, 1000, 5000, 50000]:
        blk=[M+i for i in range(k-half)]
        E=base+blk
        if len(set(E))!=k: continue
        cross,full,within=cross_block_energy(base, blk)
        ac=abs(float(corr(E)))
        ms=float(measS7(E))
        print(f"  {M:>6}{cross:>10}{ac:>10.5f}{ms:>10.5f}{iv:>10.5f}")
        last=(M,cross,ac,ms)
    print("\n  INTERPRETATION: if |corr| -> a LIMIT as M->inf governed by cross-block E_+ shrinking,")
    print("  then a bound |corr| <= f(cross-E_+) closes the multi-block case when cross-E_+ is low.")
    return None

if __name__=="__main__":
    for k in (8,9):
        run_correlation(k)
        print()
    for k in (8,9):
        run_extremality(k)
        print()
    for k in (8,9):
        run_multiblock(k)
        print()
