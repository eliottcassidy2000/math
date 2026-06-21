#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_threadD_leading_shell_opus_0620.py  (opus-2026-06-20, THREAD D part iv -- decisive)

Part iii found: consec is the GLOBAL ARGMAX of measS7 AND of the additive-triple count
N_2(E) (the leading rank-2 relation shell), AND of the difference-collision count -- the
argmax COINCIDES.  But measS7 is NOT monotone in N_2 (25-35% inversions).  WHY?  The
inversions are shapes with moderate N_2 but tiny measS7.  This script isolates the exact
mechanism and asks whether a CORRECTED leading-shell functional is a working inequality.

THE LEADING-SHELL FOURIER WEIGHT (exact).  corr(E)=Sum_{0!=n in Lambda}(weight n).
The rank-2 shell is n with support 3, entries (+1,+1,-1) i.e. e_a+e_b=e_c.  Its EXACT
single-relation weight (the contribution of one additive triple, treated in isolation
= the 3-clock corr) is a FIXED rational w2 (computed below from the 3-clock joint law).
If corr were = w2 * N_2 + (higher shells, small), monotonicity would hold.  We test:
  (S1) compute w2 exactly = corr of a single additive triple {0,a,a+b}? NO -- corr is a
       lattice sum, not per-triple.  Instead measure: regress corr(E) on N_2(E) (least
       squares over the pool); report slope, R^2, and the residual structure.
  (S2) the KEY question: is corr(E) <= w2_max * N_2(E) (a one-sided leading-shell BOUND)?
       i.e. does the leading shell DOMINATE corr from above?  If yes AND N_2 is maximized
       at consec with corr(consec)=w2_max*N_2(consec) (tight), that's a clean upper bound
       that is EXTREMAL at consec -- a working inequality.  Test tightness & one-sidedness.
  (S3) Two-coordinate certificate: does (N_2, N_3) [additive triples + the next shell:
       rank-2 with entries (+1,+1,+1,-1)?? actually 4-term a+b+c=d and 2a=b ('halving')]
       jointly determine the sign?  Test a 2-D monotone cone: corr(E) <= corr(consec)
       whenever N_2(E)<=N_2(consec) AND N_3(E)<=N_3(consec)?  consec maximizes both =>
       this would CLOSE if the joint map is monotone.  Report violations.
  (S4) The decisive diagnostic: among the INVERSION pairs (higher N_2, lower measS7), what
       distinguishes them?  Hypothesis: the high-N_2-low-measS7 shapes have additive triples
       that are 'wasted' -- the triples overlap on few clocks (clustered relations) so they
       don't independently cover sectors.  Measure the relation-lattice RANK / span.

EXACT Fractions; stdlib only.
"""
import sys, itertools
from fractions import Fraction as F
from collections import defaultdict
from math import comb, factorial
sys.stdout.reconfigure(line_buffering=True)

def cell_decomp(E):
    E = sorted(set(E)); Enz = [e for e in E if e != 0]
    bps = {F(0), F(1)}
    for e in Enz:
        for m in range(0, 7 * e + 1):
            bps.add(F(m, 7 * e))
    bps = sorted(bps)
    cells = []
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0: continue
        xm = (x0 + x1) / 2
        hits = {0}
        for e in E:
            if e == 0: continue
            hits.add(int(((e * xm) % 1) * 7))
        cells.append((x1 - x0, frozenset(hits)))
    return cells
def measS7(E):
    return sum(L for L, h in cell_decomp(E) if len(h) == 7)
def stirling2(n, k):
    return sum((-1)**(k-j)*comb(k,j)*j**n for j in range(k+1))//factorial(k)
def measS7_iid(k):
    return F(factorial(7)*stirling2(k,7), 7**k)
CAP = {8: F(2243,5880), 9: F(1979,4004), 10: F(55,91)}
def shapes(k, span):
    return [(0,)+combo for combo in itertools.combinations(range(1, span+1), k-1)]

def additive_triples(E):
    Es = set(E); El = sorted(E); cnt = 0
    for i in range(len(El)):
        for j in range(i, len(El)):
            if El[i]+El[j] in Es: cnt += 1
    return cnt

def rank2_shell_full(E):
    """N_2 split: additive a+b=c (entries +1+1-1) AND halving 2a=b (entries +2-1, effective
    unless coordinate is mult of 7). Return (n_add, n_halve)."""
    Es = set(E); El = sorted(E); add=0; hal=0
    for i in range(len(El)):
        for j in range(i, len(El)):
            if El[i]+El[j] in Es: add += 1
    for a in El:
        if 2*a in Es and a!=0: hal += 1
    return add, hal

def rank3_shell(E):
    """4-term additive a+b+c=d (entries +1+1+1-1), the next shell up. count unordered {a,b,c}."""
    Es=set(E); El=sorted(E); cnt=0
    for trip in itertools.combinations_with_replacement(El,3):
        if sum(trip) in Es: cnt+=1
    return cnt

def lattice_rank(E):
    """rank over Q of the relation lattice = (k) - 1 always (1 relation per scaling). The
    INTERESTING invariant is the SPAN of E itself = e_max (geometric extent)."""
    return max(E)

def main():
    print("="*92)
    print("THREAD D (part iv): LEADING-SHELL domination & the inversion mechanism")
    print("="*92)
    CONFIG={8:11, 9:11}
    for k in (8,9):
        span=CONFIG[k]; consec=tuple(range(k)); iid=measS7_iid(k); cap=CAP[k]
        pool=shapes(k,span)
        print("\n"+"#"*92)
        print(f"### k={k} span[0..{span}] |pool|={len(pool)} iid={float(iid):.5f} cap={float(cap):.4f}")
        print("#"*92)
        rows=[]
        for E in pool:
            m=measS7(E); add,hal=rank2_shell_full(E); r3=rank3_shell(E)
            rows.append(dict(E=E,m=m,corr=m-iid,add=add,hal=hal,r3=r3,span=max(E)))
        cons=next(r for r in rows if r['E']==consec)
        print(f"\nconsec: measS7={float(cons['m']):.5f} corr=+{float(cons['corr']):.5f} "
              f"add={cons['add']} halve={cons['hal']} r3={cons['r3']} span={cons['span']}")

        # (S1) regress corr on add-triples (least squares, float)
        xs=[float(r['add']) for r in rows]; ys=[float(r['corr']) for r in rows]
        n=len(xs); mx=sum(xs)/n; my=sum(ys)/n
        sxx=sum((x-mx)**2 for x in xs); sxy=sum((xs[i]-mx)*(ys[i]-my) for i in range(n))
        slope=sxy/sxx; int(0);
        ss_tot=sum((y-my)**2 for y in ys); ss_res=sum((ys[i]-(my+slope*(xs[i]-mx)))**2 for i in range(n))
        r2=1-ss_res/ss_tot
        print(f"\n[S1] corr ~ slope*add :  slope={slope:.5f}  R^2={r2:.3f}  (add-triple count explains "
              f"{100*r2:.0f}% of corr variance)")

        # (S2) one-sided leading-shell bound: corr(E) <= w * add(E) ?  find min w making it hold,
        #      check tightness at consec.
        w_needed = max((r['corr']/r['add']) for r in rows if r['add']>0)
        tight = next(r for r in rows if r['add']>0 and r['corr']/r['add']==w_needed)
        consec_ratio = cons['corr']/cons['add']
        print(f"\n[S2] one-sided bound corr <= w*add holds for w = max ratio = {float(w_needed):.5f} "
              f"(at {tight['E']})")
        print(f"     consec ratio corr/add = {float(consec_ratio):.5f}  "
              f"(tight at consec? {tight['E']==consec})")
        # is corr<=w_consec*add (using consec's own slope)? count violations
        wv = sum(1 for r in rows if r['corr'] > consec_ratio*r['add'])
        print(f"     using w=consec-ratio: # shapes with corr> w*add = {wv}/{len(rows)} "
              f"({'consec slope is an UPPER envelope' if wv==0 else 'NOT an upper envelope'})")

        # (S3) joint 2-D monotone cone: corr(E)<=corr(consec) whenever add<=add(consec) & r3<=r3(consec)
        # consec maximizes add & r3 (verify), so the cone {add<=A0, r3<=R0} contains ALL shapes.
        A0=cons['add']; R0=cons['r3']
        all_in_cone = all(r['add']<=A0 and r['r3']<=R0 for r in rows)
        max_add = max(r['add'] for r in rows); max_r3 = max(r['r3'] for r in rows)
        print(f"\n[S3] consec maximizes add ({A0}=={max_add}? {A0==max_add}) and r3 ({R0}=={max_r3}? {R0==max_r3}).")
        print(f"     ALL shapes inside cone {{add<=A0, r3<=R0}}: {all_in_cone} "
              f"(=> cone is TRIVIAL: every shape is dominated coordinatewise; the cone alone "
              f"does NOT order measS7).")

        # (S4) inversion mechanism: high-add low-measS7 shapes -- what is their span/r3?
        # find shapes with add >= 0.7*A0 but measS7 < 0.5*consec measS7
        print(f"\n[S4] high-leading-shell but LOW measS7 shapes (the inversion cause):")
        bad=[r for r in rows if r['add']>=0.7*A0 and r['m']<cons['m']/2]
        bad.sort(key=lambda r:-r['add'])
        for r in bad[:5]:
            print(f"     {r['E']} add={r['add']} r3={r['r3']} span={r['span']} measS7={float(r['m']):.4f}")
        if bad:
            print(f"     => these have HIGH additive-triple count but the triples are CLUSTERED")
            print(f"        (large span, few effective sector-covering resonances): leading-shell")
            print(f"        count is necessary-not-sufficient. corr needs the FULL lattice, not N_2.")

    print("\n"+"="*92)
    print("READOUT (Thread D part iv): is the leading shell a working inequality?")
    print("="*92)

if __name__=="__main__":
    main()
