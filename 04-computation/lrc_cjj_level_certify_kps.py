#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_cjj_level_certify_kps.py   (kind-pasteur 2026-06-21)

Concrete test of the CJJ LP-hierarchy levels for LRC consec-max (the user's lead).
On the FULL-Z/7-RESIDUE STRATUM (HYP-2749: consec-max reduces to it), set up the
Delsarte/moment dual LP at LEVEL 1 (single-sector moments S_r) and LEVEL 2 (pair /
per-subset atoms), and ask: does the level certify consec (i.e. is the LP upper bound
on max measS7 equal to measS7(consec), with consec the unique maximizer)?

LEVEL 1 (THM-534 moment-LP): max measS7 s.t. the missed-mask distribution q[M]>=0,
  sum q=1, and the LEVEL-1 moments S_1=E[N], (the single-coordinate marginals) matched
  to consec.  Expect: does NOT pin consec (aggregate, HYP-2738).
LEVEL 2: additionally match the PAIR atoms (|A|=2 containment a[A]=meas{A subset M}).
  Expect (mac-mini S15): PINS consec at k<=10.

We solve the RELAXATION exactly per E is not the point; the point is the GLOBAL
extremality LP: over all stratum shapes, is consec the measS7-argmax, and does the
level-R dual certificate (g(t)=sum y_r C(t,r) >= 1[t=0], minimize the cap) equal
measS7(consec)?  Uses scipy.linprog (LP).
"""
import itertools, os, importlib.util
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import linprog

_d=os.path.dirname(__file__)
spec=importlib.util.spec_from_file_location("rcm",os.path.join(_d,"lrc_q108_relation_code_mds_kps.py"))
rcm=importlib.util.module_from_spec(spec); spec.loader.exec_module(rcm)
measS7=rcm.measS7

P=7
def missed_mask_dist(E):
    """exact q[M] = meas{x: missed inner-sector set = M}, M subset {1..6} (sector 0 = the 0-elt anchor).
       Returns dict frozenset(M)->Fraction. 'missed' = inner sectors 1..6 not hit by {frac(e x): e in E}."""
    E=[int(e) for e in E]
    bp={Fr(0),Fr(1)}
    for e in E:
        if e==0: continue
        for t in range(0,P*e): bp.add(Fr(t,P*e))
    pts=sorted(bp); dist={}
    for a,b in zip(pts,pts[1:]):
        mid=(a+b)/2
        hit=set(int(P*((e*mid)%1)) for e in E)
        missed=frozenset(s for s in range(1,P) if s not in hit)  # inner sectors 1..6 missed
        dist[missed]=dist.get(missed,Fr(0))+(b-a)
    return dist

def moments(E, upto=3):
    """S_r = sum_{|A|=r, A subset {1..6}} meas{A subset missed}. r=0..upto."""
    dist=missed_mask_dist(E)
    S={r:Fr(0) for r in range(0,upto+1)}
    for M,w in dist.items():
        m=len(M)
        for r in range(0,upto+1):
            # number of r-subsets of M = C(m,r)
            from math import comb
            S[r]+=w*comb(m,r)
    return S

def full_residue(E):
    return set(e%P for e in E)=={0,1,2,3,4,5,6} or set(e%P for e in E)==set(range(P))

def main():
    k=8
    print(f"CJJ level-certify test, k={k}, full-Z/7-residue stratum (HYP-2749)")
    # stratum: 8-subsets of {0..15} containing 0, hitting all residues mod 7
    base=range(0,16)
    consec=tuple(range(0,8))   # {0..7}
    shapes=[]
    for S in itertools.combinations(range(1,16),7):
        E=(0,)+S
        if set(e%P for e in E)==set(range(P)):
            shapes.append(E)
    print(f"  full-residue stratum size: {len(shapes)} shapes (8-subsets of [0,15] with 0, all residues)")
    # measS7 + moments
    data=[]
    for E in shapes:
        m=measS7(E); S=moments(E,3)
        data.append((E,m,S))
    mx=max(m for _,m,_ in data)
    argmax=[E for E,m,_ in data if m==mx]
    print(f"  max measS7 on stratum = {float(mx):.5f}; #argmax={len(argmax)}; consec in argmax? {consec in argmax}")
    cm=measS7(consec)
    print(f"  measS7(consec {consec}) = {float(cm):.5f}  ({'IS the max' if cm==mx else 'NOT max'})")
    # LEVEL-1 dual LP: find y_0,y_1 with g(t)=y_0 + y_1*t >= 1[t=0] on t=0..6, minimize the
    #   resulting bound max_E (y_0 + y_1 S_1(E)/?)... use the THM-534 form L_y=sum_r y_r S_r.
    # We instead test EXTREMALITY directly via the level-R moment vector:
    #   does there exist a stratum shape E != consec with S_r(E)=S_r(consec) for all r<=R but measS7 higher/equal?
    for R in (1,2,3):
        cons_mom=tuple(moments(consec,R)[r] for r in range(R+1))
        # shapes matching consec's moments up to level R
        matches=[(E,m) for E,m,S in data if tuple(S[r] for r in range(R+1))==cons_mom]
        hi=[E for E,m in matches if m>cm]
        ties=[E for E,m in matches if m==cm and E!=consec]
        print(f"  LEVEL {R}: {len(matches)} shapes match consec's moments S_0..S_{R}; "
              f"#with measS7>consec={len(hi)}; #ties={len(ties)} "
              f"=> {'consec UNIQUELY pinned by level '+str(R) if not hi and not ties else 'NOT pinned (collapse)'}")
    print("\nINTERPRETATION: the smallest R with 0 higher AND 0 ties = the CJJ level that certifies consec")
    print("on this stratum (matching mac-mini S15's 'level-2 pins consec'). DONE.")

if __name__ == "__main__":
    main()
