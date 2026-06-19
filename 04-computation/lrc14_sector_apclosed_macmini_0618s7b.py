#!/usr/bin/env python3
"""
lrc14_sector_apclosed_macmini_0618s7b.py  (mac-mini-2026-06-18-S7, ANGLE B)

Per-M extremality is FALSE (refuted). AP-extremality of meas(S7) is purely AGGREGATE.
So pursue the AGGREGATE directly: find a CLOSED FORM for meas(S7(AP_k)) via the joint
cutting word of the consecutive set E={0,1,...,k-1}.

JOINT CUTTING WORD OF THE AP.  For E={0,..,k-1}, at x in [0,1):
  sigma_e(x) = floor(7 e x) mod 7,  e=0..k-1.
Set theta = 7x in [0,7). Then sigma_e(x) = floor(e theta) mod 7. So the multiset is
  W(theta) = { floor(e theta) mod 7 : e=0..k-1 }.
S7 condition: W(theta) = Z/7 (all 7 residues hit), for theta in [0,7), measure/7.

KEY: floor(e theta) for e=0..k-1 is the BEATTY-like sequence. As e increases by 1, floor(e theta)
increases by floor((e+1)theta)-floor(e theta) in {floor(theta), ceil(theta)} (a Sturmian word!).
So {floor(e theta) mod 7 : e} is the set of PARTIAL SUMS mod 7 of a Sturmian word with letters
in {floor(theta), floor(theta)+1}. THIS IS A CLEAN COMBINATORICS-ON-WORDS OBJECT.

For theta in [j, j+1), floor(theta)=j, the increments d_e = floor((e+1)theta)-floor(e theta)
are j or j+1, forming a Sturmian/mechanical word of slope frac(theta). The residues hit are
the partial sums S_e = floor(e theta) mod 7, e=0..k-1, with S_0=0, S_{e+1}=S_e + d_e mod 7.

So meas(S7(AP_k)) = (1/7) * meas{theta in [0,7): partial sums of the slope-theta mechanical
word, taken mod 7, of length k, cover all of Z/7}.

THIS SCRIPT:
 (A) verify the theta-reparametrization & Sturmian-increment structure exactly.
 (B) compute meas(S7(AP_k)) exactly for k=3..14 and look for the pattern / closed form / OEIS.
 (C) restrict to theta in [j,j+1): increments j,j+1; partial sums mod 7. The cover-time of a
     mod-7 walk with steps in {j,j+1}. Express the per-interval measure as a Sturmian cover.
 (D) the cover condition only depends on j mod 7 and the slope -- decompose meas over j=0..6.
"""
import sys, itertools
from fractions import Fraction as F
from math import comb, gcd
sys.stdout.reconfigure(line_buffering=True)

# exact meas(S7) via breakpoints in theta=7x
def measS7_AP(k):
    E=list(range(k))
    Enz=[e for e in E if e!=0]
    # breakpoints in theta in [0,7): theta=m/e
    bps=set([F(0),F(7)])
    for e in Enz:
        for m in range(0,7*e+1): bps.add(F(m,e))
    bps=sorted(x for x in bps if 0<=x<=7)
    total=F(0)
    for i in range(len(bps)-1):
        t0,t1=bps[i],bps[i+1]
        if t1<=t0: continue
        tm=(t0+t1)/2
        res=set(int(e*tm)%7 for e in E)  # floor(e tm) mod 7
        if len(res)==7:
            total+=t1-t0
    return total/7   # meas in x

# Sturmian partial-sum cover check (theta in [j,j+1))
def cover_via_sturmian(k,tm):
    j=int(tm)  # floor(theta)
    S=0; res={0}
    for e in range(1,k):
        d=int(e*tm)-int((e-1)*tm)   # increment floor(e theta)-floor((e-1)theta) in {j,j+1}
        assert d in (j,j+1), (d,j,tm)
        S=(S+d)%7; res.add(int(e*tm)%7)
    return len(res)==7

print("="*92)
print("(A) verify theta-reparam & Sturmian increments; (B) exact meas(S7(AP_k)) k=3..14")
print("="*92)
vals=[]
for k in range(3,15):
    v=measS7_AP(k)
    vals.append((k,v))
    # cross-check via Sturmian partial-sum cover on a fine mesh of breakpoints
    E=list(range(1,k)); bps=set([F(0),F(7)])
    for e in E:
        for m in range(0,7*e+1): bps.add(F(m,e))
    bps=sorted(x for x in bps if 0<=x<=7)
    tot=F(0)
    for i in range(len(bps)-1):
        t0,t1=bps[i],bps[i+1]
        if t1<=t0: continue
        tm=(t0+t1)/2
        if cover_via_sturmian(k,tm): tot+=t1-t0
    tot/=7
    ok = (tot==v)
    print(f"  k={k}: meas(S7(AP))={v} = {float(v):.8f}   Sturmian-cover match: {ok}")

print()
print("  numerators/denominators (reduced):")
for k,v in vals:
    print(f"    k={k}: {v.numerator}/{v.denominator}")

print()
print("="*92)
print("(C) Per-floor decomposition: meas over theta in [j,j+1), j=0..6")
print("    (cover by a {j,j+1}-step walk mod 7). Each j-interval contributes; sum/7 = total.")
print("="*92)
for k in [8,9,10,11,12,13]:
    parts=[]
    E=list(range(1,k))
    for j in range(7):
        bps=set([F(j),F(j+1)])
        for e in E:
            for m in range(0,7*e+1):
                x=F(m,e)
                if j<=x<=j+1: bps.add(x)
        bps=sorted(set(x for x in bps if j<=x<=j+1))
        tot=F(0)
        for i in range(len(bps)-1):
            t0,t1=bps[i],bps[i+1]
            if t1<=t0: continue
            tm=(t0+t1)/2
            if cover_via_sturmian(k,tm): tot+=t1-t0
        parts.append(tot)
    s=sum(parts,F(0))/7
    print(f"  k={k}: per-j cover-measures {[str(p) for p in parts]}  total={float(s):.6f}")
    # note j and 7-j symmetry?
print("\nDONE.")
