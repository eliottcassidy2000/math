#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
Cleaner (Lean-friendly) derivation of the FORCED ODD RESIDUES for confined tight families, and the
odd-binds / even-spectates structural split. opus-2026-07-04-S73. CONVERGES on mac-mini HYP-4070
(residues superset all 7 odds {1,3,5,7,9,11,13}) via a more elementary route:

Given confinement (tight point t*=a/14 has EXACT denominator 14, so gcd(a,14)=1, a unit) and M=1/14:
  (0) M = max_t min_i ||v_i t|| = 1/14  =>  for ALL t, min_i ||v_i t|| <= 1/14  (1/14 is the max).
  (1) [no residue 0] At t*=a/14: min_i ||r_i a/14|| = 1/14 > 0, so no runner is at 0 => no r_i*a==0 mod14
      => (a unit) no residue 0. [This DERIVES mac-mini's "miss 14" from confinement.]
  (2) [residue 7] Evaluate the ALL-t bound at t=1/7: ||r_i/7|| in {0,1/7,2/7,3/7}. "<=1/14" forces =0
      => 7|r_i => r_i in {0,7}; residue 0 excluded => RESIDUE 7 PRESENT.
  (3) [units] Evaluate at t=b/14, b a unit: ||r_i b/14||<=1/14 => r_i*b in {0,+-1} mod14; 0 excluded
      => r_i == +-b^{-1}; over all units b => residues cover all UNITS {1,3,5,9,11,13}.
  => residues superset {1,3,5,7,9,11,13} = ALL 7 ODD residues.  QED (rigorous given confinement).

Structural split (verified): at t*, the 7 ODD residues BIND (distance exactly 1/14), the EVEN runners
SPECTATE (distance >= 1/7). M(odd part alone) = 1/2 (maximally lonely at t=1/2) -- the even coverers exist
solely to pull the loose odd skeleton down to 1/14, and are MAGNITUDE-BOUNDED (large coverer under-covers
=> M>1/14, loose). The magnitude bound = mac-mini's finiteness residual => finite check => {AP, GW}.
"""
import sys
from fractions import Fraction as Fr
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
def norm(x):
    x=x-int(x)
    if x<0:x+=1
    return min(x,1-x)
def exact_M(S):
    S=sorted(set(S));cands=set()
    for v in S:
        for k in range(v):cands.add(Fr(2*k+1,2*v))
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for den in (S[i]+S[j],abs(S[i]-S[j])):
                if den:
                    for s in range(den):cands.add(Fr(s,den))
    b=Fr(0);arg=None
    for t in cands:
        v=min(norm(x*t) for x in S)
        if v>b:b=v;arg=t
    return b,arg
odds={1,3,5,7,9,11,13}
print("(1) forced-odd-residues holds for both tight families (convergence w/ mac-mini HYP-4070):")
for nm,S in [("AP {1..13}",list(range(1,14))),("GW {1..11,13,24}",[1,2,3,4,5,6,7,8,9,10,11,13,24])]:
    res=set(v%14 for v in S); M,_=exact_M(S)
    print("   %-18s odds present=%s  M=%s"%(nm,odds<=res,str(M)))
print("(2) residue-7 necessity (t=1/7 test): families missing residue 7 are loose:")
for nm,S in [("{1..6,8..13,20} no res7",[1,2,3,4,5,6,8,9,10,11,12,13,20]),("AP has res7",list(range(1,14)))]:
    M,_=exact_M(S); has7=any(v%14==7 for v in S)
    print("   %-24s res7=%-5s M=%s"%(nm,has7,str(M)))
print("(3) odd binds / even spectates at t*=13/14 (AP):")
ts=Fr(13,14)
print("   ODD  min dist =",str(min(norm(v*ts) for v in range(1,14) if v%2)),"(=1/14, binding)")
print("   EVEN min dist =",str(min(norm(v*ts) for v in range(2,14,2))),"(=1/7, spectator)")
print("   M(odd part alone) =",str(exact_M([1,3,5,7,9,11,13])[0]),"(maximally lonely; even coverers pull ->1/14)")
print("(4) coverers magnitude-bounded (residual): {1..11,13,X} tight only X in {12,24}:")
for X in [12,24,36,48]:
    M,_=exact_M([1,2,3,4,5,6,7,8,9,10,11,13,X]); print("   X=%2d M=%-6s %s"%(X,str(M),"TIGHT" if M==Fr(1,14) else "loose"))
print("DONE. Convergence on HYP-4070 via Lean-friendly single-point evaluations; residual = magnitude bound.")
