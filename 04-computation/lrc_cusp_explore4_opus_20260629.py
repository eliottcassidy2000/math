"""
BINDING STRUCTURE of the razor's edge (cyclotomic self-duality?).
At the 6 unit witnesses k/14 (k in units mod 14), which speeds are BINDING (||v k/14||=1/14)?
Claim: the 6 UNIT-speeds (v ≡ units mod 14) bind; non-unit speeds are slack (safe).
And: a mult-of-14 speed sits at the ORIGIN (danger) at every unit witness => KILLS the edge.
"""
from fractions import Fraction as F
from math import gcd
def nrm(x): f=x-(x.numerator//x.denominator); return min(f,1-f)
units=[k for k in range(1,14) if gcd(k,14)==1]
print(f"units mod 14 (phi(14)={len(units)}): {units}")
ap=list(range(1,14))
print("\nAP {1..13}: at each unit witness t=k/14, the BINDING runners (||v t||=1/14) and slack:")
allbind=set()
for k in units:
    t=F(k,14)
    bind=[v for v in ap if nrm(v*t)==F(1,14)]
    mn=min(nrm(v*t) for v in ap)
    allbind|=set(bind)
    print(f"   t={k}/14: min over runners = {mn} ; BINDING speeds = {sorted(bind)}")
print(f"\n   union of all binding speeds = {sorted(allbind)}")
print(f"   = the unit-residue speeds (v mod 14 in units): {sorted(v for v in ap if gcd(v,14) in (1,) or (v%14 in units))}")
print(f"   SLACK speeds (never binding, safe): {sorted(set(ap)-allbind)}  (these are v with 14|... or non-unit residue)")
print(f"   slack residues mod 14: {sorted(set(v%14 for v in (set(ap)-allbind)))} (non-units)\n")
# the covering KILL: a mult-of-14 speed at the unit witnesses
print("COVERING KILL: a speed w with 14|w sits at the ORIGIN at every unit witness t=k/14:")
for w in [14,28,84]:
    vals=[nrm(F(w*k,14)) for k in units]
    print(f"   w={w}: ||w*k/14|| over units = {vals}  => all 0 (danger) => KILLS the razor's edge")
print("\n=> the razor's edge (6 unit witnesses) is HELD by the 6 unit-speeds and KILLED by any mult-of-14.")
print("   Covering forces a mult-of-14 => the edge dies => M>1/14. The self-dual cyclotomic obstruction.")
# verify: which slack speeds can be doubled/replaced keeping M=1/14?
print("\nslack-substitution: replace a slack speed by a multiple safe at all unit witnesses -> still M=1/14?")
import numpy as np
def M_grid(S,Q=200000):
    t=np.arange(1,Q)/Q; f=np.ones(Q-1)
    for v in S: f=np.minimum(f,np.abs(((v*t+0.5)%1)-0.5))
    return f.max()
for drop,add in [(12,24),(8,16),(6,12),(2,16),(4,8),(7,21)]:
    S=[x for x in ap if x!=drop]+[add]
    if len(set(S))==13:
        print(f"   drop {drop}, add {add}: M={M_grid(S):.6f}  (slack {drop}: residue {drop%14}, unit={gcd(drop,14)==1})")
