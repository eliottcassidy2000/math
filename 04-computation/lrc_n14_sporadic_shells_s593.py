#!/usr/bin/env python3
"""Round 4: n=14 sporadic shell structure mod 27=3^3. Check V* tight; which antipodal
shells {a,27-a} are doubly-occupied / missed; confirm the prime-3 (multiple-of-3 shell)
mechanism. opus-2026-06-03-S593."""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
def dist(x): x%=1; return min(x,1-x)
def M(V):
    cands=set()
    for a,b in combinations(V,2):
        for D in (a+b,abs(a-b)):
            if D:
                for m in range(1,D): cands.add(F(m,D))
    for v in V: cands.add(F(1,2*v))
    return max(min(dist(v*t) for v in V) for t in cands)
def shells_mod(V,C):
    occ={}
    for v in V:
        r=v%C; s=min(r,C-r) if r else 0
        occ.setdefault(s,[]).append(v%C)
    return occ
def main():
    n=14; C=2*n-1  # 27
    Vstar=(1,2,3,4,5,6,7,8,9,10,11,13,24)
    print(f"n=14, C=2n-1={C}=3^3. Non-unit residues mod 27 (mult of 3): {[r for r in range(1,C) if gcd(r,C)>1]}")
    for V in [tuple(range(1,14)), Vstar]:
        m=M(V); occ=shells_mod(V,C)
        doubled=[s for s,lst in occ.items() if len(lst)>1]
        present=set(occ); allshells=set(range(1,n))
        missed=sorted(allshells-present)
        is_t = all(len(lst)==1 for lst in occ.values()) and len(occ)==n-1
        print(f"  V={V}: M={m}={float(m):.4f} (floor 1/14={1/14:.4f}, tight={m==F(1,14)}); transversal={is_t}")
        print(f"     doubly-occupied shells {doubled} (mult-of-3 shell? {[s%3==0 for s in doubled]}); missed shells {missed}")
        if doubled:
            for s in doubled: print(f"        shell {{{s},{C-s}}} occupied by {occ[s]}  (s mult of 3: {s%3==0})")
if __name__=='__main__': main()
