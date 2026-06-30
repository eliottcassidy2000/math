#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The primitive covering-min is a STERN-BROCOT semiconvergent [0;n-1,k]; the irregular core a(n); the n=12
transition to 1/(n-1). mac-mini-2026-06-30-S52. (ILP data from covering_min_ip_v2 / _extended / _confirm.)
"""
from __future__ import annotations
import functools
from fractions import Fraction as F
print = functools.partial(print, flush=True)

# pinned primitive covering-min (ILP, scipy.milp; n=9 cross-checked vs 2M-set exhaustive; n=12,14 at V=72)
MPRIM = {7: F(2, 13), 8: F(2, 15), 9: F(4, 33), 10: F(4, 37), 11: F(3, 31),
         12: F(1, 11), 13: F(1, 12), 14: F(1, 13), 15: F(1, 14)}


def sc(n, k):  # the Stern-Brocot ray [0; n-1, k] = k/((n-1)k+1)
    return F(k, (n-1)*k+1)


def main():
    print("THE STERN-BROCOT RAY [0; n-1, k] = k/((n-1)k+1):  floor(k=1)=1/n  ...  top(k->inf)=1/(n-1)")
    print("  the floor, the covering-min, AND the construction n/Phi_6 all live on this ONE ray:\n")
    print(f"  {'n':>3} {'floor 1/n=[0;n-1,1]':>20} {'M_prim=[0;n-1,a]':>18} {'a(n)':>5} {'constr n/Phi6=[0;n-1,n]':>24} {'margin':>9}")
    for n in sorted(MPRIM):
        M = MPRIM[n]
        # depth a: M = [0;n-1,a] => a = M.num if Farey-nbr of 1/(n-1); for n>=12 M=1/(n-1) (k->inf)
        a = M.numerator if (M.denominator - (n-1)*M.numerator == 1) else "inf(=1/(n-1))"
        constr = F(n, n*n-n+1)
        print(f"  {n:>3} {str(F(1,n)):>20} {str(M):>18} {str(a):>5} {str(constr):>24} {str(M-F(1,n)):>9}")
    print()
    print("KEY FACTS:")
    print("  * M_prim is a FAREY NEIGHBOR of 1/(n-1) (den-(n-1)num=1) for n=7..11 -> M_prim=[0;n-1,a(n)].")
    print("  * IRREGULAR CORE a(n) = 2,2,4,4,3 (n=7..11) = the smallest Stern-Brocot depth achievable by a")
    print("    primitive covering set; achievability is NON-monotonic in k (a=4 achievable at n=9 but a=2,3 NOT).")
    print("  * TRANSITION at n=12: the spread/semiconvergent family DIES (no proper [0;n-1,k] achievable);")
    print("    M_prim = 1/(n-1) (the convergent, k->inf) = the (n-1)-consecutive LRC value (an n->n-1 FOLD).")
    print("  * So M_prim CLEAN for n>=12: 1/(n-1), margin 1/(n(n-1)). LRC14 hard core = 1/13 (beats THM-523's")
    print("    empirical 7/89; > 1/14 with margin 1/182). This PINS HYP-2566 looseness = 1/(n(n-1)) for n>=12.")
    print()
    print("  Unification on the ray [0;n-1,k]:  k=1 -> 1/n (floor) ; k=a(n) -> covering-min ; k=n -> n/Phi_6")
    print("  (construction, S47/HYP-3725) ; k->inf -> 1/(n-1). Covering-min = smallest achievable k>1.")
    # verify floor and construction are on the ray
    print("\n  verify: [0;n-1,1]=1/n and [0;n-1,n]=n/Phi6 for n=14:",
          sc(14, 1), "= 1/14?", sc(14, 1) == F(1, 14), "|", sc(14, 14), "= 14/183?", sc(14, 14) == F(14, 183))


if __name__ == "__main__":
    main()
