#!/usr/bin/env python3
"""
covering_min_farey_binding_klein.py  --  klein-2026-07-01-S62

FAREY-TESSELLATION / BINDING-FIRST method for the covering-min (menu item from HYP-3779), aimed at
HYP-3779's remaining gap: speeds > n(n-1). This check is SPEED-INDEPENDENT.

The covering-min is a FAREY NEIGHBOR of 1/(n-1) (HYP-3732): M = [0;n-1,j] = j/(j(n-1)+1), rung j, binding
modulus D = j(n-1)+1 (det-1 neighbor in the Farey tessellation; the Stern-Brocot ray from 1/(n-1)).
A covering set realizing rung j must, at the binding rotation a* (coprime to D), have EVERY speed avoid
the (j-1)-neighborhood of 0 mod D -- speeds lie in the SAFE BAND B(a*) = {r : dist(r a*, D) >= j} -- AND
be a covering set (a multiple of every q in {2..n}).  A multiple of q realizes residues = multiples of
gcd(q,D) mod D.  So:

  NECESSARY condition (rung j "binding-feasible"):  EXISTS a* coprime to D  s.t. for EVERY q in {2..n},
     { multiples of gcd(q,D) mod D }  intersects  B(a*).

This is a finite check mod D (D <= (n-1)^2), INDEPENDENT of speed size. If ALL rungs j=2..n-1 are
binding-INfeasible, then NO covering set (any speeds) realizes a rung < n, so covering-min = rung n =
construction -- closing the huge-speed gap (given the covering-min is a Farey rung, HYP-3732).

VALIDATION: known covering-min rungs a(n)=2,2,4,4,3 for n=7..11 MUST be binding-feasible.
"""
from math import gcd

def dist0(x, D):
    x %= D
    return min(x, D-x)

def safe_band(a_star, D, j):
    return set(r for r in range(D) if dist0(r*a_star, D) >= j)

def rung_binding_feasible(n, j):
    D = j*(n-1)+1
    for a_star in range(1, D):
        if gcd(a_star, D) != 1:
            continue
        B = safe_band(a_star, D, j)
        ok = True
        for q in range(2, n+1):
            g = gcd(q, D)
            mults = set((g*t) % D for t in range(D//g + 1))   # {multiples of g mod D}
            if mults.isdisjoint(B):
                ok = False; break
        if ok:
            return True, a_star, D
    return False, None, j*(n-1)+1

if __name__ == "__main__":
    print("FAREY BINDING-FIRST (speed-independent): which rungs j are binding-feasible? (necessary cond.)")
    print("known covering-min rung a(n)=2,2,4,4,3 for n=7..11 (validation); rung n = construction")
    for n in range(7, 15):
        feas = [j for j in range(2, n+1) if rung_binding_feasible(n, j)[0]]
        min_feas = min(feas) if feas else None
        # known a(n) for validation
        known = {7:2,8:2,9:4,10:4,11:3}.get(n)
        tag = ""
        if known is not None:
            tag = f"  [known a(n)={known}: {'OK' if known in feas else 'MISSING -- check FAILS'}]"
        # for n>=12: is any rung < n binding-feasible? if not, covmin=construction (any speed)
        low = [j for j in feas if j < n]
        verdict = ""
        if n >= 12:
            verdict = ("  => rungs<n binding-feasible: "+str(low)) if low else "  => NO rung<n binding-feasible -> covmin=construction (ALL speeds)!"
        print(f"  n={n:2d} (D_j=j*{n-1}+1): binding-feasible rungs {feas}; smallest={min_feas}{tag}{verdict}")
    print()
    print("If (n>=12) no rung<n is binding-feasible: the covering-min is the construction for ALL speeds")
    print("(closing HYP-3779's >n(n-1) gap), since a beater's rung would have to pass this necessary check.")
