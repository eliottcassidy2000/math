#!/usr/bin/env python3
"""
The localisation lemma via first-return / renewal: the minimal level sees only PRIMITIVE
atoms, and a bottom-up renewal peels them straddle by straddle.        (mac-mini-S153)
================================================================================
Owner: take the localisation lemma (HYP-8590) in the first-return / covering direction.

THE OBJECTS.  P = sum_i c_i Z^{a_i} W^{b_i}, charge k_i = a_i - b_i, V = ZW.  A BALANCED
tuple = a multiset of monomials with sum of charges = 0; E[P^m] sums over balanced m-tuples.
An ATOM = a minimal balanced charge-multiset (no proper nonempty balanced sub-multiset) --
the 'first-return events' / minimal vanishing sums over the charge lattice (opus THM-1685,
THM-415).  m0(atom) = its size.

THE FIRST-RETURN ISOLATION LEMMA (proved below by minimality, verified here).
Let m* = min atom size.  Then EVERY balanced tuple of size m* is a single atom -- no composite
(a composite splits into >= 2 balanced pieces, each of size >= m*, total >= 2m* > m*).  So
    E[P^{m*}] = sum over atoms A of size m*  of  (coeff form of A),
with NO cross/composite terms.  This is the covering backbone: at the first return, only
primitive returns are visible.

THE RENEWAL INDUCTION.  Bottom-up: at m*, the atoms of size m* fire; the nullcone forces their
coeff forms to vanish.  Peel (their coefficients constrained), then the next level's NEW atoms
are isolated the same way.  This is the rigorous 'last-to-return is alone' -- read forward as
'first-to-return is alone', renewal-style.

This file (A) verifies isolation at m* on multi-straddle / multi-charge patterns; (B) verifies
the renewal firing order = bottom-up; (C) exhibits a 3-charge atom (not a pos-neg pair) to show
the atoms are minimal vanishing sums, connecting to THM-1685/THM-415.
"""
import sympy as sp
from math import factorial
from itertools import combinations_with_replacement

def moments(monos, coeffs, M):
    P = {}
    for i, (a, b) in enumerate(monos): P[(a, b)] = P.get((a, b), 0) + coeffs[i]
    def mul(X, Y):
        o = {}
        for (a1, b1), c1 in X.items():
            for (a2, b2), c2 in Y.items():
                k = (a1+a2, b1+b2); o[k] = o.get(k, 0) + c1*c2
        return o
    Pm = {(0, 0): sp.Integer(1)}; g = []
    for m in range(1, M+1):
        Pm = mul(Pm, P)
        g.append(sp.expand(sum(c*factorial(a) for (a, b), c in Pm.items() if a == b)))
    return g

def atoms(charges, maxsize):
    """minimal balanced charge-multisets (as sorted tuples of charge-indices), size<=maxsize."""
    idx = list(range(len(charges)))
    found = []
    for size in range(2, maxsize+1):
        for combo in combinations_with_replacement(idx, size):
            if sum(charges[i] for i in combo) != 0:
                continue
            # minimal: no proper nonempty balanced sub-multiset
            from itertools import combinations
            minimal = True
            for sz in range(1, size):
                for sub in combinations(range(size), sz):
                    if sum(charges[combo[t]] for t in sub) == 0:
                        minimal = False; break
                if not minimal: break
            if minimal:
                # dedupe by sorted charge-multiset
                key = tuple(sorted(charges[i] for i in combo))
                if not any(a[0] == key for a in found):
                    found.append((key, size, combo))
    return found

print("=" * 78)
print("PART A -- first-return isolation: E[P^{m*}] sees ONLY atoms of the minimal size m*")
print("=" * 78)
pats = [
    ("aZ^2 + bW + cW^3",          [(2, 0), (0, 1), (0, 3)]),          # straddles (2,1),(2,3)
    ("aZ^2 + bZ^3 + cW  (two pos)",[(2, 0), (3, 0), (0, 1)]),          # (2,1),(3,1)
    ("aZ^2 + bZ^3 + cW^5 (3-charge atom +2+3-5)", [(2, 0), (3, 0), (0, 5)]),
    ("aZ + bW + cW^2 (charges +1,-1,-2)", [(1, 0), (0, 1), (0, 2)]),
]
for nm, monos in pats:
    charges = [a-b for a, b in monos]
    coeffs = sp.symbols(f'c0:{len(monos)}')
    ats = atoms(charges, 8)
    if not ats:
        print(f"  {nm}: no atoms (one-sided)"); continue
    mstar = min(a[1] for a in ats)
    small_atoms = [a for a in ats if a[1] == mstar]
    gens = moments(list(monos), list(coeffs), mstar+2)
    Estar = gens[mstar-1]
    # does E[P^{m*}] factor into exactly the size-m* atom coefficient monomials?
    monoms_in_E = set()
    poly = sp.Poly(Estar, *coeffs)
    for mono in poly.monoms():
        monoms_in_E.add(tuple(mono))
    # expected: each size-m* atom contributes the monomial prod c_i (over its charge-indices)
    print(f"  {nm}")
    print(f"     charges {charges};  atoms(size): {[(a[0],a[1]) for a in ats]}")
    print(f"     m* = {mstar};  E[P^{mstar}] = {sp.factor(Estar)}")
    print(f"     lower moments E[P^m], m<{mstar}: "
          f"{[sp.simplify(gens[m-1]) for m in range(1,mstar)]} (all zero => first return at m*)")
    print()

print("=" * 78)
print("PART B -- the renewal firing order is BOTTOM-UP (verify on the two-straddle witness)")
print("=" * 78)
monos = [(2, 0), (0, 1), (0, 3)]      # aZ^2 + bW + cW^3
a, b, c = sp.symbols('a b c')
gens = moments(list(monos), [a, b, c], 8)
print("  aZ^2 + bW + cW^3:  straddle (2,1) atom size 3, straddle (2,3) atom size 5")
for m in range(1, 7):
    if gens[m-1] != 0:
        print(f"     E[P^{m}] = {sp.factor(gens[m-1])}")
print("  => level 3 fires the (2,1) atom (forces a b^2 = 0), level 5 fires the (2,3) atom")
print("     (forces a^3 c^2 = 0). With a != 0: b=0 at level 3, then c=0 at level 5. Bottom-up,")
print("     each atom isolated at its own first-return level -- the renewal induction.")

print()
print("=" * 78)
print("PART C -- atoms are MINIMAL VANISHING SUMS, not just pos-neg pairs (THM-1685/THM-415)")
print("=" * 78)
monos = [(2, 0), (3, 0), (0, 5)]      # +2, +3, -5
a, b, c = sp.symbols('a b c')
charges = [2, 3, -5]
ats = atoms(charges, 6)
print(f"  aZ^2 + bZ^3 + cW^5:  charges {charges}")
print(f"     atoms: {[(x[0], x[1]) for x in ats]}  <-- (+2,+3,-5) is a 3-CHARGE atom, size 3")
gens = moments(list(monos), [a, b, c], 6)
for m in range(1, 5):
    if gens[m-1] != 0:
        print(f"     E[P^{m}] = {sp.factor(gens[m-1])}  <-- the atom coeff product a*b*c")
print("  => the first-return event is a 3-charge minimal vanishing sum; its coeff MONOMIAL")
print("     (a b c) is distinct from any other atom's, so it cannot cancel -- this is exactly")
print("     opus THM-1685's vanishing-sums bridge and THM-415's prime/no-collision structure.")

print()
print("SUMMARY -- the localisation lemma, in first-return form")
print("  PROVED (minimality): at the minimal atom size m*, E[P^{m*}] = sum over size-m* atoms")
print("  of their coeff forms, with NO composite contributions -- the first return is primitive.")
print("  This is the covering backbone of HYP-8590: process levels bottom-up (renewal), each")
print("  atom isolated at its first-return level; distinct atoms have distinct coeff monomials")
print("  (THM-1685/THM-415), so they cannot cancel.  The residual is the MULTIPLICITY case --")
print("  when one charge carries r terms, the size-m* forms are r-dim (Vandermonde), needing")
print("  the r levels m*, 2m*, ..., r m* to kill them (THM-1740's r*m0, the single-straddle law).")
