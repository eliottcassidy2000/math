#!/usr/bin/env python3
"""
The Laplace moment engine: charge-graded (weight-graded) powers of a binary form,
and the tournament / in-transitivity reading of the atoms.       (mac-mini-S154)
================================================================================
Owner: take the atom-covering target; work the representation theory of binary forms and
how it relates to tournaments (in/transitivity itself); work the Laplace moment engine.

REP THEORY.  P(Z,W), W = Zbar, is a binary form.  The diagonal torus T = {Z->lam Z, W->lam^-1 W}
acts, and the CHARGE k = a-b of a monomial Z^a W^b is its T-WEIGHT.  The Gaussian moment
functional E[Z^A W^B] = A! delta_{AB} is T-invariant and U(1)-invariant, so
    E[P^m] = L( weight-0 projection of P^m ) = L( CT_theta(P^m) ),   L(V^k)=k!, V=ZW,
i.e. the moment is the projection of P^m onto the TRIVIAL representation (weight 0 = charge 0),
then the radial Laplace.  CT_theta = averaging over U(1) = character projection (Peter-Weyl).

THE ENGINE.  Store P as { charge k : radial polynomial in V }.  Powers convolve on the charge
lattice (weight = tensor-product grading) while multiplying radial polynomials.  Take charge 0,
apply L.  This is far faster than expanding P^m in Z,W, and it makes the weight structure --
hence the atoms -- explicit.

TOURNAMENTS / IN-TRANSITIVITY.  A charge-0 tuple is a CLOSED WALK on Z (steps = charges, start
and end 0).  An ATOM (minimal vanishing sum) is a PRIMITIVE closed walk -- one that never
returns to 0 in between: a CYCLE in the charge flow.  P is one-sided (all charges same sign)
<=> the charge order is TRANSITIVE (partial sums are monotone, never return to 0) <=> NO cycles
<=> nullcone-harmless.  Two-sided-with-cycles = INTRANSITIVE.  The atom-covering is a disjoint-
cycle-packing question -- the same shape as the repo's odd-cycle collection formula for
tournaments.  in/transitivity is the pivot.
"""
import sympy as sp
from math import factorial

# ------------------------------------------------------------------ the engine
def make_P(monos, coeffs):
    """binary form P = sum coeff_i Z^{a_i} W^{b_i}  ->  { charge : {radial_deg(V): coeff} }.
       Z^a W^b = Z^{a-b} V^b (if a>=b) or W^{b-a} V^a (if a<b); radial deg in V = min(a,b)."""
    P = {}
    for (a, b), c in zip(monos, coeffs):
        k = a - b; r = min(a, b)
        P.setdefault(k, {}).setdefault(r, 0)
        P[k][r] += c
    return P

def cmul(A, B):
    """charge-convolution + radial (V) product of two {charge:{Vdeg:coeff}} dicts."""
    out = {}
    for k1, r1 in A.items():
        for k2, r2 in B.items():
            k = k1 + k2
            tgt = out.setdefault(k, {})
            for d1, c1 in r1.items():
                for d2, c2 in r2.items():
                    # Z^{k1} V^{d1} * Z^{k2} V^{d2}: charges add; but crossing +/- charges makes V.
                    # Z^{k1} and W^{|k2|} (k2<0) combine: min cancellation.  Handle via (a,b):
                    a1, b1 = (k1 + d1, d1) if k1 >= 0 else (d1, d1 - k1)
                    a2, b2 = (k2 + d2, d2) if k2 >= 0 else (d2, d2 - k2)
                    A2, B2 = a1 + a2, b1 + b2
                    kk = A2 - B2; rr = min(A2, B2)
                    tgt2 = out.setdefault(kk, {})
                    tgt2[rr] = tgt2.get(rr, 0) + c1 * c2
                    if k in out and not out[k]:
                        pass
    # clean empties
    return {k: {d: c for d, c in r.items() if c} for k, r in out.items() if any(r.values())}

def L_radial(rpoly):
    """L(sum c_d V^d) = sum c_d d!."""
    return sum(c * factorial(d) for d, c in rpoly.items())

def moment_engine(monos, coeffs, m):
    """E[P^m] via charge-graded convolution: take charge 0, apply L."""
    P = make_P(monos, coeffs)
    Pm = {0: {0: sp.Integer(1)}}
    for _ in range(m):
        Pm = cmul(Pm, P)
    return sp.expand(L_radial(Pm.get(0, {})))

# ------------------------------------------------------------------ direct check
def moment_direct(monos, coeffs, m):
    P = {}
    for (a, b), c in zip(monos, coeffs): P[(a, b)] = P.get((a, b), 0) + c
    def mul(X, Y):
        o = {}
        for (a1, b1), c1 in X.items():
            for (a2, b2), c2 in Y.items():
                k = (a1+a2, b1+b2); o[k] = o.get(k, 0) + c1*c2
        return o
    Pm = {(0, 0): sp.Integer(1)}
    for _ in range(m): Pm = mul(Pm, P)
    return sp.expand(sum(c*factorial(a) for (a, b), c in Pm.items() if a == b))

print("=" * 78)
print("PART A -- the engine agrees with the direct moment, and is charge-graded")
print("=" * 78)
tests = [
    ("aZ^2+bW+cW^3", [(2, 0), (0, 1), (0, 3)]),
    ("aZ^2+bZ^3+cW^5 (+2,+3,-5)", [(2, 0), (3, 0), (0, 5)]),
    ("aZ+bW+cW^2+d (+1,-1,-2,0)", [(1, 0), (0, 1), (0, 2), (0, 0)]),
]
for nm, monos in tests:
    coeffs = sp.symbols(f'c0:{len(monos)}')
    ok = all(sp.expand(moment_engine(list(monos), list(coeffs), m)
                        - moment_direct(list(monos), list(coeffs), m)) == 0
             for m in range(1, 7))
    print(f"  {nm}: engine == direct for m=1..6:  {ok}")

print()
print("=" * 78)
print("PART B -- REP THEORY: E[P^m] = L(weight-0 projection), and torus-invariance")
print("=" * 78)
print("  E[P^m] = L( charge-0 part of P^m ) = projection onto the TRIVIAL rep, then radial L.")
print("  The diagonal torus Z->lam Z, W->lam^-1 W scales charge-k by lam^k and PRESERVES E,")
print("  so E[P_lam^m] = E[P^m] for all lam (P_lam = P with charges scaled).  Check:")
lam = sp.Symbol('lam')
for nm, monos in tests[:2]:
    coeffs = sp.symbols(f'c0:{len(monos)}')
    # P_lam: multiply coeff of charge-k monomial by lam^k
    coeffs_lam = [coeffs[i] * lam**(a-b) for i, (a, b) in enumerate(monos)]
    inv = all(sp.expand(moment_direct(list(monos), coeffs_lam, m)
                        - moment_direct(list(monos), list(coeffs), m)) == 0
              for m in range(1, 6))
    print(f"  {nm}: E[P_lam^m] == E[P^m] (torus-invariant weight-0 projection):  {inv}")

print()
print("=" * 78)
print("PART C -- ATOMS ARE CHARGE CYCLES: transitive = one-sided, cyclic = intransitive")
print("=" * 78)
print("  A charge-0 tuple = a CLOSED WALK on Z (steps = charges).  An atom = a PRIMITIVE cycle")
print("  (never returns to 0 in between).  One-sided P: all charges one sign => partial sums")
print("  strictly monotone => NEVER return to 0 => NO cycle => TRANSITIVE, nullcone-harmless.")
print("  Two-sided => cycles exist => the moment ideal is generated by the cycle (atom) forms.")
print()
def is_onesided(charges):
    return all(c >= 0 for c in charges if c != 0) or all(c <= 0 for c in charges if c != 0)
def has_primitive_cycle(charges, maxlen=8):
    from itertools import combinations_with_replacement
    idx = [i for i in range(len(charges)) if charges[i] != 0]
    for size in range(2, maxlen+1):
        for combo in combinations_with_replacement(idx, size):
            if sum(charges[i] for i in combo) != 0: continue
            # primitive: no proper nonempty balanced subwalk (as a multiset)
            from itertools import combinations
            prim = True
            for sz in range(1, size):
                for sub in combinations(range(size), sz):
                    if sum(charges[combo[t]] for t in sub) == 0: prim = False; break
                if not prim: break
            if prim: return True
    return False
print(f"{'charges':>18} {'one-sided?':>11} {'has cycle (atom)?':>18} {'reading':>16}")
for charges in ([1, 2, 3], [1, -1], [2, 3, -5], [1, -1, -2], [1, 2, -3, -4], [2, 4, 6]):
    os = is_onesided(charges); cyc = has_primitive_cycle(charges)
    reading = "transitive" if os else ("intransitive" if cyc else "?")
    print(f"{str(charges):>18} {str(os):>11} {str(cyc):>18} {reading:>16}")
print("  => one-sided <=> transitive (no cycle) <=> nullcone-harmless;  two-sided <=> a cycle")
print("     exists <=> the atom-covering / disjoint-cycle question -- the OCF shape for the")
print("     charge lattice.  in/transitivity is exactly the pivot between the two.")

print()
print("=" * 78)
print("PART D -- pushing the atom-covering on multi-charge atoms with the engine")
print("=" * 78)
print("  For charges {+2,+3,-5}: size-3 atom (+2,+3,-5) gives abc; the renewal must also reach")
print("  the size-8 atom (5x+3, 3x-5) = b^5 c^3.  Check the full moment ideal cuts out one-sided.")
monos = [(2, 0), (3, 0), (0, 5)]
a, b, c = sp.symbols('a b c')
gens = [moment_engine(list(monos), [a, b, c], m) for m in range(1, 9)]
nz = [(m+1, sp.factor(g)) for m, g in enumerate(gens) if g != 0]
for m, g in nz: print(f"     E[P^{m}] = {g}")
w = sp.Symbol('w')
# V(all moments) cap two-sided: two-sided = (a or b nonzero) and (c nonzero) [+2,+3 pos, -5 neg]
allpass = True
for cp in (a, b):
    G = sp.groebner(list(g for g in gens if g != 0) + [1 - w*cp*c], a, b, c, w, order='grevlex')
    if not any(x.is_number and x != 0 for x in G.exprs): allpass = False
print(f"  V(moment ideal) cap two-sided (a or b nonzero, c nonzero) = empty:  {allpass}")
print("  => even with a multi-charge (3-charge) atom, the full ideal (reaching the size-8 atom")
print("     b^5 c^3) cuts out exactly the one-sided locus.  The covering completes across levels.")
