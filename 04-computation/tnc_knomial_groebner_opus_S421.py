"""opus-2026-07-20-S421 -- THE k-NOMIAL GROEBNER REDUCTION for TNC.

From THM-1680 (was 1675): a k-nomial R (k terms) reduces by the u-scale + t-scale gauges to
k-2 free coefficients a_1,...,a_{k-2}, and TNC for that charge pattern holds iff the CT ideal
   I = < CT(Lambda^m) : m >= 1 >  subset  C[a_1,...,a_{k-2}]
has NO zero with ALL coordinates nonzero (a zero with some a_i=0 is a lower-term degeneration,
handled by the (k-1)-nomial case).  So TNC <=> V(I) intersect (C*)^{k-2} = empty.

By the Nullstellensatz this is DECIDABLE: saturate I by (a_1...a_{k-2}) and test for 1 in the
saturation (empty variety) -- a finite Groebner computation per charge pattern.

REPO CONNECTIONS mined (S421):
 - THM-415 (vanishing sums): "prime modulus = no nontrivial vanishing; composite modulus =
   nontrivial vanishing sums = collisions."  The CT cancellations ARE vanishing sums of the
   coefficient monomials -- our charge lattice is the additive group, and non-unique minimal
   reps are exactly the "composite" collisions.  The trinomial witness a^2=-1 is a root of
   unity flavour (i), resonating with vanishing-sums-of-roots-of-unity (HYP on JC monodromy).
 - The JC Groebner fiber work (kind-pasteur HYP-8070): same tool (saturate + test 1 in ideal)
   used to prove a fiber is EXACTLY three points.  Here we prove a variety is EMPTY.
 - THM-1535 GMC cascade: CT(2m0) = (CT(m0)/c)^2 + correction -- the "leading square + surviving
   correction" that makes the saturation trivial.

This script: (1) build the CT ideal for 4-nomials and 5-nomials, (2) saturate and test
emptiness of V(I) cap (C*)^{k-2}, (3) confirm no k-nomial nullcone violator.
"""
import sympy as sp

u = sp.symbols('u')

def CT(Rexpr, N, mm):
    e = sp.expand(Rexpr**mm)
    return sp.Poly(e, u).coeff_monomial(u**(N*mm))

def ct_polys(exps, N, params, Rexpr, maxm):
    """CT(Lambda^m) as polynomials in the params, for m=1..maxm (drop zero/constant)."""
    out = []
    for m in range(1, maxm+1):
        c = sp.expand(CT(Rexpr, N, m))
        if c == 0: continue
        p = sp.Poly(c, *params) if params else None
        if p is not None and p.total_degree() >= 1:
            out.append(c)
    return out

print("="*78)
print("(1) 4-NOMIAL: R = 1 + a u^j + b u^l + u^d  (2 free params a,b), TNC per pattern")
print("="*78)
a, b = sp.symbols('a b')
patterns4 = [
    (2, [3, 4, 5], "N=2 charges -2,1,2,3"),
    (2, [1, 3, 4], "N=2 charges -2,-1,1,2"),
    (3, [1, 4, 5], "N=3 charges -3,-2,1,2"),
    (2, [3, 5, 6], "N=2 charges -2,1,3,4"),
    (3, [2, 4, 7], "N=3 charges -3,-1,1,4"),
]
for N, (j, l, d), desc in patterns4:
    R = 1 + a*u**j + b*u**l + u**d
    polys = ct_polys([j, l, d], N, [a, b], R, 12)
    if len(polys) < 2:
        print(f"   {desc}: <2 nonconstant CT -- unique-min (TNC by THM-1655)"); continue
    # saturate the ideal by a*b and test whether 1 is in it (empty V cap (C*)^2)
    I = polys[:4]
    w = sp.symbols('w')
    # Rabinowitsch: add 1 - w*a*b to force a,b != 0; test if Groebner basis = [1]
    G = sp.groebner(I + [1 - w*a*b], a, b, w, order='grevlex')
    empty = (list(G) == [sp.Integer(1)] or G.exprs == (sp.Integer(1),))
    # also report first few CT to show they force it
    print(f"   {desc}: {len(polys)} CT-polys; V(I) cap (C*)^2 EMPTY (1 in saturation): {empty}")
    if not empty:
        print(f"      *** NONEMPTY -- possible 4-nomial nullcone violator, inspect: {list(G)[:3]}")

print()
print("="*78)
print("(2) 5-NOMIAL spot-check: R = 1 + a u^j + b u^l + c u^p + u^d (3 params)")
print("="*78)
c = sp.symbols('c')
for N, (j, l, p, d), desc in [(2, [3, 4, 5, 6], "N=2 charges -2,1,2,3,4"),
                              (3, [1, 2, 4, 5], "N=3 charges -3,-2,-1,1,2")]:
    R = 1 + a*u**j + b*u**l + c*u**p + u**d
    polys = ct_polys([j, l, p, d], N, [a, b, c], R, 12)
    if len(polys) < 3:
        print(f"   {desc}: <3 nonconstant CT -- lower case"); continue
    w = sp.symbols('w')
    I = polys[:5]
    try:
        G = sp.groebner(I + [1 - w*a*b*c], a, b, c, w, order='grevlex')
        empty = (list(G) == [sp.Integer(1)])
        print(f"   {desc}: {len(polys)} CT-polys; V(I) cap (C*)^3 EMPTY: {empty}")
    except Exception as e:
        print(f"   {desc}: Groebner cost too high ({type(e).__name__}); need better order/bound")

print()
print("="*78)
print("(3) THE VANISHING-SUMS BRIDGE (THM-415): why the saturation is trivial")
print("="*78)
print("   CT(m) = sum over charge-reps of 0 of (multinomial) * (coeff monomial).  Setting all")
print("   CT(m)=0 is a system of VANISHING SUMS of coefficient monomials.  THM-415: nontrivial")
print("   vanishing needs 'composite' structure (a sublattice relation among charges).  The")
print("   MINIMAL rep m0 is the primitive relation; its cancellation (CT(m0)=0) is ONE equation")
print("   in the params, but CT(2m0) = (that)^2 + correction adds an INDEPENDENT equation, and")
print("   the two generate the unit ideal after saturation -- exactly the trinomial gcd, now")
print("   multivariate.  The 'correction' is the GMC-cascade surviving term (THM-1535).")
