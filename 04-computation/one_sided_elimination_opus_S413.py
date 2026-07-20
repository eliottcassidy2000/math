"""opus-2026-07-20-S413 -- (A) THE ONE-SIDED CONJECTURE BY EXACT ELIMINATION, and
                            (B) SETTLING THE LAPLACE / GMC(1) LAYER.

(A) THE ONE-SIDED CONJECTURE, BOUNDED CHARGE SPAN.
Everything so far verified the n=2 nullcone conjecture on FINITE GRIDS of coefficients.
That is exactly the weakness I flagged in S412 ("name the domain of the sweep").  Here it is
replaced by EXACT ELIMINATION, which covers ALL complex coefficients at once.

Charge span {-1,0,+1}: every such P is
        P  =  z*A(s)  +  B(s)  +  zbar*C(s),        s = z*zbar
with A, B, C polynomials (charge +1, 0, -1 respectively).  "One-sided" means A == 0 or
C == 0.  So the conjecture, for this charge span, is exactly

        V( E[P^1], ..., E[P^M] )  =  {A = 0}  u  {C = 0}

as a variety over C.  We prove containment by RABINOWITSCH SATURATION: adjoin
1 - t*a_i*c_j and check the ideal becomes <1>, i.e. no point of the nullcone has both a
nonzero A-coefficient and a nonzero C-coefficient.

(B) THE LAPLACE / GMC(1) LAYER.
From THM-1540, the charge-0 part reduces to a ONE-DIMENSIONAL problem against the
EXPONENTIAL measure: find all g in C[s] with  int_0^inf g(s)^m e^{-s} ds = 0 for all m>=1.
Claim: only g = 0, hence the Laplace nullcone is TRIVIAL and GMC(1) for this measure holds
vacuously.  Proved here for all degrees by the saddle-point/phase-stabilisation argument,
checked numerically.
"""
import sympy as sp
from math import factorial

z, zb, t = sp.symbols('z zb t')

def E1(e):
    """E[z^a zbar^b] = a! delta_ab"""
    e = sp.expand(e)
    if e == 0: return sp.Integer(0)
    p = sp.Poly(e, z, zb); tot = 0
    for (a, b), c in zip(p.monoms(), p.coeffs()):
        if a == b: tot += c*factorial(a)
    return sp.expand(tot)

print("="*78)
print("(A) ONE-SIDED CONJECTURE, CHARGE SPAN {-1,0,+1}, BY EXACT ELIMINATION")
print("="*78)

def run(D, M):
    """A,B,C of degree <= D in s = z*zb; use M moment equations."""
    A = sp.symbols(f'a0:{D+1}'); B = sp.symbols(f'b0:{D+1}'); C = sp.symbols(f'c0:{D+1}')
    s = z*zb
    Ap = sum(A[k]*s**k for k in range(D+1))
    Bp = sum(B[k]*s**k for k in range(D+1))
    Cp = sum(C[k]*s**k for k in range(D+1))
    P = sp.expand(z*Ap + Bp + zb*Cp)
    eqs = []
    Pm = sp.Integer(1)
    for m in range(1, M+1):
        Pm = sp.expand(Pm*P)
        e = E1(Pm)
        if e != 0: eqs.append(sp.expand(e))
    gens = list(A)+list(B)+list(C)
    print(f"\n  degree D={D}: unknowns {len(gens)}, moment equations used {len(eqs)}")
    for i, e in enumerate(eqs[:3], 1):
        print(f"     E[P^{i}] = {e}")
    # saturation: is there a nullcone point with a_i != 0 AND c_j != 0 ?
    any_pt = False
    for i in range(D+1):
        for j in range(D+1):
            G = sp.groebner(eqs + [1 - t*A[i]*C[j]], *(gens+[t]), order='lex')
            trivial = (list(G.exprs) == [sp.Integer(1)])
            if not trivial:
                any_pt = True
                print(f"     a{i}*c{j}: saturation NOT trivial -> a possible point with "
                      f"A,C both nonzero (inspect)")
    if not any_pt:
        print(f"     ALL saturations give <1>  =>  every nullcone point has A == 0 or C == 0")
        print(f"     *** ONE-SIDED CONJECTURE PROVED for span {{-1,0,1}}, deg <= {D}, over ALL of C ***")
    return any_pt

run(0, 4)
run(1, 6)

print()
print("="*78)
print("(B) THE LAPLACE / GMC(1) LAYER:  int_0^inf g^m e^{-s} ds = 0 for all m => g = 0")
print("="*78)
S = sp.symbols('S')
def lap(g, m):
    e = sp.expand(g**m); p = sp.Poly(e, S)
    return sp.expand(sum(c*factorial(k) for (k,), c in zip(p.monoms(), p.coeffs())))

print("  PROOF (all degrees).  Let d = deg g.")
print("   d = 0: g = c0 and the integral is c0^m, so c0 = 0.")
print("   d >= 1: write g(s) = c_d s^d (1 + a1/s + O(1/s^2)).  The weight |g|^m e^{-s} has")
print("     log-derivative m g'/g - 1 = 0, i.e. m*d/s ~ 1, so the SADDLE is at s ~ d*m.")
print("     There m*a1/s -> a1/d, a CONSTANT, and arg g(s) -> arg(c_d): THE PHASE")
print("     STABILISES, so no oscillation remains.  Laplace's method gives")
print("        int_0^inf g^m e^{-s} ds  ~  c_d^m * (d m)! * e^{a1/d}  != 0  for large m.")
print("     Hence some moment is nonzero unless g = 0.  QED\n")
print("  numerical check of the ratio  I_m / (c_d^m (dm)!)  ->  e^{a1/d}:")
for g, d, cd, a1 in [(S - 1, 1, 1, -1), (S + 3, 1, 1, 3), (2*S**2 - S, 2, 2, sp.Rational(-1,2)),
                     (sp.I*S - 1, 1, sp.I, sp.I)]:
    print(f"   g = {g}   (d={d}, c_d={cd}, a1={a1}, predicted limit e^(a1/d) = "
          f"{sp.N(sp.exp(sp.Rational(1,1)*a1/d), 8)})")
    for m in (6, 10, 14, 18):
        I = lap(g, m)
        ratio = sp.N(I/(cd**m*factorial(d*m)), 8)
        print(f"      m={m:2d}: I_m/(c_d^m (dm)!) = {ratio}")
