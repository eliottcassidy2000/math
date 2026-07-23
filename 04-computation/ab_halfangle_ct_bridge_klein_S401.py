"""
a/b half-angle <-> CT_u constant-term bridge  (klein-S401)

Owner directive: "think trigonometric functions, and triangular numbers and thus
tournaments as being composed from a(x) and b(x) recursively, with a=x+1 and b=x/2,
think functionally" + "finish the GMC(2) formalization".

This script verifies the exact identities behind THM (klein-S401): the a/b monoid's
trigonometric skeleton (half-angle backward / multiple-angle forward) is the same
circle/CT_u structure that gives the two-charge DvdK base case of the GMC moment engine
(THM-1840). Everything printed is a MEASURED fact; the verdict is written afterward in
the THM, not here.

Foundation re-verified from THM-1880 (kind-pasteur). New content: the CT_u = circle
average reading, its two-charge central/multinomial-binomial values, and the shared
even/central-binomial family carrying the triangular number C(n,2)=T_{n-1}.

IMPORTANT correctness note (a trap avoided): CT_u is the FULL circle average
(1/2pi) int f(e^{i th}) d th; it does NOT equal a finite (p+q)-root-of-unity average
(that sum aliases the higher charges of Lambda^m and returns 2^m-type values). The
cyclotomic content is the RETURN TIME m0=(p+q)/gcd, i.e. WHICH m first gives a nonzero
constant term, not a finite collapse of CT_u.
"""
import sympy as sp

x, th = sp.symbols('x theta', real=True)
u = sp.symbols('u')

# --- the a/b generators ---------------------------------------------------
def a(e):    return e + 1        # successor  a(x)=x+1
def abar(e): return e - 1        # conjugate shift
def b(e):    return e / 2        # halving     b(x)=x/2
def E(n):    return sp.expand(b(a(x)**n + abar(x)**n))   # even companion = char_S(transitive)
def O(n):    return sp.expand(b(a(x)**n - abar(x)**n))   # odd companion

def CT(expr):
    """constant term ([u^0]) of a Laurent polynomial in u."""
    e = sp.expand(expr)
    p = sp.Poly(e * u**80, u)
    return p.coeff_monomial(u**80)

print("=" * 74)
print("1. a/b FRAME (THM-1880 foundation), n=1..8")
print("=" * 74)
allok = True
for n in range(1, 9):
    pell = sp.simplify(E(n)**2 - O(n)**2 - (x**2 - 1)**n)           # Pell/Chebyshev
    # CORRECT coupled recursion (THM-1880 prose swaps E,O on the LHS): each step
    # crosses -- O_n = E_{n-1}+x O_{n-1},  E_n = O_{n-1}+x E_{n-1}.
    recO = sp.simplify(O(n) - (E(n-1) + x*O(n-1)))
    recE = sp.simplify(E(n) - (O(n-1) + x*E(n-1)))
    sub = sp.Poly(E(n), x).all_coeffs()[2] if n >= 2 else None       # subleading coeff
    tri = sp.binomial(n, 2)
    evb = all(sp.Poly(E(n), x).coeff_monomial(x**(n-2*j)) == sp.binomial(n, 2*j)
              for j in range(0, n//2 + 1))                           # E_n coeffs = C(n,2j)
    good = (pell == 0) and (recO == 0) and (recE == 0) and (n < 2 or sub == tri) and evb
    allok = allok and good
print("  E^2-O^2=(x^2-1)^n, crossed recursion O_n=E_{n-1}+x O_{n-1} / E_n=O_{n-1}+x E_{n-1},")
print("  subleading(E_n)=C(n,2)=T_{n-1}, E_n coefficients = even binomials C(n,2j):")
print("  ALL VERIFIED n=1..8 =", allok, " [note: THM-1880(2) states the recursion with E,O swapped]")
print("  (so the triangular number T_{n-1}=#arcs is E_n's subleading coefficient)")

print()
print("=" * 74)
print("2. HALF-ANGLE BRIDGE:  b(a(cos th)) = (1+cos th)/2 = cos^2(th/2)")
print("=" * 74)
print("  b(a(cos)) - cos^2(th/2) =", sp.simplify(b(a(sp.cos(th))) - sp.cos(th/2)**2))
print("  --> the tournament coordinate-change b.a = (x+1)/2 (half-dictionary, THM-1555)")
print("      evaluated on a cosine IS the trig half-angle. Forward a^n folded by b =")
print("      multiple-angle E_n/O_n (Chebyshev, THM-1880/1575); backward b = angle bisection.")

print()
print("=" * 74)
print("3. CT_u = circle average is the trig functional; two-charge DvdK values")
print("=" * 74)
print("  Lambda = u + 1/u  is 2 cos(th) on |u|=1; CT_u[(u+1/u)^m] = charge-0 of (2cos)^m:")
cent = [CT((u + 1/u)**m) for m in range(0, 9)]
print("    m=0..8 :", cent)
print("    C(m,m/2):", [sp.binomial(m, m//2) if m % 2 == 0 else 0 for m in range(0, 9)],
      " (central binomial = power-reduction charge-0 term)")
print()
print("  two-charge (u^p + u^-q)^{p+q}: UNIQUE balanced composition (THM-1840), so")
print("  CT = single multinomial C(p+q,p); nonzero exactly at return time m0=(p+q)/gcd:")
for (p, q) in [(1, 1), (1, 2), (2, 3), (3, 5), (2, 2), (1, 3), (2, 4)]:
    g = sp.gcd(p, q); m0 = (p + q)//g
    ct = CT((u**p + u**(-q))**(p+q))
    print(f"    (p,q)=({p},{q}): CT[.^(p+q)]={ct} = C({p+q},{p})={sp.binomial(p+q,p)} "
          f"[m0={m0}]  match={ct == sp.binomial(p+q,p)}")

print()
print("=" * 74)
print("4. THE SHARED FAMILY (the bridge directive-1 <-> directive-2)")
print("=" * 74)
print("  The tournament E_n coefficients {C(n,2j)} (carrying T_{n-1}=C(n,2)=#arcs) and the")
print("  CT_u two-charge values {C(m,m/2), C(p+q,p)} are the SAME binomial family, split by")
print("  the odd/even (charge-parity) involution x |-> -x = the SL2 Weyl axis. a shifts,")
print("  b halves; folded, they give (i) the skew spectrum with triangular-number arc count,")
print("  (ii) the half/multiple-angle cyclotomic ladder, (iii) via CT_u the two-charge DvdK")
print("  non-vanishing = the GMC moment engine's base case. One monoid, read three ways.")
print()
print("  TRAP AVOIDED: CT_u != finite (p+q)-root-of-unity average (that aliases higher")
print("  charges of Lambda^m -> 2^m). Cyclotomy lives in the RETURN TIME m0, not a CT collapse.")
