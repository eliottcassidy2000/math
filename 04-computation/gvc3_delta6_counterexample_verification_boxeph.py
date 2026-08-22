#!/usr/bin/env python3
"""
gvc3_delta6_counterexample_verification_boxeph.py

Exact sympy verification of a claimed counterexample to the Generalized
Vanishing Conjecture in 3 variables, GVC(3), for Lambda = Delta^6 with

    Delta  = 4 d/dx d/dy + d^2/dt^2          (operator f -> 4 f_xy + f_tt)
    rho    = t^2 + x*y
    A      = rho + x^2
    C      = y*rho^2 - 2*x*t^2*rho - x^3*t^2
    P      = A*C^2
    Q      = x^2

Claims checked (all finite exact computations, no floats, no numpy):
  1. deg(P) = 12, P has exactly 23 terms fully expanded.
  2. x*C == rho^3 - t^2*(rho + x^2)^2.
  3. Delta^6(P) == 0, Delta^12(P^2) == 0, Delta^18(P^3) == 0.
  4. Delta^(6m+1)(Q*P^m) == 2^(8m+1)*(6m+1)!*(2m)!*(12m+3)!!/(4m+1)!!
     for m = 1, 2 (a nonzero integer constant); m = 0 edge case FAILS.
  5. Delta^(6m)(Q*P^m) != 0 for m = 1, 2.
  6. Structure probes: Delta(P), Delta(C), dual-form isotropy of grad C,
     Delta(rho^k) law, homogeneity, exact vanishing orders.

Provenance of the claim is UNCONFIRMED (cited arXiv id resolves to an
unrelated paper); everything below stands on its own arithmetic.
"""

import time
from math import factorial

import sympy as sp

x, y, t = sp.symbols('x y t')

T0 = time.time()
PASS = []
FAIL = []


def check(label, ok):
    tag = "PASS" if ok else "FAIL"
    (PASS if ok else FAIL).append(label)
    print(f"[{tag}] {label}")
    return ok


def probe(label, ok):
    """Exploratory structural probe: recorded, never counted as FAIL.
    These test MY candidate mechanisms, not the submitted claims."""
    print(f"[PROBE] {label}: {'YES' if ok else 'NO'}")
    return ok


def Delta(f):
    """Delta f = 4 f_xy + f_tt, fully expanded."""
    return sp.expand(4*sp.diff(f, x, y) + sp.diff(f, t, t))


def Delta_pow(f, k, label=None):
    """Apply Delta k times, expanded at each step; report timing."""
    g = sp.expand(f)
    t0 = time.time()
    for i in range(k):
        g = Delta(g)
    dt = time.time() - t0
    if label:
        print(f"    timing: Delta^{k}({label}) computed in {dt:.3f} s")
    return g


def dfact(n):
    """Double factorial n!! = n(n-2)(n-4)...(2 or 1); (-1)!! = 0!! = 1."""
    if n <= 0:
        return 1
    r = 1
    while n > 0:
        r *= n
        n -= 2
    return r


# ---------------------------------------------------------------- data
rho = t**2 + x*y
A = rho + x**2
C = y*rho**2 - 2*x*t**2*rho - x**3*t**2
P = sp.expand(A*C**2)
Q = x**2

print("=" * 72)
print("GVC(3) claimed counterexample: Lambda = Delta^6, P = A*C^2, Q = x^2")
print("Delta = 4 d/dx d/dy + d^2/dt^2 ; exact sympy arithmetic")
print("=" * 72)

# ---------------------------------------------------------------- claim 1
print("\n--- Claim 1: deg(P) = 12 and P has exactly 23 terms ---")
Ppoly = sp.Poly(P, x, y, t)
degP = Ppoly.total_degree()
ntermsP = len(Ppoly.terms())
print(f"    total degree of P  = {degP}")
print(f"    number of terms    = {ntermsP}")
check("deg(P) == 12", degP == 12)
check("P has exactly 23 terms", ntermsP == 23)

# homogeneity check (so Delta drops degree by exactly 2 each application)
homog = all(sum(mon) == 12 for mon, _ in Ppoly.terms())
check("P is homogeneous of degree 12", homog)
print("    P (expanded) =")
print("      ", sp.sstr(P))

# ---------------------------------------------------------------- claim 2
print("\n--- Claim 2: key identity x*C == rho^3 - t^2*(rho + x^2)^2 ---")
lhs2 = sp.expand(x*C)
rhs2 = sp.expand(rho**3 - t**2*(rho + x**2)**2)
check("x*C == rho^3 - t^2*A^2", sp.expand(lhs2 - rhs2) == 0)

# ---------------------------------------------------------------- claim 3
print("\n--- Claim 3: Delta^(6m)(P^m) == 0 for m = 1, 2, 3, 4 ---")
for m in (1, 2, 3, 4):
    f = sp.expand(P**m)
    fp = sp.Poly(f, x, y, t)
    print(f"    P^{m}: degree {fp.total_degree()}, {len(fp.terms())} terms")
    g = Delta_pow(f, 6*m, label=f"P^{m}")
    check(f"Delta^{6*m}(P^{m}) == 0", sp.simplify(g) == 0)

# ---------------------------------------------------------------- claim 4
print("\n--- Claim 4: Delta^(6m+1)(Q*P^m) == "
      "2^(8m+1)*(6m+1)!*(2m)!*(12m+3)!!/(4m+1)!! , m = 1, 2, 4 ---")
for m in (1, 2, 4):
    f = sp.expand(Q*P**m)
    fp = sp.Poly(f, x, y, t)
    d = fp.total_degree()
    print(f"    Q*P^{m}: degree {d} (predicted 12m+2 = {12*m + 2}), "
          f"{len(fp.terms())} terms")
    check(f"deg(Q*P^{m}) == 12m+2 = {12*m+2}", d == 12*m + 2)
    g = Delta_pow(f, 6*m + 1, label=f"Q*P^{m}")
    is_const = (g.free_symbols == set())
    check(f"Delta^{6*m+1}(Q*P^{m}) is a constant", is_const)
    num = 2**(8*m + 1) * factorial(6*m + 1) * factorial(2*m) * dfact(12*m + 3)
    den = dfact(4*m + 1)
    assert num % den == 0, "predicted constant is not an integer!"
    predicted = num // den
    got = sp.Integer(g)
    print(f"    computed  Delta^{6*m+1}(Q*P^{m}) = {got}")
    print(f"    predicted 2^{8*m+1}*{6*m+1}!*{2*m}!*{12*m+3}!!/{4*m+1}!! "
          f"= {predicted}")
    check(f"Delta^{6*m+1}(Q*P^{m}) == predicted constant (m={m})",
          got == predicted)
    check(f"Delta^{6*m+1}(Q*P^{m}) != 0 (m={m})", got != 0)

# m = 0 edge case: Delta(x^2) = 0, but the formula would predict 6
print("\n    m = 0 edge case:")
d_x2 = Delta(Q)
pred0 = 2**1 * factorial(1) * factorial(0) * dfact(3) // dfact(1)
print(f"    Delta(x^2) = {d_x2}; formula at m=0 would give "
      f"2^1*1!*0!*3!!/1!! = {pred0}")
check("m=0 edge case FAILS as claimed (Delta(x^2)=0 != 6)",
      d_x2 == 0 and pred0 == 6)

# ---------------------------------------------------------------- claim 5
print("\n--- Claim 5: Delta^(6m)(Q*P^m) != 0 for m = 1, 2 ---")
for m in (1, 2):
    f = sp.expand(Q*P**m)
    g = Delta_pow(f, 6*m, label=f"Q*P^{m}")
    nz = sp.expand(g) != 0
    check(f"Delta^{6*m}(Q*P^{m}) != 0 (m={m})", nz)
    if m == 1:
        print(f"    Delta^6(Q*P) = {sp.factor(g)}")

# ---------------------------------------------------------------- claim 6
print("\n--- Claim 6: structure probes (mechanism, not just verdict) ---")

DP = Delta(P)
print(f"    Delta(P) == 0 ?  {sp.expand(DP) == 0}")
if sp.expand(DP) != 0:
    DPf = sp.factor(DP)
    print(f"    Delta(P) = {DPf}")

# exact vanishing order of Delta on P and P^2
g5 = Delta_pow(P, 5)
print(f"    Delta^5(P) = {sp.factor(g5)}   (nonzero => order exactly 6)")
check("Delta^5(P) != 0 (vanishing order of P is exactly 6)",
      sp.expand(g5) != 0)
g11 = Delta_pow(sp.expand(P**2), 11)
check("Delta^11(P^2) != 0 (vanishing order of P^2 is exactly 12)",
      sp.expand(g11) != 0)

# Delta as the Laplacian of the quadratic form q = rho = t^2 + x*y:
# Gram matrix of rho in (x,y,t) is M = [[0,1/2,0],[1/2,0,0],[0,0,1]];
# its inverse gives the dual form q*(xi) = 4 xi_x xi_y + xi_t^2, which is
# precisely the symbol of Delta.  So rho is the Fourier dual of the symbol.
M = sp.Matrix([[0, sp.Rational(1, 2), 0],
               [sp.Rational(1, 2), 0, 0],
               [0, 0, 1]])
xx = sp.Matrix([x, y, t])
check("rho == v^T M v for M = Gram(t^2+xy)",
      sp.expand((xx.T*M*xx)[0] - rho) == 0)
Minv = M.inv()
xi_x, xi_y, xi_t = sp.symbols('xi_x xi_y xi_t')
xi = sp.Matrix([xi_x, xi_y, xi_t])
dual = sp.expand((xi.T*Minv*xi)[0])
print(f"    dual form q*(xi) = {dual}  (== symbol of Delta: "
      f"{sp.expand(dual - (4*xi_x*xi_y + xi_t**2)) == 0})")
check("symbol of Delta == dual quadratic form of rho",
      sp.expand(dual - (4*xi_x*xi_y + xi_t**2)) == 0)

# classical law Delta(rho^k) = 2k(2k+1) rho^(k-1) in n = 3 variables
print("    Delta(rho^k) vs 2k(2k+1)*rho^(k-1):")
for k in range(1, 5):
    lhs = Delta(rho**k)
    rhs = sp.expand(2*k*(2*k + 1)*rho**(k - 1))
    check(f"Delta(rho^{k}) == 2*{k}*(2*{k}+1)*rho^{k-1}",
          sp.expand(lhs - rhs) == 0)

# harmonicity / isotropy probes for the building blocks (candidate
# mechanisms of MINE -- negative answers are observations, not failures)
DC = Delta(C)
harmC = probe("C is Delta-harmonic (Delta(C) == 0)", sp.expand(DC) == 0)
if not harmC:
    print(f"    Delta(C) = {sp.factor(DC)}")

DA = Delta(A)
print(f"    Delta(A) = {DA}  (A = rho + x^2 is NOT harmonic)")

# isotropy of grad C w.r.t. the dual form: q*(grad C) = 4 C_x C_y + C_t^2
isoC = sp.expand(4*sp.diff(C, x)*sp.diff(C, y) + sp.diff(C, t)**2)
isoOK = probe("grad C is isotropic for the dual form (q*(grad C) == 0)",
              isoC == 0)
if not isoOK:
    print(f"    q*(grad C) = {sp.factor(isoC)}")

# would-be consequence of Hessian nilpotency (fails, consistent with above)
for k in (2, 3):
    probe(f"Delta(C^{k}) == 0 (would follow from Hessian nilpotency)",
          sp.expand(Delta(sp.expand(C**k))) == 0)

# mixed dual pairing: q*(grad A, grad C) = 2(A_x C_y + A_y C_x) + A_t C_t
mixed = sp.expand(2*(sp.diff(A, x)*sp.diff(C, y)
                     + sp.diff(A, y)*sp.diff(C, x))
                  + sp.diff(A, t)*sp.diff(C, t))
print(f"    dual pairing <grad A, grad C>_{{q*}} = {sp.factor(mixed)}")
isoA = sp.expand(4*sp.diff(A, x)*sp.diff(A, y) + sp.diff(A, t)**2)
print(f"    q*(grad A) = {sp.factor(isoA)}")
probe("q*(grad A) == 4*(A + x^2)",
      sp.expand(isoA - 4*(A + x**2)) == 0)

# Fischer-pairing reading: for homogeneous f of degree 2k,
# Delta^k(f) = <f, rho^k> up to the normalization Delta^k(rho^k).
# So claim 3 says  P^m  is Fischer-orthogonal to rho^(6m), and claim 4
# says <Q*P^m, rho^(6m+1)> != 0.  Normalization law (n = 3):
#   Delta^k(rho^k) = prod_{j=1}^{k} 2j(2j+1) = 2^k * k! * (2k+1)!!
print("    Fischer normalization Delta^k(rho^k) = 2^k*k!*(2k+1)!! :")
for k in (6, 7, 13):
    got = Delta_pow(sp.expand(rho**k), k)
    pred = 2**k * factorial(k) * dfact(2*k + 1)
    check(f"Delta^{k}(rho^{k}) == 2^{k}*{k}!*{2*k+1}!! = {pred}",
          sp.Integer(got) == pred)
# decomposition of the claim-4 constant: predicted/Delta^(6m+1)(rho^(6m+1))
for m in (1, 2):
    num = 2**(8*m + 1) * factorial(6*m + 1) * factorial(2*m) * dfact(12*m + 3)
    den = dfact(4*m + 1)
    norm = 2**(6*m + 1) * factorial(6*m + 1) * dfact(12*m + 3)
    ratio = sp.Rational(num, den) / norm
    print(f"    m={m}: claim-4 constant / Delta^(6m+1)(rho^(6m+1)) "
          f"= {ratio}  (= 2^(2m)*(2m)!/(4m+1)!! = "
          f"{sp.Rational(2**(2*m)*factorial(2*m), dfact(4*m+1))})")

# Q*P = x^2*A*C^2 = A*(x*C)^2 = A*(rho^3 - t^2*A^2)^2
qp_alt = sp.expand(A*(rho**3 - t**2*A**2)**2)
check("Q*P == A*(rho^3 - t^2*A^2)^2 (claim-2 consequence)",
      sp.expand(Q*P - qp_alt) == 0)

# ---------------------------------------------------------------- summary
print("\n" + "=" * 72)
print(f"TOTAL: {len(PASS)} PASS, {len(FAIL)} FAIL "
      f"(wall time {time.time() - T0:.2f} s)")
if FAIL:
    print("FAILED CHECKS:")
    for f in FAIL:
        print("  -", f)
    raise SystemExit(1)
print("ALL CHECKS PASSED: the checked instances (m = 1, 2, 3, 4 vanishing; "
      "m = 1, 2, 4 nonvanishing)")
print("refute GVC(3) for Lambda = Delta^6 *provided* Delta^(6m)(P^m) = 0 "
      "for ALL m (verified only m <= 4 here; all-m claim needs the paper's "
      "induction and stays CITED-UNVERIFIED).")
print("NOTE: nonvanishing for infinitely many m is what contradicts GVC's "
      "'for m >> 0' conclusion; m = 1, 2, 4 already exclude any vanishing "
      "threshold m0 <= 4.")
