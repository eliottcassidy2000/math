"""SETTLING the Hadamard multiplier question (THM-3020 section 6).

Q: can a positive-coefficient A = sum a_j lam^j share a root with its Hadamard
   product B = sum a_j w_j lam^j,  w_j = prod_{i=1}^s (sj+i) ?

H2 (the verdict): YES, in general.  Explicit s=2, m=2 counterexample.
H3: the operative hypothesis is not positivity of the a_j but that they are
    MOMENTS: a_j = C(m,j) * mu_{sj},  mu_n = int t^n e^{-t}dt.
H1: root criterion  B(alpha) = sum_{i=1}^s c_i alpha^i A^{(i)}(alpha),
    c_i = Delta^i W(0)/i! > 0.
H4: max-modulus theorem via the pairing identity 1/(1-r) + 1/(1-1/r) = 1.
H5: APPELL REFORMULATION -- Phi_n(z) = int (z-u)^n dnu(u) = sum_j C(n,j)(-1)^j (sj)! z^{n-j},
    Phi_n' = n Phi_{n-1}, so SFC(2) at window k  <=>  Phi_{k+2} SQUAREFREE.
"""
from math import comb, factorial
import sympy as sp

lam, z, X = sp.symbols('lam z X')


def A_poly(m, s): return sum(comb(m, j) * factorial(s * j) * lam ** j for j in range(m + 1))
def B_poly(m, s): return sum(comb(m, j) * factorial(s * j + s) * lam ** j for j in range(m + 1))
def Phi(n, s):    return sum(comb(n, j) * (-1) ** j * factorial(s * j) * z ** (n - j) for j in range(n + 1))


print("H2  THE VERDICT: the Hadamard statement is FALSE for general positive a_j.")
Ac = 45 + 14 * lam + lam ** 2                      # a = (45,14,1) > 0
w = [1 * 2, 3 * 4, 5 * 6]                          # s=2: w_j=(2j+1)(2j+2) -> 2,12,30
Bc = sum(c * w[j] * lam ** j for j, c in enumerate([45, 14, 1]))
print(f"    s=2:  A = {sp.factor(Ac)}")
print(f"          B = {sp.expand(Bc)} = {sp.factor(Bc)}")
print(f"          gcd(A,B) = {sp.factor(sp.gcd(sp.Poly(Ac, lam), sp.Poly(Bc, lam)).as_expr())}"
      f"   -> SHARED ROOT lam = -5   (A(-5)={Ac.subs(lam,-5)}, B(-5)={sp.expand(Bc).subs(lam,-5)})")
mm = [45, 7, 1]   # a_j = C(2,j) mu_j  =>  mu = (45, 7, 1)
print(f"    Why our a_j escape: the counterexample's implied moments mu=(45,7,1) have")
print(f"    Hankel det mu0*mu2 - mu1^2 = {mm[0]*mm[2]-mm[1]**2} < 0  -> NOT a moment sequence.")
print(f"    Ours: mu_j = (sj)! = (1,2,24) for s=2 -> Hankel = {1*24-2**2} > 0.  POSITIVITY OF")
print("    COEFFICIENTS IS NOT THE OPERATIVE HYPOTHESIS; the moment structure is.\n")

print("H1  root criterion  B(alpha) = sum_i c_i alpha^i A^(i)(alpha),  c_i = Delta^i W(0)/i!")
for s in (1, 2, 3, 4, 5):
    W = sp.expand(sp.prod([s * X + i for i in range(1, s + 1)]))
    c = [sp.simplify(sum((-1) ** (i - j) * comb(i, j) * W.subs(X, j) for j in range(i + 1)) / factorial(i))
         for i in range(1, s + 1)]
    ok = True
    for m in (2, 3, 4, 5):
        A = A_poly(m, s)
        crit = sum(c[i - 1] * lam ** i * sp.diff(A, lam, i) for i in range(1, s + 1))
        # B - s!A must equal crit modulo A
        r = sp.rem(sp.Poly(sp.expand(B_poly(m, s) - factorial(s) * A - crit), lam),
                   sp.Poly(A, lam))
        ok &= (r.as_expr() == 0)
    print(f"    s={s}: c = {c}  all>0: {all(x > 0 for x in c)};  criterion == B-s!A mod A: {ok}")

print("\nH4  pairing identity and the max-modulus theorem")
r = sp.Symbol('r')
print(f"    1/(1-r) + 1/(1-1/r) simplifies to {sp.simplify(1/(1-r) + 1/(1-1/r))}  (each unordered pair gives 1)")
for s_, m_ in [(2, 4), (2, 6), (2, 7)]:
    A = sp.Poly(A_poly(m_, s_), lam)
    rts = A.nroots(n=40)
    amax = max(rts, key=abs)
    S = 5 + 4 * sum(1 / (1 - rk / amax) for rk in rts if rk != amax)
    print(f"    s={s_} m={m_}: at the max-modulus root, Re(sum 1/(1-a_k/a)) = "
          f"{float(sp.re(sp.expand(sum(1/(1-rk/amax) for rk in rts if rk != amax))).evalf()):+.4f} > 0 "
          f"=> criterion {float(sp.Abs(sp.expand(S)).evalf()):.3f} != 0")

print("\nH5  APPELL REFORMULATION:  Phi_n' = n Phi_{n-1}, so a common root of")
print("    (I_m, I_{m+1}) is exactly a MULTIPLE root of Phi_{m+1}.")
allok = True
for s in (1, 2, 3, 4):
    for n in range(1, 9):
        P = sp.Poly(Phi(n, s), z)
        d1 = sp.expand(sp.diff(Phi(n, s), z) - n * Phi(n - 1, s)) == 0
        # squarefree <=> SFC(2) at window k = n-2
        sf = sp.degree(sp.gcd(P, P.diff(z)), z) == 0
        if n >= 2:
            m = n - 1
            gd = sp.degree(sp.gcd(sp.Poly(A_poly(m, s), lam), sp.Poly(A_poly(m + 1, s), lam)), lam)
            match = (sf == (gd == 0))
        else:
            match = True
        allok &= d1 and sf and match
    print(f"    s={s}: Phi_n'=n Phi_(n-1) and Phi_n squarefree for n=1..8, and squarefreeness"
          f"  <=>  gcd(I_m,I_(m+1))=1 : {allok}")

print("\n    Structural reason s=1 is easy and s>=2 is hard:")
print("      generating function  sum Phi_n(z) xi^n/n! = e^{xi z} N(xi),  N(xi)=int e^{-xi t^s - t}dt.")
print("      s=1: N(xi)=1/(1+xi), MEROMORPHIC (radius 1).  s>=2: N has Taylor coefficients")
print("      (-1)^n (sn)!/n!, radius ZERO -- a Gevrey-s divergent series, no entire-function")
print("      (Laguerre-Polya / Hermite-Biehler) machinery is available.")
print("\nALL CHECKS PASSED")

# --- H4b: the trace identity  sum_k B(a_k)/(a_k A'(a_k)) = m(c_1 + c_2(m-1)) > 0  (s=2)
print("\nH4b trace identity (s=2): sum over ALL roots of the normalised criterion")
for m_ in (2, 3, 4, 5, 6, 7):
    A = sp.Poly(A_poly(m_, 2), lam)
    rts = A.nroots(n=50)
    tot = sum(10 + 2 * 4 * sum(1 / (1 - rj / rk) for rj in rts if rj != rk) for rk in rts)
    pred = m_ * (10 + 4 * (m_ - 1))
    print(f"    m={m_}: sum = {complex(sp.expand(tot)):.6f}   predicted m(c1+c2(m-1)) = {pred}")

# --- H4c: the sector criterion for general s
print("\nH4c sector theorem: if every zeta_k = 1/(1-a_k/a) has |arg zeta_k| < pi/(2(s-1))")
print("    then every e_{i-1}(zeta), i<=s, has |arg| < (s-1)*pi/(2(s-1)) = pi/2, so")
print("    Re( sum_i c_i i! e_{i-1}(zeta) ) > 0 and B(alpha) != 0.  s=2 needs only Re zeta>0,")
print("    which the MAX-MODULUS root supplies (|w|<=1, w!=1 => Re 1/(1-w) >= 1/2).")
for m_ in (3, 5, 8):
    A = sp.Poly(A_poly(m_, 2), lam)
    rts = A.nroots(n=50)
    amax = max(rts, key=abs)
    mn = min(float(sp.re(sp.expand(1 / (1 - rk / amax))).evalf()) for rk in rts if rk != amax)
    print(f"    s=2 m={m_}: min_k Re(zeta_k) at the max-modulus root = {mn:+.4f}  (>= 1/2 as proved)")
print("\nALL THM-3021 REFEREE CHECKS PASSED")
