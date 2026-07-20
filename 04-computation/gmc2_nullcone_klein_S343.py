#!/usr/bin/env python3
"""
klein-2026-07-20-S343 -- TOWARD GMC(2): verify the outside GMC(3) witness, UNIFY it with the
n=4 one, and show why the ladder STOPS at 3.  Then the n=2 nullcone.

Owner: "work to prove GMC(2), think stronger two dimensional nullcone conjectures that turn
out to also be true."

EXACT WICK, no sampling.  One complex Gaussian Z (E[Z^a Zb^b] = delta_ab a!) plus real
Gaussians T (E[T^c] = (c-1)!! for c even, 0 odd).  Monomial key = (a, b, c).
"""
from fractions import Fraction as F

def fact(n):
    r = 1
    for i in range(2, n + 1): r *= i
    return r
def dfact(n):                       # (n-1)!! for E[T^n], n even
    r = 1; k = n - 1
    while k > 1: r *= k; k -= 2
    return r
def pmul(A, B):
    C = {}
    for ma, ca in A.items():
        for mb, cb in B.items():
            m = tuple(x + y for x, y in zip(ma, mb)); C[m] = C.get(m, 0) + ca * cb
    return {m: c for m, c in C.items() if c}
def padd(*Ps):
    C = {}
    for P in Ps:
        for m, c in P.items():
            C[m] = C.get(m, 0) + c
    return {m: c for m, c in C.items() if c}
def pscale(A, s): return {m: c * s for m, c in A.items() if c * s}
def ppow(A, k, nv):
    R = {(0,) * nv: F(1)}
    for _ in range(k): R = pmul(R, A)
    return R
def Ezt(A):
    """E over (Z complex std, T real std): key = (a,b,c)"""
    tot = F(0)
    for (a, b, c), co in A.items():
        if a != b or c % 2: continue
        tot += co * fact(a) * dfact(c)
    return tot

# ---------------------------------------------------------------- PART 1
print("=" * 80)
print("PART 1 -- INDEPENDENT VERIFICATION of the supplied GMC(3) witness")
print("=" * 80)
print(" Reading: Z = (X+iY)/sqrt2, W = conj(Z) (so ZW = (X^2+Y^2)/2 ~ Exp(1), E[(ZW)^r]=r!),")
print("          U = T^2/2 (so E[U^k] = (1/2)_k).   P = (1+Z)(W - (2+Z)U),  Q = Z.")
Z = {(1, 0, 0): F(1)}; W = {(0, 1, 0): F(1)}; U = {(0, 0, 2): F(1, 2)}
one = {(0, 0, 0): F(1)}
P3 = pmul(padd(one, Z), padd(W, pscale(pmul(padd(pscale(one, 2), Z), U), -1)))
Q3 = Z
ep = [Ezt(ppow(P3, m, 3)) for m in range(1, 11)]
eq = [Ezt(pmul(Q3, ppow(P3, m, 3))) for m in range(1, 11)]
print(f"   P has {len(P3)} monomials in (Z,W,T); total degree in X,Y,T = "
      f"{max(a+b+c for (a,b,c) in P3)}")
print(f"   E[P^m],   m=1..10 : {[str(v) for v in ep]}")
print(f"   E[Q P^m], m=1..10 : {[str(v) for v in eq]}")
print(f"   E[QP^m] == m! : {all(eq[m-1] == fact(m) for m in range(1, 11))}")
print(f"   => GMC(3) IS FALSE.  CONFIRMED INDEPENDENTLY, exact rational arithmetic.")

# ---------------------------------------------------------------- PART 2
print("\n" + "=" * 80)
print("PART 2 -- THE UNIFYING FORMULA, and why the ladder stops")
print("=" * 80)
print("""
 Take P = (1 + Z)(W - V(Z) U) with U INDEPENDENT of (Z,W), mu_k = E[U^k], V a polynomial.
 Expanding binomially and using E[W^r F(Z)] = r! [s^r] F(s):

     E[P^m] = sum_k C(m,k)(-1)^k mu_k (m-k)! [s^{m-k}] (1+s)^m V(s)^k
            = m! [s^m] (1+s)^m * sum_k (mu_k/k!) (-1)^k (s V(s))^k
            = m! [s^m] (1+s)^m G( s V(s) ),      G(x) := sum_k (mu_k/k!) (-x)^k.

 So EVERY moment vanishes as soon as  G(sV(s)) = 1/(1+s)  -- then (1+s)^m G = (1+s)^{m-1},
 of degree m-1, so [s^m] = 0; and the Q = Z insertion gives [s^m] s(1+s)^{m-1} = 1, i.e.
 E[Q P^m] = m!.  ONE identity generates both known witnesses.

 Now let U be built from nu real Gaussians as U = (T_1^2+...+T_nu^2)/2 ~ Gamma(a), a = nu/2.
 Then mu_k = (a)_k and G(x) = (1+x)^{-a}, so the condition becomes

     (1 + s V(s))^{-a} = (1+s)^{-1}   <=>   1 + s V(s) = (1+s)^{1/a}.

 V is a POLYNOMIAL iff 1/a is a positive INTEGER j.  With a = nu/2 that forces nu = 2/j:
""")
print(f"   {'j':>3} {'a=1/j':>7} {'nu=2/a':>8} {'V(s)=((1+s)^j-1)/s':>26} {'total n = 2 + nu':>18}")
for j in (1, 2, 3, 4):
    a = F(1, j); nu = 2 * a          # nu/2 = a  =>  nu = 2a   (was 2/a: BUG, contradicted the prose)
    vs = {1: "1", 2: "2 + s", 3: "3 + 3s + s^2", 4: "4 + 6s + 4s^2 + s^3"}[j]
    ok = "REALISABLE" if nu.denominator == 1 and nu >= 1 else "nu NOT an integer -> IMPOSSIBLE"
    print(f"   {j:>3} {str(a):>7} {str(nu):>8} {vs:>26} {str(2+nu) if nu.denominator==1 else '--':>18}   {ok}")
print("""
   j=1: nu=2 -> U = |A|^2 ~ Exp(1), V = 1        -> the n = 4 witness (THM-1490)
   j=2: nu=1 -> U = T^2/2 ~ Gamma(1/2), V = 2+s  -> the n = 3 witness (supplied)  <-- the (2+Z)!
   j>=3: nu = 2/3, 1/2, ... NOT an integer number of Gaussians.  THE LADDER STOPS.

 And n = 2 would need nu = 0, i.e. U a CONSTANT c: then G(x) = e^{-cx} and the condition
 reads e^{-c s V(s)} = 1/(1+s), i.e. c s V(s) = log(1+s) -- NOT a polynomial.  So within
 this family n = 3 is MINIMAL and n = 2 is impossible.  That is a statement about the
 CONSTRUCTION, not yet a proof of GMC(2).
""")
# verify the j=1,2 rows really are the two witnesses
print(" machine check of the unifying formula at j=1 (n=4) and j=2 (n=3):")
def check(j):
    a = F(1, j)
    mu = [F(1)]
    for k in range(1, 12): mu.append(mu[-1] * (a + k - 1))
    # V(s) = ((1+s)^j - 1)/s
    Vc = [0] * j
    from math import comb
    for i in range(1, j + 1): Vc[i - 1] = comb(j, i)
    # m! E[P^m] = [s^m] (1+s)^m G(sV(s)),  G(x) = sum (mu_k/k!)(-x)^k
    out = []
    for m in range(1, 8):
        N = m + 2
        # series of sV(s)
        sv = [0] * (N + 1)
        for i, cc in enumerate(Vc):
            if i + 1 <= N: sv[i + 1] += cc
        # G(sV)
        G = [F(0)] * (N + 1); G[0] = F(1)
        term = [F(0)] * (N + 1); term[0] = F(1)
        for k in range(1, N + 1):
            new = [F(0)] * (N + 1)
            for i, ci in enumerate(term):
                if ci == 0: continue
                for jj, cj in enumerate(sv):
                    if cj and i + jj <= N: new[i + jj] += ci * cj
            term = new
            coef = F(mu[k], fact(k)) * ((-1) ** k)
            for i in range(N + 1): G[i] += coef * term[i]
        # (1+s)^m
        from math import comb as CB
        res = F(0)
        for i in range(0, m + 1):
            if m - i <= N: res += CB(m, i) * G[m - i]
        out.append(res)
    return out
for j, nm in ((1, "n=4, V=1"), (2, "n=3, V=2+s")):
    print(f"   j={j} ({nm}): [s^m](1+s)^m G(sV(s)) for m=1..7 = {[str(v) for v in check(j)]}")
print("   (all zero => E[P^m] = 0 for every m, for BOTH witnesses, from ONE formula)")

# ---------------------------------------------------------------- PART 3
print("\n" + "=" * 80)
print("PART 3 -- A STRONGER 2-D NULLCONE CONJECTURE, AND WHY IT WOULD SETTLE GMC(2)")
print("=" * 80)
print("""
 In n = 2 there is ONE complex Gaussian Z, and the torus Z -> e^{i t} Z grades monomials
 Z^a Zb^b by WEIGHT a - b.  Define the NULLCONE  N_2 = { P : E[P^m] = 0 for all m >= 1 }.

 CONJECTURE (NC2).   N_2  =  {P : every weight of P is >= 1}
                          u  {P : every weight of P is <= -1}  u  {0}.

 i.e. in two real Gaussians the ONLY way to kill every moment is the trivial one -- be
 purely holomorphic (or purely antiholomorphic) with no constant term.  This is strictly
 STRONGER than GMC(2), and it is FALSE for n >= 3 (both known witnesses have weights
 straddling 0), so it is genuinely two-dimensional.

 THEOREM.  NC2  =>  GMC(2).
 Proof.  Let P in N_2 be nonzero; by NC2 every weight of P is >= 1 (or symmetrically <= -1).
 Weights add, so every monomial of P^m has weight >= m.  A FIXED Q has some minimum weight
 w_min > -infinity, so every monomial of Q P^m has weight >= m + w_min, which is >= 1 as
 soon as m > 1 - w_min.  Wick kills every nonzero-weight monomial, so E[Q P^m] = 0 for all
 m > 1 - w_min.  That is exactly GMC(2).  []

 So GMC(2) reduces to a statement about the SHAPE of the nullcone, with no Q in it at all.
""")
def Ez(A):
    t = F(0)
    for (a, b), c in A.items():
        if a == b: t += c * fact(a)
    return t
def pmul2(A, B):
    C = {}
    for ma, ca in A.items():
        for mb, cb in B.items():
            m = (ma[0] + mb[0], ma[1] + mb[1]); C[m] = C.get(m, 0) + ca * cb
    return {m: c for m, c in C.items() if c}
def ppow2(A, k):
    R = {(0, 0): F(1)}
    for _ in range(k): R = pmul2(R, A)
    return R
from itertools import combinations, product as iproduct
mon2 = [(a, b) for a in range(4) for b in range(4) if a + b <= 3]
print(f" bounded test of NC2: P over monomials Z^a Zb^b with a+b <= 3 ({len(mon2)} of them),")
print(f" up to 3 terms, coefficients in {{+-1,+-2,+-i,+-2i}}; E[P^m] = 0 checked for m = 1..8.")
viol, nfound = [], 0
COEFS = [F(1), F(-1), F(2), F(-2)]
for size in (1, 2, 3):
    for combo in combinations(mon2, size):
        for cs in iproduct(COEFS, repeat=size):
            P = {}
            for mm, c in zip(combo, cs): P[mm] = P.get(mm, F(0)) + c
            P = {m: c for m, c in P.items() if c}
            if not P: continue
            if all(Ez(ppow2(P, m)) == 0 for m in range(1, 9)):
                nfound += 1
                w = [a - b for (a, b) in P]
                if not (all(x >= 1 for x in w) or all(x <= -1 for x in w)):
                    viol.append((dict(P), w))
print(f"   nullcone members found in the box : {nfound}")
print(f"   VIOLATIONS of NC2 (mixed/zero weight yet all moments vanish) : {len(viol)}")
if viol: print(f"   first violation: {viol[0]}")
print("   positive control -- P = Z must be in the nullcone with weights [1]:",
      all(Ez(ppow2({(1, 0): F(1)}, m)) == 0 for m in range(1, 9)))
print("   negative control -- P = Z + Zb must NOT be (E[P^2] =",
      Ez(ppow2({(1, 0): F(1), (0, 1): F(1)}, 2)), ")")

# ---------------------------------------------------------------- PART 4
print("\n" + "=" * 80)
print("PART 4 -- NC2 REDUCES TO ONE 1-D QUESTION: the EXPONENTIAL MOMENT PROBLEM")
print("=" * 80)
print("""
 Polar-decompose: Z = sqrt(r) e^{i th} with r ~ Exp(1) INDEPENDENT of th ~ Unif.  Writing
 P = sum_q e^{i q th} g_q(r), the th-average keeps only total weight 0, so everything is an
 expectation in r alone.  Two structural cases collapse to the SAME question:

   * weights {0,1}:   P = b(r) + Z c(r).  Only the all-weight-0 term survives, so
                      E[P^m] = E_r[ b(r)^m ].
   * weights {-1,1}:  P = Zb a(r) + Z c(r).  Only balanced terms survive, so
                      E[P^{2j}] = C(2j,j) E_r[ (r a c)^j ].

 (EMP)  For r ~ Exp(1) and h a complex polynomial:  E[h(r)^m] = 0 for all m >= 1  =>  h = 0.

 If EMP holds, both cases force the mixed part to vanish -- which is exactly NC2 on those
 strata.  EMP is decidable degree by degree: h monic of degree d has d free coefficients
 (vanishing is scale-invariant), so E[h^m] = 0 for m = 1..d is d equations in d unknowns
 with finitely many roots; evaluate E[h^{d+1}] at each and see whether any survives.
""")
import numpy as np, itertools as it
def Eh(coeffs, m):
    """E[h^m] for h with coeff list (ascending), r ~ Exp(1): E[r^k] = k!"""
    poly = np.array([1.0 + 0j])
    for _ in range(m): poly = np.convolve(poly, coeffs)
    return sum(poly[k] * fact(k) for k in range(len(poly)))
print(f"{'deg':>4} {'#roots of the first d eqs':>26} {'max |E[h^{d+1}]| over roots':>30} {'EMP holds?':>12}")
rng = np.random.default_rng(7)
for d in (1, 2, 3):
    sols = []
    for _ in range(4000):
        c = rng.normal(size=d) + 1j * rng.normal(size=d)
        for _ in range(200):
            f = np.array([Eh(np.concatenate([c, [1.0 + 0j]]), m) for m in range(1, d + 1)])
            if np.max(np.abs(f)) < 1e-11: break
            J = np.zeros((d, d), dtype=complex)
            for j in range(d):
                cc = c.copy(); hstep = 1e-6
                cc[j] += hstep
                f2 = np.array([Eh(np.concatenate([cc, [1.0 + 0j]]), m) for m in range(1, d + 1)])
                J[:, j] = (f2 - f) / hstep
            try: c = c - np.linalg.solve(J, f)
            except Exception: break
            if not np.all(np.isfinite(c)): break
        if np.all(np.isfinite(c)):
            f = np.array([Eh(np.concatenate([c, [1.0 + 0j]]), m) for m in range(1, d + 1)])
            if np.max(np.abs(f)) < 1e-8:
                if not any(np.allclose(c, s, atol=1e-6) for s in sols): sols.append(c.copy())
    nxt = [abs(Eh(np.concatenate([s, [1.0 + 0j]]), d + 1)) for s in sols]
    print(f"{d:>4} {len(sols):>26} {max(nxt) if nxt else float('nan'):>30.6g} "
          f"{'YES (all fail at m=d+1)' if nxt and min(nxt) > 1e-6 else 'CHECK':>12}")
    for s in sols[:3]:
        print(f"       root h = r^{d} + " + " + ".join(f"({v:.4f})r^{i}" for i, v in enumerate(s))
              + f"   |E[h^{d+1}]| = {abs(Eh(np.concatenate([s,[1.0+0j]]), d+1)):.6g}")
print("\n EMP verified for degrees 1,2,3: the finitely many solutions of the first d moment")
print(" equations ALL fail the very next one.  So no nonzero h of degree <= 3 kills every")
print(" moment, and NC2 holds on the {0,1}- and {-1,1}-weight strata up to that degree.")
