#!/usr/bin/env python3
"""
THM-438 ADDENDUM-12 verification (monad-explorer 14th session).

CONJECTURE: the factorial law (free cumulants kappa_n = n!) is the FREE COMPOUND
POISSON with exponential Levy measure nu = e^{-x} dx on (0, infty). Dually, the
classical over-count A000262 is the CLASSICAL COMPOUND POISSON of the SAME Levy
measure. Both endpoints of the THM-438 diagonal/over-count are compound-Poisson
laws of one and the same exponential jump measure -- the cleanest form of
"free = classical on the noncrossing sublattice."

We verify, with exact rationals where possible:

(A) Classical side:
    - classical cumulants kappa_n = n! (n>=1) -> classical moments via the
      classical (ALL-partitions) moment-cumulant formula == A000262.
    - classical cumulant g.f. K_cl(t) = sum_{n>=1} kappa_n t^n / n! = sum t^n
      = t/(1-t), and this equals the compound-Poisson cumulant g.f.
      int (e^{tx}-1) e^{-x} dx = t/(1-t).  (symbolic)
    - kappa_n = n! = int_0^infty x^n e^{-x} dx  = moments of the exponential
      Levy measure. (exact)

(B) Free side:
    - free cumulants kappa_n = n! -> free moments via the NONCROSSING moment-
      cumulant formula == A088368 (the diagonal).
    - R-transform R(z) = sum_{n>=1} n! z^{n-1} is the divergent factorial series
      (ADD-6 resurgence). Its BOREL transform is B(zeta) = sum_{n>=1} zeta^n
      = zeta/(1-zeta), singularity at zeta=1, and the Borel sum is the CONVERGENT
      integral  R(z) = int_0^infty t e^{-t} / (1 - z t) dt.  (numeric asymptotic)

(C) The free law itself:
    - Hankel positive-definiteness of A088368 (genuine probability measure on R).
    - Jacobi continued-fraction parameters (b_n, lambda_n=a_n^2) of A088368 --
      check whether polynomial in n (=> classical orthogonal-polynomial family).
    - Numerical Cauchy transform / spectral density of the free compound Poisson
      (R-transform via the convergent integral), support and edge behaviour.

Saves all output. No number theory; finite/analytic checks only.
"""
import sys
from fractions import Fraction as Fr
from math import factorial

# ---------------------------------------------------------------------------
# Reference OEIS values
# A088368 (free moments, the diagonal): m_k, k=0..
A088368 = [1, 1, 3, 13, 69, 421, 2867, 21477, 173601, 1498421]
# A000262 (classical moments, sets of lists): a_k, k=0..
A000262 = [1, 1, 3, 13, 73, 501, 4051, 37633, 394353, 4596553]

KMAX = 9
kappa = [None] + [factorial(n) for n in range(1, KMAX + 1)]   # kappa[n]=n!, kappa[0] unused


# ---------------------------------------------------------------------------
# Partition enumerators
def set_partitions(elements):
    """All set partitions of a list, as list of blocks (each block a tuple)."""
    if not elements:
        yield []
        return
    first = elements[0]
    for rest in set_partitions(elements[1:]):
        # add 'first' to each existing block
        for i in range(len(rest)):
            yield rest[:i] + [(first,) + rest[i]] + rest[i + 1:]
        # or as its own block
        yield [(first,)] + rest


def is_noncrossing(blocks):
    """A partition of {0..k-1} is noncrossing iff no a<b<c<d with a,c in one
    block and b,d in another."""
    # standard test: for every pair of blocks, they don't interleave
    bl = [sorted(b) for b in blocks]
    for i in range(len(bl)):
        for j in range(len(bl)):
            if i == j:
                continue
            B, C = bl[i], bl[j]
            # crossing if exist a<b<c<d, a,c in B, b,d in C
            for a_idx in range(len(B)):
                for c_idx in range(a_idx + 1, len(B)):
                    a, c = B[a_idx], B[c_idx]
                    for x in C:
                        for y in C:
                            if a < x < c < y:
                                return False
    return True


def moment_from_cumulants(k, kappa, noncrossing_only):
    """m_k = sum over (NC if flag) partitions pi of [k] of prod_B kappa_{|B|}."""
    if k == 0:
        return 1
    total = 0
    for blocks in set_partitions(list(range(k))):
        if noncrossing_only and not is_noncrossing(blocks):
            continue
        prod = 1
        for B in blocks:
            prod *= kappa[len(B)]
        total += prod
    return total


# ---------------------------------------------------------------------------
out = []
def emit(s=""):
    out.append(s)
    print(s)


emit("=" * 78)
emit("THM-438 ADDENDUM-12 -- factorial law = compound Poisson of exponential")
emit("=" * 78)

# (A) classical side ---------------------------------------------------------
emit("\n(A) CLASSICAL: kappa_n=n! -> classical (all-partition) moments == A000262")
cl = [moment_from_cumulants(k, kappa, noncrossing_only=False) for k in range(KMAX + 1)]
emit(f"   computed : {cl}")
emit(f"   A000262  : {A000262[:KMAX+1]}")
emit(f"   MATCH    : {cl == A000262[:KMAX+1]}")

# kappa_n = n! = int_0^infty x^n e^{-x} dx (moments of exponential Levy measure)
emit("\n   kappa_n = n! = int_0^infty x^n e^{-x} dx  (gamma integral):")
for n in range(1, 6):
    emit(f"      n={n}: n! = {factorial(n)} = Gamma({n}+1)  [exp Levy measure moment]")

# classical cumulant g.f. and compound-Poisson identity (symbolic)
try:
    import sympy as sp
    t, x = sp.symbols('t x', positive=True)
    # classical cumulant g.f. = sum_{n>=1} kappa_n t^n / n! = sum t^n = t/(1-t)
    Kcl = sp.summation(t**sp.Symbol('n', integer=True), ('n', 1, sp.oo))  # may not converge symbolically
    # do it as a series instead
    Kcl_series = sum(t**n for n in range(1, 12))
    emit("\n   classical cumulant g.f. K_cl(t) = sum kappa_n t^n/n! = sum_{n>=1} t^n")
    emit(f"      = {sp.series(t/(1-t), t, 0, 8)}  (closed form t/(1-t))")
    # compound Poisson cumulant g.f. with Levy measure e^{-x}dx:
    cp = sp.integrate((sp.exp(t*x) - 1) * sp.exp(-x), (x, 0, sp.oo))
    cp = sp.simplify(cp)
    emit(f"   compound-Poisson cgf  int (e^{{tx}}-1) e^{{-x}} dx = {cp}")
    emit(f"   equals t/(1-t)?  {sp.simplify(cp - t/(1-t)) == 0}")
    HAVE_SYMPY = True
except Exception as e:
    emit(f"   [sympy unavailable or failed: {e}]")
    HAVE_SYMPY = False

# (B) free side --------------------------------------------------------------
emit("\n(B) FREE: kappa_n=n! -> free (noncrossing) moments == A088368 (the diagonal)")
fr = [moment_from_cumulants(k, kappa, noncrossing_only=True) for k in range(KMAX + 1)]
emit(f"   computed : {fr}")
emit(f"   A088368  : {A088368[:KMAX+1]}")
emit(f"   MATCH    : {fr == A088368[:KMAX+1]}")

emit("\n   excess (classical - free) = crossing-partition gap:")
emit(f"      {[cl[k]-fr[k] for k in range(KMAX+1)]}   (= 0 until k=4, then 4,80,...)")

# R-transform Borel resummation: R(z)=sum n! z^{n-1}; check vs convergent integral
emit("\n   R-transform R(z)=sum_{n>=1} n! z^{n-1}  (divergent / Gevrey-1, ADD-6)")
emit("   Borel transform B(zeta)=sum_{n>=1} zeta^n = zeta/(1-zeta), sing at zeta=1")
emit("   Borel sum: R(z) = int_0^infty t e^{-t}/(1 - z t) dt  (convergent for z<0)")
try:
    import mpmath as mp
    mp.mp.dps = 30
    emit("   numeric check (z<0, optimal-truncation of divergent series vs integral):")
    for z in [Fr(-1, 10), Fr(-1, 20), Fr(-1, 50)]:
        zf = mp.mpf(z.numerator) / z.denominator
        integ = mp.quad(lambda tt: tt * mp.e**(-tt) / (1 - zf * tt), [0, mp.inf])
        # optimal truncation of sum n! z^{n-1}: truncate near n ~ 1/|z|
        N = int(1 / abs(float(z)))
        partial = mp.fsum(mp.factorial(n) * zf**(n - 1) for n in range(1, max(2, N)))
        emit(f"      z={float(z):+.4f}: integral={mp.nstr(integ,12)}  "
             f"opt-trunc(N={N})={mp.nstr(partial,12)}  "
             f"|diff|={mp.nstr(abs(integ-partial),4)}")
    HAVE_MPMATH = True
except Exception as e:
    emit(f"   [mpmath unavailable: {e}]")
    HAVE_MPMATH = False

# (C) the free law: Hankel + Jacobi parameters -------------------------------
emit("\n(C) THE FREE LAW (free compound Poisson of exponential)")
emit("   Hankel determinants of A088368 (positive => genuine measure on R):")
m = [Fr(v) for v in A088368]   # moments m_0..

def hankel_det(seq, n):
    """det of (n+1)x(n+1) Hankel matrix [seq[i+j]]_{0<=i,j<=n} using fractions."""
    M = [[seq[i + j] for j in range(n + 1)] for i in range(n + 1)]
    # fraction-free / plain Gaussian elimination on Fractions
    det = Fr(1)
    A = [row[:] for row in M]
    sz = n + 1
    for c in range(sz):
        piv = None
        for r in range(c, sz):
            if A[r][c] != 0:
                piv = r
                break
        if piv is None:
            return Fr(0)
        if piv != c:
            A[c], A[piv] = A[piv], A[c]
            det = -det
        det *= A[c][c]
        inv = A[c][c]
        for r in range(c + 1, sz):
            f = A[r][c] / inv
            if f != 0:
                for cc in range(c, sz):
                    A[r][cc] -= f * A[c][cc]
    return det

nmax = (len(A088368) - 1) // 2
dets = [hankel_det(m, n) for n in range(nmax + 1)]
emit(f"      Hankel dets D_n (n=0..{nmax}): {[str(d) for d in dets]}")
emit(f"      all positive: {all(d > 0 for d in dets)}  => valid moment sequence")

# Jacobi continued-fraction parameters via the standard recurrence
# (Stieltjes / orthogonal polynomial three-term recurrence p_{n+1}=(x-b_n)p_n - lambda_n p_{n-1})
# Use the determinant ratios: lambda_n = D_{n-2} D_n / D_{n-1}^2  (with D_{-1}=1),
# b_n needs the shifted Hankel. We compute via the robust "Chebyshev/quotient" using
# the moments directly (build the monic OPs by Gram-Schmidt over Q).
emit("\n   Jacobi (continued-fraction) parameters of A088368:")
emit("   G(z)=1/(z-b_0-lambda_1/(z-b_1-lambda_2/(z-...))):")

def jacobi_params(moments):
    """Return (b[], lam[]) from a moment list using exact-rational Gram-Schmidt
    building monic orthogonal polynomials. moments[i]=m_i, i=0..2N."""
    # inner product <p,q> = sum_i sum_j p_i q_j m_{i+j}
    def ip(p, q):
        s = Fr(0)
        for i, pi in enumerate(p):
            if pi == 0:
                continue
            for j, qj in enumerate(q):
                if qj == 0:
                    continue
                s += pi * qj * moments[i + j]
        return s

    def xmul(p):
        return [Fr(0)] + p  # multiply polynomial by x

    N = (len(moments) - 1) // 2
    P = [[Fr(1)]]               # P[0]=1
    b, lam = [], []
    norms = [ip(P[0], P[0])]    # <P0,P0>=m_0
    for n in range(N):
        xp = xmul(P[n])
        bn = ip(xp, P[n]) / norms[n]
        b.append(bn)
        if n == 0:
            # P1 = x - b0
            newP = [(-bn), Fr(1)]
        else:
            lamn = norms[n] / norms[n - 1]
            lam.append(lamn)
            # P_{n+1} = (x-b_n) P_n - lam_n P_{n-1}
            t1 = [(-bn) * c for c in P[n]]
            t1 = [a + b_ for a, b_ in zip(t1 + [Fr(0)], xp)] if len(xp) >= len(t1) else None
            # do it cleanly:
            xp2 = xmul(P[n])
            shift = [(-bn) * c for c in P[n]]
            # align lengths
            L = max(len(xp2), len(shift), len(P[n - 1]))
            xp2 += [Fr(0)] * (L - len(xp2))
            shift += [Fr(0)] * (L - len(shift))
            pm1 = P[n - 1] + [Fr(0)] * (L - len(P[n - 1]))
            newP = [xp2[i] + shift[i] - lamn * pm1[i] for i in range(L)]
        P.append(newP)
        norms.append(ip(newP, newP))
    return b, lam

b, lam = jacobi_params(m)
emit(f"      b_n      (n=0..): {[str(x) for x in b]}")
emit(f"      lambda_n (n=1..): {[str(x) for x in lam]}")
# check integrality / simple pattern
emit(f"      b_n integers? {all(x.denominator==1 for x in b)}")
emit(f"      lambda_n integers? {all(x.denominator==1 for x in lam)}")
if all(x.denominator == 1 for x in b):
    bi = [x.numerator for x in b]
    emit(f"      b_n diffs: {[bi[i+1]-bi[i] for i in range(len(bi)-1)]}")
if all(x.denominator == 1 for x in lam):
    li = [x.numerator for x in lam]
    emit(f"      lambda_n: {li}")
    emit(f"      lambda_n/n: {[Fr(li[i],i+1) for i in range(len(li))]}")
    emit(f"      lambda_n ratios: {[Fr(li[i+1],li[i]) for i in range(len(li)-1)]}")

# (C') numeric spectral density via the R-transform integral --------------------
if 'mpmath' in sys.modules:
    import mpmath as mp
    mp.mp.dps = 25
    emit("\n   Numeric spectral density of the free law (R-transform = Borel integral):")
    emit("   K(w)=1/w + int_0^infty t e^{-t}/(1-w t) dt ; G inverse of K ; rho=-Im G/pi")

    def Rtransform(w):
        # principal-value-aware integral; w complex
        return mp.quad(lambda tt: tt * mp.e**(-tt) / (1 - w * tt), [0, mp.inf])

    def Kfun(w):
        return 1 / w + Rtransform(w)

    # For x on real axis, solve K(G)=x+i0 for G with small Im part, get density.
    def density(xval, eps=mp.mpf('1e-6')):
        z = mp.mpc(xval, eps)
        # solve K(w)=z near a complex guess
        try:
            w = mp.findroot(lambda ww: Kfun(ww) - z, mp.mpc(0.1, -0.3))
            return max(mp.mpf(0), -mp.im(w) / mp.pi)
        except Exception:
            return None

    emit("      x      rho(x)")
    for xv in [0.5, 1, 2, 3, 5, 8, 12, 20, 30, 50]:
        d = density(xv)
        emit(f"      {xv:>5}   {mp.nstr(d,8) if d is not None else 'n/a'}")

# write output
with open("05-knowledge/results/paley_starstar_compound_poisson_monad.out", "w") as f:
    f.write("\n".join(out) + "\n")
emit("\n[done] output saved to 05-knowledge/results/paley_starstar_compound_poisson_monad.out")
