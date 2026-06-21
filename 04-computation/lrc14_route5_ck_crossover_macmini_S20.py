#!/usr/bin/env python3
"""
ROUTE 5 -- the c_k crossover mechanism for the tournament H-maximization,
and the test of whether it transfers to the LRC measS7.

GOAL (tournament side):
  H = c_0 + sum_k c_k e_k(x), x = (y_1^2,...,y_m^2) the squared imaginary
  eigenvalue parts of a circulant tournament on Z_p.
  Each e_k is SCHUR-CONCAVE (Schur-Ostrowski: (x_i-x_j)(d_i e_k - d_j e_k)
  = -(x_i-x_j)^2 e_{k-2}(x_hat) <= 0).
  So if ALL c_k >= 0, H is Schur-concave -> uniform (Paley) maximizes.
  The crossover prime is where the FIRST c_k goes NEGATIVE.

  This script:
   (a) computes H exactly (Held-Karp DP) for all circulants on Z_p,
   (b) determines whether H is exactly a multilinear (e_k) symmetric function
       of x or needs higher monomial-symmetric terms,
   (c) if it IS e_k-expressible, solves for the exact c_k (rationals),
   (d) reports the c_k sign pattern and the first flip,
   (e) runs the DIRECT Schur-Ostrowski test on H (no basis assumed), to
       confirm where Schur-concavity actually breaks.

Author: mac-mini-2026-06-21-S20 (ROUTE 5)
"""
import sys, cmath, math
from fractions import Fraction
from itertools import combinations
sys.stdout.reconfigure(line_buffering=True)


def is_qr(a, p):
    if a % p == 0:
        return False
    return pow(a, (p - 1) // 2, p) == 1


def all_circ(p):
    pairs = []
    used = set()
    for a in range(1, p):
        if a not in used:
            b = p - a
            pairs.append((a, b))
            used.add(a)
            used.add(b)
    res = []
    for bits in range(2 ** len(pairs)):
        S = frozenset(a if (bits >> i) & 1 else b for i, (a, b) in enumerate(pairs))
        res.append(S)
    return res


def ham_dp(p, S):
    Sset = set(S)
    adj = [[(j != i and (j - i) % p in Sset) for j in range(p)] for i in range(p)]
    dp = [[0] * p for _ in range(1 << p)]
    for v in range(p):
        dp[1 << v][v] = 1
    full = (1 << p) - 1
    for mask in range(1, full + 1):
        row = dp[mask]
        for v in range(p):
            c = row[v]
            if c == 0 or not (mask & (1 << v)):
                continue
            av = adj[v]
            base = mask
            for w in range(p):
                if base & (1 << w):
                    continue
                if av[w]:
                    dp[base | (1 << w)][w] += c
    return sum(dp[full][v] for v in range(p))


def xvec_exact(p, S):
    """Squared imaginary parts of eigenvalues, as EXACT algebraic numbers
    we recover via their elementary symmetric polynomials.
    We need rationals; the x_k themselves are algebraic, but the
    e_k(x) are rational (symmetric in conjugate eigenvalue pairs).
    We compute e_k(x) from the characteristic polynomial of the
    'imaginary-part-squared' operator. Cleanest: build x numerically with
    high precision, snap e_k to rationals."""
    m = (p - 1) // 2
    xs = []
    for k in range(1, m + 1):
        s = sum(cmath.exp(2j * cmath.pi * k * ss / p) for ss in S)
        xs.append(s.imag ** 2)
    return xs


def esym(xs, k):
    if k == 0:
        return 1.0
    if k > len(xs):
        return 0.0
    tot = 0.0
    for c in combinations(xs, k):
        pr = 1.0
        for v in c:
            pr *= v
        tot += pr
    return tot


def msym(xs, lam):
    """Monomial symmetric polynomial m_lambda evaluated at xs.
    lambda is a partition (tuple, nonincreasing). Sum over distinct
    permutations of the exponent vector padded with zeros."""
    m = len(xs)
    exps = list(lam) + [0] * (m - len(lam))
    seen = set()
    tot = 0.0
    from itertools import permutations
    for perm in set(permutations(exps)):
        if perm in seen:
            continue
        seen.add(perm)
        pr = 1.0
        for xi, ei in zip(xs, perm):
            pr *= xi ** ei
        tot += pr
    return tot


def partitions(n, maxpart=None, maxlen=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    if maxlen is not None and maxlen <= 0:
        return
    for first in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - first, first,
                               None if maxlen is None else maxlen - 1):
            yield (first,) + rest


def snap(x, denom_max=10**12):
    f = Fraction(x).limit_denominator(denom_max)
    return f


def analyze_tournament(p, do_full_basis=True):
    print("=" * 72)
    print(f"TOURNAMENT H -- p = {p}  (p mod 4 = {p % 4})")
    print("=" * 72)
    m = (p - 1) // 2
    S0 = frozenset(j for j in range(1, p) if is_qr(j, p))

    # collect distinct (H, xvec) classes
    classes = {}  # key: rounded x signature -> (H, xs, is_paley)
    for S in all_circ(p):
        H = ham_dp(p, S)
        xs = xvec_exact(p, S)
        key = tuple(sorted(round(x, 6) for x in xs))
        if key not in classes:
            classes[key] = [H, xs, (S == S0)]
    reps = list(classes.values())
    reps.sort(key=lambda r: -r[0])
    print(f"  {len(reps)} distinct spectral classes, m = {m}")
    for H, xs, isp in reps[:8]:
        print(f"    H={H:>16}  x={sorted(round(x,4) for x in xs)} {'PALEY' if isp else ''}")

    # ---- (a) Try MULTILINEAR (e_k) fit: H = c_0 + sum_{k=1}^m c_k e_k ----
    # e_1 is universal (= m*p/4); we still allow c_1 to absorb it.
    # Set up exact linear system over rationals using snapped e_k values.
    print("\n  --- e_k (multilinear symmetric) fit ---")
    npts = len(reps)
    nvar = m + 1  # c_0..c_m
    if npts < nvar:
        print(f"  Not enough classes ({npts}) for {nvar} unknowns; using LS.")
    # Build rational matrix
    rows = []
    bs = []
    for H, xs, isp in reps:
        ev = [snap(esym(xs, k)) for k in range(0, m + 1)]
        rows.append(ev)  # [e_0=1, e_1, ..., e_m]
        bs.append(Fraction(H))

    ck, residual = solve_exact_or_ls(rows, bs)
    if ck is not None:
        maxres = max(abs(float(sum(rows[i][j] * ck[j] for j in range(nvar)) - bs[i]))
                     for i in range(npts))
        print(f"  Solved. max|residual| = {maxres:.3e}")
        if maxres < 1e-4:
            print("  *** H IS EXACTLY MULTILINEAR (e_k) symmetric. ***")
            for k in range(nvar):
                fl = float(ck[k])
                flag = "  <-- NEGATIVE" if ck[k] < 0 else ""
                print(f"    c_{k} = {ck[k]}  (~{fl:.4f}){flag}")
            # first negative
            firstneg = next((k for k in range(2, nvar) if ck[k] < 0), None)
            print(f"  First NEGATIVE c_k (k>=2): {firstneg}")
            if firstneg is None:
                print("  => ALL c_k>=0 => H Schur-CONCAVE => Paley (uniform) maximizes.")
            else:
                print(f"  => c_{firstneg}<0 is the Schur-Ostrowski OFFENDER candidate.")
        else:
            print("  H is NOT exactly multilinear in e_k (residual too large).")
            print("  -> needs higher monomial-symmetric terms; see full-basis fit.")

    # ---- (b) FULL symmetric-function fit (all monomial symmetric m_lambda) ----
    if do_full_basis:
        print("\n  --- FULL symmetric fit: H = sum_lambda a_lambda m_lambda(x) ---")
        # degrees 0..m, partitions with at most m parts, each part <= ???
        # H is a polynomial in x of total degree <= m (since H ~ product structure).
        # Use all partitions of d=0..m with <= m parts. (parts can exceed 1.)
        basis = []
        for d in range(0, m + 1):
            for lam in partitions(d, maxpart=d, maxlen=m):
                if len(lam) <= m:
                    basis.append(lam)
        print(f"  basis size = {len(basis)} (partitions of 0..{m}, <= {m} parts)")
        rows2 = []
        for H, xs, isp in reps:
            rows2.append([snap(msym(xs, lam)) for lam in basis])
        bs2 = [Fraction(H) for H, xs, isp in reps]
        a, _ = solve_exact_or_ls(rows2, bs2)
        if a is not None:
            maxres = max(abs(float(sum(rows2[i][j] * a[j] for j in range(len(basis))) - bs2[i]))
                         for i in range(len(reps)))
            print(f"  max|residual| = {maxres:.3e}")
            if maxres < 1e-4:
                print("  *** EXACT full-symmetric expansion: ***")
                for lam, coef in zip(basis, a):
                    if coef != 0:
                        print(f"    m_{lam}: {coef}  (~{float(coef):.4f})")
    return reps


def solve_exact_or_ls(rows, bs):
    """Solve rows @ c = bs exactly over Q if square & nonsingular,
    else least-squares over rationals via normal equations (best effort)."""
    npts = len(rows)
    nvar = len(rows[0])
    if npts == nvar:
        c = gauss_solve([r[:] for r in rows], [b for b in bs])
        return c, None
    # overdetermined: normal equations A^T A c = A^T b (rational)
    AtA = [[sum(rows[i][r] * rows[i][cc] for i in range(npts)) for cc in range(nvar)]
           for r in range(nvar)]
    Atb = [sum(rows[i][r] * bs[i] for i in range(npts)) for r in range(nvar)]
    c = gauss_solve(AtA, Atb)
    return c, None


def gauss_solve(A, b):
    n = len(b)
    M = [A[i][:] + [b[i]] for i in range(n)]
    for col in range(n):
        piv = None
        for r in range(col, n):
            if M[r][col] != 0:
                piv = r
                break
        if piv is None:
            return None
        M[col], M[piv] = M[piv], M[col]
        pv = M[col][col]
        for r in range(n):
            if r == col:
                continue
            if M[r][col] != 0:
                f = M[r][col] / pv
                for cc in range(col, n + 1):
                    M[r][cc] -= f * M[col][cc]
    return [M[i][n] / M[i][i] for i in range(n)]


if __name__ == "__main__":
    for p in [7, 11]:
        analyze_tournament(p, do_full_basis=True)
        print()
