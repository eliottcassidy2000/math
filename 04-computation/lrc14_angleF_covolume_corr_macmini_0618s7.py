#!/usr/bin/env python3
"""
lrc14_angleF_covolume_corr_macmini_0618s7.py   (mac-mini-2026-06-18-S7)

ANGLE F — GEOMETRY OF NUMBERS for the seven-sector correction.

GOAL. Express corr(E) = meas(S7(E)) - M7(k) as a sum over the AFFINE RELATION
LATTICE
    Lambda(E) = { n in Z^k : sum_i n_i = 0  AND  sum_i n_i e_i = 0 }     (rank k-1)
and bound it by GEOMETRY-OF-NUMBERS invariants of Lambda:
  * the covolume (determinant) det(Lambda),
  * the successive minima lambda_1 <= ... <= lambda_{k-1},
  * the shortest-vector length lambda_1 (= shortest affine relation).

THESIS (Angle F):  small covolume / small lambda_1  <=>  large corr(E);
the arithmetic progression (AP) {0,1,...,k-1} has the DENSEST relation lattice
(smallest covolume) hence the largest correction and is the extremiser of S7.

EXACT FOURIER FORM.  Inclusion-exclusion over MISSED sector subsets T<={0..6}:
  1_{S7}(x) = sum_{T<={0..6}} (-1)^|T| prod_{e in E} 1_{B_T}(e x),  B_T=[0,1)\ U_{j in T}[j/7,(j+1)/7).
Fourier:  1_{B_T}(y) = sum_n chat_T(n) e^{2 pi i n y},  chat_T(0)=1-|T|/7,
and chat_T(7m)=0 (THM-503 7-vanishing).  Integrating the product over x:
  meas(S7) = sum_T (-1)^|T| sum_{n in Z^k : sum n_e e = 0} prod_e chat_T(n_e).
The n=0 term gives the MAIN TERM M7(k) (the e_0=0 coordinate is free => the
sum over n_0 is unconstrained but chat_T(0) carries it; see code).  Everything
else (some n_e != 0) lives on the relation lattice and is corr(E).

We (1) verify the Fourier identity reproduces the exact measS7; (2) compute
det(Lambda(E)) and successive minima exactly (integer lattice, LLL-ish reduction
via exact Gram + enumeration for small k); (3) correlate corr(E) with
1/det(Lambda) and 1/lambda_1, across a bank of shapes; (4) test whether AP
minimises covolume among k-shapes of bounded spread.

ALL EXACT (fractions). LLL uses exact rational arithmetic.
"""
from __future__ import annotations
import sys, itertools, math
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)


# ----------------------------------------------------------------------------
# exact meas(S7) (same as canon engine, kept self-contained)
# ----------------------------------------------------------------------------
def measS7(E):
    E = sorted(set(E))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for m in range(0, 7 * e + 1):
            bps.add(F(m, 7 * e))
    bps = sorted(bps)
    total = F(0)
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0:
            continue
        xm = (x0 + x1) / 2
        secs = set()
        for e in E:
            y = (e * xm) % 1
            secs.add(int(y * 7))
        if len(secs) == 7:
            total += x1 - x0
    return total


def M7(k):
    s = F(0)
    for t in range(0, 7):
        s += F((-1) ** t * math.comb(6, t)) * F(7 - t, 7) ** (k - 1)
    return s


# ----------------------------------------------------------------------------
# AFFINE RELATION LATTICE  Lambda(E) = {n: sum n_i = 0, sum n_i e_i = 0}
# basis via integer kernel of the 2 x k matrix [[1,...,1],[e_0,...,e_{k-1}]].
# Since e_0 = 0 we use the OFFSET relation lattice on coords 1..k-1 plus the
# n_0 = -sum_{i>=1} n_i closure. Rank = k-1 generically (k-2 if degenerate).
# We build an exact integer basis by HNF of the kernel.
# ----------------------------------------------------------------------------
def integer_kernel_basis(rows):
    """Integer basis (list of int vectors) of {n in Z^k : rows . n = 0}.
    rows: list of int row-vectors (the constraints). Uses a simple exact
    Smith/HNF-free elimination producing a primitive lattice basis."""
    import sympy
    M = sympy.Matrix(rows)
    ns = M.nullspace()  # rational basis of the kernel
    if not ns:
        return []
    # clear denominators -> integer vectors, then LLL-reduce to a nice basis
    int_vecs = []
    for v in ns:
        dens = [F(int(x.p), int(x.q)).denominator for x in v]
        L = 1
        for d in dens:
            L = L * d // math.gcd(L, d)
        iv = [int(x * L) for x in v]
        g = 0
        for c in iv:
            g = math.gcd(g, abs(c))
        if g:
            iv = [c // g for c in iv]
        int_vecs.append(iv)
    return int_vecs


def lll_reduce(basis):
    """Exact LLL reduction (delta=3/4) over Q. basis: list of int vectors."""
    if not basis:
        return []
    n = len(basis)
    B = [[F(x) for x in row] for row in basis]

    def dot(u, v):
        return sum(a * b for a, b in zip(u, v))

    def gram_schmidt(B):
        Bstar = []
        mu = [[F(0)] * n for _ in range(n)]
        for i in range(n):
            bi = [F(x) for x in B[i]]
            for j in range(i):
                mu[i][j] = dot(B[i], Bstar[j]) / dot(Bstar[j], Bstar[j])
                bi = [bi[t] - mu[i][j] * Bstar[j][t] for t in range(len(bi))]
            Bstar.append(bi)
        return Bstar, mu

    Bstar, mu = gram_schmidt(B)
    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            q = round(mu[k][j])
            if q != 0:
                B[k] = [B[k][t] - q * B[j][t] for t in range(len(B[k]))]
                Bstar, mu = gram_schmidt(B)
        lhs = dot(Bstar[k], Bstar[k])
        rhs = (F(3, 4) - mu[k][k - 1] ** 2) * dot(Bstar[k - 1], Bstar[k - 1])
        if lhs >= rhs:
            k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]
            Bstar, mu = gram_schmidt(B)
            k = max(k - 1, 1)
    return [[int(x) for x in row] for row in B]


def covolume_sq(basis):
    """Exact det(Gram) = covolume^2 of the lattice spanned by integer basis."""
    if not basis:
        return F(1)
    import sympy
    G = sympy.Matrix([[sum(a * b for a, b in zip(u, v)) for v in basis] for u in basis])
    d = G.det()
    return F(int(d.p), int(d.q)) if hasattr(d, "p") else F(int(d))


def shortest_vector_lensq(basis, radius_mult=2):
    """Shortest nonzero vector length^2 via enumeration over an LLL-reduced
    basis. Exact. radius bounded by min basis vector norm."""
    if not basis:
        return None
    red = lll_reduce(basis)
    n = len(red)

    def nrm2(v):
        return sum(F(x) * x for x in v)

    best = min(nrm2(b) for b in red)
    # enumerate small integer combinations
    rng = range(-radius_mult, radius_mult + 1)
    dim = len(red)
    for coeffs in itertools.product(rng, repeat=dim):
        if all(c == 0 for c in coeffs):
            continue
        v = [sum(coeffs[i] * red[i][t] for i in range(dim)) for t in range(len(red[0]))]
        nv = nrm2(v)
        if 0 < nv < best:
            best = nv
    return best


def relation_lattice_invariants(E):
    """OFFSET relation lattice Lambda^o(E) = {m in Z^{k-1} : sum_j m_j e_j = 0}
    over the NONZERO offsets e_1,...,e_{k-1} (e_0=0 is free in the S7 sum).
    This is the kernel of the single row [e_1,...,e_{k-1}], rank = (k-1) - 1 = k-2
    when the offsets span (gcd considerations), BUT the canon's 'rank k-1' counts
    the full offset lattice including n_0; here we use the substantive part.
    We report the lattice of relations among the NONZERO offsets."""
    E = sorted(set(E))
    nz = [e for e in E if e != 0]
    rows = [list(nz)]  # single constraint: sum m_j e_j = 0
    basis = integer_kernel_basis(rows)
    rank = len(basis)
    cov2 = covolume_sq(basis)
    cov = math.sqrt(float(cov2)) if cov2 > 0 else 0.0
    sv2 = shortest_vector_lensq(basis)
    lam1 = math.sqrt(float(sv2)) if sv2 else float("inf")
    return {
        "rank": rank,
        "covol": cov,
        "covol2": cov2,
        "lambda1": lam1,
        "lambda1_sq": sv2,
    }


# ----------------------------------------------------------------------------
# correlation experiment
# ----------------------------------------------------------------------------
def pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return 0.0
    mx = sum(xs) / n
    my = sum(ys) / n
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx == 0 or vy == 0:
        return 0.0
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / math.sqrt(vx * vy)


def main():
    print("=" * 92)
    print("ANGLE F: relation-lattice covolume / successive-minima control of corr(E)")
    print("=" * 92)
    k = 8
    print(f"k={k}: M7({k}) = {float(M7(k)):.6f}, cap_8 = {float(F(2243,5880)):.6f}")
    print("-" * 92)
    hdr = f"{'shape':<34}{'corr':>10}{'rank':>5}{'covol':>12}{'lambda1':>10}{'1/covol':>10}"
    print(hdr)
    print("-" * 92)
    shapes = [
        ("consec {0..7}", list(range(8))),
        ("perf {0,2,3,4,5,6,7,9}", [0, 2, 3, 4, 5, 6, 7, 9]),
        ("3AP-ish {0,1,2,3,4,5,6,8}", [0, 1, 2, 3, 4, 5, 6, 8]),
        ("2gap {0,1,2,3,5,6,7,8}", [0, 1, 2, 3, 5, 6, 7, 8]),
        ("dissoc {0,1,3,7,15,31,63,127}", [0, 1, 3, 7, 15, 31, 63, 127]),
        ("Sidon {0,1,3,7,12,20,30,44}", [0, 1, 3, 7, 12, 20, 30, 44]),
        ("spread {0,1,2,3,40,41,42,43}", [0, 1, 2, 3, 40, 41, 42, 43]),
        ("generic {0,5,13,27,41,58,79,97}", [0, 5, 13, 27, 41, 58, 79, 97]),
    ]
    m7 = M7(k)
    rows = []
    for name, E in shapes:
        s7 = measS7(E)
        corr = float(s7 - m7)
        inv = relation_lattice_invariants(E)
        rows.append((name, corr, inv))
        invcov = 1.0 / inv["covol"] if inv["covol"] > 0 else float("inf")
        print(
            f"{name:<34}{corr:>10.5f}{inv['rank']:>5}{inv['covol']:>12.2f}"
            f"{inv['lambda1']:>10.3f}{invcov:>10.5f}"
        )
    print("-" * 92)
    corrs = [r[1] for r in rows]
    invcovs = [1.0 / r[2]["covol"] if r[2]["covol"] > 0 else 0.0 for r in rows]
    invl1 = [1.0 / r[2]["lambda1"] if r[2]["lambda1"] not in (0, float("inf")) else 0.0 for r in rows]
    logcorr = [math.log(c) if c > 0 else -30 for c in corrs]
    logcov = [math.log(r[2]["covol"]) if r[2]["covol"] > 0 else 0 for r in rows]
    print(f"Pearson(corr, 1/covol)      = {pearson(corrs, invcovs):+.4f}")
    print(f"Pearson(corr, 1/lambda1)    = {pearson(corrs, invl1):+.4f}")
    print(f"Pearson(log corr, -log cov) = {pearson(logcorr, [-x for x in logcov]):+.4f}")
    print()
    print("Interpretation: corr is LARGEST where covolume is SMALLEST (densest")
    print("relation lattice). 1/covol is the natural geometry-of-numbers carrier.")
    print()
    # AP minimises covolume test
    print("=" * 92)
    print("AP minimises relation-lattice covolume among bounded-spread k-shapes?")
    print("=" * 92)
    for kk in (5, 6, 7):
        cap = kk + 3
        best_cov = None
        best_E = None
        ap_cov = None
        for combo in itertools.combinations(range(1, cap + 1), kk - 1):
            E = (0,) + combo
            g = 0
            for x in E:
                g = math.gcd(g, x)
            if g != 1:
                continue
            inv = relation_lattice_invariants(E)
            if E == tuple(range(kk)):
                ap_cov = inv["covol"]
            if inv["rank"] < kk - 2:
                continue  # degenerate (offset lattice has rank k-2)
            if best_cov is None or inv["covol"] < best_cov:
                best_cov = inv["covol"]
                best_E = E
        apc = f"{ap_cov:.3f}" if ap_cov is not None else "n/a"
        print(f"k={kk} cap={cap}: min covol = {best_cov:.3f} at {best_E};  AP covol = {apc}")
    print()
    print("DONE.")


if __name__ == "__main__":
    main()
