#!/usr/bin/env python3
"""
collatz_procgen_20260923_cube_theta_lattice.py

Cube-theta lane (HYP-9127): 2-adic lattice experiments ("hidden structure" hunt).

  T1  rational reconstruction: no rational a/b with |a|,|b| <= 2^((N-1)/2) equals X
      (independent Python half-Euclid at N = 10^5; PARI bestappr at N = 10^6 and 10^7 if gp exists)
  T2  2-adic approximation profile: for N on a grid, the shortest (a,b) with a == bX mod 2^N;
      ratio log2 max(|a|,|b|) / (N/2) (generic ~ 1; an anomalously good approximation dips below 1)
      for X, Theta_2(rho), random 2-adic integers, and a lacunary control sum 2^(3^k) (dips to 2/3)
  T3  LLL relation hunts among 1, X, X^2, X^3, ..., the auxiliary series Y_j = sum k^j rho^(k^3),
      Theta_2(rho) and products, at N = 20000 and 40000 bits; achieved exponent vs the
      Dirichlet/pigeonhole baseline, with random and algebraic controls.

Usage: python3 collatz_procgen_20260923_cube_theta_lattice.py [--quick]
"""
import math
import os
import random
import shutil
import subprocess
import sys
import tempfile
import time
from fractions import Fraction

import flint
import gmpy2
from gmpy2 import mpz

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from collatz_procgen_20260923_cube_theta_core import (  # noqa: E402
    L_Y3, M0_Y3, X_mod, aux_mod, theta2_mod, theta_mod, v2, log2abs)

QUICK = "--quick" in sys.argv
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
SCR = os.path.join(ROOT, "scratch", "procgen_cube")


def hdr(s):
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)
    sys.stdout.flush()


def _half_euclid(x, N, A):
    """Euclid on (2^N, x): pairs (r, t) with r == t x (mod 2^N); stop at the first r < A.
    Returns the last two pairs (r0, t0), (r1, t1): they form a basis of the lattice
    {(a, b) : a == b x (mod 2^N)} (determinant +-2^N)."""
    r0, r1 = mpz(1) << N, mpz(x) % (mpz(1) << N)
    t0, t1 = mpz(0), mpz(1)
    while r1 != 0 and r1 >= A:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        t0, t1 = t1, t0 - q * t1
    return (r0, t0), (r1, t1)


def ratrecon(x, N):
    """Wang's rational reconstruction with A = B = 2^(N/2 - 1) (so 2AB < 2^N): returns the unique
    a/b (b odd, |a| < A, 0 < |b| < B, a == b x mod 2^N) if it exists, else None.
    None is a certificate: no rational of height < 2^(N/2 - 1) is congruent to x mod 2^N."""
    A = mpz(1) << (N // 2 - 1)
    (_, _), (r1, t1) = _half_euclid(x, N, A)
    if t1 != 0 and abs(t1) < A and gmpy2.gcd(r1, t1) == 1:
        return (r1, t1)
    return None


def short_2d(x, N):
    """Lagrange-Gauss reduced shortest vector (a, b) of {(a,b): a == b x mod 2^N}."""
    A = mpz(1) << (N // 2)
    (r0, t0), (r1, t1) = _half_euclid(x, N, A)
    u, v = [r0, t0], [r1, t1]
    n2 = lambda w: w[0] * w[0] + w[1] * w[1]
    while True:
        if n2(u) < n2(v):
            u, v = v, u
        # now |v| <= |u| ; reduce u by v
        if n2(v) == 0:
            break
        mu = (u[0] * v[0] + u[1] * v[1] + n2(v) // 2) // n2(v)
        if mu == 0:
            break
        u = [u[0] - mu * v[0], u[1] - mu * v[1]]
    cands = [v, u, [u[0] + v[0], u[1] + v[1]], [u[0] - v[0], u[1] - v[1]]]
    cands = [w for w in cands if w[0] or w[1]]
    w = min(cands, key=lambda w: max(abs(w[0]), abs(w[1])))
    return w[0], w[1]


# ----------------------------------------------------------------------------- T1
def section_T1():
    hdr("T1  rational reconstruction certificates for X = sum rho^(k^3), rho = 2^10/3^9")
    N = 20000 if QUICK else 100000
    t0 = time.time()
    X = X_mod(N)
    a, b = short_2d(X, N)
    ok = ((a - b * X) % (mpz(1) << N)) == 0
    rr = ratrecon(X, N)
    print(f"[python] N = {N}: Lagrange-shortest (a,b) with a == bX mod 2^N: log2 max(|a|,|b|) = "
          f"{max(log2abs(a), log2abs(b)):.2f} (N/2 = {N/2:.0f}); congruence holds: {ok}  ({time.time()-t0:.1f}s)")
    print(f"         Wang reconstruction (A = B = 2^(N/2-1)): {'NONE' if rr is None else rr}"
          + (f"  => FINITE-EXACT: X is not a rational of height < 2^{N//2-1}" if rr is None else ""))
    gp = shutil.which("gp")
    if gp is None:
        print("[pari] gp not found: skipped")
        return
    Ns = [10 ** 5, 10 ** 6] if QUICK else [10 ** 5, 10 ** 6, 10 ** 7]
    os.makedirs(SCR, exist_ok=True)
    for Nn in Ns:
        t0 = time.time()
        Xn = X_mod(Nn)
        fn = os.path.join(SCR, f"X_mod_2e{int(math.log10(Nn))}.txt")
        with open(fn, "w") as f:
            f.write("0x" + Xn.digits(16))   # GMP conversion (subquadratic); GP reads 0x literals
        script = (f'x = eval(read("{fn}")); N = {Nn}; '
                  f'r = bestappr(Mod(x, 2^N)); '
                  f'print(if(type(r)=="t_VEC" && #r==0, "NONE", Str("FOUND ", r))); quit;')
        with tempfile.NamedTemporaryFile("w", suffix=".gp", delete=False, dir=SCR) as g:
            g.write(script)
            gname = g.name
        out = subprocess.run([gp, "-q", "--stacksize=600000000", gname], stdin=subprocess.DEVNULL,
                             capture_output=True, text=True, timeout=3600).stdout.strip()
        os.unlink(gname)
        os.unlink(fn)
        print(f"[pari] N = {Nn}: bestappr(Mod(X, 2^N)) -> {out[:80]}  ({time.time()-t0:.1f}s)")
        if out == "NONE":
            print(f"        => FINITE-EXACT: X is not a rational a/b with |a|, |b| <= 2^{(Nn-1)//2}"
                  " (PARI: unique reconstruction bound sqrt(N/2))")


# ----------------------------------------------------------------------------- T2
def section_T2():
    hdr("T2  2-adic approximation profile: h(N) = log2 of the shortest (a,b), ratio h/(N/2)")
    Nmax = 6000 if QUICK else 20000
    grid = sorted(set([int(64 * 1.03 ** i) for i in range(0, 400) if int(64 * 1.03 ** i) <= Nmax]))
    rng = random.Random(20260923)
    numbers = [
        ("X (cubes, rho)", X_mod(Nmax)),
        ("Theta_2(rho) (squares)", theta2_mod(Nmax)),
        ("Y_1 = sum k rho^(k^3)", aux_mod(Nmax, 1)),
        ("random 2-adic #1", mpz(rng.getrandbits(Nmax))),
        ("random 2-adic #2", mpz(rng.getrandbits(Nmax))),
        ("control sum 2^(3^k)", sum((mpz(1) << (3 ** k)) for k in range(0, 12) if 3 ** k < Nmax)),
        ("control sum 2^(k^3)", sum((mpz(1) << (k ** 3)) for k in range(0, 40) if k ** 3 < Nmax)),
    ]
    print(f"grid: {len(grid)} values of N in [64, {Nmax}] (geometric, ratio 1.03)")
    print(" number                      | all N: min ratio (at N) | N >= 2000: min ratio (at N) | mean | #N>=2000 with ratio < 0.95")
    for name, x in numbers:
        rs = []
        for N in grid:
            a, b = short_2d(x % (mpz(1) << N), N)
            h = max(log2abs(a), log2abs(b))
            rs.append((h / (N / 2), N))
        mn = min(rs)
        big = [(r, N) for r, N in rs if N >= 2000]
        mnb = min(big)
        mean = sum(r for r, _ in big) / len(big)
        low = sum(1 for r, _ in big if r < 0.95)
        print(f" {name:27s} | {mn[0]:.4f} ({mn[1]:6d})        | {mnb[0]:.4f} ({mnb[1]:6d})            | {mean:.4f} | {low}")
    print("Small-N dips for X and Y_1 are the partial sums S_K (N <= 10(K+1)^3, height 14.26 K^3 bits),")
    print("which beat the Dirichlet scale only for K <= 3; beyond N ~ 2000 nothing dips.")
    print("Reading: ratio < 1 by a fixed factor on a sequence of N = an approximation beating Dirichlet")
    print("(irrationality exponent > 2). The lacunary controls show the detector works: sum 2^(3^k) dips")
    print("to 2/3 at N = 3^(k+1) (Hadamard gaps) and the unit cubic control sum 2^(k^3) dips only by")
    print("O(1/k) (its partial sums have exponent (k+1)^3/k^3).")


# ----------------------------------------------------------------------------- T3
def relation_lattice(xis, N):
    """Basis of {c in Z^m : sum c_i xi_i == 0 mod 2^N}, xi_0 = 1."""
    m = len(xis)
    mod = mpz(1) << N
    rows = []
    rows.append([int(mod)] + [0] * (m - 1))
    for i in range(1, m):
        r = [0] * m
        r[0] = int((-xis[i]) % mod)
        r[i] = 1
        rows.append(r)
    return flint.fmpz_mat(rows)


def hunt(name, xis, N, Ncheck):
    t0 = time.time()
    m = len(xis)
    B = relation_lattice([x % (mpz(1) << N) for x in xis], N).lll()
    best = None
    for i in range(B.nrows()):
        v = [int(B[i, j]) for j in range(m)]
        nm = max(abs(t) for t in v)
        if nm and (best is None or nm < best[0]):
            best = (nm, v)
    nm, v = best
    lg = math.log2(nm)
    # certificate: every nonzero lattice vector (in particular every EXACT relation, which lies in the
    # lattice at every precision) has Euclidean norm >= min_i |b_i*| (Gram-Schmidt of the LLL basis)
    rows = [[Fraction(int(B[i, j])) for j in range(m)] for i in range(B.nrows())]
    gs = []
    for r in rows:
        w = r[:]
        for g, gg in gs:
            mu = sum(a * b for a, b in zip(r, g)) / gg
            w = [a - mu * b for a, b in zip(w, g)]
        gs.append((w, sum(a * a for a in w)))
    min_gs = min(gg for _, gg in gs)
    lg_gs = 0.5 * (math.log2(min_gs.numerator) - math.log2(min_gs.denominator))
    # check at higher precision
    modc = mpz(1) << Ncheck
    s = sum(mpz(v[i]) * (xis[i] % modc) for i in range(m)) % modc
    val = v2(s)
    tau = N / lg
    print(f" {name:34s} | m={m:2d} | log2|c|max = {lg:9.1f} | N/m = {N/m:8.1f} | ratio {lg/(N/m):.4f} |"
          f" tau/m = {tau/m:.4f} | v2 at 2^{Ncheck}: {val if val != math.inf else 'EXACT'} |"
          f" no exact relation with |c|_2 < 2^{lg_gs:.1f}  [{time.time()-t0:.1f}s]")
    sys.stdout.flush()
    return lg / (N / m)


def section_T3():
    hdr("T3  LLL relation hunts (Dirichlet baseline: shortest relation vector ~ 2^(N/m))")
    Ns = [20000] if QUICK else [20000, 40000]
    for N in Ns:
        Nc = 2 * N
        X = X_mod(Nc)
        Y1, Y2, Y3, Y4 = (aux_mod(Nc, j) for j in (1, 2, 3, 4))
        T2 = theta2_mod(Nc)
        # values of G(u,v) = sum u^k v^(k^2) rho^(k^3) at non-tail points (u,v) = (rho^a, rho^b)
        Gk = lambda a, b: theta_mod(Nc, lambda k: k ** 3 + b * k * k + a * k)
        mod = mpz(1) << Nc
        rng = random.Random(N)
        R1, R2, R3 = (mpz(rng.getrandbits(Nc)) for _ in range(3))
        # algebraic control: the 2-adic cube root of 17 (x^3 = 17; Hensel since 3x^2 is odd)
        cr = mpz(1)
        for _ in range(int(math.log2(Nc)) + 3):
            # Newton: x <- x - (x^3 - 17)/(3x^2)
            cr = (cr - (cr ** 3 - 17) * gmpy2.invert(3 * cr * cr, mod)) % mod
        assert (cr ** 3 - 17) % mod == 0
        print(f"-- N = {N} bits (checked at 2^{Nc}) --")
        powX = [mpz(1)]
        for _ in range(7):
            powX.append((powX[-1] * X) % mod)
        tests = [
            ("1, X, X^2, X^3 (degree <= 3)", powX[:4]),
            ("1, X, ..., X^6 (degree <= 6)", powX[:7]),
            ("1, X, Y1, Y2 (G and u-derivs at (1,1))", [mpz(1), X, Y1, Y2]),
            ("1, X, Y1, Y2, Y3, Y4", [mpz(1), X, Y1, Y2, Y3, Y4]),
            ("1, X, Theta2, X*Theta2", [mpz(1), X, T2, (X * T2) % mod]),
            ("1, X, G(rho,1), G(1,rho), G(rho,rho)", [mpz(1), X, Gk(1, 0), Gk(0, 1), Gk(1, 1)]),
            ("1, X, X^2, Y1, X*Y1, Y1^2", [mpz(1), X, (X * X) % mod, Y1, (X * Y1) % mod, (Y1 * Y1) % mod]),
            ("CONTROL random: 1, r, r^2, r^3", [mpz(1), R1, (R1 * R1) % mod, (R1 ** 3) % mod]),
            ("CONTROL random: 1, r1, r2, r3", [mpz(1), R1, R2, R3]),
            ("CONTROL algebraic: 1, c, c^2, c^3", [mpz(1), cr, (cr * cr) % mod, (cr ** 3) % mod]),
        ]
        for name, xis in tests:
            hunt(name, xis, N, Nc)
    print("Reading: ratio ~ 1 (tau/m ~ 1) is the Dirichlet/pigeonhole baseline for a generic vector;")
    print("a relation would show ratio << 1 and survive the check at double precision (EXACT).")


if __name__ == "__main__":
    t0 = time.time()
    print("collatz_procgen_20260923_cube_theta_lattice.py" + (" --quick" if QUICK else ""))
    section_T1()
    section_T2()
    section_T3()
    print(f"\n[lattice done in {time.time()-t0:.1f}s]")
