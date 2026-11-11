#!/usr/bin/env python3
"""S91 helper (g20_local), part A: the class law of a gluing of finite seam rank r, checked on random distributions.

Setting: p prime (2 included), V = (Z/p)^r, B(u, v) = u^T M v, D1, D2 probability distributions on V with rational
weights, u ~ D1 and v ~ D2 independent, delta_i = D_i - p^(-r), pi_i = D_i(0), E_i = max_{xi != 0} |hat D_i(xi)|.

THEOREM (nondegenerate B).  P(B(u,v) = c) = 1/p + (pi1 + pi2 - p^(-r)) ([c = 0] - 1/p) + Rhat(c),
    Rhat(c) = delta1^T K_c delta2,  K_c(u, v) = [B(u,v) = c] - 1/p,   sum_c Rhat(c) = 0,
    Rhat = 0 if D1 or D2 is uniform, and with ||.|| the l2 norm on V
        |Rhat(c)| <= sigma_c ||delta1|| ||delta2||,   ||delta_i||^2 <= (1 - p^(-r)) E_i^2,
        sigma_c = p^((r-1)/2)  (c != 0),   sigma_0 = (p-1) p^((r-2)/2)  (r >= 2),   sigma_0 = (p-1)/p  (r = 1),
    and sigma_c is the exact norm of the bilinear form on sum-zero vectors (p odd; for p = 2, sigma_1 = sigma_0).
    In particular |Rhat(c)| <= ((p-1)/p) p^(r/2) E1 E2 (the bound asked for) holds, and it can be improved to
        (1 - p^(-r)) p^((r-1)/2) E1 E2  for c != 0,   (1 - p^(-r)) (p-1) p^((r-2)/2) E1 E2  for c = 0, r >= 2,
        ((p-1)/p)^2 E1 E2  for c = 0, r = 1.

What is checked here:
  1. the identity, sum Rhat = 0, Rhat = delta1^T K_c delta2, Rhat = 0 for a uniform factor: EXACT rationals;
  2. the l2 bound, squared, EXACT rationals (sigma_c^2 is rational);
  3. the E-bounds (floating point Fourier coefficients, relative tolerance 1e-9), with the worst observed ratios;
  4. the norms sigma_c: numpy singular values of the kernels restricted to sum-zero vectors (floating point);
  5. degenerate B of rank r' < r: the law with pi_i replaced by D1(left kernel), D2(right kernel), r by r';
  6. affine readout B(u,v) + l1(u) + l2(v) + k: the biased class is c0 = k - B(a, b) and the divisible states are the
     blind points u* = -a, v* = -b (l1 = B(., b), l2 = B(a, .));
  7. the general form: readout B(U,V) + z1 + z2 with joint laws of (U, z1), (V, z2):
        P(X = c) = 1/p + [(mu1*nu2)(c) - pi2/p] + [(nu1*mu2)(c) - pi1/p] - p^(-r) [(mu1*mu2)(c) - 1/p] + Rhat(c),
     mu_i = law of z_i, nu_i = law of z_i on {state_i = 0} (mass pi_i): the "one biased class" becomes a profile;
     Rhat = 0 when U is uniform and independent of z1; |Rhat| <= ((p-1)/p) p^(r/2) (1 - p^(-r)) E1 E2 with E_i the
     largest |E e((xi.U + t z_i)/p)| over xi != 0, t != 0.
Universe: p in {2, 3, 5, 7}, r in {1, 2, 3} (p^r <= 343), several random nondegenerate M (not symmetric) and random
distributions per case: dense random, sparse random (many zero weights), near-uniform, and one with an atom at 0.
"""
import itertools
import random
import sys
from fractions import Fraction

import numpy as np

rng = random.Random(20261110)
ok = True
worst = {}


def check(name, cond):
    global ok
    if not cond:
        print("FAIL  " + name)
    ok &= bool(cond)


def vecs(p, r):
    return list(itertools.product(range(p), repeat=r))


def det_mod(Mx, p):
    A = [row[:] for row in Mx]
    n = len(A)
    det = 1
    for c in range(n):
        piv = None
        for i in range(c, n):
            if A[i][c] % p:
                piv = i
                break
        if piv is None:
            return 0
        if piv != c:
            A[c], A[piv] = A[piv], A[c]
            det = -det
        det = det * A[c][c] % p
        inv = pow(A[c][c], p - 2, p)
        for i in range(c + 1, n):
            f = A[i][c] * inv % p
            if f:
                A[i] = [(x - f * y) % p for x, y in zip(A[i], A[c])]
    return det % p


def rand_invertible(p, r):
    while True:
        Mx = [[rng.randrange(p) for _ in range(r)] for _ in range(r)]
        if det_mod(Mx, p):
            return Mx


def bil(Mx, u, v, p):
    return sum(u[i] * Mx[i][j] * v[j] for i in range(len(u)) for j in range(len(v))) % p


def rand_dist(p, r, kind):
    N = p ** r
    if kind == "dense":
        w = [rng.randint(1, 50) for _ in range(N)]
    elif kind == "sparse":
        w = [rng.randint(1, 50) if rng.random() < 0.3 else 0 for _ in range(N)]
        if sum(w) == 0:
            w[rng.randrange(N)] = 1
    elif kind == "near":
        w = [1000 + rng.randint(-20, 20) for _ in range(N)]
    elif kind == "atom":
        w = [rng.randint(1, 10) for _ in range(N)]
        w[0] += 5 * N
    elif kind == "uniform":
        w = [1] * N
    else:
        raise ValueError
    return w


def fourier_max(w, p, r, mask=None):
    """largest |hat D(xi)| over xi != 0 (optionally only xi with mask[xi] true), floating point"""
    arr = np.array(w, dtype=np.float64).reshape((p,) * r)
    arr = arr / arr.sum()
    F = np.abs(np.fft.fftn(arr)).reshape(-1)
    F[0] = 0.0
    if mask is not None:
        F = F * np.array(mask, dtype=np.float64)
    return float(F.max())


def law_counts(w1, w2, table, p):
    """integer numerators of P(readout = c): sum of w1[u] w2[v] over pairs with table[u][v] = c"""
    out = [0] * p
    for iu, a in enumerate(w1):
        if a == 0:
            continue
        row = table[iu]
        for iv, b in enumerate(w2):
            if b:
                out[row[iv]] += a * b
    return out


def sigma2(p, r, c):
    """square of the exact norm of delta1^T K_c delta2 on sum-zero vectors (rational)"""
    if c % p != 0 and p > 2:
        return Fraction(p) ** (r - 1)
    if r >= 2:
        return Fraction((p - 1) ** 2) * Fraction(p) ** (r - 2)
    return Fraction(p - 1, p) ** 2


def e_bound_improved(p, r, c):
    s = float(sigma2(p, r, c)) ** 0.5
    return s * (1 - p ** (-r))


def run_nondegenerate(p, r, trials):
    V = vecs(p, r)
    N = p ** r
    kinds = ["dense", "sparse", "near", "atom"]
    for t in range(trials):
        Mx = rand_invertible(p, r)
        table = [[bil(Mx, u, v, p) for v in V] for u in V]
        for k1 in kinds:
            for k2 in (kinds if t == 0 else [rng.choice(kinds)]):
                w1, w2 = rand_dist(p, r, k1), rand_dist(p, r, k2)
                W1, W2 = sum(w1), sum(w2)
                cnt = law_counts(w1, w2, table, p)
                check("total", sum(cnt) == W1 * W2)
                pi1, pi2 = Fraction(w1[0], W1), Fraction(w2[0], W2)
                d1 = [Fraction(x, W1) - Fraction(1, N) for x in w1]
                d2 = [Fraction(x, W2) - Fraction(1, N) for x in w2]
                n1 = sum(x * x for x in d1)
                n2 = sum(x * x for x in d2)
                E1, E2 = fourier_max(w1, p, r), fourier_max(w2, p, r)
                rh = []
                for c in range(p):
                    P = Fraction(cnt[c], W1 * W2)
                    first = Fraction(1, p) + (pi1 + pi2 - Fraction(1, N)) * ((1 if c == 0 else 0) - Fraction(1, p))
                    rh.append(P - first)
                check(f"sum Rhat = 0 (p={p}, r={r})", sum(rh) == 0)
                # kernel form (exact) on one class per trial to save time
                c = rng.randrange(p)
                kern = sum(d1[iu] * d2[iv] * ((1 if table[iu][iv] == c else 0) - Fraction(1, p))
                           for iu in range(N) for iv in range(N)) if N <= 343 else None
                if kern is not None:
                    check(f"Rhat = delta1 K_c delta2 (p={p}, r={r})", kern == rh[c])
                for c in range(p):
                    s2 = sigma2(p, r, c)
                    check(f"l2 bound exact (p={p}, r={r}, c={c}, {k1}/{k2})", rh[c] * rh[c] <= s2 * n1 * n2)
                    check(f"||delta||^2 <= (1 - p^-r) E^2 (p={p}, r={r})",
                          float(n1) <= (1 - p ** (-r)) * E1 * E1 * (1 + 1e-9) + 1e-15)
                    asked = (p - 1) / p * p ** (r / 2) * E1 * E2
                    impr = e_bound_improved(p, r, c) * E1 * E2
                    a = abs(float(rh[c]))
                    check(f"asked E-bound (p={p}, r={r}, c={c}, {k1}/{k2})", a <= asked * (1 + 1e-9) + 1e-15)
                    check(f"improved E-bound (p={p}, r={r}, c={c}, {k1}/{k2})", a <= impr * (1 + 1e-9) + 1e-15)
                    key = (p, r, "c=0" if c == 0 else "c!=0")
                    if asked > 0:
                        wa, wi, wl = worst.get(key, (0.0, 0.0, 0.0))
                        l2 = (float(s2) * float(n1) * float(n2)) ** 0.5
                        worst[key] = (max(wa, a / asked), max(wi, a / impr), max(wl, a / l2 if l2 > 0 else 0.0))
        # uniform factor: Rhat = 0 exactly
        w1 = rand_dist(p, r, "dense")
        wu = rand_dist(p, r, "uniform")
        for (a, b) in ((w1, wu), (wu, w1)):
            cnt = law_counts(a, b, table, p)
            Wa, Wb = sum(a), sum(b)
            pa, pb = Fraction(a[0], Wa), Fraction(b[0], Wb)
            for c in range(p):
                first = Fraction(1, p) + (pa + pb - Fraction(1, N)) * ((1 if c == 0 else 0) - Fraction(1, p))
                check(f"Rhat = 0 for a uniform factor (p={p}, r={r})", Fraction(cnt[c], Wa * Wb) == first)


def run_norms(p, r):
    """singular values of Pi K_c Pi (Pi = projection on sum-zero vectors) against sigma_c, floating point"""
    V = vecs(p, r)
    N = p ** r
    Mx = rand_invertible(p, r)
    tab = np.array([[bil(Mx, u, v, p) for v in V] for u in V])
    Pi = np.eye(N) - np.ones((N, N)) / N
    out = []
    for c in range(p):
        K = (tab == c).astype(np.float64) - 1.0 / p
        sv = np.linalg.svd(Pi @ K @ Pi, compute_uv=False)
        top = float(sv[0])
        want = float(sigma2(p, r, c)) ** 0.5
        check(f"norm of K_{c} on sum-zero vectors = sigma_c (p={p}, r={r}): {top:.6f} vs {want:.6f}",
              abs(top - want) < 1e-8)
        out.append(top)
    return out


def run_degenerate(p, r, rk, trials):
    V = vecs(p, r)
    N = p ** r
    idx = {v: i for i, v in enumerate(V)}
    for _ in range(trials):
        X, Y = rand_invertible(p, r), rand_invertible(p, r)
        Dg = [[1 if (i == j and i < rk) else 0 for j in range(r)] for i in range(r)]
        Mx = [[sum(X[i][a] * Dg[a][b] * Y[b][j] for a in range(r) for b in range(r)) % p for j in range(r)]
              for i in range(r)]
        table = [[bil(Mx, u, v, p) for v in V] for u in V]
        Lk = [iu for iu in range(N) if all(x == 0 for x in table[iu])]                 # left kernel
        Rk = [iv for iv in range(N) if all(table[iu][iv] == 0 for iu in range(N))]     # right kernel
        check("kernel sizes", len(Lk) == p ** (r - rk) and len(Rk) == p ** (r - rk))
        # annihilators: xi in L-perp iff xi.u = 0 for all u in L
        maskL = [all(sum(a * b for a, b in zip(xi, V[iu])) % p == 0 for iu in Lk) for xi in V]
        maskR = [all(sum(a * b for a, b in zip(xi, V[iv])) % p == 0 for iv in Rk) for xi in V]
        for kind in ("dense", "sparse", "atom"):
            w1, w2 = rand_dist(p, r, kind), rand_dist(p, r, rng.choice(["dense", "sparse", "near"]))
            W1, W2 = sum(w1), sum(w2)
            cnt = law_counts(w1, w2, table, p)
            pi1 = Fraction(sum(w1[i] for i in Lk), W1)
            pi2 = Fraction(sum(w2[i] for i in Rk), W2)
            E1 = fourier_max(w1, p, r, maskL)
            E2 = fourier_max(w2, p, r, maskR)
            rh = []
            for c in range(p):
                first = Fraction(1, p) + (pi1 + pi2 - Fraction(1, p ** rk)) * ((1 if c == 0 else 0) - Fraction(1, p))
                rh.append(Fraction(cnt[c], W1 * W2) - first)
            check(f"degenerate: sum Rhat = 0 (p={p}, r={r}, rank={rk})", sum(rh) == 0)
            for c in range(p):
                impr = e_bound_improved(p, rk, c) * E1 * E2
                check(f"degenerate: improved E-bound with r' (p={p}, r={r}, rank={rk}, c={c})",
                      abs(float(rh[c])) <= impr * (1 + 1e-9) + 1e-15)
        # uniform factor
        w1, wu = rand_dist(p, r, "dense"), rand_dist(p, r, "uniform")
        cnt = law_counts(w1, wu, table, p)
        pi1 = Fraction(sum(w1[i] for i in Lk), sum(w1))
        for c in range(p):
            first = Fraction(1, p) + pi1 * ((1 if c == 0 else 0) - Fraction(1, p))
            check("degenerate: exact law for a uniform second factor", Fraction(cnt[c], sum(w1) * N) == first)


def run_affine(p, r, trials):
    V = vecs(p, r)
    N = p ** r
    idx = {v: i for i, v in enumerate(V)}
    for _ in range(trials):
        Mx = rand_invertible(p, r)
        a = tuple(rng.randrange(p) for _ in range(r))
        b = tuple(rng.randrange(p) for _ in range(r))
        k = rng.randrange(p)
        # l1(u) = B(u, b), l2(v) = B(a, v)
        table = [[(bil(Mx, u, v, p) + bil(Mx, u, b, p) + bil(Mx, a, v, p) + k) % p for v in V] for u in V]
        c0 = (k - bil(Mx, a, b, p)) % p
        ustar = idx[tuple((-x) % p for x in a)]
        vstar = idx[tuple((-x) % p for x in b)]
        check("blind points: the readout is constant = c0 on u = u* and on v = v*",
              all(x == c0 for x in table[ustar]) and all(table[iu][vstar] == c0 for iu in range(N)))
        for kind in ("dense", "sparse", "atom"):
            w1, w2 = rand_dist(p, r, kind), rand_dist(p, r, rng.choice(["dense", "sparse", "near"]))
            W1, W2 = sum(w1), sum(w2)
            cnt = law_counts(w1, w2, table, p)
            pi1, pi2 = Fraction(w1[ustar], W1), Fraction(w2[vstar], W2)
            E1, E2 = fourier_max(w1, p, r), fourier_max(w2, p, r)
            d1 = [Fraction(x, W1) - Fraction(1, N) for x in w1]
            d2 = [Fraction(x, W2) - Fraction(1, N) for x in w2]
            n1, n2 = sum(x * x for x in d1), sum(x * x for x in d2)
            rh = []
            for c in range(p):
                first = Fraction(1, p) + (pi1 + pi2 - Fraction(1, N)) * ((1 if c == c0 else 0) - Fraction(1, p))
                rh.append(Fraction(cnt[c], W1 * W2) - first)
            check(f"affine: sum Rhat = 0 (p={p}, r={r})", sum(rh) == 0)
            for c in range(p):
                s2 = sigma2(p, r, (c - c0) % p)
                check(f"affine: l2 bound exact around c0 (p={p}, r={r})", rh[c] * rh[c] <= s2 * n1 * n2)
                check(f"affine: improved E-bound (p={p}, r={r})",
                      abs(float(rh[c])) <= e_bound_improved(p, r, (c - c0) % p) * E1 * E2 * (1 + 1e-9) + 1e-15)


def run_profile(p, r, trials):
    """readout B(U, V) + z1 + z2, joint laws on (Z/p)^r x Z/p"""
    V = vecs(p, r)
    N = p ** r
    for _ in range(trials):
        Mx = rand_invertible(p, r)
        tab = [[bil(Mx, u, v, p) for v in V] for u in V]
        for kind in ("dense", "sparse", "product"):
            # w[iu][z]
            if kind == "product":       # U uniform and independent of z1: Rhat must vanish
                mz = [rng.randint(0, 9) for _ in range(p)]
                if sum(mz) == 0:
                    mz[0] = 1
                w1 = [[mz[z] for z in range(p)] for _ in range(N)]
            else:
                w1 = [[(rng.randint(1, 30) if (kind == "dense" or rng.random() < 0.4) else 0) for _ in range(p)]
                      for _ in range(N)]
                if sum(map(sum, w1)) == 0:
                    w1[0][0] = 1
            w2 = [[rng.randint(0, 30) for _ in range(p)] for _ in range(N)]
            W1, W2 = sum(map(sum, w1)), sum(map(sum, w2))
            cnt = [0] * p
            for iu in range(N):
                for z1 in range(p):
                    a = w1[iu][z1]
                    if a == 0:
                        continue
                    for iv in range(N):
                        base = tab[iu][iv] + z1
                        row = w2[iv]
                        for z2 in range(p):
                            if row[z2]:
                                cnt[(base + z2) % p] += a * row[z2]
            mu1 = [Fraction(sum(w1[iu][z] for iu in range(N)), W1) for z in range(p)]
            mu2 = [Fraction(sum(w2[iv][z] for iv in range(N)), W2) for z in range(p)]
            nu1 = [Fraction(w1[0][z], W1) for z in range(p)]
            nu2 = [Fraction(w2[0][z], W2) for z in range(p)]
            pi1, pi2 = sum(nu1), sum(nu2)

            def conv(f, h):
                return [sum(f[z] * h[(c - z) % p] for z in range(p)) for c in range(p)]

            m1n2, n1m2, m1m2 = conv(mu1, nu2), conv(nu1, mu2), conv(mu1, mu2)
            rh = []
            for c in range(p):
                first = (Fraction(1, p) + (m1n2[c] - pi2 / p) + (n1m2[c] - pi1 / p)
                         - Fraction(1, N) * (m1m2[c] - Fraction(1, p)))
                rh.append(Fraction(cnt[c], W1 * W2) - first)
            check(f"profile: sum Rhat = 0 (p={p}, r={r})", sum(rh) == 0)
            if kind == "product":
                check(f"profile: Rhat = 0 when U is uniform and independent of z1 (p={p}, r={r})",
                      all(x == 0 for x in rh))
            # E_i: joint Fourier coefficients at (xi != 0, t != 0)
            def emax(w, Wt):
                arr = np.array(w, dtype=np.float64).reshape((p,) * r + (p,)) / Wt
                F = np.abs(np.fft.fftn(arr)).reshape(N, p)
                return float(F[1:, 1:].max()) if p > 1 else 0.0
            E1, E2 = emax(w1, W1), emax(w2, W2)
            bound = (p - 1) / p * p ** (r / 2) * (1 - p ** (-r)) * E1 * E2
            for c in range(p):
                check(f"profile: E-bound (p={p}, r={r}, {kind})", abs(float(rh[c])) <= bound * (1 + 1e-9) + 1e-15)


def main():
    cases = [(2, 1), (2, 2), (2, 3), (2, 4), (3, 1), (3, 2), (3, 3), (5, 1), (5, 2), (5, 3), (7, 1), (7, 2), (7, 3),
             (11, 1), (11, 2), (13, 1), (13, 2)]
    print("1-3. nondegenerate law: exact identity, exact l2 bound, E-bounds")
    for (p, r) in cases:
        N = p ** r
        trials = 3 if N <= 49 else (2 if N <= 169 else 1)
        run_nondegenerate(p, r, trials)
    print(f"   {'p':>3} {'r':>2} {'class':>6} | worst |Rhat| / asked bound | / improved E-bound | / l2 bound")
    for key in sorted(worst):
        wa, wi, wl = worst[key]
        print(f"   {key[0]:3d} {key[1]:2d} {key[2]:>6} | {wa:26.4f} | {wi:18.4f} | {wl:10.4f}")
    print("4. exact norms sigma_c of K_c on sum-zero vectors (numpy singular values)")
    for (p, r) in [(2, 1), (2, 2), (2, 3), (3, 1), (3, 2), (3, 3), (5, 1), (5, 2), (5, 3), (7, 1), (7, 2), (11, 2)]:
        tops = run_norms(p, r)
        print(f"   p = {p}, r = {r}: sigma_0 = {tops[0]:.6f}, sigma_c (c != 0) = {tops[1]:.6f}")
    print("5. degenerate B (rank r' < r)")
    for (p, r, rk) in [(3, 2, 1), (3, 3, 1), (3, 3, 2), (5, 2, 1), (5, 3, 2), (7, 2, 1), (2, 3, 2), (2, 3, 1)]:
        run_degenerate(p, r, rk, 2)
    print("6. affine readout B(u,v) + l1(u) + l2(v) + k")
    for (p, r) in [(2, 2), (3, 1), (3, 2), (3, 3), (5, 1), (5, 2), (7, 2), (11, 1)]:
        run_affine(p, r, 2)
    print("7. profile law for B(U,V) + z1 + z2")
    for (p, r) in [(2, 1), (2, 2), (3, 1), (3, 2), (5, 1), (5, 2), (7, 1)]:
        run_profile(p, r, 2)
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
