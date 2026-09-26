#!/usr/bin/env python3
"""Wave-13 orchestrator check: for positive drift (q = 5, 7) the arbitrary
fixed-horizon edit price is exponentially small, although the undecided set
has positive density.

Setting (THM-4478): T_q(n) = n/2 or (q n + 1)/2.  An L-step descent
modification is any G with: for every n >= 2 some 1 <= j <= L has G^j(n) < n.
eps_L(q) = inf over G of the upper density of E(G) = {v : G(v) != T_q(v)}.

Catch-high construction (parameter W > 1).  B = sources n >= 2 with
T^j(n) >= n for all 1 <= j <= L.  For n in B let tau(n) = least j in [1, L-1]
with T^j(n) >= W n.  E = {T^tau(n)(n) : n in B, tau(n) finite}
                        U {n in B : tau(n) infinite},  G = 1 on E, T elsewhere.
Every E_high point v has a source n <= v/W, so E_high has upper density <= 1/W;
E_low lies (up to finitely many n) in the residue classes whose parity word
of length L-1 has final slope q^e/2^(L-1) < W.
"""
import math
import numpy as np


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)
    print("  ok:", msg)


def H2(x):
    return -x * math.log2(x) - (1 - x) * math.log2(1 - x)


def balance_beta(q):
    c = math.log(2) / math.log(q)
    lo, hi = 0.0, 1 / c - 1 - 1e-12
    hi = min(hi, 0.5)
    f = lambda b: b - (1 - H2(min(0.4999999, c * (1 + b))))
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(mid) > 0:
            hi = mid
        else:
            lo = mid
    return lo


def construct_and_verify(q, L, W, N):
    """Build E from sources <= N; return (|E cap [1,N]|/N, |B cap [2,N]|/N, #E_high, #E_low)."""
    n = np.arange(0, N + 1, dtype=np.int64)
    x = n.copy()
    desc = np.zeros(N + 1, dtype=bool)
    tau = np.full(N + 1, -1, dtype=np.int64)
    catch = np.zeros(N + 1, dtype=np.int64)
    for j in range(1, L + 1):
        x = np.where(x % 2 == 1, (q * x + 1) // 2, x // 2)
        desc |= x < n
        if j <= L - 1:
            newly = (tau < 0) & (x >= W * n)
            tau[newly] = j
            catch[newly] = x[newly]
    bad = ~desc
    bad[:2] = False
    Ehigh = set(int(v) for v in catch[bad & (tau > 0)])
    Elow = set(int(v) for v in np.nonzero(bad & (tau < 0))[0])
    E = Ehigh | Elow
    # verify descent of G for every 2 <= n <= N (partial E suffices: extra E points only help)
    xs = np.arange(2, N + 1, dtype=np.int64)
    start = xs.copy()
    done = np.zeros(len(xs), dtype=bool)
    Earr = np.array(sorted(E), dtype=np.int64)
    for j in range(1, L + 1):
        inE = np.isin(xs, Earr)
        xs = np.where(inE, 1, np.where(xs % 2 == 1, (q * xs + 1) // 2, xs // 2))
        done |= xs < start
    nE = sum(1 for v in E if v <= N)
    return nE / N, bad.sum() / N, done.all(), len([v for v in Ehigh if v <= N]), len([v for v in Elow if v <= N])


def word_band_density(q, Lm1, W):
    """fraction of words of length Lm1 whose final slope q^e/2^Lm1 < W (superset of E_low's classes)."""
    tot = 0
    for e in range(Lm1 + 1):
        if q ** e < W * 2 ** Lm1:
            tot += math.comb(Lm1, e)
    return tot / 2 ** Lm1


def band_population(q, L, K):
    """rho_L(K): fraction of words of length L with all prefix slopes q^(e_j)/2^j in [1, K] (THM-4478 Theorem A)."""
    cur = {0: 1}
    for j in range(1, L + 1):
        nxt = {}
        for e, cnt in cur.items():
            for b in (0, 1):
                e2 = e + b
                if 2 ** j <= q ** e2 <= K * 2 ** j:
                    nxt[e2] = nxt.get(e2, 0) + cnt
        cur = nxt
    return sum(cur.values()) / 2 ** L


def thm4478_lower(q, L, K):
    R = sum(k // q + 1 for k in range(L))
    M = (int(math.floor(math.log(K, q) + 1e-12)) + 1) * R
    return band_population(q, L, K) / (K * M)


print("1. exponents")
for q in (3, 5, 7, 9):
    c = math.log(2) / math.log(q)
    print(f"   q={q}: c_q = log_q 2 = {c:.6f}, 1 - H(c_q) = {1 - H2(c):.6f}" +
          (f", catch-high balance exponent theta_q = {balance_beta(q):.6f}" if c < 0.5 else " (negative drift: THM-4478 applies)"))
b5 = balance_beta(5)
check(0.0118 < b5 < 0.0120 and abs((1 - H2(math.log(2) / math.log(5))) - 0.013916) < 5e-6,
      "q=5: theta_5 = 0.0119 (upper construction) and 1 - H(log_5 2) = 0.01392 (lower bound exponent)")

print("2. undecided density of 5n+1 stays positive")
for L in (10, 14, 18, 22):
    # exact: fraction of words of length L with all prefix slopes > 1
    cur = {0: 1}
    for j in range(1, L + 1):
        nxt = {}
        for e, cnt in cur.items():
            for b in (0, 1):
                if 5 ** (e + b) > 2 ** j:
                    nxt[e + b] = nxt.get(e + b, 0) + cnt
        cur = nxt
    print(f"   L={L}: |Bad_L(5)|/2^L = {sum(cur.values())/2**L:.4f}")

print("3. the catch-high construction, built and verified on sources <= N")
N = 400000
for q in (5, 7):
    for L in (8, 12, 16):
        beta = balance_beta(q)
        for W in (2.0, 4.0, 8.0, 2 ** (beta * L) if 2 ** (beta * L) > 1.5 else 1.5):
            dens, bad, ok, nh, nl = construct_and_verify(q, L, W, N)
            bound = 1 / W + word_band_density(q, L - 1, W)
            assert ok, (q, L, W)
            assert dens <= bound + 0.01, (q, L, W, dens, bound)
            print(f"   q={q} L={L:2d} W={W:7.3f}: undecided {bad:.4f}  edit density {dens:.4f} "
                  f"(high {nh/N:.4f}, low {nl/N:.4f})  bound 1/W + band = {bound:.4f}  all n<=N descend: {ok}")
check(True, "q=5,7: G = 1 on E makes every 2 <= n <= 4e5 descend within L; the edit density respects 1/W + band(W)")

print("4. the bound as L grows (exact arithmetic on words), q = 5, W = 2^(theta L)")
th = balance_beta(5)
prev = None
for L in (50, 100, 200, 400, 800):
    W = 2 ** (th * L)
    ub = 1 / W + word_band_density(5, L - 1, W)
    print(f"   L={L:4d}: upper bound {ub:.3e}   (1/L) log2 = {math.log2(ub)/L:.5f}")
    if prev is not None:
        assert ub < prev
    prev = ub
check(math.log2(ub) / 800 < -0.0110, "the upper bound decays like 2^(-0.0119 L + o(L)) for q = 5")

print("5. THM-4478 Theorem A transferred to q = 5 (R_L with floor(k/q), window log_q K)")
for L in (8, 12, 16, 20, 24):
    best = max(thm4478_lower(5, L, K) for K in (1.5, 2, 3, 5, 8, 16, 32))
    ub = min(1 / W + word_band_density(5, L - 1, W) for W in (1.5, 2, 3, 4, 8, 16, 64, 256))
    print(f"   L={L:2d}: lower bound {best:.3e} <= eps_L(5) <= {ub:.3e}")
    assert best <= ub
check(True, "lower (capacity) <= upper (catch-high) at every tested L")
