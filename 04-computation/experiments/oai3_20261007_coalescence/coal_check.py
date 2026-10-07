#!/usr/bin/env python3
"""Independent check of 'Theorem H' (pick lane): Haar coalescence of affinely related
Collatz (Terras) orbits u = 3^k v + e.  State (k, N) with e = N / 3^max(0,-k).
Written from scratch (no lane code reused)."""
import random, sys, math
from fractions import Fraction

def step(k, N, beta):
    """one Terras step of the pair chain; beta = parity of v (fresh fair bit)."""
    if k >= 0:
        if N % 2 == 0:
            if beta == 0: return k, N // 2
            return k, (3*N + 1 - 3**k) // 2
        if beta == 0: return k + 1, (3*N + 1) // 2
        if k >= 1: return k - 1, (N - 3**(k-1)) // 2
        return -1, (3*N - 1) // 2            # k = 0 -> -1, e' = (3N-1)/6, stored N' = (3N-1)/2
    a = -k
    if N % 2 == 0:
        if beta == 0: return k, N // 2
        return k, (3*N + 3**a - 1) // 2
    if beta == 0: return k + 1, (N + 3**(a-1)) // 2
    return k - 1, (3*N - 1) // 2

def e_of(k, N):
    return Fraction(N, 3**max(0, -k))

def check_exact(trials=300, steps=1500, K=1700, seed=1):
    """chain vs direct 2-adic orbits mod 2^K (K - steps bits of slack)."""
    rng = random.Random(seed); bad = 0; total = 0
    for _ in range(trials):
        k0 = rng.randint(-6, 6); N0 = rng.randint(-10**6, 10**6)
        y = rng.getrandbits(K); mod = 1 << K
        inv3 = pow(3, -1, mod)
        e0num = N0 * pow(inv3, max(0, -k0), mod)
        u = (pow(3, k0, mod) if k0 >= 0 else pow(inv3, -k0, mod)) * y + e0num
        u %= mod; v = y; k, N = k0, N0; m = mod
        for n in range(steps):
            beta = v & 1
            # direct steps mod m
            v = (v // 2) if beta == 0 else (3*v + 1) // 2
            ub = u & 1
            u = (u // 2) if ub == 0 else (3*u + 1) // 2
            m >>= 1; v %= m; u %= m
            # predicted parity of u: beta xor (e mod 2)
            if ub != (beta ^ (N & 1)): bad += 1
            k, N = step(k, N, beta)
            # check u = 3^k v + e mod m
            t3 = pow(3, k, m) if k >= 0 else pow(pow(3, -1, m), -k, m)
            en = N * (pow(pow(3, -1, m), max(0, -k), m)) % m
            if (t3 * v + en - u) % m != 0: bad += 1
            total += 1
    return bad, total

def m_of(h):
    if h == 0: return None
    return 1 if h % 2 else 2 + ((h & -h).bit_length() - 1)

def run_length_stats(paths=20000, T=4000, seed=2):
    """empirical run lengths (consecutive even-e steps) by level |k|; compare with m_h + 1."""
    rng = random.Random(seed); stats = {}
    for _ in range(paths):
        k, N = 0, rng.choice([1, 3, 5, -7, 11, 1001])
        cur = None
        for t in range(T):
            if k == 0 and N == 0: break
            h = abs(k)
            if N % 2 == 0:
                if cur is None: cur = [h, 0]
                cur[1] += 1
            else:
                if cur is not None and cur[0] >= 1:
                    s = stats.setdefault(cur[0], [0, 0, 0]); s[0] += 1; s[1] += cur[1]; s[2] = max(s[2], cur[1])
                cur = None
            k, N = step(k, N, rng.getrandbits(1))
    return stats

def mean_excursion_ratio(e0, theta=0.5, reps=4000, Tcap=200000, seed=3):
    """E|e at first return to k=0|^theta / |e0|^theta from (0, e0) (one excursion)."""
    rng = random.Random(seed); acc = 0.0; n = 0; capped = 0
    for _ in range(reps):
        k, N = 0, e0; left = False
        for t in range(Tcap):
            k, N = step(k, N, rng.getrandbits(1))
            if k != 0: left = True
            if left and k == 0: break
        else:
            capped += 1; continue
        acc += abs(N) ** theta; n += 1
    return acc / n / abs(e0) ** theta, n, capped

def survival(start, paths, Ts, seed):
    rng = random.Random(seed); Tmax = max(Ts); alive = [0]*len(Ts)
    for _ in range(paths):
        k, N = start; t = 0
        while t < Tmax and not (k == 0 and N == 0):
            k, N = step(k, N, rng.getrandbits(1)); t += 1
        tt = t if (k == 0 and N == 0) else math.inf
        for i, T in enumerate(Ts):
            if tt > T: alive[i] += 1
    return [a / paths for a in alive]

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("exact", "all"):
        bad, tot = check_exact()
        print(f"[exact] chain vs direct 2-adic orbits: {tot} steps, {bad} mismatches"); sys.stdout.flush()
    if which in ("runs", "all"):
        st = run_length_stats()
        print("[runs] level h: count, mean run length, max, bound m_h+1")
        ok = True
        for h in sorted(st)[:24]:
            c, s, mx = st[h]; mean = s / c; b = m_of(h) + 1
            ok &= mean <= b + 0.05
            print(f"  h={h:3d} n={c:7d} mean={mean:.3f} max={mx:3d} bound={b}")
        print("[runs] mean <= m_h + 1 at all levels:", ok); sys.stdout.flush()
    if which in ("drift", "all"):
        th = 0.5; s = 2**th * (1 - math.sqrt(1 - 0.75**th)); rho = s * 2**-th
        print(f"[drift] theta={th}: s={s:.4f}, rho = s 2^-theta = {rho:.4f}")
        for e0 in [10**3 + 1, 10**6 + 1, 10**9 + 7, 10**12 + 39, 10**15 + 37, -(10**9 + 7)]:
            r, n, cap = mean_excursion_ratio(e0)
            print(f"  e0={e0:>18d}: E|e_ret|^1/2 / |e0|^1/2 = {r:.4f}  (n={n}, capped={cap})")
        sys.stdout.flush()
    if which in ("surv", "all"):
        Ts = [100, 400, 1600, 6400, 25600]
        q = survival((0, 1), 6000, Ts, 4)
        print("[surv] start (0,1) [y vs y+1]:", " ".join(f"T={T}: q={x:.4f} sqrtT*q={math.sqrt(T)*x:.2f}" for T, x in zip(Ts, q)))
