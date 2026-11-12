# verify_structure.py  (lane energy, S94)  -- exact integer checks of the lemmas in PROOFS.txt
#  (L1) eta(h) = 1 + 2h - sum_{i<j} kappa(j-i) h_i h_j,  kappa(1) = 2, kappa(d) = 2^(d-1) (d >= 2)
#  (L3) 2^n (eta(h) - 1) = h (2^(n+1) - rev_n h) + 2^n rho(h),
#       rho(h) = sum_i h_i [1/2 - (3/4) h_{i-1} + sum_{j <= i-2} h_j 2^(j-i-1)],
#       b/2 - N1/4 <= rho(h) < 3b/4  (b = number of blocks of ones, N1 = number of adjacent pairs of ones)
#  (L4) eta(1|r) = eta(1|rc(r))  (reverse-complement of the n-1 low bits)
#  (L2) eta(h) - eta(h') = L0(eps) - sum_{j not in S} K_j lambda_j(eps)   (random pairs)
#       lambda_j = X_{>j} + Xi_{<j} + eps_{j+1} + eps_{j-1}
#  (NET) number of unordered pairs h != h' of n-bit words with |h c(h) - h' c(h')| < T 2^n, for several T
import sys, random, time
import numpy as np
from numba import njit

def eta_list(N):
    e = [0] * N
    e[0] = 1
    for h in range(1, N):
        p = h >> 1
        e[h] = e[p] + (p + 2 - (p & 1) if h & 1 else 2 * p)
    return e

def kappa(d):
    return 2 if d == 1 else 1 << (d - 1)

def rev(x, n):
    r = 0
    for i in range(n):
        if (x >> i) & 1:
            r |= 1 << (n - 1 - i)
    return r

def check_small(nmax):
    N = 1 << nmax
    e = eta_list(N)
    ok1 = ok3 = ok4 = True
    maxrho = [-1e9] * (nmax + 1)
    minrho = [1e9] * (nmax + 1)
    for n in range(1, nmax + 1):
        for h in range(1 << (n - 1), 1 << n):
            bits = [(h >> i) & 1 for i in range(n)]
            ones = [i for i in range(n) if bits[i]]
            # L1
            q2 = 0
            for a in range(len(ones)):
                for b in range(a + 1, len(ones)):
                    q2 += kappa(ones[b] - ones[a])
            if e[h] != 1 + 2 * h - q2:
                ok1 = False
            # L3
            R = 0
            R4m = 0
            for i in ones:
                if i == 0 or bits[i - 1] == 0:
                    R += 4 << n >> 3          # 2^(n-1), kept as a multiple of 2^n/8 * ... (exact for n >= 1)
                else:
                    R4m += 1                  # contributes -2^(n-2) = -2^n/4
                for j in ones:
                    if j <= i - 2:
                        R += 1 << (n + j - i - 1)
            # exact identity, multiplied by 4 to keep integers:  4*2^n*(eta-1) = 4*h*c(h) + 4*2^n*rho
            R4 = 4 * R - R4m * (1 << n)
            if 4 * ((e[h] - 1) << n) != 4 * h * ((1 << (n + 1)) - rev(h, n)) + R4:
                ok3 = False
            N1 = R4m
            b = len(ones) - N1
            # bounds  b/2 - N1/4 <= rho < 3b/4   <=>   (2b - N1) 2^n <= R4 < 3 b 2^n
            if not ((2 * b - N1) * (1 << n) <= R4 < 3 * b * (1 << n)):
                ok3 = False
            maxrho[n] = max(maxrho[n], R4 / (4 << n))
            minrho[n] = min(minrho[n], R4 / (4 << n))
            # L4
            r = h - (1 << (n - 1))
            rc = rev(((1 << (n - 1)) - 1) ^ r, n - 1) if n > 1 else 0
            if e[h] != e[(1 << (n - 1)) + rc]:
                ok4 = False
    print(f"L1 quadratic form, all h < 2^{nmax}: {ok1}")
    print(f"L3 product formula with remainder (exact) and b/2 - N1/4 <= rho < 3b/4, all n <= {nmax}: {ok3}")
    print("   min rho(h) over n-bit words / n:", " ".join(f"{n}:{minrho[n]/n:.3f}" for n in range(4, nmax + 1, 2)))
    print("   max rho(h) over n-bit words / n:", " ".join(f"{n}:{maxrho[n]/n:.3f}" for n in range(4, nmax + 1, 2)))
    print(f"L4 involution eta(1|r) = eta(1|rc r), all n <= {nmax}: {ok4}")
    return e

def check_pairs(e, nmax, trials=200000):
    rnd = random.Random(12345)
    ok2 = oklam = True
    for _ in range(trials):
        n = rnd.randint(3, nmax)
        h = rnd.randrange(1 << (n - 1), 1 << n)
        hp = rnd.randrange(1 << (n - 1), 1 << n)
        if h == hp:
            continue
        eps = [((h >> i) & 1) - ((hp >> i) & 1) for i in range(n)]
        K = [((h >> i) & 1) & ((hp >> i) & 1) for i in range(n)]
        S = [i for i in range(n) if eps[i] != 0]
        L0 = 2 * sum(eps[i] << i for i in S)
        tw = 0
        for a in range(len(S)):
            for b in range(a + 1, len(S)):
                tw += kappa(S[b] - S[a]) * (eps[S[a]] + eps[S[b]])
        assert tw % 2 == 0
        L0 -= tw // 2
        tot = 0
        for j in range(n):
            if eps[j] == 0:
                lam = sum(kappa(abs(i - j)) * eps[i] for i in S)
                if K[j]:
                    tot += lam
                # lambda_j = X_{>j} + Xi_{<j} + eps_{j+1} + eps_{j-1}
                Xg = sum(eps[i] << (i - j - 1) for i in S if i > j)
                Xl = sum(eps[i] << (j - 1 - i) for i in S if i < j)
                adj = (eps[j + 1] if j + 1 < n else 0) + (eps[j - 1] if j >= 1 else 0)
                if lam != Xg + Xl + adj:
                    oklam = False
        if e[h] - e[hp] != L0 - tot:
            ok2 = False
    print(f"L2 collision equation linear in the common bits ({trials} random pairs, n <= {nmax}): {ok2}; lambda_j formula: {oklam}")

@njit(cache=True)
def vtable(n):
    lo = 1 << (n - 1)
    N = 1 << n
    v = np.empty(N - lo, dtype=np.int64)
    for h in range(lo, N):
        r = 0
        x = h
        for i in range(n):
            r = (r << 1) | (x & 1)
            x >>= 1
        v[h - lo] = h * ((1 << (n + 1)) - r)
    return v

@njit(cache=True)
def count_close(vs, tol):
    # number of unordered pairs i < j with vs[j] - vs[i] < tol (vs sorted)
    m = vs.shape[0]
    cnt = 0
    i = 0
    for j in range(m):
        while vs[j] - vs[i] >= tol:
            i += 1
        cnt += j - i
    return cnt

def net_counts(n0, n1):
    print("NET: unordered pairs h != h' of n-bit words with |Phi(h) - Phi(h')| < T, Phi(h) = h (2^(n+1) - rev_n h)/2^n")
    print("  n  words  [T: pairs/words, pairs/(words*T)] for T = 1, 4, 16, 5n/8, 3n/4;  exact collisions/word for comparison in cells_*.out")
    for n in range(n0, n1 + 1):
        t0 = time.time()
        v = vtable(n)
        v.sort()
        words = v.shape[0]
        row = []
        for T4 in (4, 16, 64, (5 * n + 1) // 2, 3 * n):          # T = T4/4 ; 5n/8 and 3n/4
            tol = (T4 << n) // 4
            c = count_close(v, tol)
            row.append(f"T={T4/4:g}: {c/words:.4f} {c/words/(T4/4):.4f}")
        print(f"  {n:2d} {words:9d}  " + " | ".join(row) + f"   ({time.time()-t0:.1f}s)")
        sys.stdout.flush()

if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 14
    n0 = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    n1 = int(sys.argv[3]) if len(sys.argv) > 3 else 24
    e = check_small(nmax)
    check_pairs(e, nmax)
    net_counts(n0, n1)
