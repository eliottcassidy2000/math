# Auditor recompute-5128, extras (numerical content of THM-5128 (B3), (E), and one README number):
#  (1) (B3): 0 <= q <= h/2 and |rho(h) - (1 - phi(h)/2)| < n/(2h+1) for h with n binary digits;
#  (2) (E): D_n = #{h < 2^n : x = 1 mod 4} - #{x = 3 mod 4}, S_n = sum (-1)^q(h); the quadratic form
#      q(h) = sum_i h_i (h_(i+1) + h_(i+2)) mod 2; corank of the band matrix B_n over F_2; T_k T_k^* = 2I modulo 2^j;
#      the bound |sum_{h < 2^n} e(k q(h)/2^j)| <= 2^((n+j+1)/2);
#  (3) the rate at m = 1039 (README section 3: 0.764), by sparse eigenvalues.
import sys
import math
import numpy as np
from numba import njit
from a0_defs import psi_table32


@njit(cache=True)
def b3_check(psi, nb):
    worst = 0.0  # max of |rho - (1 - phi/2)| * (2h+1) / n
    bad = 0
    for h in range(1, 1 << nb):
        n = 0
        t = h
        phi = 0.0
        i = 0
        while t > 0:
            if t & 1:
                phi += 2.0 ** (-(i + 1))
            t >>= 1
            i += 1
            n += 1
        q = h - psi[h]
        if q < 0 or 2 * q > h:
            bad += 1
        x = 2.0 * psi[h] + 1.0
        y = 2.0 * h + 1.0
        dev = abs(x / y - (1.0 - phi / 2.0)) * y / n
        if dev > worst:
            worst = dev
    return worst, bad


@njit(cache=True)
def two_adic(psi, nb, jmax):
    D = np.zeros(nb + 1, dtype=np.int64)   # block contributions
    S = np.zeros(nb + 1, dtype=np.int64)
    badq = 0
    hist = np.zeros((nb + 1, 1 << jmax), dtype=np.int64)  # hist[n, q mod 2^jmax] over the block n
    hist[0, 0] = 1
    D[0] = 1   # h = 0: x = 1
    S[0] = 1
    for n in range(1, nb + 1):
        for h in range(1 << (n - 1), 1 << n):
            x = 2 * np.int64(psi[h]) + 1
            q = h - np.int64(psi[h])
            if x % 4 == 1:
                D[n] += 1
            else:
                D[n] -= 1
            if q & 1:
                S[n] -= 1
            else:
                S[n] += 1
            # quadratic form sum_i h_i (h_(i+1) + h_(i+2)) mod 2
            par = 0
            t = h
            while t > 0:
                if t & 1:
                    par ^= ((t >> 1) & 1) ^ ((t >> 2) & 1)
                t >>= 1
            if par != (q & 1):
                badq += 1
            hist[n, q & ((1 << jmax) - 1)] += 1
    return D, S, badq, hist


def corank_band(n):
    """corank over F_2 of the n x n matrix with ones at distances 1 and 2 from the diagonal (rows as python ints)."""
    rows = []
    for i in range(n):
        r = 0
        for d in (1, 2):
            if i - d >= 0:
                r |= 1 << (i - d)
            if i + d < n:
                r |= 1 << (i + d)
        rows.append(r)
    rank = 0
    for col in range(n):
        piv = None
        for i in range(rank, n):
            if (rows[i] >> col) & 1:
                piv = i
                break
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        for i in range(n):
            if i != rank and (rows[i] >> col) & 1:
                rows[i] ^= rows[rank]
        rank += 1
    return n - rank


def T2(j, k):
    """twisted matrix modulo 2^j on states w = h mod 2^(j+1): T[w -> 2w] = 1, T[w -> 2w+1] = e(k c(w)/2^j), c(w) = (w + (w mod 2))/2."""
    size = 1 << (j + 1)
    T = np.zeros((size, size), dtype=np.complex128)
    for w in range(size):
        T[w, (2 * w) % size] += 1.0
        c = (w + (w & 1)) // 2
        T[w, (2 * w + 1) % size] += np.exp(2j * np.pi * k * c / (1 << j))
    return T


def main():
    nb = int(sys.argv[1]) if len(sys.argv) > 1 else 26
    psi = psi_table32(1 << nb)
    worst, bad = b3_check(psi, min(nb, 24))
    print(f"== (B3), h < 2^{min(nb, 24)} ==")
    print(f"   violations of 0 <= q <= h/2: {bad};  max over h of |rho - (1 - phi/2)| (2h+1)/n = {worst:.4f}  (claim: < 1)")
    jmax = 6
    D, S, badq, hist = two_adic(psi, nb, jmax)
    Dc = np.cumsum(D)
    Sc = np.cumsum(S)
    print(f"== (E), h < 2^n, n <= {nb} ==")
    print(f"   q(h) = sum_i h_i (h_(i+1) + h_(i+2)) (mod 2): mismatches for 1 <= h < 2^{nb}: {badq}")
    print("   n:   " + " ".join(f"{n:7d}" for n in range(1, nb + 1)))
    print("   D_n: " + " ".join(f"{int(Dc[n]):7d}" for n in range(1, nb + 1)))
    print("   S_n: " + " ".join(f"{int(Sc[n]):7d}" for n in range(1, nb + 1)))
    zero_n = [n for n in range(1, nb + 1) if Dc[n] == 0]
    print(f"   D_n = 0 exactly for n in {zero_n};  n mod 12: {sorted(set(n % 12 for n in zero_n))}   (claim: n = 1, 4, 5, 7, 9, 10 mod 12)")
    nz = sorted(set(n % 12 for n in range(1, nb + 1) if Dc[n] != 0))
    print(f"   residues mod 12 with D_n != 0: {nz}")
    for n, v in [(18, -512), (20, 1024), (27, 16384)]:
        if n <= nb:
            print(f"   D_{n} = {int(Dc[n])}  (claimed {v})")
    cor = [corank_band(n) for n in range(1, 49)]
    print("   corank of the band matrix B_n over F_2, n = 1..48: " + " ".join(str(c) for c in cor))
    pattern = {1: 1, 2: 0, 3: 1, 4: 2, 5: 1, 0: 0}
    print(f"   pattern 1, 0, 1, 2, 1, 0 for n = 1, 2, 3, 4, 5, 0 (mod 6): {all(cor[n - 1] == pattern[n % 6] for n in range(1, 49))}")
    # |S_n| against 2^(n - rank/2), rank = n - corank
    ok = all(Sc[n] == 0 or abs(int(Sc[n])) == 2 ** (n - (n - cor[n - 1]) // 2) for n in range(1, nb + 1))
    print(f"   S_n is 0 or +-2^(n - rank/2) for n <= {nb}: {ok}")
    ok = all(Dc[n] == 0 or abs(int(Dc[n])) == 2 ** (n - (n - cor[n - 1]) // 2) for n in range(1, nb + 1))
    print(f"   D_n is 0 or +-2^(n - rank/2) for n <= {nb}: {ok}")
    # unitarity
    worst_u = 0.0
    for j in range(1, 7):
        for k in range(1, 1 << j, 2):
            T = T2(j, k)
            worst_u = max(worst_u, np.abs(T @ T.conj().T - 2 * np.eye(T.shape[0])).max())
    print(f"   max |T_k T_k^* - 2I| over j <= 6, odd k: {worst_u:.2e}")
    # bound on character sums modulo 2^j
    ch = np.cumsum(hist, axis=0)
    worst_ratio = 0.0
    for j in range(1, jmax + 1):
        for k in range(1, 1 << j, 2):
            ph = np.exp(2j * np.pi * k * (np.arange(1 << jmax) % (1 << j)) / (1 << j))
            for n in range(1, nb + 1):
                val = abs(np.sum(ch[n] * ph))
                worst_ratio = max(worst_ratio, val / 2 ** ((n + j + 1) / 2))
    print(f"   max over j <= {jmax}, odd k, n <= {nb} of |sum_(h < 2^n) e(k q(h)/2^j)| / 2^((n+j+1)/2) = {worst_ratio:.4f}  (claim: <= 1)")
    # rate at m = 1039
    if len(sys.argv) > 2 and sys.argv[2] == "1039":
        import scipy.sparse as sp
        from scipy.sparse.linalg import eigs
        m = 1039
        inv2 = (m + 1) // 2
        eta = np.repeat(np.arange(m), 2)
        beta = np.tile(np.arange(2), m)
        idx = 2 * eta + beta
        j0 = 2 * ((2 * eta) % m)
        j1 = 2 * ((2 * eta + 1) % m) + 1
        c = ((eta + beta) * inv2) % m
        best = (0.0, 0)
        vals = []
        rng = np.random.default_rng(1)
        for k in range(1, (m - 1) // 2 + 1):
            ph = np.exp(2j * np.pi * ((k * c) % m) / m)
            T = sp.csr_matrix((np.concatenate([np.ones(2 * m), ph]), (np.concatenate([idx, idx]), np.concatenate([j0, j1]))), shape=(2 * m, 2 * m))
            try:
                ev = eigs(T, k=4, which="LM", return_eigenvectors=False, v0=rng.standard_normal(2 * m) + 0j, maxiter=20000, tol=1e-9)
                r = float(np.max(np.abs(ev))) / 2
            except Exception as e:  # fall back to power iteration
                v = rng.standard_normal(2 * m) + 1j * rng.standard_normal(2 * m)
                lg = 0.0
                for it in range(6000):
                    v = T @ v
                    nv = np.linalg.norm(v)
                    v /= nv
                    if it >= 1000:
                        lg += math.log(nv)
                r = math.exp(lg / 5000) / 2
            vals.append(r)
            if r > best[0]:
                best = (r, k)
        print(f"== rate at m = 1039: max_k rho(T_k)/2 = {best[0]:.6f} at k = {best[1]} (and m - k); min over k = {min(vals):.6f}   (README section 3: 0.764)")


if __name__ == "__main__":
    main()
