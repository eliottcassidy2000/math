# Auditor recompute-5128, task 1: mixing rates max_k rho(T_k)/2 of the twisted digit operators.
# Operator built from the definition in THM-5128 (B):
#   states (eta, beta) = (h mod m, h mod 2);  T_k[(eta,beta) -> (2 eta, 0)] = 1,
#   T_k[(eta,beta) -> (2 eta + 1, 1)] = zeta^(k (eta + beta)/2)   (division by 2 modulo m), zeta = e(1/m).
# Checks:
#   (a) identity  sum_{h < 2^n} zeta^(k q(h) + l h) = e_(0,0)^T T_k^n v_l  against brute force (all k, l; several n);
#   (b) brute-force decay rate of the exponential sums against rho(T_k);
#   (c) the same sums written with x, y:  sum e((k x + l y)/m)  (character (-2k) of q, 2(k + l) of h, plus a constant phase).
import sys
import numpy as np
from numba import njit
from a0_defs import psi_table

np.set_printoptions(linewidth=200)


def build_T(m, k):
    inv2 = (m + 1) // 2
    T = np.zeros((2 * m, 2 * m), dtype=np.complex128)
    for eta in range(m):
        for beta in range(2):
            i = 2 * eta + beta
            j0 = 2 * ((2 * eta) % m) + 0
            j1 = 2 * ((2 * eta + 1) % m) + 1
            T[i, j0] += 1.0
            c = ((eta + beta) * inv2) % m
            T[i, j1] += np.exp(2j * np.pi * ((k * c) % m) / m)
    return T


def rates(m):
    """rho(T_k)/2 for k = 1..m-1 (and k = 0 for control)."""
    out = {}
    for k in range(0, m):
        ev = np.linalg.eigvals(build_T(m, k))
        mods = np.sort(np.abs(ev))[::-1]
        out[k] = (mods[0] / 2.0, mods[1] / 2.0)
    return out


@njit(cache=True)
def block_hist(psi, m, nmax):
    """H[n, a, b] = #{h < 2^n : q(h) = a, h = b (mod m)} for n = 0..nmax (cumulative)."""
    H = np.zeros((nmax + 1, m, m), dtype=np.int64)
    # block n holds h with exactly n binary digits; h = 0 is put in n = 0
    H[0, 0, 0] = 1
    for n in range(1, nmax + 1):
        lo = 1 << (n - 1)
        hi = 1 << n
        for h in range(lo, hi):
            qv = h - psi[h]
            H[n, qv % m, h % m] += 1
    for n in range(1, nmax + 1):
        H[n] += H[n - 1]
    return H


def brute_sums(H, m):
    """S[n, k, l] = sum_{h < 2^n} zeta^(k q + l h)."""
    nmax = H.shape[0] - 1
    a = np.arange(m)
    E = np.exp(2j * np.pi * np.outer(a, a) / m)  # E[k, a] = zeta^(k a)
    S = np.zeros((nmax + 1, m, m), dtype=np.complex128)
    for n in range(nmax + 1):
        S[n] = E @ H[n].astype(np.float64) @ E.T
    return S


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    psi = psi_table(1 << nmax)
    primes_main = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 43, 73, 127]
    claimed = {3: 0.8660, 5: 0.7418, 7: 0.9320, 11: 0.8587, 13: 0.8272, 17: 0.8289, 19: 0.8053, 23: 0.7931,
               29: 0.7894, 31: 0.9002, 43: 0.8342, 73: 0.8516, 127: 0.8754}
    print("== rates max_{k != 0} rho(T_k)/2 (operator built from the definition) ==")
    print("m  rate(mine)   claimed   diff      argmax k   second modulus/2 at argmax   rate at k=0 (control, =1)   theta_m/2")
    allrates = {}
    for m in primes_main:
        r = rates(m)
        allrates[m] = r
        kbest = max(range(1, m), key=lambda k: r[k][0])
        rate = r[kbest][0]
        theta = np.log(1 / rate) / np.log(2 / rate ** 2)
        print(f"{m:4d} {rate:.6f}   {claimed[m]:.4f}   {rate - claimed[m]:+.5f}   {kbest:4d}   {r[kbest][1]:.6f}   {r[0][0]:.6f}   {theta / 2:.4f}")
    print()
    print("== all primes up to 89 (claim: the primes not listed are between 0.776 and 0.806) ==")
    listed = {3, 5, 7, 11, 13, 17, 19, 23, 31, 43, 73, 127}
    from sympy import primerange
    others = []
    for m in primerange(3, 90):
        r = allrates[m] if m in allrates else rates(m)
        rate = max(r[k][0] for k in range(1, m))
        tag = "listed" if m in listed else "other"
        if m not in listed:
            others.append((rate, m))
        print(f"  m = {m:3d}  rate = {rate:.6f}  ({tag})")
    print(f"  other primes up to 89: min = {min(others)[0]:.6f} (m = {min(others)[1]}), max = {max(others)[0]:.6f} (m = {max(others)[1]})")
    print()
    print("== extra: odd composite moduli mentioned in README section 3 ==")
    for m in [9, 15, 21, 27, 63]:
        r = rates(m)
        rate = max(r[k][0] for k in range(1, m))
        kbest = max(range(1, m), key=lambda k: r[k][0])
        print(f"  m = {m:3d}  rate = {rate:.6f}  argmax k = {kbest}")
    print()
    print(f"== brute force, h < 2^n, n <= {nmax} ==")
    for m in [3, 5, 7, 11, 13, 17, 19, 23, 31]:
        H = block_hist(psi, m, nmax)
        S = brute_sums(H, m)
        # (a) identity against the operator
        worst = 0.0
        a = np.arange(m)
        for k in range(m):
            T = build_T(m, k)
            row = np.zeros(2 * m, dtype=np.complex128)
            row[0] = 1.0
            for n in range(0, nmax + 1):
                if n > 0:
                    row = row @ T
                # v_l(eta, beta) = zeta^(l eta)
                for l in range(m):
                    v = np.exp(2j * np.pi * ((l * (np.arange(2 * m) // 2)) % m) / m)
                    val = row @ v
                    err = abs(val - S[n, k, l]) / max(1.0, abs(S[n, k, l]))
                    worst = max(worst, err)
        # (b) decay: for each k != 0, A_n(k) = max_l |S_n(k, l)|; slope of log A_n over the last 10 n; also l2 version
        r = allrates[m] if m in allrates else rates(m)
        n0 = nmax - 10
        best_fit = 0.0
        lines = []
        for k in range(1, m):
            A = np.max(np.abs(S[:, k, :]), axis=1)
            L2 = np.sqrt(np.sum(np.abs(S[:, k, :]) ** 2, axis=1))
            ns = np.arange(n0, nmax + 1)
            slope_max = np.polyfit(ns, np.log(A[n0:nmax + 1]), 1)[0]
            slope_l2 = np.polyfit(ns, np.log(L2[n0:nmax + 1]), 1)[0]
            fit = np.exp(slope_l2) / 2
            best_fit = max(best_fit, fit)
            lines.append((k, r[k][0], np.exp(slope_max) / 2, fit, (L2[nmax] / 2 ** nmax) ** (1 / nmax)))
        # worst-k summary
        kbest = max(range(1, m), key=lambda k: r[k][0])
        print(f"m = {m}: identity sum = e^T T_k^n v_l, max relative error over all k, l, n <= {nmax}: {worst:.2e}")
        print(f"   spectral rate max_k = {r[kbest][0]:.6f} (k = {kbest}); brute-force fitted decay (l2 over l, n in [{n0},{nmax}]) max_k = {best_fit:.6f}")
        if m <= 13:
            for (k, sr, fm, fl, root) in lines:
                print(f"     k = {k:2d}: rho(T_k)/2 = {sr:.6f}   fitted (max over l) = {fm:.6f}   fitted (l2 over l) = {fl:.6f}   (|S_n|_2/2^n)^(1/n) at n = {nmax}: {root:.6f}")
        # (c) the sums in x, y
        psi_n = psi[: 1 << nmax]
        h = np.arange(1 << nmax, dtype=np.int64)
        x = 2 * psi_n + 1
        y = 2 * h + 1
        worst_c = 0.0
        for (k, l) in [(1, 0), (1, 1), (2, 1), (m - 1, 2)]:
            ph = (k * (x % m) + l * (y % m)) % m
            cnt = np.bincount(ph, minlength=m)
            val = np.sum(cnt * np.exp(2j * np.pi * np.arange(m) / m))
            # k x + l y = (k + l)(2h + 1) - 2k q
            kq = (-2 * k) % m
            lh = (2 * (k + l)) % m
            pred = np.exp(2j * np.pi * ((k + l) % m) / m) * S[nmax, kq, lh]
            worst_c = max(worst_c, abs(val - pred) / max(1.0, abs(val)))
        print(f"   x,y-character sums agree with the (q,h) sums (4 pairs (k,l), n = {nmax}): max rel err {worst_c:.2e}")
        # normalised size of the largest non-trivial sum at nmax
        Sn = np.abs(S[nmax]).copy()
        Sn[0, 0] = 0
        kk, ll = np.unravel_index(np.argmax(Sn), Sn.shape)
        print(f"   largest non-trivial |S_n(k,l)|/2^n at n = {nmax}: {Sn[kk, ll] / 2 ** nmax:.3e} at (k,l) = ({kk},{ll});  rate^n = {r[kbest][0] ** nmax:.3e}")
        # k = 0, l != 0 control: bounded
        print(f"   k = 0: max_l!=0 |S_n(0,l)| at n = {nmax}: {np.max(np.abs(S[nmax, 0, 1:])):.3f} (claimed bounded by m/2 = {m / 2})")


if __name__ == "__main__":
    main()
