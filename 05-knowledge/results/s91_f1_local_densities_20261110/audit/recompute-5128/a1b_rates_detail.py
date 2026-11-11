# Auditor recompute-5128, task 1 (detail):
#  (1) does the starting vector e_(0,0) see the dominant eigenvalue of T_k (is the spectral rate the true decay rate of the sums)?
#  (2) independent confirmation of rho(T_k) by a second method (norm growth of T_k^n, Gelfand) for the moduli whose
#      value matters for a finding (79, and the brief's 3, 5, 7, 11, 13);
#  (3) brute-force decay over a longer range for m = 7, 31, 79.
import sys
import numpy as np
from numba import njit
from a0_defs import psi_table
from a1_rates import build_T, block_hist, brute_sums


def gelfand(T, n=4096):
    """||T^n||_2^(1/n) by repeated squaring with renormalisation (upper estimate converging to rho)."""
    A = T.copy()
    logs = 0.0
    k = 1
    while k < n:
        nrm = np.linalg.norm(A, 2)
        A = A / nrm
        logs = 2 * (logs + np.log(nrm))
        A = A @ A
        k *= 2
    logs += np.log(np.linalg.norm(A, 2))
    return np.exp(logs / k)


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 26
    print("== (1) component of e_(0,0) and v_l on the dominant eigenvalues of T_k ==")
    for m in [3, 5, 7, 11, 13, 17, 19, 23, 31, 79]:
        best = None
        for k in range(1, m):
            T = build_T(m, k)
            lam, V = np.linalg.eig(T)
            if best is None or np.max(np.abs(lam)) > best[0]:
                best = (np.max(np.abs(lam)), k, lam, V, T)
        rho, k, lam, V, T = best
        Vinv = np.linalg.inv(V)
        e0 = np.zeros(2 * m, dtype=np.complex128)
        e0[0] = 1
        a = e0 @ V  # coefficients of the row vector
        order = np.argsort(-np.abs(lam))
        out = []
        for i in order[:4]:
            bmax = 0.0
            for l in range(m):
                v = np.exp(2j * np.pi * ((l * (np.arange(2 * m) // 2)) % m) / m)
                b = Vinv[i] @ v
                bmax = max(bmax, abs(a[i] * b))
            out.append((abs(lam[i]) / 2, np.angle(lam[i]) / (2 * np.pi), bmax))
        g = gelfand(T) / 2
        print(f"m = {m:3d}, k = {k:3d}: rho/2 (eig) = {rho / 2:.6f}; Gelfand ||T^4096||^(1/4096)/2 = {g:.6f}; cond(V) = {np.linalg.cond(V):.2e}")
        for (r, ang, c) in out:
            print(f"      |lambda|/2 = {r:.6f}, arg/2pi = {ang:+.4f}, max_l |coefficient in S_n(k,l)| = {c:.4e}")
    print()
    print(f"== (3) brute force decay, n up to {nmax} ==")
    psi = psi_table(1 << nmax)
    for m in [7, 31, 79]:
        H = block_hist(psi, m, nmax)
        S = brute_sums(H, m)
        # for each k the l2 norm over l of S_n(k, .) ; print the worst k (largest at nmax)
        L2 = np.sqrt(np.sum(np.abs(S) ** 2, axis=2))  # [n, k]
        kb = 1 + int(np.argmax(L2[nmax, 1:]))
        print(f"m = {m}: k with the largest sum at n = {nmax}: {kb}")
        print("   n   |S_n(k,.)|_2 / 2^n     ratio to previous / 2")
        for n in range(nmax - 12, nmax + 1):
            print(f"   {n:2d}  {L2[n, kb] / 2 ** n:.6e}   {L2[n, kb] / L2[n - 1, kb] / 2:.6f}")
        ns = np.arange(nmax - 12, nmax + 1)
        slope = np.polyfit(ns, np.log(L2[nmax - 12:nmax + 1, kb]), 1)[0]
        print(f"   fitted rate over the last 13 n: {np.exp(slope) / 2:.6f}")
        # max over k of fitted rates
        fr = []
        for k in range(1, m):
            sl = np.polyfit(ns, np.log(L2[nmax - 12:nmax + 1, k]), 1)[0]
            fr.append((np.exp(sl) / 2, k))
        print(f"   max over k of fitted rates: {max(fr)[0]:.6f} at k = {max(fr)[1]}")


if __name__ == "__main__":
    main()
