#!/usr/bin/env python3
"""Lyapunov symmetric-maximizer conjecture at order six (FINITE-NUMERICAL; mac-mini six-seven, 2026-10-06).

L_A(X) = AX + XA^T.  Conjecture: ||L_A|Sym_n|| >= ||L_A|Skew_n|| (Frobenius-induced operator norms).
Known for n <= 5; Kressner-Vandereycken (arXiv:2608.20875) give an integer counterexample at n = 7,
padding gives every n >= 7; n = 6 is OPEN (repo: 07-reflections/octonions-at-the-center-..., section 8).
We maximise h(A) = log(sigma_skew / sigma_sym) (scale invariant; counterexample iff h > 0) with the exact
gradient U V^T + U^T V of a simple top singular value.

modes
  search  n trials seed : structured random starts (lower block-triangular layer patterns, sparse diagonals,
                          random permutations) + L-BFGS.  CONTROL: run it at n = 7, where counterexamples exist.
  kv                    : the KV matrix, its local maximum, and the 7 coordinate deletions (6x6) re-optimised.
  strip                 : continuation along the KV branch with mu(A) = lambda_min(A^T A + A A^T)/|A|_F^2 <= m
                          (m -> 0 means A -> A_6 (+) 0, a 6x6 matrix padded by a zero row and column, which is a
                          counterexample iff A_6 is); at each m the near-null direction is stripped and the 6x6
                          block is re-optimised from 7 starts.
Run: python3 sixseven_20261006_lyapunov_six.py kv | strip | search 6 2000 1
"""
import sys, time, itertools
import numpy as np
from scipy.optimize import minimize

A_KV = np.array([[0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0], [6, -5, -13, 0, 0, 0, 0], [-14, -18, 0, 0, 0, 0, 0],
                 [12, -12, 11, 0, 0, 0, 0], [0, 0, 0, -6, 0, -14, 0], [0, 0, 0, 6, -18, 0, 0]], float)


def bases(n):
    S, K = [], []
    for i in range(n):
        E = np.zeros((n, n)); E[i, i] = 1; S.append(E.ravel())
    for i in range(n):
        for j in range(i + 1, n):
            E = np.zeros((n, n)); E[i, j] = E[j, i] = 2 ** -.5; S.append(E.ravel())
            F = np.zeros((n, n)); F[i, j] = 2 ** -.5; F[j, i] = -2 ** -.5; K.append(F.ravel())
    return np.array(S).T, np.array(K).T


class Lyap:
    def __init__(self, n):
        self.n = n
        self.BS, self.BK = bases(n)
        self.I = np.eye(n)

    def top(self, A, B):
        M = B.T @ (np.kron(self.I, A) + np.kron(A, self.I)) @ B
        U, s, Vt = np.linalg.svd(M)
        return s[0], (B @ U[:, 0]).reshape(self.n, self.n), (B @ Vt[0]).reshape(self.n, self.n)

    def fg(self, a):
        """returns (-h, -grad h)"""
        A = a.reshape(self.n, self.n)
        ss, Us, Vs = self.top(A, self.BS)
        sk, Uk, Vk = self.top(A, self.BK)
        g = (Uk @ Vk.T + Uk.T @ Vk) / sk - (Us @ Vs.T + Us.T @ Vs) / ss
        return -(np.log(sk) - np.log(ss)), -g.ravel()

    def h(self, A):
        return -self.fg(np.asarray(A).ravel())[0]

    def climb(self, A, maxiter=20000):
        r = minimize(self.fg, np.asarray(A).ravel(), jac=True, method="L-BFGS-B", options={"maxiter": maxiter})
        return -r.fun, r.x.reshape(self.n, self.n)


def check_kron():
    L = Lyap(5); A = np.random.randn(5, 5); X = np.random.randn(5, 5)
    return np.abs(((np.kron(L.I, A) + np.kron(A, L.I)) @ X.ravel()).reshape(5, 5) - (A @ X + X @ A.T)).max()


def mode_kv():
    L7 = Lyap(7)
    s_sym = L7.top(A_KV, L7.BS)[0] ** 2; s_skew = L7.top(A_KV, L7.BK)[0] ** 2
    print(f"KV matrix: sigma_sym^2 = {s_sym:.6f} < 1196 < sigma_skew^2 = {s_skew:.6f}; h = {np.log(s_skew / s_sym) / 2:.4e}")
    hmax, Aopt = L7.climb(A_KV)
    print(f"local maximum of h from KV: {hmax:.4e}  (sigma_skew/sigma_sym = {np.exp(hmax):.6f})")
    L6 = Lyap(6)
    for k in range(7):
        B = np.delete(np.delete(Aopt, k, 0), k, 1)
        print(f"  delete coordinate {k}: h(6x6) = {L6.h(B):+.3e} -> local max {L6.climb(B)[0]:+.3e}")
    np.save("sixseven_20261006_lyapunov_kv_localmax.npy", Aopt)


def mode_strip():
    L7, L6 = Lyap(7), Lyap(6)

    def mu(a):
        A = a.reshape(7, 7); return np.linalg.eigvalsh(A.T @ A + A @ A.T)[0] / (A * A).sum()

    def strip(A):
        w, V = np.linalg.eigh(A.T @ A + A @ A.T)
        Q = np.column_stack([V[:, 1:], V[:, 0]]); return (Q.T @ A @ Q)[:6, :6]

    a = L7.climb(A_KV)[1].ravel()
    best6 = -1
    rows = []
    for m in np.geomspace(0.15, 3.2e-3, 26):
        r = minimize(lambda x: L7.fg(x)[0], a, jac=lambda x: L7.fg(x)[1], method="SLSQP",
                     constraints=[{"type": "ineq", "fun": lambda x: m - mu(x)}], options={"maxiter": 5000, "ftol": 1e-15})
        a = r.x; B6 = strip(a.reshape(7, 7)); h6 = L6.h(B6)
        hb = L6.climb(B6)[0]
        for t in range(6):
            x0 = B6 + 1e-2 * np.abs(B6).max() * np.random.default_rng(t).standard_normal((6, 6))
            hb = max(hb, L6.climb(x0)[0])
        best6 = max(best6, hb); rows.append((m, -r.fun))
        print(f"m = {m:.3e}: 7x7 h = {-r.fun:+.4e} | stripped 6x6 h = {h6:+.4e}, best 6x6 local max = {hb:+.4e}", flush=True)
    ms, hs = np.array(rows).T
    sel = (hs > 1e-6) & (ms < 0.02)
    p = np.polyfit(np.log(ms[sel]), np.log(hs[sel]), 1)
    print(f"branch fit for m < 0.02: h ~ {np.exp(p[1]):.3f} * m^{p[0]:.3f};  best 6x6 value reached: {best6:+.4e}")


def mode_search(n, trials, seed):
    rng = np.random.default_rng(seed); L = Lyap(n)
    comps = [c for k in range(2, 5) for c in itertools.product(range(1, n), repeat=k) if sum(c) == n]
    best, t0 = -9, time.time()
    for t in range(trials):
        c = comps[rng.integers(len(comps))]
        lay = np.repeat(np.arange(len(c)), c)
        perm = rng.permutation(n) if rng.random() < 0.2 else np.arange(n)
        A = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                if lay[i] == lay[j] + 1 or (lay[i] > lay[j] + 1 and rng.random() < 0.3):
                    A[i, j] = rng.standard_normal()
                if i == j and rng.random() < 0.35:
                    A[i, j] = rng.standard_normal()
        A = A[np.ix_(perm, perm)] + 1e-3 * rng.standard_normal((n, n))
        best = max(best, L.climb(A, 4000)[0])
    print(f"search n={n} seed={seed}: {trials} structured starts, best h = {best:.4e} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    print("kron identity check:", check_kron())
    mode = sys.argv[1] if len(sys.argv) > 1 else "kv"
    if mode == "kv":
        mode_kv()
    elif mode == "strip":
        mode_strip()
    else:
        mode_search(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
