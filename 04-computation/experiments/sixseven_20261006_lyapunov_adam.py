#!/usr/bin/env python3
"""Adam search for the Lyapunov symmetric-maximizer gap at orders 6 and 7 (FINITE-NUMERICAL; six-seven session).

Kressner-Vandereycken (arXiv:2608.20875, section 4) report that L-BFGS from random starts stays on the
"equality ridge" g(A) = 0, while Adam (first-moment momentum + RMS scaling) escapes it.  We reproduce their
recipe: maximise g(A) = ||L_A|Skew|| - ||L_A|Sym|| on the Frobenius unit sphere, gradient from the top singular
vectors projected to the tangent space, Adam with beta1 = 0.9, beta2 = 0.999, step 0.02 for 800 iterations,
then 0.006 with quadratic decay to 15% (800 more), renormalising after every step.  POSITIVE CONTROL: n = 7
(and 9), where counterexamples exist; then the same budget at n = 6.
Run: python3 sixseven_20261006_lyapunov_adam.py n starts seed [iterations_per_phase=800] [lr=0.02]
"""
import sys, time
import numpy as np


def bases(n):
    S, K = [], []
    for i in range(n):
        E = np.zeros((n, n)); E[i, i] = 1; S.append(E.ravel())
    for i in range(n):
        for j in range(i + 1, n):
            E = np.zeros((n, n)); E[i, j] = E[j, i] = 2 ** -.5; S.append(E.ravel())
            F = np.zeros((n, n)); F[i, j] = 2 ** -.5; F[j, i] = -2 ** -.5; K.append(F.ravel())
    return np.array(S).T, np.array(K).T


def make(n):
    BS, BK = bases(n); I = np.eye(n)

    def top(A, B):
        M = B.T @ (np.kron(I, A) + np.kron(A, I)) @ B
        U, s, Vt = np.linalg.svd(M)
        return s[0], (B @ U[:, 0]).reshape(n, n), (B @ Vt[0]).reshape(n, n)

    def g_and_grad(A):
        ss, Us, Vs = top(A, BS); sk, Uk, Vk = top(A, BK)
        return sk - ss, (Uk @ Vk.T + Uk.T @ Vk) - (Us @ Vs.T + Us.T @ Vs), ss, sk
    return g_and_grad


def adam_run(n, gg, rng, it1=800, it2=800, lr1=0.02, lr2=0.006):
    A = rng.standard_normal((n, n)); A /= np.linalg.norm(A)
    m = np.zeros_like(A); v = np.zeros_like(A); b1, b2, eps = 0.9, 0.999, 1e-12
    best, bestA = -1.0, A.copy()
    for t in range(1, it1 + it2 + 1):
        lr = lr1 if t <= it1 else lr2 * (1 - 0.85 * ((t - it1) / it2) ** 2)
        g, G, _, _ = gg(A)
        if g > best:
            best, bestA = g, A.copy()
        G = G - np.sum(G * A) * A            # tangent projection (ascent direction)
        m = b1 * m + (1 - b1) * G; v = b2 * v + (1 - b2) * G * G
        mh = m / (1 - b1 ** t); vh = v / (1 - b2 ** t)
        A = A + lr * mh / (np.sqrt(vh) + eps)
        A /= np.linalg.norm(A)
    g, _, ss, sk = gg(A)
    if g > best:
        best, bestA = g, A.copy()
    return best, bestA


if __name__ == "__main__":
    n, starts, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    it1 = int(sys.argv[4]) if len(sys.argv) > 4 else 800
    lr1 = float(sys.argv[5]) if len(sys.argv) > 5 else 0.02
    rng = np.random.default_rng(seed); gg = make(n)
    t0 = time.time(); hits = 0; best = -1; bestA = None; vals = []
    for s in range(starts):
        b, A = adam_run(n, gg, rng, it1, it1, lr1, 0.3 * lr1)
        vals.append(b)
        if b > 1e-9:
            hits += 1
        if b > best:
            best, bestA = b, A
    vals = np.array(vals)
    print(f"n={n} seed={seed} it={it1}+{it1} lr={lr1}: {starts} Adam runs, {hits} with gap > 1e-9; best gap {best:.4e} "
          f"(ratio sigma_skew/sigma_sym - 1 = {best / max(1e-300, gg(bestA)[2]):.3e}); "
          f"quantiles of best gap: {np.quantile(vals, [0.5, 0.9, 0.99])}; {time.time() - t0:.0f}s", flush=True)
    np.save(f"lyap_adam_best_n{n}_s{seed}.npy", bestA)
