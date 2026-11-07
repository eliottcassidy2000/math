# mac-mini-2026-11-06-microcosm. THM-4900 (2): closed-form spectrum of the nested-gate tower over one vertex and closed-form characteristic polynomial for any core.
import numpy as np, math, random
def gate(M):
    n = M.shape[0]; G = np.zeros((n + 2, n + 2)); G[1:n + 1, 1:n + 1] = M
    G[0, 1:n + 1] = 1; G[1:n + 1, 0] = -1; G[1:n + 1, n + 1] = 1; G[n + 1, 1:n + 1] = -1; G[n + 1, 0] = 1; G[0, n + 1] = -1
    return G
ok = True
for k in range(1, 11):
    M = np.zeros((1, 1))
    for _ in range(k): M = gate(M)
    lam = np.sort(np.abs(np.linalg.eigvals(M).imag))
    pred = [0.0, math.sqrt(3), math.sqrt(3)]
    for j in range(k - 1):
        t = (2 * j + 1) * math.pi / (2 * (2 * k - 1))
        v = math.sqrt(4 + 1 / math.tan(t) ** 2); pred += [v, v]
    pred = np.sort(np.array(pred))
    err = np.max(np.abs(lam - pred))
    ok &= err < 1e-7
    if k in (2, 3, 6, 10): print('k = %2d: |eigenvalues| = %s ; predicted max error %.1e' % (k, np.round(lam[::2], 4), err))
print('closed-form spectrum of the nested-gate tower Gamma^k(K1), k = 1..10:', ok)
# general core: P_{Gamma^k(T)} = [(s+1)^(2k) Psi+(T) - (s-1)^(2k) Psi-(T)]/(2s), Psi = z N + (-2 +- s) P
random.seed(7); ok2 = True
for trial in range(5):
    n = random.randint(1, 5); M = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            sgn = 1 if random.random() < 0.5 else -1; M[i, j] = sgn; M[j, i] = -sgn
    for z in (0.7, 1.3, 2.9):
        s = math.sqrt(z * z + 4); u = np.ones(n)
        P = np.linalg.det(z * np.eye(n) - M); N = P * (u @ np.linalg.solve(z * np.eye(n) - M, u))
        Pp, Pm = z * N + (-2 + s) * P, z * N + (-2 - s) * P
        G = M.copy()
        for k in range(1, 5):
            G = gate(G)
            lhs = np.linalg.det(z * np.eye(G.shape[0]) - G)
            rhs = ((s + 1) ** (2 * k) * Pp - (s - 1) ** (2 * k) * Pm) / (2 * s)
            ok2 &= abs(lhs - rhs) <= 1e-8 * max(1, abs(lhs))
print('P(Gamma^k T) = [(s+1)^2k Psi+ - (s-1)^2k Psi-]/(2s) for random cores, k <= 4:', ok2)
