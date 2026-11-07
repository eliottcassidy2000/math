# mac-mini-2026-11-06-microcosm. THM-4900 (2),(3),(5): numerical tower limits (gate -> (2+sqrt(z^2+4))/z, doubling -> z phi, padding -> 1), commutation at z = 1, D(Gamma T) gate-free.
import numpy as np, itertools, random, math
def F(M, z):
    n = M.shape[0]; u = np.ones(n); return float(u @ np.linalg.solve(z * np.eye(n) - M, u))
def gate(M):
    n = M.shape[0]; G = np.zeros((n + 2, n + 2)); G[1:n + 1, 1:n + 1] = M
    G[0, 1:n + 1] = 1; G[1:n + 1, 0] = -1; G[1:n + 1, n + 1] = 1; G[n + 1, 1:n + 1] = -1; G[n + 1, 0] = 1; G[0, n + 1] = -1
    return G
def pad(M):
    G = gate(M); n = M.shape[0]; G[n + 1, 0] = -1; G[0, n + 1] = 1; return G
def double(M):
    n = M.shape[0]; I = np.eye(n); return np.block([[M, M + I], [M - I, -M]])
def rand_dom(n):
    M = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            s = 1 if random.random() < 0.5 else -1; M[i, j] = s; M[j, i] = -s
    return M
random.seed(4)
for z in (0.5, 1.0, 2.0, 3.0):
    lim = (2 + math.sqrt(z * z + 4)) / z
    rows = []
    for trial in range(3):
        M = rand_dom(random.randint(1, 5))
        for k in range(12): M = gate(M)
        rows.append(F(M, z))
    print('z = %.1f: F after 12 gates on random cores:' % z, ['%.10f' % r for r in rows], ' limit (2+sqrt(z^2+4))/z = %.10f' % lim)
# doubling tower -> z * phi(T); padding tower -> 1
for z in (0.5, 2.0):
    M0 = rand_dom(4); phi = F(M0, 1.0); M = M0.copy()
    for k in range(6): M = double(M)
    Mp = M0.copy()
    for k in range(10): Mp = pad(Mp)
    print('z = %.1f: F(D^6 T) = %.6f vs z*phi(T) = %.6f ; F(E^10 T) = %.6f' % (z, F(M, z), z * phi, F(Mp, z)))
# at z = 1: D and Gamma commute in phi; D(Gamma T) has no gate
M0 = rand_dom(5)
print('phi(D Gamma T) = %.10f, phi(Gamma D T) = %.10f, phi(Gamma Gamma D T) = %.10f' % (F(double(gate(M0)), 1), F(gate(double(M0)), 1), F(gate(gate(double(M0))), 1)))
def has_gate(M):
    n = M.shape[0]; A = (M > 0)
    for s in range(n):
        outs = [y for y in range(n) if A[s, y]]
        if len(outs) == 1:
            t = outs[0]
            if [x for x in range(n) if A[x, t]] == [s]: return True
    return False
cnt = 0
for trial in range(50):
    M0 = rand_dom(random.randint(2, 5)); cnt += has_gate(double(gate(M0)))
print('D(Gamma(T)) has a gate in', cnt, 'of 50 random cases')
