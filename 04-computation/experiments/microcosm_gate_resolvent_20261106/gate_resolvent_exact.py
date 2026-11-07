# mac-mini-2026-11-06-microcosm. THM-4900 (1),(3),(4): exact (sympy) marked-resolvent identities for the gate, padding and skew-Sylvester doubling; phi = F(1) invariants; Fibonacci ratios of the one-vertex tower.
import itertools, random
import sympy as sp
z = sp.symbols('z')
def dom(n, arc):
    return sp.Matrix(n, n, lambda i, j: 0 if i == j else (1 if arc[i][j] else -1))
def F_of(M):
    n = M.shape[0]; u = sp.ones(n, 1)
    R = (z * sp.eye(n) - M).inv()
    return sp.cancel((u.T * R * u)[0, 0])
def gate(M):   # order: b (source-like), core, a (sink-like); b -> core -> a -> b
    n = M.shape[0]; G = sp.zeros(n + 2, n + 2)
    G[1:n + 1, 1:n + 1] = M
    for i in range(1, n + 1): G[0, i] = 1; G[i, 0] = -1; G[i, n + 1] = 1; G[n + 1, i] = -1
    G[n + 1, 0] = 1; G[0, n + 1] = -1
    return G
def pad(M):    # E: s -> core -> t, s -> t
    G = gate(M); n = M.shape[0]; G[n + 1, 0] = -1; G[0, n + 1] = 1; return G
def double(M):
    n = M.shape[0]; I = sp.eye(n)
    return sp.Matrix(sp.BlockMatrix([[M, M + I], [M - I, -M]]))
random.seed(3); ok = True
for n in range(1, 5):
    for trial in range(4):
        arc = [[False] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                if random.random() < 0.5: arc[i][j] = True
                else: arc[j][i] = True
        M = dom(n, arc); F = F_of(M); p = M.charpoly(z).as_expr()
        G = gate(M); FG = F_of(G); pG = G.charpoly(z).as_expr()
        ok &= sp.simplify(FG - ((z**2 + 9) * F + 2 * z) / (z**2 + 1 + 2 * z * F)) == 0
        ok &= sp.simplify(pG - p * (z**2 + 1 + 2 * z * F)) == 0
        E = pad(M); ok &= sp.simplify(E.charpoly(z).as_expr() - pG) == 0
        D = double(M); FD1 = F_of(D).subs(z, 1)
        ok &= sp.simplify(FD1 - F.subs(z, 1)) == 0                                   # doubling fixes phi = F(1)
        ok &= sp.simplify(F_of(E).subs(z, 1) - 1) == 0                                 # padding sends phi to 1
        ok &= sp.simplify(FG.subs(z, 1) - (5 * F.subs(z, 1) + 1) / (F.subs(z, 1) + 1)) == 0   # golden Mobius map
print('Gamma marked resolvent, p_Gamma = p_E, phi(D T) = phi(T), phi(E T) = 1, phi(Gamma T) = (5 phi + 1)/(phi + 1):', ok)
# nested gate tower over one vertex: Fibonacci ratios
M = sp.zeros(1, 1); fib = [0, 1]
for _ in range(40): fib.append(fib[-1] + fib[-2])
F_ = lambda k: fib[k] if k >= 0 else (-1) ** (k + 1) * fib[-k]
vals = []
for k in range(0, 6):
    v = F_of(M).subs(z, 1); vals.append(v)
    ok_k = v == sp.Rational(F_(2 * k + 2), F_(2 * k - 1))
    print('k = %d gates: phi = %s ; Fib(2k+2)/Fib(2k-1) = %s/%s : %s' % (k, v, F_(2 * k + 2), F_(2 * k - 1), ok_k))
    M = gate(M)
print('phi^3 = 2 + sqrt5 =', float(2 + sp.sqrt(5)))
