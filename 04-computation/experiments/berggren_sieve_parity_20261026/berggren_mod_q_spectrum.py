# The Berggren tree modulo a prime q: the three branches act on the projective line P1(F_q) (slopes t = m/n),
# U: t -> 2 - 1/t,  A: t -> 2 + 1/t,  D: t -> t + 2.  Two operators on functions f on P1(F_q):
#   T f(x)   = (f(Ux) + f(Ax) + f(Dx)) / 3          (uniform random branch: equidistribution of residues by depth),
#   T_chi f  = (f(Ux) - f(Ax) + f(Dx)) / 3          (signed by the branch character).
# Reported: the second-largest eigenvalue modulus of T (a uniform gap = the engine of a growing-modulus sieve on
# the tree), and the spectral radius of T_chi (the rate at which the branch parity decorrelates from the residue
# class: the signed mass of depth-k nodes over a fixed residue class decays like rho_chi^k, relative to 1).
import numpy as np
from sympy import primerange
INF = 'inf'
def act(q, t, g):
    if g == 'D': return INF if t == INF else (t + 2) % q
    if t == INF: return 2 % q                          # 2 -+ 1/inf = 2
    if t == 0: return INF
    inv = pow(t, -1, q)
    return (2 - inv) % q if g == 'U' else (2 + inv) % q
rows = []
for q in primerange(3, 400):
    pts = list(range(q)) + [INF]; idx = {p: i for i, p in enumerate(pts)}; N = len(pts)
    T = np.zeros((N, N)); Tc = np.zeros((N, N))
    for p in pts:
        for g, s in (('U', 1), ('A', -1), ('D', 1)):
            T[idx[p], idx[act(q, p, g)]] += 1 / 3; Tc[idx[p], idx[act(q, p, g)]] += s / 3
    ev = np.sort(np.abs(np.linalg.eigvals(T)))[::-1]
    evc = np.sort(np.abs(np.linalg.eigvals(Tc)))[::-1]
    # norm of T on mean-zero functions (singular value), a stronger mixing measure for non-normal T
    P = np.eye(N) - np.ones((N, N)) / N
    sv = np.linalg.svd(P @ T @ P, compute_uv=False)[0]
    rows.append((q, ev[1], sv, evc[0]))
print("q  |lambda_2(T)|  ||T restricted to mean-zero||  rho(T_chi)")
for r in rows:
    if r[0] < 40 or r[0] % 50 in (1, 3, 7, 9) and r[0] > 100 or r[0] > 380:
        print(f"{r[0]:4d}  {r[1]:.4f}  {r[2]:.4f}  {r[3]:.4f}")
print(f"max over 3 <= q < 400: |lambda_2| {max(r[1] for r in rows):.4f}, mean-zero norm {max(r[2] for r in rows):.4f}, "
      f"rho(T_chi) {max(r[3] for r in rows):.4f}")
