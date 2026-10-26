# Corrected Hashimoto structure with multiplicities.  For a dart t of P1(F_q) let M(t) be the multiset
# [Ut, At, Dt, tau t] (tau = D o iota).  Then (K1) M(t) depends only on the vertex head(t) := D(V t) and, as a set,
# equals D(V t); (K2) B = P_U + P_A + P_D is the non-backtracking operator: B[e, .] = M(e) minus one copy of tau(e);
# (K3) A, the 4-regular vertex adjacency with multiplicities, (A g)(v) = sum_{f in M(v)} g(head(f)), satisfies the
# Ihara-Bass correspondence: every eigenvalue lambda != +-1 of B solves lambda^2 - mu lambda + 3 = 0 for an eigenvalue
# mu of A; (K4) Ramanujan: every eigenvalue mu != 4 of A has |mu| <= 2 sqrt 3.
import sys
import numpy as np
from collections import Counter
from sympy import primerange
INF = 'inf'
def inv(q, t): return 0 if t == INF else (INF if t == 0 else pow(t, -1, q))
def neg(q, t): return t if t == INF else (-t) % q
def Dm(q, t): return INF if t == INF else (t + 2) % q
def S(q, t): return neg(q, inv(q, t))
def U(q, t): return Dm(q, S(q, t))
def A_(q, t): return Dm(q, S(q, neg(q, t)))
def tau(q, t): return Dm(q, neg(q, t))
QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 120
worst = 0; allok = True
for q in primerange(5, QMAX):
    pts = list(range(q)) + [INF]; idx = {p: i for i, p in enumerate(pts)}; N = len(pts)
    M = {t: Counter([U(q, t), A_(q, t), Dm(q, t), tau(q, t)]) for t in pts}
    headkey = {t: frozenset(M[t].items()) for t in pts}
    # vertices = distinct head multisets; K1: darts with the same head set have the same multiset
    sets = {}
    for t in pts: sets.setdefault(frozenset(M[t]), set()).add(headkey[t])
    K1 = all(len(v) == 1 for v in sets.values())
    vkeys = list({headkey[t] for t in pts}); vid = {k: i for i, k in enumerate(vkeys)}
    head = {t: vid[headkey[t]] for t in pts}
    B = np.zeros((N, N)); Bnb = np.zeros((N, N))
    for t in pts:
        for g in (U, A_, Dm): B[idx[t], idx[g(q, t)]] += 1
        c = Counter(M[t]); c[tau(q, t)] -= 1
        for f, k in c.items(): Bnb[idx[t], idx[f]] += k
    K2 = np.array_equal(B, Bnb)
    nv = len(vkeys); Ad = np.zeros((nv, nv))
    for k in vkeys:
        for f, mult in dict(k).items(): Ad[vid[k], head[f]] += mult
    muA = np.linalg.eigvals(Ad)
    lamB = np.linalg.eigvals(B)
    # Ihara: for each mu, roots of l^2 - mu l + 3 should appear among eigenvalues of B (numerically; tolerance loose
    # because B is non-normal)
    roots = np.concatenate([np.roots([1, -m, 3]) for m in muA])
    def near(z, arr, tol=1e-4): return np.min(np.abs(arr - z)) < tol
    K3 = all(near(r, lamB, 1e-3) or abs(abs(r) - 1) < 1e-6 for r in roots)
    nontriv = [m for m in muA if abs(m - 4) > 1e-7]
    mx = max(abs(m) for m in nontriv)
    worst = max(worst, mx / (2 * np.sqrt(3)))
    sym = np.allclose(Ad, Ad.T)
    ok = K1 and K2 and mx <= 2 * np.sqrt(3) + 1e-9
    allok &= ok
    if q < 40 or not ok:
        print(f"q={q}: vertices {nv}, A symmetric {sym}, (K1) {K1}, (K2) B == NB {K2}, (K3) Ihara roots found {K3}, "
              f"max nontrivial |mu| = {mx:.4f} (2 sqrt3 = {2*np.sqrt(3):.4f})")
print(f"all primes 5 <= q < {QMAX}: K1, K2 and Ramanujan hold: {allok}; worst max|mu| / 2sqrt3 = {worst:.4f}")
