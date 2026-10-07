"""Audit A: kappa(q) = chi(Cay(F_{q^2}, mu_{q+1})) by SAT (glucose4 / cadical195), independent encoding.
Also: knight torus G_7 = Cay(F_49, mu_8) check by explicit isomorphism; mod-7 reduction claim of note §7.1."""
import sys, time, itertools
from auditA_gf import GF
from pysat.solvers import Solver

def prime_power(q):
    for p in range(2, q + 1):
        if q % p == 0:
            m, r = 0, q
            while r % p == 0: r //= p; m += 1
            assert r == 1
            return p, m

def cay(q):
    p, m = prime_power(q)
    F = GF(p, 2 * m)
    mu = [x for x in range(1, F.q) if F.pw(x, q + 1) == 1]
    assert len(mu) == q + 1
    E = set()
    for x in range(F.q):
        for u in mu:
            y = F.add[x][u]
            E.add((min(x, y), max(x, y)))
    return F, mu, sorted(E)

def colourable(n, E, k, fix, solver='glucose4', budget=None):
    v = lambda i, c: i * k + c + 1
    s = Solver(name=solver)
    for i in range(n):
        s.add_clause([v(i, c) for c in range(k)])
    for a, b in E:
        for c in range(k):
            s.add_clause([-v(a, c), -v(b, c)])
    for i, c in fix:
        s.add_clause([v(i, c)])
    if budget:
        s.conf_budget(budget)
        r = s.solve_limited()
    else:
        r = s.solve()
    if r:
        mdl = set(l for l in s.get_model() if l > 0)
        col = [next(c for c in range(k) if v(i, c) in mdl) for i in range(n)]
        assert all(col[a] != col[b] for a, b in E)
    s.delete()
    return r

if __name__ == "__main__":
  qs = [int(a) for a in sys.argv[1:]] or [2, 3, 4, 5, 7, 8, 9, 11, 16]
  for q in qs:
      F, mu, E = cay(q)
      n = F.q
      adj = [set() for _ in range(n)]
      for a, b in E: adj[a].add(b); adj[b].add(a)
      # symmetry breaking: vertex 0 colour 0, a neighbour colour 1, a common neighbour (if any) colour 2
      u = mu[0]
      common = sorted(adj[0] & adj[u])
      t0 = time.time()
      res = {}
      for k in range(2, 9):
          fix = [(0, 0), (u, 1)] + ([(common[0], 2)] if common and k >= 3 else [])
          if common and k < 3:
              res[k] = False; continue
          r = colourable(n, E, k, fix, solver='glucose4' if q <= 11 else 'cadical195', budget=(None if q <= 11 else 3_000_000))
          res[k] = r
          if r: break
      print(f"q={q}: |V|={n}, deg={len(mu)}, triangle={'yes' if common else 'no'}, results {res}, {time.time()-t0:.1f}s", flush=True)
