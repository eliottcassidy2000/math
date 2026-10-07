import sys, time
from auditA_kappa import cay, colourable
q, k, solver = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
F, mu, E = cay(q)
n = F.q
adj = [set() for _ in range(n)]
for a, b in E: adj[a].add(b); adj[b].add(a)
u = mu[0]
common = sorted(adj[0] & adj[u])
fix = [(0, 0), (u, 1)] + ([(common[0], 2)] if common and k >= 3 else [])
t0 = time.time()
r = colourable(n, E, k, fix, solver=solver)
print(f"q={q} k={k} solver={solver}: {'SAT' if r else 'UNSAT'} ({time.time()-t0:.1f}s)", flush=True)
