"""Randomized-restart witness hunt for gamma=1/2 D0=0, M=11, 9-bias envelope.
Reuses laneC2 structures; random value ordering per restart + greedy-guided
first choice. Feasible verdict only (no exhaustion claim)."""
import importlib.util, sys, random
from math import comb

spec = importlib.util.spec_from_file_location(
    "laneC2", __file__.replace("amm12592_m11_witness_hunt", "amm12592_laneC2_finiteM_dfs"))
mod = importlib.util.module_from_spec(spec)
sys.modules["laneC2"] = mod
spec.loader.exec_module(mod)

BIASES = [(1, 2), (1, 3), (2, 3), (2, 5), (3, 5), (1, 4), (3, 4),
          (1285, 2181), (8847357, 11821757)]
g1, g2, D0, Mtry = 1, 2, 0, 11
d = lambda m: (g1 * m) // g2 + D0
A = lambda m: m + d(m) + 1
Lmax = A(Mtry)
nP = len(BIASES)

rows = []
for m in range(1, Mtry + 1):
    dm = d(m)
    cells = []
    for k in range(dm + 1):
        cap = comb(dm, k)
        for (z, o) in ((m + dm - k, k + 1), (k + 1, m + dm - k)):
            w = tuple(n**z * (dd - n)**o * dd**(Lmax - z - o) for (n, dd) in BIASES)
            cells.append((cap, cap % 2, w))
    rows.append(cells)

R = [tuple(sum(c[0] * c[2][i] for c in rows[m]) for i in range(nP)) for m in range(Mtry)]
Env = [tuple((n**(m + 2) + (dd - n)**(m + 2)) * dd**(Lmax - m - 2) for (n, dd) in BIASES)
       for m in range(Mtry)]
allowed = []
for m in range(Mtry):
    best = None
    for mp in range(m, Mtry):
        tot = tuple(Env[mp][i] + sum(R[r][i] for r in range(m + 1, mp + 1)) for i in range(nP))
        best = tot if best is None else tuple(min(a, b) for a, b in zip(best, tot))
    allowed.append(best)
suffix = []
for m in range(Mtry):
    cells = rows[m]
    sw = [tuple(0 for _ in range(nP))]
    for c in reversed(cells):
        sw.append(tuple(sw[-1][i] + c[0] * c[2][i] for i in range(nP)))
    sw.reverse()
    suffix.append(sw)

sys.setrecursionlimit(100000)

def hunt(seed, node_budget):
    rng = random.Random(seed)
    nodes = [0]
    def dfs(m, ci, S, choice_list):
        nodes[0] += 1
        if nodes[0] > node_budget:
            raise TimeoutError
        cells = rows[m]
        if ci == len(cells):
            for i in range(nP):
                if abs(S[i]) > Env[m][i] or abs(S[i]) > allowed[m][i]:
                    return None
            if m + 1 == Mtry:
                return list(choice_list)
            return dfs(m + 1, 0, S, choice_list)
        rest = suffix[m][ci]
        for i in range(nP):
            if abs(S[i]) > rest[i] + allowed[m][i]:
                return None
        cap, par, w = cells[ci]
        vals = [v for v in range(-cap, cap + 1) if (v - par) % 2 == 0]
        # bias-guided soft ordering with noise
        vals.sort(key=lambda v: abs(S[0] + v * w[0]) * (1 + rng.random()))
        for v in vals:
            S2 = tuple(S[i] + v * w[i] for i in range(nP)) if v else S
            choice_list.append(v)
            r = dfs(m, ci + 1, S2, choice_list)
            if r is not None:
                return r
            choice_list.pop()
        return None
    try:
        return dfs(0, 0, tuple(0 for _ in range(nP)), [])
    except TimeoutError:
        return None

for seed in range(60):
    w = hunt(seed, 3_000_000)
    print(f"seed {seed}: {'WITNESS ' + str(w) if w else 'miss'}", flush=True)
    if w:
        # independent re-verify
        S = [0]*nP; idx = 0
        for m in range(Mtry):
            for (cap, par, wv) in rows[m]:
                v = w[idx]; idx += 1
                assert abs(v) <= cap and (v - par) % 2 == 0
                for i in range(nP):
                    S[i] += v * wv[i]
            for i in range(nP):
                assert abs(S[i]) <= Env[m][i], (m, i)
        print("WITNESS RE-VERIFIED EXACTLY: gamma=1/2 D0=0 M=11 FEASIBLE", flush=True)
        break
else:
    print("NO WITNESS in 60 restarts x 3M nodes (not an infeasibility proof)", flush=True)
