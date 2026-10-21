#!/usr/bin/env python3
"""Finite-horizon increments of the Nyquist terminal energy across one 64-scale step, against the exact new-annulus source
budget 153N (native_cancellation_105_20261018g.md, equations (3)-(4)). Imports the inherited checker unchanged: source rows
[s, s%15, s%35, terminal u, clock t, depth d, a, b, norm]; signed vector v_s = (-1)^t (a,b); W_u = sum of v_s over sources with
terminal u; E = sum_u Q(W_u) (Eisenstein norm). For a P-depth horizon D only sources with route depth <= D are counted.
E(N,D) uses sources s < 4N (s odd, s != 5 mod 8); E(64N,D) uses s < 256N. Reports, for each D: E(N,D), E(64N,D), the
increment, its decomposition sum_u [Q(H_u) + cross(W_u,H_u)] (H = new sources 4N <= s < 256N), and c(D) = increment/(153N)."""
import sys, json, time
from collections import defaultdict
sys.path.insert(0, '/Users/e/Documents/GitHub/math-wt-mersenne-561/04-computation/experiments')
import native_cancellation_105_20261018g as NC
inh = NC.inherited
ranks = {int(l.split()[0]): int(l.split()[1]) for l in NC.FRONTIER.read_text().splitlines() if l.strip() and not l.startswith('#')}
routes = NC.Routes(frozenset(ranks))
def rows_for(n, routes):
    return [NC.source_row(s, n, routes) for s in range(1, 4 * n, 2) if s % 8 != 5]
def energy_by_depth(rows, maxd):
    # W_u(D) for each D: accumulate sources sorted by depth
    by_depth = defaultdict(list)
    for r in rows: by_depth[r[5]].append(r)
    W = defaultdict(lambda: (0, 0)); out = {}
    for D in range(0, maxd + 1):
        for r in by_depth.get(D, []):
            W[r[3]] = NC.add(W[r[3]], NC.signed(r))
        out[D] = sum(inh.norm(v) for v in W.values() if v != (0, 0))
    return out
for N in [int(a) for a in sys.argv[1:]] or [4096]:
    t0 = time.time()
    old = rows_for(N, routes); print(f"N = {N}: {len(old)} original rows, max depth {max(r[5] for r in old)}, {time.time()-t0:.0f}s", flush=True)
    new_scale = rows_for(64 * N, routes); print(f"64N = {64*N}: {len(new_scale)} rows, max depth {max(r[5] for r in new_scale)}, {time.time()-t0:.0f}s", flush=True)
    # transport check: original sources keep their charges at 64N (equation (2))
    bys = {r[0]: r for r in new_scale}
    assert all((bys[r[0]][6], bys[r[0]][7]) == (r[6], r[7]) and bys[r[0]][3] == r[3] and bys[r[0]][4] == r[4] for r in old), 'scale-64 transport'
    newsrc = [r for r in new_scale if r[0] >= 4 * N]
    S_new = sum(r[8] for r in newsrc); print(f"new-annulus source norm {S_new} = 153N = {153*N}: {S_new == 153*N}; new sources {len(newsrc)}", flush=True)
    maxd = max(max(r[5] for r in old), max(r[5] for r in new_scale))
    Eold = energy_by_depth(old, maxd); Enew = energy_by_depth(new_scale, maxd)
    # decomposition at each horizon
    byd_old = defaultdict(list); byd_new = defaultdict(list)
    for r in old: byd_old[r[5]].append(r)
    for r in newsrc: byd_new[r[5]].append(r)
    W = defaultdict(lambda: (0, 0)); H = defaultdict(lambda: (0, 0)); table = []
    worst = (0, None)
    for D in range(0, maxd + 1):
        for r in byd_old.get(D, []): W[r[3]] = NC.add(W[r[3]], NC.signed(r))
        for r in byd_new.get(D, []): H[r[3]] = NC.add(H[r[3]], NC.signed(r))
        qH = sum(inh.norm(v) for v in H.values()); cr = sum(NC.cross(W[u], H[u]) for u in H if W[u] != (0, 0))
        inc = Enew[D] - Eold[D]; assert inc == qH + cr, (D, inc, qH, cr)
        c = inc / (153 * N)
        table.append((D, Eold[D], Enew[D], inc, qH, cr, round(c, 4)))
        if c > worst[0]: worst = (c, D)
    print("D | E(N,D) | E(64N,D) | increment | sum Q(H_u) | sum cross(W_u,H_u) | c = increment/153N")
    for row in table:
        if row[0] <= 12 or row[0] % 25 == 0 or row[0] == maxd: print('  ', row)
    print(f"max c over horizons: {worst[0]:.4f} at D = {worst[1]}; full arrival c = {table[-1][6]}; E(N)/S(N) = {Eold[maxd]/sum(r[8] for r in old):.4f}, E(64N)/S(64N) = {Enew[maxd]/sum(r[8] for r in new_scale):.4f}", flush=True)
    json.dump({'N': N, 'rows': table, 'max_c': worst}, open(f'annulus_horizon_N{N}.json', 'w'))
