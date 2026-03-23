#!/usr/bin/env python3
"""
three_view_s20cq.py — Three-View Invariant Engine for G_n/Z_2
kind-pasteur-2026-03-23-S20cq

MANDATE: Every invariant computed in TRIPLICATE:
  (B) BLUE subgraph (SC-preserving edges)
  (K) BLACK subgraph (SC-changing edges)
  (C) COMBINED (all edges)

Then for each invariant, test relationships:
  C = B + K ?  (additive)
  C = B * K ?  (multiplicative)
  C = f(B, K)? (other formula)
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations, combinations
from collections import defaultdict, Counter, deque
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  THREE-VIEW INVARIANT ENGINE FOR G_n/Z_2")
print("  Blue | Black | Combined — always all three")
print("  kind-pasteur-2026-03-23-S20cq")
print("=" * 80)

# ============================================================================
# TOURNAMENT HELPERS (compact)
# ============================================================================

def tadj(n, bits):
    a = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): a[i][j] = 1
            else: a[j][i] = 1
            idx += 1
    return a

def Hdp(a, n):
    dp = [0]*((1<<n)*n)
    for v in range(n): dp[(1<<v)*n+v] = 1
    for S in range(1, 1<<n):
        if bin(S).count('1') >= n: continue
        for v in range(n):
            if not(S&(1<<v)): continue
            val = dp[S*n+v]
            if val == 0: continue
            for u in range(n):
                if S&(1<<u): continue
                if a[v][u]: dp[(S|(1<<u))*n+u] += val
    return sum(dp[((1<<n)-1)*n+v] for v in range(n))

def canon(a, n):
    sc = [sum(a[i][j] for j in range(n)) for i in range(n)]
    sg = defaultdict(list)
    for v in range(n): sg[sc[v]].append(v)
    gs = [sg[s] for s in sorted(set(sc))]
    if all(len(g)==1 for g in gs):
        p = [g[0] for g in gs]
        return tuple(a[p[i]][p[j]] for i in range(n) for j in range(n))
    best = None
    def gp(gs):
        if not gs: yield []; return
        for p in permutations(gs[0]):
            for r in gp(gs[1:]): yield list(p)+r
    for p in gp(gs):
        f = tuple(a[p[i]][p[j]] for i in range(n) for j in range(n))
        if best is None or f < best: best = f
    return best

def comp(a, n):
    return [[1-a[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]

def ssq(a, n):
    return tuple(sorted(sum(a[i][j] for j in range(n)) for i in range(n)))

def c3c(a, n):
    c = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if a[i][j] and a[j][k] and a[k][i]: c += 1
                if a[i][k] and a[k][j] and a[j][i]: c += 1
    return c

# ============================================================================
# BUILD G_n/Z_2 WITH BLUE/BLACK DECOMPOSITION
# ============================================================================

def build(n):
    m = comb(n,2); total = 1<<m; t0 = time.time()
    iso = defaultdict(list)
    for b in range(total):
        a = tadj(n, b); iso[canon(a, n)].append(b)

    cl = []
    c2c = {}
    for idx, (cn, mem) in enumerate(sorted(iso.items())):
        a = tadj(n, mem[0])
        cl.append({'cid':idx,'canon':cn,'mem':mem,'adj':a,
                   'H':Hdp(a,n),'sc':canon(comp(a,n),n)==cn,
                   'score':ssq(a,n),'c3':c3c(a,n),
                   'aut':factorial(n)//len(mem),'size':len(mem)})
        c2c[cn] = idx

    for d in cl:
        cc = canon(comp(d['adj'],n),n)
        d['comp_cid'] = c2c.get(cc,-1)

    edges = set(); ecol = {}
    for d in cl:
        ci = d['cid']
        for b in d['mem']:
            a2 = tadj(n, b)
            for ai in range(m):
                f = b^(1<<ai); fa = tadj(n, f); fc = canon(fa, n)
                nb = c2c.get(fc)
                if nb is not None and nb != ci:
                    e = (min(ci,nb),max(ci,nb)); edges.add(e)
                    if e not in ecol:
                        ecol[e] = 'blue' if d['sc']==cl[nb]['sc'] else 'black'

    mi = {}; mid = 0
    for d in cl:
        ci = d['cid']
        if ci in mi: continue
        cp = d['comp_cid']; mi[ci] = mid
        if cp != ci: mi[cp] = mid
        mid += 1

    Vm = mid
    me = set(); mec = {}; col_count = 0
    for (a,b) in edges:
        ma,mb = mi[a],mi[b]
        if ma != mb:
            e = (min(ma,mb),max(ma,mb)); me.add(e)
            c = ecol.get((min(a,b),max(a,b)),'blue')
            if e not in mec: mec[e] = c
            elif mec[e] != c: mec[e] = 'mixed'

    mc = []; seen = set()
    for d in cl:
        mv = mi[d['cid']]
        if mv not in seen:
            seen.add(mv)
            mc.append({'mid':mv,'H':d['H'],'sc':d['sc'],'score':d['score'],
                       'c3':d['c3'],'aut':d['aut'],'size':d['size']})
    mc.sort(key=lambda x: x['mid'])

    # Build three adjacency matrices
    Ab = np.zeros((Vm,Vm),dtype=int)  # blue
    Ak = np.zeros((Vm,Vm),dtype=int)  # black
    Ac = np.zeros((Vm,Vm),dtype=int)  # combined
    for e, c in mec.items():
        a,b = e; Ac[a][b]=1; Ac[b][a]=1
        if c in ('blue','mixed'): Ab[a][b]=1; Ab[b][a]=1
        if c in ('black','mixed'): Ak[a][b]=1; Ak[b][a]=1

    blue_e = sum(1 for c in mec.values() if c=='blue')
    black_e = sum(1 for c in mec.values() if c=='black')
    mixed_e = sum(1 for c in mec.values() if c=='mixed')
    print(f"  G_{n}/Z_2: V={Vm} E={len(me)} (blue={blue_e} black={black_e} mixed={mixed_e}) [{time.time()-t0:.1f}s]")

    return {'n':n,'V':Vm,'classes':mc,'Ab':Ab,'Ak':Ak,'Ac':Ac,
            'edges':me,'edge_colors':mec}


# ============================================================================
# THREE-VIEW INVARIANT COMPUTATION
# ============================================================================

def spectrum(A):
    return np.sort(np.linalg.eigvalsh(A.astype(float)))[::-1]

def components(A, V):
    vis = [False]*V; comps = []
    for s in range(V):
        if vis[s]: continue
        c = []; q = deque([s]); vis[s] = True
        while q:
            v = q.popleft(); c.append(v)
            for u in range(V):
                if A[v][u]>0 and not vis[u]:
                    vis[u] = True; q.append(u)
        comps.append(c)
    return comps

def count_triangles(A, V):
    A2 = A.astype(float) @ A.astype(float)
    return int(round(np.trace(A.astype(float) @ A2) / 6))

def walks(A, V, k_max=8):
    P = np.eye(V); w = [V]
    for _ in range(k_max):
        P = P @ A.astype(float)
        w.append(int(round(np.trace(P))))
    return w

def clique_number(A, V, max_size=10):
    omega = 1
    for s in range(2, min(V+1, max_size+1)):
        if comb(V,s) > 200000: break
        found = False
        for sub in combinations(range(V), s):
            if all(A[sub[i]][sub[j]]>0 for i in range(len(sub)) for j in range(i+1,len(sub))):
                omega = s; found = True; break
        if not found: break
    return omega

def independence_number(A, V, max_size=10):
    alpha = 1
    for s in range(2, min(V+1, max_size+1)):
        if comb(V,s) > 200000: break
        found = False
        for sub in combinations(range(V), s):
            if all(A[sub[i]][sub[j]]==0 for i in range(len(sub)) for j in range(i+1,len(sub))):
                alpha = s; found = True; break
        if not found: break
    return alpha

def independence_poly(A, V):
    alpha = [0]*(V+1); alpha[0] = 1
    for s in range(1, V+1):
        if comb(V,s) > 500000: break
        cnt = 0
        for sub in combinations(range(V), s):
            if all(A[sub[i]][sub[j]]==0 for i in range(len(sub)) for j in range(i+1,len(sub))):
                cnt += 1
        alpha[s] = cnt
        if cnt == 0: break
    while alpha and alpha[-1]==0: alpha.pop()
    return alpha

def forman_ricci_stats(A, V):
    deg = np.sum(A, axis=1).astype(int)
    curvs = []
    for i in range(V):
        for j in range(i+1, V):
            if A[i][j] > 0:
                t = sum(1 for k in range(V) if k!=i and k!=j and A[i][k]>0 and A[j][k]>0)
                curvs.append(4 - deg[i] - deg[j] + 3*t)
    if not curvs: return {'min':0,'max':0,'avg':0,'total':0,'count':0}
    return {'min':min(curvs),'max':max(curvs),'avg':float(np.mean(curvs)),
            'total':sum(curvs),'count':len(curvs)}


def three_view(data):
    """Compute all invariants in triplicate."""
    n = data['n']; V = data['V']
    Ab = data['Ab']; Ak = data['Ak']; Ac = data['Ac']

    views = {'Blue': Ab, 'Black': Ak, 'Combined': Ac}
    results = {}

    for name, A in views.items():
        E = int(np.sum(A)//2)
        deg = np.sum(A, axis=1).astype(int)
        r = {'E': E}

        # Degrees
        r['deg_min'] = int(np.min(deg))
        r['deg_max'] = int(np.max(deg))
        r['deg_avg'] = float(np.mean(deg))
        r['deg_sum'] = int(np.sum(deg))  # = 2E

        # Components
        comps = components(A, V)
        r['num_components'] = len(comps)
        r['largest_component'] = max(len(c) for c in comps)
        r['component_sizes'] = sorted([len(c) for c in comps], reverse=True)

        # Spectrum
        eigs = spectrum(A)
        r['spectral_radius'] = float(eigs[0])
        r['spectral_gap'] = float(eigs[0]-eigs[1]) if V > 1 else 0
        r['graph_energy'] = float(np.sum(np.abs(eigs)))

        # Laplacian
        D = np.diag(deg.astype(float))
        L = D - A.astype(float)
        leigs = np.sort(np.linalg.eigvalsh(L))
        r['alg_connectivity'] = float(leigs[1]) if V > 1 else 0
        r['lap_radius'] = float(leigs[-1])

        # Triangles
        r['triangles'] = count_triangles(A, V)

        # Walks
        w = walks(A, V, k_max=6)
        r['walks'] = w

        # Clique/independence
        r['omega'] = clique_number(A, V)
        r['alpha'] = independence_number(A, V)

        # Independence polynomial
        if V <= 40:
            ip = independence_poly(A, V)
            r['indep_poly'] = ip
            r['I_at_2'] = sum(a*(2**k) for k,a in enumerate(ip))

        # Forman curvature
        fc = forman_ricci_stats(A, V)
        r['forman_min'] = fc['min']
        r['forman_max'] = fc['max']
        r['forman_avg'] = fc['avg']
        r['forman_total'] = fc['total']

        # Kirchhoff
        nz = leigs[leigs > 1e-10]
        if len(nz) > 0:
            r['kirchhoff'] = float(V * np.sum(1.0/nz))

        # Spanning trees
        if len(nz) == V-1 and len(nz) > 0:
            log_st = float(np.sum(np.log(nz)) - np.log(V))
            r['log_spanning_trees'] = log_st
            if log_st < 50:
                r['spanning_trees'] = round(np.exp(log_st))

        # Diameter
        dist = np.full((V,V), V+1)
        for s in range(V):
            dist[s][s] = 0; q = deque([s]); vis = {s}
            while q:
                v = q.popleft()
                for u in range(V):
                    if A[v][u]>0 and u not in vis:
                        vis.add(u); dist[s][u] = dist[s][v]+1; q.append(u)
        reach = dist[dist <= V]
        r['diameter'] = int(np.max(reach)) if len(reach) > 0 else -1
        r['wiener'] = int(np.sum(dist[dist <= V])//2)
        nonzero = dist[(dist>0)&(dist<=V)]
        r['avg_distance'] = float(np.mean(nonzero)) if len(nonzero)>0 else 0

        # Clustering
        cc_sum = 0; cc_cnt = 0
        for v in range(V):
            nb = [u for u in range(V) if A[v][u]>0]
            k = len(nb)
            if k >= 2:
                t = sum(1 for i in range(len(nb)) for j in range(i+1,len(nb)) if A[nb[i]][nb[j]]>0)
                cc_sum += 2*t/(k*(k-1)); cc_cnt += 1
        r['avg_clustering'] = cc_sum/cc_cnt if cc_cnt > 0 else 0

        results[name] = r

    return results


# ============================================================================
# FORMULA DISCOVERY
# ============================================================================

def discover_formulas(results, n):
    """Test relationships between blue (B), black (K), and combined (C) invariants."""
    B = results['Blue']
    K = results['Black']
    C = results['Combined']

    print(f"\n  FORMULA DISCOVERY (n={n}):")
    print(f"  {'Invariant':<25} {'Blue':>10} {'Black':>10} {'Combined':>10} | {'B+K':>10} {'B*K':>10} | {'Relation':>15}")
    print(f"  {'-'*95}")

    tests = [
        ('E (edges)', B['E'], K['E'], C['E']),
        ('deg_sum', B['deg_sum'], K['deg_sum'], C['deg_sum']),
        ('deg_max', B['deg_max'], K['deg_max'], C['deg_max']),
        ('deg_avg', B['deg_avg'], K['deg_avg'], C['deg_avg']),
        ('triangles', B['triangles'], K['triangles'], C['triangles']),
        ('spectral_radius', B['spectral_radius'], K['spectral_radius'], C['spectral_radius']),
        ('graph_energy', B['graph_energy'], K['graph_energy'], C['graph_energy']),
        ('alg_connectivity', B['alg_connectivity'], K['alg_connectivity'], C['alg_connectivity']),
        ('omega (clique)', B['omega'], K['omega'], C['omega']),
        ('alpha (indep)', B['alpha'], K['alpha'], C['alpha']),
        ('num_components', B['num_components'], K['num_components'], C['num_components']),
        ('diameter', B['diameter'], K['diameter'], C['diameter']),
        ('wiener', B['wiener'], K['wiener'], C['wiener']),
        ('avg_clustering', B['avg_clustering'], K['avg_clustering'], C['avg_clustering']),
        ('forman_total', B['forman_total'], K['forman_total'], C['forman_total']),
    ]

    for name, b, k, c in tests:
        bk_sum = b + k
        bk_prod = b * k if isinstance(b, int) and isinstance(k, int) else b * k
        rel = ""
        if isinstance(c, int) and isinstance(bk_sum, int):
            if c == bk_sum: rel = "C = B+K"
            elif c == bk_prod: rel = "C = B*K"
            elif bk_sum != 0 and c != 0:
                ratio = c / bk_sum
                if abs(ratio - round(ratio)) < 0.001:
                    rel = f"C = {round(ratio)}*(B+K)"
        elif isinstance(c, float):
            if abs(c - (b+k)) < 0.001: rel = "C = B+K"
            elif abs(b+k) > 0.001:
                ratio = c / (b+k)
                if abs(ratio - round(ratio, 2)) < 0.01:
                    rel = f"C ~ {ratio:.2f}*(B+K)"

        if isinstance(c, float):
            print(f"  {name:<25} {b:>10.4f} {k:>10.4f} {c:>10.4f} | {bk_sum:>10.4f} {bk_prod:>10.4f} | {rel:>15}")
        else:
            print(f"  {name:<25} {b:>10} {k:>10} {c:>10} | {bk_sum:>10} {'':>10} | {rel:>15}")


# ============================================================================
# MAIN
# ============================================================================

all_data = {}

for n in [3, 4, 5, 6]:
    print(f"\n{'#'*80}")
    print(f"  n = {n}")
    print(f"{'#'*80}")

    data = build(n)
    results = three_view(data)
    all_data[n] = (data, results)

    # Print results
    for view in ['Blue', 'Black', 'Combined']:
        r = results[view]
        print(f"\n  --- {view} view ---")
        print(f"  E={r['E']}  deg=[{r['deg_min']},{r['deg_max']}]  avg_deg={r['deg_avg']:.2f}")
        print(f"  components={r['num_components']}  sizes={r['component_sizes'][:5]}")
        print(f"  rho={r['spectral_radius']:.4f}  energy={r['graph_energy']:.4f}  lam2={r['alg_connectivity']:.4f}")
        print(f"  triangles={r['triangles']}  omega={r['omega']}  alpha={r['alpha']}")
        print(f"  diameter={r['diameter']}  wiener={r['wiener']}  avg_dist={r['avg_distance']:.4f}")
        print(f"  clustering={r['avg_clustering']:.4f}  Forman_avg={r['forman_avg']:.3f}  Forman_total={r['forman_total']}")
        if 'I_at_2' in r:
            print(f"  I(2)={r['I_at_2']}  I(x)={r.get('indep_poly','?')}")
        if 'spanning_trees' in r:
            print(f"  spanning_trees={r['spanning_trees']}")
        if 'kirchhoff' in r:
            print(f"  Kirchhoff={r['kirchhoff']:.4f}")

    # Formula discovery
    discover_formulas(results, n)


# ============================================================================
# MASTER SEQUENCE TABLE
# ============================================================================

print(f"\n\n{'='*80}")
print("  MASTER THREE-VIEW SEQUENCE TABLE")
print(f"{'='*80}")

inv_names = ['E', 'triangles', 'omega', 'alpha', 'diameter', 'num_components',
             'spectral_radius', 'graph_energy', 'alg_connectivity', 'wiener',
             'forman_total', 'avg_clustering']

for inv in inv_names:
    print(f"\n  {inv}:")
    for view in ['Blue', 'Black', 'Combined']:
        vals = []
        for n in [3,4,5,6]:
            v = all_data[n][1][view].get(inv, '?')
            if isinstance(v, float):
                vals.append(f"{v:.3f}")
            else:
                vals.append(str(v))
        print(f"    {view:>10}: {', '.join(vals)}")

# Independence polynomial sequences
print(f"\n  I(2) (meta-H):")
for view in ['Blue', 'Black', 'Combined']:
    vals = []
    for n in [3,4,5,6]:
        v = all_data[n][1][view].get('I_at_2', '?')
        vals.append(str(v))
    print(f"    {view:>10}: {', '.join(vals)}")

# Walk sequences
print(f"\n  Closed walks w(k) for k=0..6:")
for view in ['Blue', 'Black', 'Combined']:
    for n in [3,4,5,6]:
        w = all_data[n][1][view].get('walks', [])
        print(f"    {view:>10} n={n}: {w}")


# ============================================================================
# CROSS-VIEW RELATIONSHIPS
# ============================================================================

print(f"\n\n{'='*80}")
print("  UNIVERSAL CROSS-VIEW RELATIONSHIPS")
print(f"{'='*80}")

print(f"\n  Testing: C = B + K ? (edge additivity)")
for n in [3,4,5,6]:
    B = all_data[n][1]['Blue']['E']
    K = all_data[n][1]['Black']['E']
    C = all_data[n][1]['Combined']['E']
    match = "YES" if C == B+K else f"NO (B+K={B+K}, C={C})"
    print(f"    n={n}: B={B}, K={K}, C={C} -> C=B+K? {match}")

print(f"\n  Testing: C = B + K ? (triangle additivity)")
for n in [3,4,5,6]:
    B = all_data[n][1]['Blue']['triangles']
    K = all_data[n][1]['Black']['triangles']
    C = all_data[n][1]['Combined']['triangles']
    mixed = C - B - K
    print(f"    n={n}: B={B}, K={K}, C={C}, mixed_tri={mixed}")

print(f"\n  Testing: deg_combined = deg_blue + deg_black ? (vertex by vertex)")
for n in [3,4,5,6]:
    d = all_data[n][0]
    V = d['V']
    deg_b = np.sum(d['Ab'], axis=1)
    deg_k = np.sum(d['Ak'], axis=1)
    deg_c = np.sum(d['Ac'], axis=1)
    match = np.all(deg_c == deg_b + deg_k)
    print(f"    n={n}: deg_C = deg_B + deg_K for ALL vertices? {match}")

print(f"\n  Testing: A_combined = A_blue + A_black ? (matrix additivity)")
for n in [3,4,5,6]:
    d = all_data[n][0]
    match = np.all(d['Ac'] == d['Ab'] + d['Ak'])
    print(f"    n={n}: A_C = A_B + A_K? {match}")

print(f"\n  Spectral decomposition: spec(C) vs spec(B) + spec(K)?")
for n in [3,4,5,6]:
    d = all_data[n][0]
    eB = np.sort(np.linalg.eigvalsh(d['Ab'].astype(float)))[::-1]
    eK = np.sort(np.linalg.eigvalsh(d['Ak'].astype(float)))[::-1]
    eC = np.sort(np.linalg.eigvalsh(d['Ac'].astype(float)))[::-1]
    # Sum of eigenvalues (trace)
    tr_B = np.sum(eB)
    tr_K = np.sum(eK)
    tr_C = np.sum(eC)
    # Energy
    en_B = np.sum(np.abs(eB))
    en_K = np.sum(np.abs(eK))
    en_C = np.sum(np.abs(eC))
    print(f"    n={n}: tr: B={tr_B:.1f} K={tr_K:.1f} C={tr_C:.1f} (all ~0)")
    print(f"           energy: B={en_B:.2f} K={en_K:.2f} C={en_C:.2f} B+K={en_B+en_K:.2f} "
          f"{'SUPER' if en_C < en_B+en_K else 'SUB'}-additive")

print(f"\n  Blue/Black edge ratio:")
for n in [3,4,5,6]:
    B = all_data[n][1]['Blue']['E']
    K = all_data[n][1]['Black']['E']
    ratio = B/K if K > 0 else float('inf')
    print(f"    n={n}: B/K = {B}/{K} = {ratio:.3f}")

print(f"\n  Kirchhoff index relationship:")
for n in [3,4,5,6]:
    kB = all_data[n][1]['Blue'].get('kirchhoff', '?')
    kK = all_data[n][1]['Black'].get('kirchhoff', '?')
    kC = all_data[n][1]['Combined'].get('kirchhoff', '?')
    print(f"    n={n}: B={kB}  K={kK}  C={kC}")

print(f"\n\n  DONE.")
print("=" * 80)
