#!/usr/bin/env python3
"""
gn_merged_metrics_s218.py — Comprehensive metrics for G_n and G_n/Z_2
opus-2026-03-23-S218

Compute EVERY conceivable metric for both the original and merged metagraphs
at n=3..8, hunting for new OEIS-worthy sequences and patterns.
Uses nauty for class enumeration at n>=7.
"""

import sys
import subprocess
from math import comb, factorial, gcd, log2
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  COMPREHENSIVE METRICS FOR G_n AND G_n/Z_2")
print("  opus-2026-03-23-S218")
print("=" * 80)

# ============================================================================
# HELPERS
# ============================================================================

def get_tournaments(n):
    """Get tournament class representatives via nauty (n>=3) or brute force."""
    m = comb(n, 2)
    PAIRS = [(i,j) for i in range(n) for j in range(i+1,n)]

    if n <= 6:
        # Brute force canonical form
        from itertools import permutations as perms
        ALL_PERMS = list(perms(range(n)))

        def bits_to_adj(bits):
            adj = [[0]*n for _ in range(n)]
            for k,(i,j) in enumerate(PAIRS):
                if bits & (1<<k): adj[i][j]=1
                else: adj[j][i]=1
            return adj

        def canonical(adj):
            best = None
            for perm in ALL_PERMS:
                form = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
                if best is None or form < best: best = form
            return best

        classes = {}
        for bits in range(1<<m):
            adj = bits_to_adj(bits)
            c = canonical(adj)
            if c not in classes:
                classes[c] = {'bits': bits, 'adj': adj, 'canon': c, 'size': 0}
            classes[c]['size'] += 1

        result = list(classes.values())
        for i, r in enumerate(result):
            r['cid'] = i
        return result, PAIRS
    else:
        # Use nauty
        res = subprocess.run(['gentourng', str(n)], capture_output=True, text=True)
        lines = [l.strip() for l in (res.stdout or '').split('\n')
                 if len(l.strip())==m and all(c in '01' for c in l.strip())]

        def bits_to_adj_nauty(bits):
            adj = [[0]*n for _ in range(n)]
            for k,(i,j) in enumerate(PAIRS):
                if bits & (1<<(m-1-k)): adj[i][j]=1
                else: adj[j][i]=1
            return adj

        result = []
        for i, line in enumerate(lines):
            bits = int(line, 2)
            adj = bits_to_adj_nauty(bits)
            result.append({'bits': bits, 'adj': adj, 'cid': i, 'nauty': True})
        return result, PAIRS

def score_seq(adj, n):
    return tuple(sorted(sum(adj[i][j] for j in range(n)) for i in range(n)))

def c3_count(adj, n):
    c3 = 0
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                if adj[i][j] and adj[j][k] and adj[k][i]: c3+=1
                if adj[i][k] and adj[k][j] and adj[j][i]: c3+=1
    return c3

def H_dp(adj, n):
    dp = {}
    for v in range(n): dp[(1<<v,v)] = 1
    for S in range(1, 1<<n):
        for v in range(n):
            if not(S&(1<<v)): continue
            val = dp.get((S,v),0)
            if val==0: continue
            for u in range(n):
                if S&(1<<u): continue
                if adj[v][u]:
                    key=(S|(1<<u),u)
                    dp[key] = dp.get(key,0) + val
    return sum(dp.get(((1<<n)-1,v),0) for v in range(n))

def canonical_form_pruned(adj, n):
    """Degree-pruned canonical form."""
    scores = [sum(adj[i][j] for j in range(n)) for i in range(n)]
    score_groups = defaultdict(list)
    for v in range(n): score_groups[scores[v]].append(v)
    groups = [score_groups[s] for s in sorted(set(scores))]
    if all(len(g)==1 for g in groups):
        perm = [g[0] for g in groups]
        return tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
    best = None
    def gen(gs):
        if not gs: yield []; return
        for p in permutations(gs[0]):
            for rest in gen(gs[1:]): yield list(p)+rest
    for perm in gen(groups):
        form = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best: best = form
    return best

def complement_adj(adj, n):
    return [[1-adj[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]

def aut_count(adj, n):
    count = 0
    for perm in permutations(range(n)):
        ok = True
        for i in range(n):
            for j in range(n):
                if adj[i][j] != adj[perm[i]][perm[j]]:
                    ok = False; break
            if not ok: break
        if ok: count += 1
    return count

# ============================================================================
# MAIN LOOP OVER n
# ============================================================================

all_results = {}

for n in range(3, 9):
    m = comb(n, 2)
    nfact = factorial(n)

    print(f"\n{'='*80}")
    print(f"  n = {n}")
    print(f"{'='*80}")
    t0 = time.time()

    # Get class representatives
    classes, PAIRS = get_tournaments(n)
    num_classes = len(classes)
    print(f"  {num_classes} iso classes loaded ({time.time()-t0:.1f}s)")

    # Compute invariants for each class
    print(f"  Computing invariants...")
    t1 = time.time()

    for ci, cl in enumerate(classes):
        adj = cl['adj']
        cl['score'] = score_seq(adj, n)
        cl['c3'] = c3_count(adj, n)
        cl['H'] = H_dp(adj, n)
        if n <= 7:
            cl['aut'] = aut_count(adj, n) if n <= 7 else 1
        else:
            cl['aut'] = nfact // (cl.get('size', nfact))  # estimate from class size

        # SC status
        comp = complement_adj(adj, n)
        comp_canon = canonical_form_pruned(comp, n)
        own_canon = canonical_form_pruned(adj, n)
        cl['sc'] = (comp_canon == own_canon)
        cl['comp_canon'] = comp_canon
        cl['own_canon'] = own_canon

        if (ci+1) % 2000 == 0:
            print(f"    {ci+1}/{num_classes} ({time.time()-t1:.1f}s)")

    # Build hash lookup
    hash_to_cids = defaultdict(list)
    canon_to_cid = {}
    for cl in classes:
        hash_to_cids[(cl['score'], cl['c3'], cl['H'])].append(cl['cid'])
        canon_to_cid[cl['own_canon']] = cl['cid']

    # Find complement pairs
    for cl in classes:
        comp_cid = canon_to_cid.get(cl['comp_canon'], -1)
        cl['comp_cid'] = comp_cid

    print(f"  Invariants done ({time.time()-t1:.1f}s)")

    # Compute edges
    print(f"  Computing edges...")
    t2 = time.time()

    def lookup(adj):
        ss = score_seq(adj, n)
        c3 = c3_count(adj, n)
        h = H_dp(adj, n)
        cids = hash_to_cids.get((ss, c3, h))
        if cids is None: return -1
        if len(cids) == 1: return cids[0]
        canon = canonical_form_pruned(adj, n)
        return canon_to_cid.get(canon, -1)

    edges = set()
    self_loops = set()
    adj_list = defaultdict(set)  # for graph analysis
    degrees = {}
    edge_weights = Counter()

    for ci, cl in enumerate(classes):
        bits = cl['bits']
        nbrs = set()
        for arc in range(m):
            if n <= 6:
                fb = bits ^ (1 << arc)
                fa_adj = [[0]*n for _ in range(n)]
                for k, (i,j) in enumerate(PAIRS):
                    if fb & (1<<k): fa_adj[i][j]=1
                    else: fa_adj[j][i]=1
            else:
                fb = bits ^ (1 << (m-1-arc))
                fa_adj = [[0]*n for _ in range(n)]
                for k, (i,j) in enumerate(PAIRS):
                    if fb & (1<<(m-1-k)): fa_adj[i][j]=1
                    else: fa_adj[j][i]=1

            nb = lookup(fa_adj)
            if nb == cl['cid']:
                self_loops.add(cl['cid'])
            elif nb >= 0:
                nbrs.add(nb)
                edge = (min(cl['cid'], nb), max(cl['cid'], nb))
                edge_weights[edge] += 1

        for nb in nbrs:
            edge = (min(cl['cid'], nb), max(cl['cid'], nb))
            edges.add(edge)
            adj_list[cl['cid']].add(nb)
            adj_list[nb].add(cl['cid'])
        degrees[cl['cid']] = len(nbrs)

        if (ci+1) % 2000 == 0:
            elapsed = time.time()-t2
            print(f"    {ci+1}/{num_classes} ({elapsed:.1f}s, {len(edges)} edges)")

    print(f"  Edges done: {len(edges)} ({time.time()-t2:.1f}s)")

    # ================================================================
    # COMPUTE ALL METRICS
    # ================================================================

    sc_ids = [cl['cid'] for cl in classes if cl['sc']]
    ns_ids = [cl['cid'] for cl in classes if not cl['sc']]

    # Basic
    V = num_classes
    E = len(edges)
    density = 2*E/(V*(V-1)) if V > 1 else 0
    deg_vals = list(degrees.values())
    avg_deg = sum(deg_vals)/V if V > 0 else 0

    # Blue/black
    blue_e = sum(1 for (a,b) in edges if classes[a]['sc'] == classes[b]['sc'])
    black_e = E - blue_e
    scsc = sum(1 for (a,b) in edges if classes[a]['sc'] and classes[b]['sc'])
    nsns = sum(1 for (a,b) in edges if not classes[a]['sc'] and not classes[b]['sc'])
    scns = black_e

    # Merged graph
    merged_node = {}
    mid = 0
    for cl in classes:
        ci = cl['cid']
        if ci in merged_node: continue
        if cl['sc']:
            merged_node[ci] = mid; mid += 1
        else:
            merged_node[ci] = mid
            comp = cl['comp_cid']
            if comp >= 0 and comp != ci:
                merged_node[comp] = mid
            mid += 1

    V_merged = mid
    merged_edges = set()
    collapsed = 0
    for (a,b) in edges:
        ma, mb = merged_node[a], merged_node[b]
        if ma == mb: collapsed += 1
        else: merged_edges.add((min(ma,mb), max(ma,mb)))

    E_merged = len(merged_edges)
    twin = E - E_merged - collapsed
    merged_density = 2*E_merged/(V_merged*(V_merged-1)) if V_merged > 1 else 0

    # Merged degree sequence
    merged_deg = Counter()
    for (a,b) in merged_edges:
        merged_deg[a] += 1
        merged_deg[b] += 1
    merged_deg_vals = [merged_deg.get(i, 0) for i in range(V_merged)]

    # H-level analysis
    h_levels = defaultdict(list)
    for cl in classes:
        h_levels[cl['H']].append(cl['cid'])
    num_h_levels = len(h_levels)
    width = max(len(v) for v in h_levels.values())
    level_edges = sum(1 for (a,b) in edges if classes[a]['H'] == classes[b]['H'])

    # H range
    h_min = min(h_levels.keys())
    h_max = max(h_levels.keys())

    # |Aut| distribution
    aut_dist = Counter(cl.get('aut', 1) for cl in classes)

    # Degree stats
    deg_min = min(deg_vals)
    deg_max = max(deg_vals)
    deg_median = sorted(deg_vals)[V//2]

    # BFS diameter (only for small graphs)
    diameter = -1
    avg_path = -1
    if V <= 1000:
        from collections import deque
        max_d = 0
        total_d = 0
        count_d = 0
        for start in range(V):
            dist = [-1]*V
            dist[start] = 0
            q = deque([start])
            while q:
                u = q.popleft()
                for v in adj_list[u]:
                    if dist[v] == -1:
                        dist[v] = dist[u] + 1
                        q.append(v)
            for d in dist:
                if d > 0:
                    max_d = max(max_d, d)
                    total_d += d
                    count_d += 1
        diameter = max_d
        avg_path = total_d / count_d if count_d > 0 else 0

    # Independence number (greedy estimate for large graphs)
    if V <= 500:
        # Greedy independent set
        remaining = set(range(V))
        indep = set()
        for v in sorted(remaining, key=lambda x: degrees.get(x, 0)):
            if v in remaining:
                indep.add(v)
                remaining -= adj_list[v]
                remaining.discard(v)
        alpha_greedy = len(indep)
    else:
        alpha_greedy = -1

    # Clique number (greedy)
    if V <= 500:
        max_clique = 1
        for v in range(V):
            clique = {v}
            candidates = adj_list[v].copy()
            for u in sorted(candidates, key=lambda x: -degrees.get(x, 0)):
                if all(u in adj_list[w] for w in clique):
                    clique.add(u)
            max_clique = max(max_clique, len(clique))
        omega_greedy = max_clique
    else:
        omega_greedy = -1

    # Score sequence distribution
    score_dist = Counter(cl['score'] for cl in classes)
    num_score_classes = len(score_dist)

    # c3 range
    c3_min = min(cl['c3'] for cl in classes)
    c3_max = max(cl['c3'] for cl in classes)

    # Twin fraction and collapse fraction
    twin_frac = twin / E if E > 0 else 0
    collapse_frac = collapsed / E if E > 0 else 0
    merged_frac = E_merged / E if E > 0 else 0

    # Self-loop fraction
    sl_frac = len(self_loops) / V

    # Store results
    R = {
        'n': n, 'V': V, 'E': E, 'density': density,
        'avg_deg': avg_deg, 'deg_min': deg_min, 'deg_max': deg_max, 'deg_median': deg_median,
        'blue': blue_e, 'black': black_e,
        'scsc': scsc, 'nsns': nsns, 'scns': scns,
        'SC': len(sc_ids), 'NS': len(ns_ids),
        'V_merged': V_merged, 'E_merged': E_merged,
        'collapsed': collapsed, 'twin': twin,
        'merged_density': merged_density,
        'merged_deg_min': min(merged_deg_vals) if merged_deg_vals else 0,
        'merged_deg_max': max(merged_deg_vals) if merged_deg_vals else 0,
        'merged_avg_deg': sum(merged_deg_vals)/V_merged if V_merged > 0 else 0,
        'num_h_levels': num_h_levels, 'width': width, 'level_edges': level_edges,
        'h_min': h_min, 'h_max': h_max,
        'aut_dist': dict(aut_dist),
        'self_loops': len(self_loops), 'sl_frac': sl_frac,
        'diameter': diameter, 'avg_path': avg_path,
        'alpha_greedy': alpha_greedy, 'omega_greedy': omega_greedy,
        'num_score_classes': num_score_classes,
        'c3_min': c3_min, 'c3_max': c3_max,
        'twin_frac': twin_frac, 'collapse_frac': collapse_frac,
        'merged_frac': merged_frac,
        'm': m,
    }
    all_results[n] = R

    elapsed = time.time() - t0
    print(f"\n  n={n} COMPLETE in {elapsed:.1f}s")
    print(f"    V={V}, E={E}, SC={len(sc_ids)}, V_m={V_merged}, E_m={E_merged}")

# ============================================================================
# PRINT ALL METRICS
# ============================================================================

print(f"\n{'='*80}")
print(f"  COMPLETE METRICS TABLE")
print(f"{'='*80}")

# Print header row
metrics = [
    ('V', 'V'), ('E', 'E'), ('m', 'm'), ('SC', 'SC'), ('NS', 'NS'),
    ('density', 'dens'), ('avg_deg', 'avg_d'), ('deg_min', 'd_min'), ('deg_max', 'd_max'),
    ('blue', 'blue'), ('black', 'black'), ('scsc', 'SC²'), ('nsns', 'NS²'), ('scns', 'SC×NS'),
    ('V_merged', 'Vm'), ('E_merged', 'Em'), ('collapsed', 'coll'), ('twin', 'twin'),
    ('merged_density', 'm_dens'), ('merged_avg_deg', 'm_avg'),
    ('merged_deg_min', 'm_dmin'), ('merged_deg_max', 'm_dmax'),
    ('num_h_levels', '#H'), ('width', 'width'), ('level_edges', 'lvl_E'),
    ('h_min', 'Hmin'), ('h_max', 'Hmax'),
    ('self_loops', 'SL'), ('sl_frac', 'SLf'),
    ('diameter', 'diam'), ('avg_path', 'avg_p'),
    ('alpha_greedy', 'α'), ('omega_greedy', 'ω'),
    ('num_score_classes', '#scores'),
    ('twin_frac', 'tw_f'), ('collapse_frac', 'co_f'),
]

for key, label in metrics:
    vals = []
    for n in sorted(all_results):
        v = all_results[n].get(key, '?')
        if isinstance(v, float):
            vals.append(f"{v:.4f}")
        else:
            vals.append(str(v))
    print(f"  {label:>8s}: {', '.join(vals)}")

# ============================================================================
# EXTRACT NEW OEIS SEQUENCES
# ============================================================================

print(f"\n{'='*80}")
print(f"  NEW OEIS-WORTHY SEQUENCES (from G_n and G_n/Z_2)")
print(f"{'='*80}")

seqs = {
    'E_orig': [all_results[n]['E'] for n in sorted(all_results)],
    'E_merged': [all_results[n]['E_merged'] for n in sorted(all_results)],
    'collapsed': [all_results[n]['collapsed'] for n in sorted(all_results)],
    'twin': [all_results[n]['twin'] for n in sorted(all_results)],
    'blue': [all_results[n]['blue'] for n in sorted(all_results)],
    'black': [all_results[n]['black'] for n in sorted(all_results)],
    'SC-SC': [all_results[n]['scsc'] for n in sorted(all_results)],
    'NS-NS': [all_results[n]['nsns'] for n in sorted(all_results)],
    'SC-NS': [all_results[n]['scns'] for n in sorted(all_results)],
    'SC': [all_results[n]['SC'] for n in sorted(all_results)],
    'V_merged': [all_results[n]['V_merged'] for n in sorted(all_results)],
    '#H_levels': [all_results[n]['num_h_levels'] for n in sorted(all_results)],
    'width': [all_results[n]['width'] for n in sorted(all_results)],
    'level_edges': [all_results[n]['level_edges'] for n in sorted(all_results)],
    'self_loops': [all_results[n]['self_loops'] for n in sorted(all_results)],
    'H_max': [all_results[n]['h_max'] for n in sorted(all_results)],
    'diameter': [all_results[n]['diameter'] for n in sorted(all_results)],
    'alpha': [all_results[n]['alpha_greedy'] for n in sorted(all_results)],
    'omega': [all_results[n]['omega_greedy'] for n in sorted(all_results)],
    '#scores': [all_results[n]['num_score_classes'] for n in sorted(all_results)],
    'deg_max': [all_results[n]['deg_max'] for n in sorted(all_results)],
    'merged_deg_max': [all_results[n]['merged_deg_max'] for n in sorted(all_results)],
}

for name, seq in sorted(seqs.items()):
    print(f"  {name:>15s}: {seq}")

# ============================================================================
# RATIO ANALYSIS
# ============================================================================

print(f"\n{'='*80}")
print(f"  RATIO ANALYSIS")
print(f"{'='*80}")

for n in sorted(all_results):
    R = all_results[n]
    print(f"\n  n={n}:")
    print(f"    twin/E = {R['twin_frac']:.6f}")
    print(f"    collapsed/E = {R['collapse_frac']:.6f}")
    print(f"    E_merged/E = {R['merged_frac']:.6f}")
    print(f"    blue/E = {R['blue']/R['E']:.6f}" if R['E'] > 0 else "")
    print(f"    level_edges/E = {R['level_edges']/R['E']:.6f}" if R['E'] > 0 else "")
    print(f"    SC-SC/blue = {R['scsc']/R['blue']:.6f}" if R['blue'] > 0 else "")
    print(f"    avg_deg/m = {R['avg_deg']/R['m']:.6f}")
    print(f"    merged_avg_deg/m = {R['merged_avg_deg']/R['m']:.6f}")
    if R['diameter'] >= 0:
        print(f"    diameter = {R['diameter']}")

total_time = sum(time.time() - t0 for _ in [0])  # rough
print(f"\n{'='*80}")
print(f"  DONE")
print(f"{'='*80}")
