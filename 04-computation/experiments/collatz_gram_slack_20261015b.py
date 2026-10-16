#!/usr/bin/env python3
"""Exact fixed-source Gram certificates from the union of checked native prefixes.
The corpus and graph are frozen before independent ROOT labels are computed.
No floating-point quantity enters an acceptance decision.
"""
import hashlib
import itertools
import json
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from collatz_clock_mesh_audit_20261015a import envelope, need

ROOT = Path(__file__).resolve().parents[2]

def T(n):
    return (3*n+1)//2 if n & 1 else n//2

def digest(value):
    return hashlib.sha256(json.dumps(value, separators=(',', ':'), sort_keys=True).encode()).hexdigest()

def observe(N, depth):
    edges = {}
    for n in range(1, N):
        x = n
        for _ in range(depth):
            if x == 1:
                break
            y = T(x)
            edges[x] = y
            x = y
    return edges

def close_graph(edges, sources):
    """Compile exact finite paths. A witnessed non-ROOT cycle has ROOT label zero."""
    cache = {1: (1, 1, 0)}
    cycles = []
    for n in sources:
        x = n
        path, where = [], {}
        while x not in cache:
            if x not in edges:
                cache[x] = (x, 1, 0)
                break
            if x in where:
                cut = where[x]
                cycle = path[cut:]
                cycles.append(cycle)
                for y in cycle:
                    cache[y] = (0, 0, 0)
                path = path[:cut]
                break
            where[x] = len(path)
            path.append(x)
            x = edges[x]
        e, s, d = cache[x]
        for x in reversed(path):
            s, d = -s, d+1
            cache[x] = (e, s, d)
    # Local audit is a linear-size proof of every compiled source identity.
    cyclic = {n for cycle in cycles for n in cycle}
    for cycle in cycles:
        need(1 not in cycle, 'non-ROOT cycle')
        need(all(edges[x] == cycle[(i+1) % len(cycle)] for i, x in enumerate(cycle)), 'literal cycle')
    for x, (e, s, d) in cache.items():
        if x == 1:
            need((e, s, d) == (1, 1, 0), 'marked ROOT')
        elif x in cyclic:
            need((e, s, d) == (0, 0, 0), 'cycle support zero')
        elif x not in edges:
            need((e, s, d) == (x, 1, 0), 'unresolved sink retained')
        else:
            ee, ss, dd = cache[edges[x]]
            need((e, s, d) == (ee, -ss, dd+1), 'signed edge and integer step receipt')
    return cache, cycles

def rows_from_cache(k, cache):
    X, b = 2**k, 2**(k//2)
    prefix = defaultdict(int)
    next_source = 1
    rows = []
    for a in range(X, 2*X, b):
        while next_source < a:
            n = next_source
            e, s, _ = cache[n]
            prefix[e] += (2 if n % 3 == 2 else -1)*s
            next_source += 1
        row = {e: -b*v for e, v in prefix.items() if v}
        for n in range(a, a+b):
            e, s, _ = cache[n]
            w = 2 if n % 3 == 2 else -1
            row[e] = row.get(e, 0) + (-(a+b-n)*w + (0, -1, 1)[n % 3])*s
        rows.append({e: v for e, v in row.items() if v})
    return rows

def gram(rows, cols):
    values = [[row.get(e, 0) for e in cols] for row in rows]
    top = max(abs(x) for row in values for x in row)
    need(len(values)*top*top < 2**63, 'proved int64 matrix product budget')
    data = np.array(values, dtype=np.int64)
    G = (data.T @ data).tolist()
    # Independent direct dot products: diagonal plus all entries for small matrices.
    pairs = itertools.product(range(len(cols)), repeat=2) if len(cols) <= 50 else ((i, i) for i in range(len(cols)))
    for i, j in pairs:
        need(G[i][j] == sum(row[i]*row[j] for row in values), 'integer Gram dot product')
    return G

def transport_slack(G, cols, controller):
    """controller[e]=(target, sign); a sign zero deletes a certified zero label."""
    target_cols = sorted({controller[e][0] for e in cols if controller[e][1]})
    idx = {e: i for i, e in enumerate(target_cols)}
    P = [[0]*len(idx) for _ in idx]
    N = [[0]*len(idx) for _ in idx]
    deleted = 0
    for i, e in enumerate(cols):
        ee, se = controller[e]
        for j, f in enumerate(cols):
            ff, sf = controller[f]
            value = se*sf*G[i][j]
            if not se or not sf:
                deleted += abs(G[i][j])
                continue
            u, v = idx[ee], idx[ff]
            if value >= 0:
                P[u][v] += value
            else:
                N[u][v] -= value
    new = [[p-n for p, n in zip(rp, rn)] for rp, rn in zip(P, N)]
    loss = 2*sum(min(p, n) for rp, rn in zip(P, N) for p, n in zip(rp, rn))
    variance = sum((F(4*p*n, p+n) for rp, rn in zip(P, N) for p, n in zip(rp, rn) if p+n), F(0))
    old_norm = sum(abs(v) for row in G for v in row)
    new_norm = sum(abs(v) for row in new for v in row)
    need(old_norm-new_norm == deleted+loss, 'exact signed pair-fibre loss')
    need(variance/2 <= loss <= variance, 'collision variance comparison')
    return target_cols, new, {'old_norm': old_norm, 'new_norm': new_norm,
        'deleted_mass': deleted, 'cancelled_mass': loss,
        'pair_variance_exact': [variance.numerator, variance.denominator]}

def independent_labels(N):
    # Separate first-descent grounding, used ONLY after all certificates are frozen.
    labels = [0]*N
    labels[1] = 1
    for n in range(2, N):
        x, steps = n, 0
        while x >= n:
            x, steps = T(x), steps+1
            need(steps <= 10000, 'finite audit orbit cap')
        labels[n] = (-1 if steps & 1 else 1)*labels[x]
    return labels

def case(k, depth):
    X, b, N = 2**k, 2**(k//2), 2**(k+1)
    edges = observe(N, depth)
    cache, cycles = close_graph(edges, range(1, N))
    need(all(T(x) == y and x != 1 for x, y in edges.items()), 'authenticated native edges')
    rows = rows_from_cache(k, cache)
    cols = sorted({1} | {e for row in rows for e in row})
    G = gram(rows, cols)
    norm = sum(abs(v) for row in G for v in row)
    row_norm = sum(sum(abs(v) for v in row.values())**2 for row in rows)
    fine = sum((envelope(a, b) for a in range(X, 2*X, b)), F(0))
    coarse = F(norm, 4*b)
    Z = coarse+fine
    A = 2*X*X-4*Z-X
    passed = A >= 0 and A*A >= 16*X*Z
    # One fresh common extension, independently recompile its rows and compare.
    controller = {e: (1, 1) if e == 1 else (T(e), -1) for e in cols}
    new_cols, new_G, slack = transport_slack(G, cols, controller)
    new_rows = []
    for row in rows:
        out = defaultdict(int)
        for e, v in row.items():
            ee, s = controller[e]
            out[ee] += s*v
        new_rows.append(dict(out))
    need(gram(new_rows, new_cols) == new_G, 'independent row transport vs pair transport')
    # Freeze source receipts and coefficients before computing the audit labels.
    frozen = digest([sorted(edges.items()), [cache[n] for n in range(1, N)], rows])
    labels = independent_labels(N)
    J, currents = 0, [0]*N
    for n in range(1, N-1):
        currents[n] = J
        J += (0, -1, 1)[n % 3]*(labels[n+1]-labels[n])//2
    currents[N-1] = J
    actual_coarse = F(0)
    actual_energy = sum(currents[n]**2 for n in range(X, N))
    unknown = 0
    for n in range(1, N):
        e, s, d = cache[n]
        unknown += e != 1
        if e == 1:
            need(labels[n] == s, 'independent first-descent root sign')
        # Direct finite path replay on the complete small-source universe.
        if k <= 12:
            x = n
            for _ in range(d):
                need(x in edges and x != 1, 'receipt step belongs to observed graph')
                x = T(x)
            need(x == e, 'literal signed endpoint replay')
    for a, row in zip(range(X, N, b), rows):
        M = sum(currents[a:a+b])
        actual_coarse += F(M*M, b)
        need(abs(2*M-row.get(1, 0)) <= sum(abs(v) for e, v in row.items() if e != 1), 'actual source mean enclosed')
    need(actual_coarse <= coarse and actual_energy <= Z, 'joint and fine-scale audit')
    need(frozen == digest([sorted(edges.items()), [cache[n] for n in range(1, N)], rows]), 'no oracle mutation of certificate')
    result = {'k': k, 'prefix_depth_per_source': depth, 'source_count': N-1,
        'distinct_authenticated_edges': len(edges), 'max_compiled_steps': max(cache[n][2] for n in range(1, N)),
        'unresolved_source_count': unknown, 'terminal_columns': cols, 'nonroot_cycles': cycles,
        'gram_integer_numerator': G, 'gram_denominator': 4*b,
        'coarse_upper_exact': [coarse.numerator, coarse.denominator],
        'coarse_upper_over_X2': float(coarse/X**2),
        'rowwise_S_upper_over_X2': float((F(row_norm, 4*b)+fine)/X**2),
        'joint_S_upper_over_X2': float(Z/X**2), 'joint_E_test_passes': passed,
        'coarse_4_over_9_passes': 9*norm <= 16*b*X*X,
        'actual_S_over_X2': actual_energy/X**2, 'actual_coarse_over_X2': float(actual_coarse/X**2),
        'one_step_exact_slack': slack, 'frozen_graph_source_rows_sha256': frozen}
    print({key: result[key] for key in ('k', 'prefix_depth_per_source', 'distinct_authenticated_edges',
           'max_compiled_steps', 'unresolved_source_count', 'coarse_upper_over_X2', 'joint_E_test_passes')}, flush=True)
    return result

def hostile_controls():
    # Exact native tiny universe, including asynchronous 3->5->8 and 8 already a source.
    cache, cycles = close_graph(observe(9, 1), range(1, 9))
    need(cache[3] == (1, 1, 5) if False else cache[3] == (1, -1, 5), 'native cross-clock joining')
    # A free terminal preserves norm forever under a collision-free renaming.
    _, _, s = transport_slack([[7]], [9], {9: (14, -1)})
    need(s['cancelled_mass'] == 0, 'no automatic strict contraction')
    # Toy maps, explicitly not asserted native Collatz edges: both nonroot cycle parities have ROOT support zero.
    for edges in ({2: 3, 3: 4, 4: 2}, {2: 3, 3: 2}):
        cache, cycles = close_graph(edges, [2, 3])
        need(all(cache[n][1] == 0 for n in (2, 3)) and len(cycles) == 1, 'nonroot-cycle support rule')
    # Signed permutation A=-P has A^3=-I, A^6=I, yet A^T I A=I.
    A = np.array([[0, -1, 0], [0, 0, -1], [-1, 0, 0]], dtype=np.int64)
    I = np.eye(3, dtype=np.int64)
    need(np.array_equal(np.linalg.matrix_power(A, 3), -I), 'linear period six')
    need(np.array_equal(A.T @ I @ A, I), 'quadratic stationary hostile')
    need(round(np.linalg.det(I-A)) == 2, 'no nonzero fixed vector')
    return {'native_cross_clock_source_3_steps': 5, 'zero_loss_free_column': True,
        'both_nonroot_cycle_parities_grounded_to_zero': True, 'quadratic_forgets_sign_hostile': True}

def main():
    controls = hostile_controls()
    # Declared finite universe; no adaptive removal of failures.
    cases = [case(k, depth) for k, depth in ((10,8),(10,16),(12,16),(14,16),(14,24),(14,32),(16,32),(18,36))]
    out = {'status': 'FINITE-EXACT; uniform joint-terminal Gram bound OPEN',
           'controller': 'union of all depth-D native prefixes from every source 1<=n<2X, ROOT stopped; then exact graph closure',
           'controls': controls, 'cases': cases}
    (ROOT/'05-knowledge/results/collatz_gram_slack_20261015b.json').write_text(json.dumps(out, indent=2)+'\n')

if __name__ == '__main__':
    main()
