# Cores of the linear-forest loop expansion (THM-4757 / HYP-9367): connected multigraphs, loops allowed, minimum
# degree >= 3 (a loop counts 2), E - V = excess, up to isomorphism; hence V <= 2 * excess.  Faster than
# loop_engine.cores (which runs over all multisets of E pairs): labelled multigraphs are generated per degree sequence
# by assigning multiplicities pair by pair with degree budgets, and canonised by the minimum over the permutations
# that preserve the (non-increasing) degree sequence.  Reproduces loop_engine.cores for excess 1, 2 (3 and 15 cores).
# Usage: python3 cores_enum.py EXCESS  -> writes cores_e<EXCESS>.json.  mac-mini-2026-10-22-loop4.
import sys, json, itertools, time
from collections import defaultdict

def degree_sequences(V, total):
    """non-increasing sequences of V integers >= 3 with the given sum"""
    def rec(i, rem, hi):
        if i == V:
            if rem == 0: yield ()
            return
        left = V - i - 1
        for d in range(min(hi, rem - 3 * left), 2, -1):
            for rest in rec(i + 1, rem - d, d): yield (d,) + rest
    yield from rec(0, total, total)

def labelled_multigraphs(degs):
    V = len(degs)
    pairs = [(a, b) for a in range(V) for b in range(a, V)]
    last = {a: pairs.index((a, V - 1)) if a < V - 1 else pairs.index((V - 1, V - 1)) for a in range(V)}
    closes = defaultdict(list)
    for a, k in last.items(): closes[k].append(a)
    r = list(degs); mult = [0] * len(pairs)
    def rec(k):
        if k == len(pairs):
            yield tuple(mult); return
        a, b = pairs[k]
        mmax = r[a] // 2 if a == b else min(r[a], r[b])
        for m in range(mmax + 1):
            if a == b: r[a] -= 2 * m
            else: r[a] -= m; r[b] -= m
            mult[k] = m
            if all(r[x] == 0 for x in closes[k]):
                yield from rec(k + 1)
            if a == b: r[a] += 2 * m
            else: r[a] += m; r[b] += m
        mult[k] = 0
    for mlt in rec(0):
        edges = []
        for (a, b), m in zip(pairs, mlt): edges += [(a, b)] * m
        yield edges

def connected(V, edges):
    par = list(range(V))
    def f(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for a, b in edges: par[f(a)] = f(b)
    return len({f(x) for x in range(V)}) == 1

def degree_preserving_perms(degs):
    groups = defaultdict(list)
    for v, d in enumerate(degs): groups[d].append(v)
    blocks = [groups[d] for d in sorted(groups, reverse=True)]
    for choice in itertools.product(*[list(itertools.permutations(b)) for b in blocks]):
        p = [0] * len(degs)
        for b, img in zip(blocks, choice):
            for v, w in zip(b, img): p[v] = w
        yield p

def canon(edges, perms):
    return min(tuple(sorted((min(p[a], p[b]), max(p[a], p[b])) for a, b in edges)) for p in perms)

def refine(V, M, colors):
    """colour refinement (1-dim Weisfeiler-Leman) of a multigraph given by its multiplicity matrix M (M[v][v] =
    loops); colours are relabelled canonically by the sorted order of their signatures"""
    k = len(set(colors))
    while True:
        sig = [(colors[v], tuple(sorted((colors[w], M[v][w]) for w in range(V) if w != v and M[v][w])))
               for v in range(V)]
        uniq = sorted(set(sig))
        idx = {s: i for i, s in enumerate(uniq)}
        colors = [idx[s] for s in sig]
        if len(uniq) == k: return colors
        k = len(uniq)

def canon_ir(V, M):
    """canonical form by individualisation-refinement: the lexicographically smallest upper-triangular
    multiplicity vector over all leaves of the search tree (no automorphism pruning; fine for V <= 8)"""
    deg = [2 * M[v][v] + sum(M[v][w] for w in range(V) if w != v) for v in range(V)]
    init = sorted(set((deg[v], M[v][v]) for v in range(V)))
    colors = [init.index((deg[v], M[v][v])) for v in range(V)]
    best = [None]
    def search(colors):
        colors = refine(V, M, colors)
        if len(set(colors)) == V:
            order = sorted(range(V), key=lambda v: colors[v])
            vec = tuple(M[order[i]][order[j]] for i in range(V) for j in range(i, V))
            if best[0] is None or vec < best[0]: best[0] = vec
            return
        cnt = defaultdict(int)
        for c in colors: cnt[c] += 1
        c0 = min(c for c, k in cnt.items() if k > 1)
        for w in range(V):
            if colors[w] == c0:
                search([2 * c if v == w else 2 * c + 1 for v, c in enumerate(colors)])
    search(colors)
    return best[0]

def vec_to_edges(V, vec):
    edges, k = [], 0
    for i in range(V):
        for j in range(i, V):
            edges += [[i, j]] * vec[k]; k += 1
    return edges

def cores_ir(excess, vmax=None):
    """as cores(), with individualisation-refinement canonical forms (needed for V = 7, 8 at excess 4)"""
    found = []
    for V in range(1, (vmax or 2 * excess) + 1):
        E = V + excess
        for degs in degree_sequences(V, 2 * E):
            seen = set()
            for edges in labelled_multigraphs(degs):
                if not connected(V, edges): continue
                M = [[0] * V for _ in range(V)]
                for a, b in edges:
                    if a == b: M[a][a] += 1
                    else: M[a][b] += 1; M[b][a] += 1
                c = canon_ir(V, M)
                if c in seen: continue
                seen.add(c); found.append((V, vec_to_edges(V, c)))
    return found

def cores(excess):
    found = []
    for V in range(1, 2 * excess + 1):
        E = V + excess
        for degs in degree_sequences(V, 2 * E):
            perms = list(degree_preserving_perms(degs))
            seen = set()
            for edges in labelled_multigraphs(degs):
                if not connected(V, edges): continue
                c = canon(edges, perms)
                if c in seen: continue
                seen.add(c); found.append((V, [list(x) for x in c]))
    return found

def configuration_check(max_excess, top_file=None, vmax=None):
    """Completeness and automorphism check: exp(sum over cores of x^V y^E/|G|), truncated at excess <= max_excess,
    must equal the weighted count of all (possibly disconnected) multigraphs with loops and minimum degree >= 3,
    (1/V!) sum over labelled degree sequences of (2E-1)!!/prod d_v! (configuration model), at every (V, E)."""
    from fractions import Fraction
    from math import factorial
    def dfact(m):
        r = 1
        for k in range(m, 0, -2): r *= k
        return r
    def total(V, E):
        s = Fraction(0)
        def rec(i, rem, acc):
            nonlocal s
            if i == V:
                if rem == 0: s += acc
                return
            for d in range(3, rem - 3 * (V - i - 1) + 1): rec(i + 1, rem - d, acc / factorial(d))
        rec(0, 2 * E, Fraction(1))
        return s * dfact(2 * E - 1) / factorial(V)
    C = {}
    for e in range(1, max_excess + 1):
        fn = top_file if (e == max_excess and top_file) else f'cores_e{e}.json'
        for c in json.load(open(fn)):
            key = (c['V'], c['V'] + e); C[key] = C.get(key, 0) + Fraction(1, c['aut'])
    def mul(P, Q):
        R = {}
        for (v1, e1), a in P.items():
            for (v2, e2), b in Q.items():
                if (e1 - v1) + (e2 - v2) <= max_excess: R[(v1 + v2, e1 + e2)] = R.get((v1 + v2, e1 + e2), 0) + a * b
        return R
    A = {(0, 0): Fraction(1)}; term = {(0, 0): Fraction(1)}
    for k in range(1, max_excess + 1):
        term = {key: v / k for key, v in mul(term, C).items()}
        for key, v in term.items(): A[key] = A.get(key, 0) + v
    pts = [(V, E) for (V, E) in A if (V, E) != (0, 0) and (vmax is None or E - V < max_excess or V <= vmax)]
    return all(total(V, E) == A[(V, E)] for (V, E) in pts), len(pts)

if __name__ == '__main__':
    sys.path.insert(0, '.')
    from loop_engine import automorphisms
    exc = int(sys.argv[1])
    vmax = int(sys.argv[2]) if len(sys.argv) > 2 else None          # optional: only cores with V <= vmax
    t0 = time.time()
    cs = cores(exc) if exc <= 3 and vmax is None else cores_ir(exc, vmax)   # excess <= 3 files keep their labelling
    out = []
    for V, edges in cs:
        out.append({'V': V, 'edges': edges, 'aut': automorphisms(V, [tuple(e) for e in edges])})
    by_v = defaultdict(int)
    for c in out: by_v[c['V']] += 1
    # cross-check: sum over cores of 1/|Aut| (labelled count / (V! ...)) is reported for comparison with references
    from fractions import Fraction
    s = sum(Fraction(1, c['aut']) for c in out)
    print(f"excess {exc}: {len(out)} cores ({dict(sorted(by_v.items()))} by V), sum 1/|Aut| = {s} = {float(s):.6f}, {time.time() - t0:.1f}s")
    fname = f'cores_e{exc}.json' if vmax is None else f'cores_e{exc}_v{vmax}.json'
    json.dump(out, open(fname, 'w'))
    import os
    if vmax is None and all(os.path.exists(f'cores_e{e}.json') for e in range(1, exc + 1)):
        ok, npts = configuration_check(exc)
        print(f"configuration-model check, excess <= {exc}, {npts} (V, E) points: {'exact agreement' if ok else 'MISMATCH'}")
    if vmax is not None and all(os.path.exists(f'cores_e{e}.json') for e in range(1, exc)):
        ok, npts = configuration_check(exc, fname, vmax)
        print(f"configuration-model check, excess <= {exc} (top level V <= {vmax}), {npts} (V, E) points: {'exact agreement' if ok else 'MISMATCH'}")
