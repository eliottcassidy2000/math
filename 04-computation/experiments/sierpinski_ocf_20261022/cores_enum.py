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

if __name__ == '__main__':
    sys.path.insert(0, '.')
    from loop_engine import automorphisms
    exc = int(sys.argv[1])
    t0 = time.time()
    cs = cores(exc)
    out = []
    for V, edges in cs:
        out.append({'V': V, 'edges': edges, 'aut': automorphisms(V, [tuple(e) for e in edges])})
    by_v = defaultdict(int)
    for c in out: by_v[c['V']] += 1
    # cross-check: sum over cores of 1/|Aut| (labelled count / (V! ...)) is reported for comparison with references
    from fractions import Fraction
    s = sum(Fraction(1, c['aut']) for c in out)
    print(f"excess {exc}: {len(out)} cores ({dict(sorted(by_v.items()))} by V), sum 1/|Aut| = {s} = {float(s):.6f}, {time.time() - t0:.1f}s")
    json.dump(out, open(f'cores_e{exc}.json', 'w'))
