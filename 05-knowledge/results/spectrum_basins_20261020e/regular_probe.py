"""Complete labelled regular tournaments n=1,3,5,7; exact triangle reversal controls."""
from collections import Counter, defaultdict, deque
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
import networkx as nx


def check(ok, detail):
    if not ok:
        raise RuntimeError(detail)


def regular_tournaments(n):
    r = (n-1)//2
    out = [0]*n
    deg = [0]*n
    def rec(i):
        if i == n:
            yield tuple(out)
            return
        remaining = list(range(i+1, n))
        need = r-deg[i]
        if not 0 <= need <= len(remaining):
            return
        for wins in combinations(remaining, need):
            wins = set(wins)
            for j in remaining:
                if j in wins:
                    out[i] |= 1 << j
                    deg[i] += 1
                else:
                    out[j] |= 1 << i
                    deg[j] += 1
            if all(deg[j] <= r <= deg[j]+n-i-2 for j in remaining):
                yield from rec(i+1)
            for j in remaining:
                if j in wins:
                    out[i] ^= 1 << j
                    deg[i] -= 1
                else:
                    out[j] ^= 1 << i
                    deg[j] -= 1
    yield from rec(0)


def hpaths(out):
    n = len(out)
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            count = dp[mask][v]
            if count:
                available = out[v] & ~mask
                while available:
                    bit = available & -available
                    available -= bit
                    dp[mask | bit][bit.bit_length()-1] += count
    return sum(dp[-1])


def is_directed(out, cycle):
    return all(out[a] >> b & 1 for a, b in zip(cycle, cycle[1:]+cycle[:1]))


def reverse(out, edges):
    target = list(out)
    for a, b in edges:
        target[a] ^= 1 << b
        target[b] ^= 1 << a
    return tuple(target)


def triangle_neighbors(out):
    for triple in combinations(range(len(out)), 3):
        a, b, c = triple
        if not is_directed(out, triple):
            b, c = c, b
        if is_directed(out, (a, b, c)):
            yield (a, b, c), reverse(out, [(a,b), (b,c), (c,a)])


def graph(out):
    g = nx.DiGraph()
    g.add_nodes_from(range(len(out)))
    g.add_edges_from((i,j) for i, mask in enumerate(out) for j in range(len(out)) if mask >> j & 1)
    return g


def triangle_compiler(out, cyc):
    """Return literal successive directed triangles realizing one cycle reversal."""
    check(is_directed(out, cyc), ('not directed', out, cyc))
    if len(cyc) == 3:
        return [cyc]
    a, b, c = cyc[:3]
    short = (a,) + cyc[2:]
    tri = (a,b,c)
    if out[c] >> a & 1:
        first = reverse(out, [(a,b), (b,c), (c,a)])
        return [tri] + triangle_compiler(first, short)
    seq = triangle_compiler(out, short)
    return seq + [tri]


def cycle_controls():
    total = 0
    for n in range(3, 6):
        pairs = list(combinations(range(n), 2))
        for code in range(1 << len(pairs)):
            out = [0]*n
            for bit, (i,j) in enumerate(pairs):
                a,b = (i,j) if code >> bit & 1 else (j,i)
                out[a] |= 1 << b
            out = tuple(out)
            for size in range(3, n+1):
                for subset in combinations(range(n), size):
                    for tail in permutations(subset[1:]):
                        cyc = (subset[0],) + tail
                        if not is_directed(out, cyc):
                            continue
                        seq = triangle_compiler(out, cyc)
                        check(len(seq) == size-2, ('length', cyc, seq))
                        cur = out
                        for tri in seq:
                            check(is_directed(cur, tri), ('compiler directedness', cyc, seq))
                            cur = reverse(cur, list(zip(tri, tri[1:]+tri[:1])))
                        check(cur == reverse(out, list(zip(cyc, cyc[1:]+cyc[:1]))), ('compiler target', cyc, seq))
                        total += 1
    return total


def main():
    rows = []
    for n in (1,3,5,7):
        ts = sorted(regular_tournaments(n))
        universe = set(ts)
        check(len(universe) == len(ts), ('duplicates', n))
        hs = {t: hpaths(t) for t in ts}
        reps, repgraphs, classes = [], [], {}
        buckets = defaultdict(list)
        for t in ts:
            g = graph(t)
            for cid in buckets[hs[t]]:
                if nx.is_isomorphic(g, repgraphs[cid]):
                    break
            else:
                cid = len(reps)
                reps.append(t)
                repgraphs.append(g)
                buckets[hs[t]].append(cid)
            classes[t] = cid
        # Independent permutation count at each ordinary isomorphism representative.
        path_reversal_controls = 0
        for t in reps:
            paths = [p for p in permutations(range(n)) if all(t[a] >> b & 1 for a,b in zip(p,p[1:]))]
            brute = len(paths)
            check(brute == hs[t], ('independent H', n, t))
            check(hs[t] % 2 == 1, ('odd H', n, t))
            if n >= 3:
                for p in paths:
                    u = reverse(t, list(zip(p,p[1:])))
                    degree = [(n-1)//2]*n
                    degree[p[0]] -= 1
                    degree[p[-1]] += 1
                    check([x.bit_count() for x in u] == degree, ('R score', n, t, p))
                    k = tuple(((1 << n)-1) ^ (1 << v) ^ u[v] for v in range(n))
                    check([x.bit_count() for x in k] == [n-1-d for d in degree], ('K score', n, t, p))
                    changed = sum(((t[i] ^ u[i]) >> (i+1)).bit_count() for i in range(n))
                    check(changed == n-1 and changed % 2 == 0, ('R parity', n, p))
                    path_reversal_controls += 1
        adj, differences, quotient, witnesses = {}, Counter(), Counter(), {}
        triangle_count = n*(n-1)*(n-2)//6 - n*((n-1)//2)*((n-3)//2)//2
        for t in ts:
            neighbors = list(triangle_neighbors(t))
            check(len(neighbors) == triangle_count, ('triangle count', n, t))
            adj[t] = []
            for tri, u in neighbors:
                check(u in universe, ('regularity', n, t, tri))
                adj[t].append(u)
                delta = hs[u]-hs[t]
                check(delta % 2 == 0, ('H parity', t, u))
                differences[delta] += 1
                quotient[(classes[t], classes[u])] += 1
                key = (classes[t], classes[u], delta)
                witnesses.setdefault(key, {'source_out_masks': t, 'triangle': tri, 'target_out_masks': u,
                                           'source_H': hs[t], 'target_H': hs[u]})
            for i,j in combinations(range(n),2):
                u = reverse(t, [(i,j)])
                check([x.bit_count() for x in u] != [(n-1)//2]*n, ('single arc hostile', n, t))
        seen, q = {ts[0]}, deque([ts[0]])
        while q:
            for u in adj[q.popleft()]:
                if u not in seen:
                    seen.add(u)
                    q.append(u)
        check(seen == universe, ('labelled connectivity', n))
        rows.append({'n': n, 'labelled_count': len(ts), 'ordinary_isomorphism_classes': len(reps),
                     'H_distribution_labelled': dict(sorted(Counter(hs.values()).items())),
                     'triangle_count_per_tournament': triangle_count,
                     'all_class_representative_HP_R_and_K_controls': path_reversal_controls,
                     'directed_triangle_move_H_deltas': dict(sorted(differences.items())),
                     'connected': True,
                     'classes': [{'id': i, 'out_masks': t, 'H': hs[t],
                                  'labelled_orbit_size': sum(v == i for v in classes.values())}
                                 for i,t in enumerate(reps)],
                     'directed_quotient_move_counts_including_loops': [{'from': i, 'to': j, 'count': count}
                                                                     for (i,j), count in sorted(quotient.items())],
                     'move_witnesses': [{'from': i, 'to': j, 'delta_H': d, **w}
                                        for (i,j,d),w in sorted(witnesses.items())]})
    payload = {'status': 'FINITE-EXACT', 'universe': 'All labelled regular tournaments on vertex sets 0..n-1 for n=1,3,5,7; every directed triangle reversal',
               'source_sha256': sha256(Path(__file__).read_bytes()).hexdigest(), 'rows': rows,
               'independent_controls': ['H counted by permutations at every ordinary class representative',
                                        'every single arc flip destroys regularity', 'all directed cycles in all tournaments n=3,4,5 compiled to k-2 triangle reversals'],
               'cycle_compiler_cases': cycle_controls(),
               'claims_not_made': ['triangle moves are Hamiltonian-path metagraph edges', 'H invariant under degree-preserving moves',
                                   'unlabelled regular tournament classified by H at every order']}
    output = Path(__file__).with_suffix('.json')
    output.write_text(json.dumps(payload, indent=2, sort_keys=True)+'\n', encoding='utf-8')
    print(json.dumps({'all_checks': 'PASS', 'summary': [{k:r[k] for k in ('n','labelled_count','ordinary_isomorphism_classes','H_distribution_labelled','directed_triangle_move_H_deltas')}
                                                       for r in rows], 'cycle_compiler_cases': payload['cycle_compiler_cases'],
                      'receipt_sha256': sha256(output.read_bytes()).hexdigest()}, sort_keys=True))


if __name__ == '__main__':
    main()
