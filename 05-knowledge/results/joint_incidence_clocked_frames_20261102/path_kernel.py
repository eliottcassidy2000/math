"""Joint path intersections vs a separate Shannon compiler, exact integers.

The counts are labelled guarded tournament responses, NOT A000568.
"""
from collections import Counter
from functools import lru_cache
from itertools import combinations, permutations
from math import comb, isqrt
from pathlib import Path
import json


def require(ok, why):
    if not ok:
        raise RuntimeError(why)


def graph_paths(n, square):
    edges = [(i, j) for i, j in combinations(range(n), 2)
             if not square or isqrt(i+j+2)**2 == i+j+2]
    if not square:
        return edges, list(permutations(range(n)))
    adj = [set() for _ in range(n)]
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
    leaves = [v for v in range(n) if len(adj[v]) == 1]
    if any(not s for s in adj) or len(leaves) > 2:
        return edges, []
    paths = []
    def dfs(path, unseen):
        if not unseen:
            paths.append(tuple(path))
            return
        for v in sorted(adj[path[-1]] & unseen):
            dfs(path+[v], unseen-{v})
    for start in leaves[:1] if leaves else range(n):
        dfs([start], set(range(n))-{start})
    if leaves:
        paths += [p[::-1] for p in paths[:]]
    require(len(paths) == len(set(paths)), "path enumeration uniqueness")
    return edges, paths


def patterns(edges, paths):
    index = {e: i for i, e in enumerate(edges)}
    result = []
    for path in paths:
        pos = neg = 0
        for u, v in zip(path, path[1:]):
            bit = 1 << index[min(u, v), max(u, v)]
            if u < v:
                pos |= bit
            else:
                neg |= bit
        result.append((pos, neg))
    require(len(result) == len(set(result)), "distinct path arc sets")
    return result


def moment_spectrum(pats, edge_count, max_j=None):
    """Enumerate compatible subsets as cliques, retaining union incidence."""
    compatibility = []
    for i, (p, n) in enumerate(pats):
        compatibility.append(sum(1 << j for j, (q, m) in enumerate(pats)
                                 if j > i and not (p & m or n & q)))
    moments = Counter()
    visits = 0
    def visit(candidates, pos, neg, k):
        nonlocal visits
        visits += 1
        moments[k] += 1 << (edge_count-(pos | neg).bit_count())
        if max_j is not None and k >= max_j:
            return
        while candidates:
            bit = candidates & -candidates
            candidates -= bit
            i = bit.bit_length()-1
            p, n = pats[i]
            visit(candidates & compatibility[i], pos | p, neg | n, k+1)
    visit((1 << len(pats))-1, 0, 0, 0)
    if max_j is not None:
        return dict(moments), visits
    top = max(moments)
    spectrum = {k: sum((-1)**(j-k)*comb(j, k)*moments[j] for j in range(k, top+1))
                for k in range(top+1)}
    require(all(v >= 0 for v in spectrum.values()) and sum(spectrum.values()) == 1 << edge_count,
            "moment inversion positivity and total")
    return {k: v for k, v in spectrum.items() if v}, visits


def shannon_spectrum(pats, edge_count):
    """Separate decision recursion. Equal residual patterns keep multiplicity."""
    def canonical(items):
        c = Counter()
        for p, n, m in items:
            c[p, n] += m
        return tuple((p, n, m) for (p, n), m in sorted(c.items()))
    @lru_cache(None)
    def solve(state):
        shift = sum(m for p, n, m in state if p | n == 0)
        state = tuple((p, n, m) for p, n, m in state if p | n)
        active = 0
        for p, n, _ in state:
            active |= p | n
        if not active:
            return ((shift, 1),)
        bit = active & -active
        total = Counter()
        for yes in (False, True):
            child = canonical((p & ~bit, n & ~bit, m) for p, n, m in state
                              if not ((n if yes else p) & bit))
            child_active = 0
            for p, n, _ in child:
                child_active |= p | n
            factor = 1 << (active.bit_count()-1-child_active.bit_count())
            for count, mass in solve(child):
                total[count+shift] += factor*mass
        return tuple(sorted(total.items()))
    state = canonical((p, n, 1) for p, n in pats)
    union = 0
    for p, n in pats:
        union |= p | n
    factor = 1 << (edge_count-union.bit_count())
    result = {k: v*factor for k, v in solve(state)}
    return result, solve.cache_info().misses


def main():
    rows = []
    for square, orders in ((False, (3, 4, 5)), (True, (15, 16, 17, 23, 25))):
        for n in orders:
            edges, paths = graph_paths(n, square)
            pats = patterns(edges, paths)
            expected_square = {15: 2, 16: 2, 17: 2, 23: 6, 25: 20}
            if square:
                require(len(paths) == expected_square[n], "inherited square-sum control")
            shannon, nodes = shannon_spectrum(pats, len(edges))
            if square:
                intersection, visits = moment_spectrum(pats, len(edges))
                require(intersection == shannon, "independent spectrum agreement")
            else:
                literal = Counter(sum((bits & p) == p and (bits & m) == 0 for p, m in pats)
                                  for bits in range(1 << len(edges)))
                require(dict(literal) == shannon, "complete tournament exhaustive control")
                moments, visits = moment_spectrum(pats, len(edges), max_j=2)
                require(all(moments[j] == sum(comb(k, j)*v for k, v in literal.items())
                            for j in range(3)), "two-path moments")
            require(sum(shannon.values()) == 1 << len(edges), "orientation total")
            require(sum(k*v for k, v in shannon.items()) == len(paths)*(1 << (len(edges)-n+1)),
                    "first moment identity")
            row = {"guard": "square-sum" if square else "complete", "n": n,
                         "guard_edges": len(edges), "directed_guard_paths": len(paths),
                         "guard_orientation_histogram": dict(sorted(shannon.items())),
                         "outside_free_bits": comb(n, 2)-len(edges),
                         "intersection_nodes": visits, "shannon_states": nodes}
            if square:
                compatible = [(i, j) for i, (p, a) in enumerate(pats) for j, (q, b) in enumerate(pats)
                              if i < j and not (p & b or a & q)]
                row["compatible_pairs"] = compatible
                if n == 25:
                    require(len(compatible) == 4 and len({i for pair in compatible for i in pair}) == 8,
                            "four disjoint compatibility edges")
                    unions = [(pats[i][0] | pats[i][1] | pats[j][0] | pats[j][1]).bit_count()
                              for i, j in compatible]
                    require(unions == [26]*4, "joint arc cardinalities")
                    row["compatible_union_sizes"] = unions
                    row["compatible_path_witness"] = [[v+1 for v in paths[i]] for i in compatible[0]]
                else:
                    require(not compatible, "disjoint early path cubes")
            rows.append(row)
    Path(__file__).with_name("path_kernel_verification.json").write_text(json.dumps(rows, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(rows, sort_keys=True))


if __name__ == "__main__":
    main()
