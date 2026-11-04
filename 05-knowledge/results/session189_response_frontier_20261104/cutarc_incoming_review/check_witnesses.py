#!/usr/bin/env python3
"""Independent permutation enumeration, no imports from incoming checker."""
import itertools
import json

def require(test, message):
    if not test:
        raise RuntimeError(message)

def paths(a, vertices):
    return [p for p in itertools.permutations(vertices)
            if all(a[u] >> v & 1 for u, v in zip(p, p[1:]))]

def strong(a, vertices):
    for u in vertices:
        reached = {u}
        while True:
            more = reached | {v for w in reached for v in vertices if a[w] >> v & 1}
            if more == reached:
                break
            reached = more
        if reached != set(vertices):
            return False
    return True

def phi(a, x):
    return tuple(sum((1 << v) for v in range(len(a)) if v != u and
                     ((a[u] >> v & 1) if x in (u, v) else (a[v] >> u & 1)))
                 for u in range(len(a)))

def responses(a, z, q=None):
    vertices = tuple(range(len(a)))
    total = 0
    for bits in range(1 << len(a)):
        A = tuple(v for v in vertices if bits >> v & 1)
        B = tuple(v for v in vertices if not bits >> v & 1)
        if z not in B or (q is not None and q not in A):
            continue
        left = paths(a, A)
        if q is not None:
            left = [p for p in left if p[-1] == q]
        right = [p for p in paths(a, B) if p[0] == z]
        total += len(left) * len(right)
    return total

out = []
for a, expected_r, expected_h in [
    ((6,60,8,1,13,29), (60,5,8,50,38,6), 33),
    ((22,12,8,1,46,15), (44,53,8,50,5,20), 29),
]:
    V = tuple(range(6))
    require(all(((a[u] >> v) ^ (a[v] >> u)) & 1 for u in V for v in V if u != v), 'tournament')
    hp = paths(a, V)
    cycles = [p for p in hp if p[0] == 0 and a[p[-1]] & 1]
    nc = [v for v in V if strong(a, tuple(w for w in V if w != v))]
    exits = [v for v in V if a[v].bit_count() == 1]
    entries = [v for v in V if a[v].bit_count() == 4]
    require(len(cycles) == 1 and len(hp) == 15 and len(nc) == 3, 'UH invariants')
    require(len(entries) == len(exits) == 2 and set(exits) & set(nc) == {2}, 'canonical mark')
    cuts = 0
    for u in V:
        for v in V:
            if a[u] >> v & 1:
                b = list(a)
                b[u] ^= 1 << v
                cuts += not strong(b, V)
    require(cuts == 4, 'cut arcs')
    r = phi(a, 2)
    require(r == expected_r and phi(r, 2) == a, 'named Phi')
    rhp = paths(r, V)
    sigma = sum(p[-1] == 3 for p in paths(r, (0,1,3,4,5)))
    require(len(rhp) == expected_h and sigma == 1 and strong(r,V), 'ready counts')
    require(responses(r,2) == expected_h + 1, 'F_x identity')
    # Test the general update even for non-admissible p.
    for p in V:
        if p == 2:
            continue
        child = tuple(row | ((1 << 6) if u not in (2,p) else 0) for u,row in enumerate(r)) + ((1 << 2) | (1 << p),)
        actual = len(paths(child, tuple(range(7))))
        predicted = (expected_h - sum(t[-1] in (2,p) for t in rhp)
                     + responses(r,2) + responses(r,p)
                     - responses(r,p,2) - responses(r,2,p))
        require(actual == predicted, 'path law')
    out.append({'UH_rows':a, 'ready_rows':r, 'unique_cycle':cycles[0], 'H_UH':len(hp),
                'H_ready':len(rhp), 'sigma':sigma, 'cut_arcs':cuts, 'noncritical':nc,
                'all_five_distinct_p_updates_checked':True})
print(json.dumps(out, indent=2))
