#!/usr/bin/env python3
"""Independent exact guard-core decoder. Standard library; no producer imports."""
from collections import Counter, deque
from itertools import permutations
import json
import random
import sys

CHECKS = 0

def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)

def arc(t, a, b):
    return (t[a] >> b) & 1

def make(n, mask, fixed_path=False):
    rows = [0] * n
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            forward = (fixed_path and j == i + 1) or ((mask >> k) & 1)
            a, b = (i, j) if forward else (j, i)
            rows[a] |= 1 << b
            if not fixed_path or j != i + 1:
                k += 1
    return tuple(rows)

def relabel(t, p):
    return tuple(sum(arc(t, p[i], p[j]) << j for j in range(len(p))) for i in range(len(p)))

def legal(t, p):
    return all(arc(t, a, b) for a, b in zip(p, p[1:]))

def closure(rows, deleted=None):
    n = len(rows)
    available = ((1 << n) - 1) ^ (0 if deleted is None else 1 << deleted)
    reach = [((rows[i] | (1 << i)) & available) if available >> i & 1 else 0 for i in range(n)]
    for k in range(n):
        if available >> k & 1:
            for i in range(n):
                if reach[i] >> k & 1:
                    reach[i] |= reach[k]
    return tuple(reach)

def components(rows):
    n = len(rows)
    reach = closure(rows)
    unseen = set(range(n))
    cs = []
    while unseen:
        a = min(unseen)
        c = tuple(v for v in sorted(unseen) if reach[a] >> v & 1 and reach[v] >> a & 1)
        unseen.difference_update(c)
        cs.append(c)
    # This order is independently derived from reachability of the core.
    snapshot = tuple(cs)
    cs.sort(key=lambda c: -sum(bool(reach[c[0]] >> d[0] & 1) for d in snapshot))
    for i, c in enumerate(cs):
        for d in cs[i + 1:]:
            need(bool(reach[c[0]] >> d[0] & 1) and not (reach[d[0]] >> c[0] & 1), 'core SCC order not total')
    return tuple(cs)

def traceable_core(t):
    """DFS of directed Hamilton words; retain only adjacency support and count."""
    n = len(t)
    full = (1 << n) - 1
    result = [0] * n
    count = 0
    path = []
    def visit(used):
        nonlocal count
        if used == full:
            count += 1
            for a, b in zip(path, path[1:]):
                result[a] |= 1 << b
            return
        choices = (full ^ used) if not path else t[path[-1]] & (full ^ used)
        while choices:
            bit = choices & -choices
            choices ^= bit
            path.append(bit.bit_length() - 1)
            visit(used | bit)
            path.pop()
    visit(0)
    return tuple(result), count

def core_by_subset_dp(t):
    """Independent membership: concatenate a prefix ending a and suffix starting b."""
    n = len(t)
    full = (1 << n) - 1
    rev = tuple(sum(arc(t, j, i) << j for j in range(n)) for i in range(n))
    def endings(rows):
        dp = [0] * (1 << n)
        for i in range(n):
            dp[1 << i] = 1 << i
        for mask in range(1, full + 1):
            ends = dp[mask]
            while ends:
                bit = ends & -ends
                ends ^= bit
                a = bit.bit_length() - 1
                available = rows[a] & (full ^ mask)
                while available:
                    b = available & -available
                    available ^= b
                    dp[mask | b] |= b
        return dp
    end, start = endings(t), endings(rev)
    core = [0] * n
    for mask in range(1, full):
        for a in range(n):
            if end[mask] >> a & 1:
                core[a] |= t[a] & start[full ^ mask]
    return tuple(core)

def decode(core):
    """Reconstruct the tournament from its labelled traceable-arc core."""
    n = len(core)
    out = list(core)
    cs = components(core)
    witnesses = []
    for ci, c in enumerate(cs):
        for d in cs[ci + 1:]:
            for a in c:
                for b in d:
                    out[a] |= 1 << b
        if len(c) == 1:
            continue
        local = tuple(sum(arc(core, a, b) << j for j, b in enumerate(c)) for a in c)
        reach = [closure(local, z) for z in range(len(c))]
        for i, a in enumerate(c):
            for j in range(i + 1, len(c)):
                b = c[j]
                if arc(core, a, b) or arc(core, b, a):
                    continue
                directions = []
                for z in range(len(c)):
                    if z in (i, j):
                        continue
                    ij = bool(reach[z][i] >> j & 1)
                    ji = bool(reach[z][j] >> i & 1)
                    if ij != ji:
                        directions.append((ij, c[z]))
                need(bool(directions), 'missing pair has no single-vertex separator')
                need(len({v for v, _ in directions}) == 1, 'inconsistent separator directions')
                forward, z = directions[0]
                u, v = (a, b) if forward else (b, a)
                out[u] |= 1 << v
                witnesses.append((u, v, z))
    for a in range(n):
        need(not (out[a] >> a & 1), 'decoder loop')
        for b in range(a + 1, n):
            need(arc(out, a, b) + arc(out, b, a) == 1, 'decoder not a tournament')
    return tuple(out), witnesses

def local_moves(n):
    ps = set()
    base = list(range(n))
    for a in range(n):
        for b in range(a + 2, n):
            ps.add(tuple(base[:a] + base[a + 1:b + 1] + [a] + base[b + 1:]))
            ps.add(tuple(base[:a] + [b] + base[a:b] + base[b + 1:]))
    return tuple(sorted(ps))

def verify_one(t, independent=False, mutate=False):
    core, h = traceable_core(t)
    recovered, witnesses = decode(core)
    if mutate and len(t) > 1:
        changed = list(recovered)
        changed[0] ^= 1 << 1
        changed[1] ^= 1
        recovered = tuple(changed)
    need(recovered == t, 'decoded relation differs from original')
    need(h % 2 == 1, 'Redei odd parity control')
    if independent:
        need(core_by_subset_dp(t) == core, 'DFS/subset-DP core disagreement')
        need(traceable_core(core)[0] == core, 'core not idempotent')
        need(traceable_core(core)[1] == h, 'core changes full Hamilton language count')
    cs = components(t)
    need(cs == components(core), 'SCC partition/order mismatch')
    strong = len(cs) == 1
    if strong:
        need(closure(core) == closure(t), 'strong reachability mismatch')
        for z in range(len(t)):
            need(closure(core, z) == closure(t, z), 'strong one-deletion reachability mismatch')
    return h, len(witnesses), strong, sum(r.bit_count() for r in t) - sum(r.bit_count() for r in core)

def main():
    if '--mutate' in sys.argv:
        verify_one(make(3, 0, True), mutate=True)
        raise RuntimeError('negative control was not detected')
    summary = []
    for n in range(1, 8):
        count = 1 << max(0, (n - 1) * (n - 2) // 2)
        maxima = [0, 0, 0]
        strong = 0
        for mask in range(count):
            h, restored, is_strong, absent = verify_one(make(n, mask, True), independent=n <= 5)
            strong += is_strong
            maxima = [max(maxima[0], h), max(maxima[1], restored), max(maxima[2], absent)]
        summary.append({'n': n, 'fixed_path_presentations': count, 'strong_presentations': strong,
                        'maximum_H': maxima[0], 'maximum_missing_internal_pairs': maxima[1],
                        'maximum_all_missing_pairs': maxima[2], 'decoder': 'PASS'})
    labelled = []
    for n in range(1, 6):
        total = 1 << (n * (n - 1) // 2)
        for mask in range(total):
            verify_one(make(n, mask), independent=True)
        labelled.append({'n': n, 'all_labelled_tournaments': total, 'decoder': 'PASS'})
    rng = random.Random(700314159)
    random_controls = []
    for n in range(8, 12):
        total = 40
        strong = 0
        for _ in range(total):
            t = make(n, rng.getrandbits(n * (n - 1) // 2))
            h, _, s, _ = verify_one(t, independent=True)
            strong += s
        random_controls.append({'n': n, 'seeded_samples': total, 'strong_samples': strong})
    local = []
    for n in range(1, 5):
        ps = local_moves(n)
        cells = Counter(tuple(legal(make(n, mask, True), p) for p in ps)
                        for mask in range(1 << max(0, (n - 1) * (n - 2) // 2)))
        local.append({'n': n, 'states': sum(cells.values()), 'one_step_elementary_guard_cells': len(cells)})
    x, y = make(4, 0, True), make(4, 5, True)
    ps = local_moves(4)
    need(tuple(legal(x, p) for p in ps) == tuple(legal(y, p) for p in ps), 'local hostile guards differ')
    p = (1, 3, 0, 2)
    need((legal(x, p), legal(y, p)) == (False, True), 'whole-permutation hostile not separated')
    trace = ((1, 2, 3, 0), (0, 2, 3, 1))
    need(legal(x, trace[0]) and legal(y, trace[0]), 'first local move not common')
    a, b = relabel(x, trace[0]), relabel(y, trace[0])
    need((legal(a, trace[1]), legal(b, trace[1])) == (False, True), 'depth-two local trace not separating')
    # Strong traceability is genuinely incomplete: ordinary mask10 at n5.
    hostile = make(5, 10)
    core, h = traceable_core(hostile)
    need(arc(hostile, 4, 1) and not arc(core, 4, 1), 'nontraceable arc hostile lost')
    decoded, witnesses = decode(core)
    need(decoded == hostile and (4, 1) in {(u, v) for u, v, _ in witnesses}, 'hidden arc not recovered')
    receipt = {'status': 'PASS', 'checks': CHECKS, 'fixed_path_exhaustion': summary,
               'labelled_exhaustion': labelled, 'seed': 700314159, 'seeded_random_controls': random_controls,
               'local_guard_partition': local,
               'minimal_local_hostile': {'n': 4, 'masks': [0, 5], 'rows': [x, y],
                   'full_permutation': p, 'legal': [False, True], 'shortest_local_trace': trace,
                   'local_guard_vector': [legal(x, q) for q in ps], 'named_local_moves': ps},
               'nontraceable_strong_hostile': {'n': 5, 'ordinary_mask': 10, 'rows': hostile,
                   'core': core, 'H': h, 'restoration_witnesses': witnesses}}
    print(json.dumps(receipt, sort_keys=True, indent=2))

if __name__ == '__main__':
    main()
