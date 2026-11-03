#!/usr/bin/env python3
"""Exact independent controls for Cayley–Dickson tournament bit circuits."""
import itertools
import json
import math
from collections import deque
from functools import lru_cache
from pathlib import Path


def check(value, label):
    if not value:
        raise RuntimeError(label)


@lru_cache(None)
def sign(n, a, b):
    if n == 0:
        return 0
    half = 1 << (n-1)
    u, r, v, t = a % half, a // half, b % half, b // half
    if r == 0 and t == 0:
        return sign(n-1, u, v)
    if r == 0 and t == 1:
        return sign(n-1, v, u)
    if r == 1 and t == 0:
        return sign(n-1, u, v) ^ bool(v)
    return sign(n-1, v, u) ^ bool(v) ^ 1


def conjugate(v):
    return (v[0],) + tuple(-x for x in v[1:])


def add(v, w):
    return tuple(x+y for x, y in zip(v, w))


def neg(v):
    return tuple(-x for x in v)


# Independent dense integer-vector doubling, with no use of sign().
@lru_cache(None)
def multiply(v, w):
    if len(v) == 1:
        return (v[0]*w[0],)
    h = len(v)//2
    a, b, c, d = v[:h], v[h:], w[:h], w[h:]
    return add(multiply(a, c), neg(multiply(conjugate(d), b))) + add(multiply(d, a), multiply(b, conjugate(c)))


def rank(rows):
    basis = {}
    for value in rows:
        while value:
            p = value.bit_length()-1
            if p in basis:
                value ^= basis[p]
            else:
                basis[p] = value
                break
    return len(basis), basis


def anf(values):
    out = list(values)
    n = (len(out)-1).bit_length()
    for j in range(n):
        for x in range(1 << n):
            if x & (1 << j):
                out[x] ^= out[x ^ (1 << j)]
    return [x for x in range(1 << n) if out[x]]


def full_permutation_group(n):
    N = 1 << n
    generators = [tuple((x ^ a) + N*(s ^ sign(n, a, x)) for s in range(2) for x in range(N)) for a in range(1, N)]
    identity = tuple(range(2*N))
    found, queue = {identity}, deque([identity])
    while queue:
        p = queue.popleft()
        for g in generators:
            q = tuple(g[p[i]] for i in range(2*N))
            if q not in found:
                found.add(q)
                queue.append(q)
    return len(found)


def table_controls(n):
    N = 1 << n
    d, half = n-1, N//2
    rows = [[sign(n, a, x) for x in range(N)] for a in range(N)]
    max_terms, degrees = 0, set()
    for a in range(N):
        check(rows[a][0] == 0, 'unit right identity')
        if a:
            check(all(rows[a][x] ^ rows[a][x ^ a] == 1 for x in range(N)), 'signed basis inverse')
        u, r = a % half, a // half
        for x in range(N):
            y, t = x % half, x // half
            c = 0 if not u else 1 ^ (y == 0) ^ (y == u)
            check(rows[a][x] == (sign(n-1, u, y) ^ (r and bool(y)) ^ (t and (c ^ r))), 'doubling row recurrence')
        terms = anf(rows[a])
        if a:
            degree = max(m.bit_count() for m in terms)
            check(degree == (1 if n == 1 else n-1), 'exact degree')
            degrees.add(degree)
            max_terms = max(max_terms, len(terms))
        for x in range(N):
            circuit_sign = 0
            for mask in terms:
                circuit_sign ^= (x & mask) == mask
            check(circuit_sign == rows[a][x], 'normal form circuit all inputs')
    if n < 3:
        return {'n': n, 'dimension': N, 'enumerated_group_order': full_permutation_group(n), 'nonzero_row_degrees': sorted(degrees), 'max_anf_terms': max_terms}
    defect_rows = []
    for a in range(N):
        for b in range(N):
            value = sum((rows[b][x] ^ rows[a][x ^ b] ^ rows[a ^ b][x]) << x for x in range(N))
            lo = value & ((1 << half)-1)
            hi = value >> half
            check(lo.bit_count() % 2 == 0, 'defect bottom parity')
            check((lo ^ hi) in (0, (1 << half)-1), 'defect constant last derivative')
            defect_rows.append(value)
    kernel_rank, basis = rank(defect_rows)
    check(kernel_rank == half, 'exact defect span')
    for value in basis.values():
        for j in range(n):
            shifted = sum(((value >> (x ^ (1 << j))) & 1) << x for x in range(N))
            check(rank(list(basis.values()) + [shifted])[0] == half, 'kernel translation invariance')
    for u in range(1, half):
        a, b = u, half
        got = [rows[b][x] ^ rows[a][x ^ b] ^ rows[a ^ b][x] for x in range(N)]
        check(got == [1 ^ (x % half == 0) ^ (x % half == u) for x in range(N)], 'kernel pair basis')
    edges = set()
    for s in range(2):
        for x in range(N):
            source = x + s*N
            for a in range(1, N):
                target = (x ^ a) + N*(s ^ rows[a][x])
                edges.add(tuple(sorted((source, target))))
    expected = {(i, j) for i in range(2*N) for j in range(i+1, 2*N) if i % N != j % N}
    check(edges == expected, 'cocktail-party graph')
    result = {'n': n, 'dimension': N, 'kernel_rank': kernel_rank, 'group_order_from_kernel_proof': 1 << (half+n), 'transport_memory_bits': half+n, 'nonzero_row_degrees': sorted(degrees), 'max_anf_terms': max_terms, 'cocktail_vertices': 2*N, 'cocktail_edges': len(edges)}
    if n <= 4:
        result['enumerated_group_order'] = full_permutation_group(n)
        check(result['enumerated_group_order'] == result['group_order_from_kernel_proof'], 'independent group enumeration')
    return result


def adjacency(n, mask):
    out = [0]*n
    for bit, (i, j) in enumerate(itertools.combinations(range(n), 2)):
        if mask >> bit & 1:
            out[i] |= 1 << j
        else:
            out[j] |= 1 << i
    return out


def strong(out):
    n = len(out)
    for root in range(n):
        seen, frontier = 1 << root, 1 << root
        while frontier:
            v = (frontier & -frontier).bit_length()-1
            frontier ^= 1 << v
            new = out[v] & ~seen
            seen |= new
            frontier |= new
        if seen != (1 << n)-1:
            return False
    return True


def profile(out):
    n, N = len(out), 1 << len(out)
    paths = [[0]*n for _ in range(N)]
    H = [0]*N
    for v in range(n):
        paths[1 << v][v] = 1
    for mask in range(1, N):
        for v, count in enumerate(paths[mask]):
            if count:
                choices = out[v] & ~mask
                while choices:
                    bit = choices & -choices
                    choices ^= bit
                    paths[mask | bit][bit.bit_length()-1] += count
        H[mask] = sum(paths[mask])
    partition = [[0]*(n+1) for _ in range(N)]
    partition[0][0] = 1
    for mask in range(1, N):
        root = mask & -mask
        sub = mask
        while sub:
            if sub & root:
                other = mask ^ sub
                for j in range(1, mask.bit_count()+1):
                    partition[mask][j] += H[sub] * partition[other][j-1]
            sub = (sub-1) & mask
    return partition[-1][1:]


def profile_by_cuts(out):
    n = len(out)
    totals = [0]*(n+1)
    for perm in itertools.permutations(range(n)):
        bad = 0
        for i in range(n-1):
            if not (out[perm[i]] >> perm[i+1] & 1):
                bad |= 1 << i
        for cut in range(1 << (n-1)):
            if cut & bad == bad:
                totals[cut.bit_count()+1] += 1
    for j in range(1, n+1):
        check(totals[j] % math.factorial(j) == 0, 'unordered paths divisibility')
    return [totals[j] // math.factorial(j) for j in range(1, n+1)]


def relabel(mask, perm):
    old = adjacency(len(perm), mask)
    new = [0]*len(perm)
    for i in range(len(perm)):
        for j in range(len(perm)):
            if old[i] >> j & 1:
                new[perm[i]] |= 1 << perm[j]
    return sum(((new[i] >> j) & 1) << bit for bit, (i, j) in enumerate(itertools.combinations(range(len(perm)), 2)))


def tournament_bits(bits):
    blocks = [adjacency(6, [153, 185][b]) for b in bits]
    n = 6*len(bits)
    out = [0]*n
    for j, block in enumerate(blocks):
        later = ((1 << n)-1) ^ ((1 << (6*(j+1)))-1)
        for i, value in enumerate(block):
            out[6*j+i] = (value << (6*j)) | later
    mask = sum(((out[i] >> j) & 1) << bit for bit, (i, j) in enumerate(itertools.combinations(range(n), 2)))
    return out, mask


def physical_controls():
    p = [29, 131, 162, 77, 15, 1]
    for mask in [153, 185]:
        out = adjacency(6, mask)
        check(strong(out), 'strong binary seed')
        check(profile(out) == p, 'seed DP profile')
        check(profile_by_cuts(out) == p, 'seed cuts profile')
    check(153 ^ 185 == 1 << 5, 'literal seed arc12')
    check(185 not in {relabel(153, perm) for perm in itertools.permutations(range(6))}, 'nonisomorphic binary seed')
    profiles12 = [profile(tournament_bits(bits)[0]) for bits in itertools.product(range(2), repeat=2)]
    check(all(p == profiles12[0] for p in profiles12), 'all four physical two-block full profiles')
    states = {}
    for s in range(2):
        for x in range(8):
            bits = [s] + [(x >> j) & 1 for j in range(3)]
            states[(s, x)] = tournament_bits(bits)[1]
    check(len(set(states.values())) == 16, 'sixteen literal states')
    total_macro_arcs = 0
    for state in states:
        s, x = state
        for a in range(1, 8):
            now, itinerary = (s, x), [state]
            if sign(3, a, x):
                now = (now[0] ^ 1, now[1]); itinerary.append(now)
            for j in range(3):
                if a >> j & 1:
                    now = (now[0], now[1] ^ (1 << j)); itinerary.append(now)
            check(now == (s ^ sign(3, a, x), x ^ a), 'macro endpoint')
            check(len(itinerary)-1 == a.bit_count()+sign(3, a, x), 'short macro exact length')
            for left, right in zip(itinerary, itinerary[1:]):
                check((states[left] ^ states[right]).bit_count() == 1, 'literal one-arc update')
                total_macro_arcs += 1
    associator = []
    for state in states:
        s, x = state
        for a, inverse in [(2, False), (1, False), (3, True)]:
            s ^= sign(3, a, x) ^ inverse
            x ^= a
        check((s, x) == (state[0] ^ bool(state[1] & 4), state[1]), 'octonion defect controlled sign')
        associator.append({'input': list(state), 'output': [int(s), x]})
    return {'fixed_masks': [153, 185], 'changed_arc': [1, 2], 'seed_profile': p, 'four_two_block_profiles': profiles12, 'octonion_states': 16, 'tournament_order': 24, 'common_hamilton_count': 29**4, 'macro_arc_checks': total_macro_arcs, 'associator': associator}


def main():
    checked_products = 0
    for n in range(1, 6):
        N = 1 << n
        units = [tuple(int(i == j) for i in range(N)) for j in range(N)]
        for a in range(N):
            for b in range(N):
                got = multiply(units[a], units[b])
                expected = tuple((1-2*sign(n, a, b)) if i == (a ^ b) else 0 for i in range(N))
                check(got == expected, 'independent integer-vector product')
                checked_products += 1
    result = {'status': 'all explicit checks passed', 'independent_basis_products': checked_products, 'levels': [table_controls(n) for n in range(1, 8)], 'physical_tournament_controls': physical_controls()}
    destination = Path(__file__).with_name('results.json')
    destination.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'status': result['status'], 'independent_basis_products': checked_products, 'transport_bits_n3_to_n7': [row['transport_memory_bits'] for row in result['levels'][2:]], 'physical_macro_arc_checks': result['physical_tournament_controls']['macro_arc_checks']}, sort_keys=True))


if __name__ == '__main__':
    main()
