"""Exact finite L/K quotients, with independent algebraic and physical controls.

Run from any directory; output is deterministic JSON on stdout. No assertion
statement is used, so every check survives python -O. The universe is all
three-coordinate residues, with triangle seed (2,0,1) and both legal letters.
"""
from array import array
from collections import Counter, deque
from itertools import permutations, product
import json


def require(test, message):
    if not test:
        raise RuntimeError(message)


def step(v, letter, modulus=None):
    x, b, d = v
    result = (x + 2*b, d - b, x + d) if letter == 'L' else (2*x + 2*b, 0, d)
    return tuple(t % modulus for t in result) if modulus else result


def lift(word):
    v = (2, 0, 1)
    for c in word:
        v = step(v, c)
    return v


def encode(v, m):
    x, b, d = v
    return (x*m + b)*m + d


def decode(k, m):
    xb, d = divmod(k, m)
    x, b = divmod(xb, m)
    return x, b, d


def reach(m):
    size = m**3
    parent = array('i', [-1])*size
    distance = array('H', [0])*size
    letter = bytearray(size)
    seed = encode((2 % m, 0, 1 % m), m)
    parent[seed] = seed
    queue = deque([seed])
    count = 1
    while queue:
        k = queue.popleft()
        v = decode(k, m)
        for c in 'LK':
            u = encode(step(v, c, m), m)
            if parent[u] == -1:
                parent[u] = k
                distance[u] = distance[k] + 1
                letter[u] = ord(c)
                queue.append(u)
                count += 1

    def word(v):
        k = encode(v, m)
        require(parent[k] >= 0, 'word requested for unreachable state')
        letters = []
        while k != seed:
            letters.append(chr(letter[k]))
            k = parent[k]
        return ''.join(reversed(letters))

    return parent, distance, count, word


def phase7(v):
    """Bit 0 = even L tail from b=0; bit 1 = odd L tail."""
    x, b, d = v
    a = (6*x + 3*b + 5*d) % 7
    u = (b - a) % 7
    w = (5*(d-b) - a) % 7
    norm = (6*u*u + 3*u*w + w*w) % 7
    allowed = (0, 3, 5, 6)
    return int((a*a-4*norm) % 7 in allowed) + 2*int((a*a+4*norm) % 7 in allowed)


def phase11(v):
    x, b, d = (z % 11 for z in v)
    if (2*x + 6*b + d) % 11:
        return 3
    a, c = (4*x-b) % 11, (b-3*x) % 11
    if a == c == 0:
        return 3
    if a == 0 or c == 0:
        return 0
    return 1 if (6*c*pow(a, -1, 11)) % 11 in (1, 3, 4, 5, 9) else 2


I = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
LM = ((1, 2, 0), (0, -1, 1), (1, 0, 1))


def matmul(a, b, m):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(3)) % m for j in range(3)) for i in range(3))


def matrix_order(m):
    a, n = I, 0
    while True:
        a, n = matmul(LM, a, m), n+1
        if a == I:
            return n


def explicit_phases(m, period):
    result = [0]*(m**3)
    for x, d in product(range(m), repeat=2):
        v = (x, 0, d)
        for j in range(period):
            result[encode(v, m)] |= 1 << (j % 2)
            v = step(v, 'L', m)
    return result


def count_paths(rows):
    n = len(rows)
    size = 1 << n
    incoming = [sum(1 << u for u in range(n) if rows[u] >> v & 1) for v in range(n)]
    dp = [[0]*n for _ in range(size)]
    totals = [1]+[0]*(size-1)
    for mask in range(1, size):
        for v in range(n):
            if not mask >> v & 1:
                continue
            rest = mask ^ (1 << v)
            if rest == 0:
                dp[mask][v] = 1
            else:
                dp[mask][v] = sum(dp[rest][u] for u in range(n) if rest & incoming[v] & (1 << u))
        totals[mask] = sum(dp[mask])
    return totals, dp[-1]


def physical(word):
    rows, root = [2, 4, 1], 0
    for letter in word:
        n = len(rows)
        rows[root] |= 1 << n
        rows.append(((1 << n)-1) ^ (1 << root))
        if letter == 'L':
            root = n
    return tuple(rows), root


def physical_state(rows, root):
    counts, _ = count_paths(rows)
    n = len(rows)
    rev = tuple(sum(1 << u for u in range(n) if rows[u] >> v & 1) for v in range(n))
    _, reverse_ends = count_paths(rev)
    h = counts[-1]
    d = counts[((1 << n)-1) ^ (1 << root)]
    return h-d, d-reverse_ends[root], d


def validate_physical(rows, root):
    n = len(rows)
    for u in range(n):
        require(not rows[u] >> u & 1, 'loop in physical tournament')
        for v in range(u):
            require(((rows[u] >> v) & 1) + ((rows[v] >> u) & 1) == 1, 'physical pair orientation failed')
        seen, todo = {u}, [u]
        while todo:
            a = todo.pop()
            for b in range(n):
                if rows[a] >> b & 1 and b not in seen:
                    seen.add(b)
                    todo.append(b)
        require(len(seen) == n, 'physical tournament lost strongness')
    require(sum(rows[u] >> root & 1 for u in range(n)) == 1, 'marked root lost indegree one')


def main():
    receipt = {'scope': 'Exact closure of all positive L/K words from the rooted triangle; residues retain x,b,d.'}
    phases = {}
    for m, classifier in ((7, phase7), (11, phase11)):
        period = matrix_order(m)
        actual = explicit_phases(m, period)
        expected = [classifier(decode(k, m)) for k in range(m**3)]
        require(actual == expected, f'algebraic parity classifier mismatch modulo {m}')
        phases[m] = actual
        receipt[f'local_{m}'] = {'L_order': period, 'phase_histogram': dict(sorted(Counter(actual).items()))}

    results = {}
    for m in (7, 11, 14, 77):
        parent, distance, count, word = reach(m)
        planes, h_ready, h_all = [], set(), set()
        for k in range(m**3):
            x, b, d = v = decode(k, m)
            if m in (7, 11):
                expected = bool(phases[m][k])
            elif m == 14:
                expected = x % 2 == 0 and d % 2 == 1 and bool(phase7(v) & (1 << (b % 2)))
            else:
                expected = bool(phase7(v) & phase11(v))
            require((parent[k] >= 0) == expected, f'closure classifier mismatch modulo {m} at {v}')
            if parent[k] >= 0:
                h_all.add((x+d) % m)
                if b == 1:
                    h_ready.add((x+d) % m)
            if b == 0 and (m != 14 or (x % 2 == 0 and d % 2 == 1)):
                require(parent[k] >= 0, f'plane reachability failure at {m}:{v}')
                planes.append(k)
        farthest_plane = max(planes, key=lambda k: (distance[k], k))
        pw = word(decode(farthest_plane, m))
        require(tuple(z % m for z in lift(pw)) == decode(farthest_plane, m), 'integer lift failure')
        results[m] = parent, word
        receipt[f'closure_{m}'] = {
            'states': count, 'universe': m**3, 'maximum_shortest_word': max(distance),
            'reachable_b_zero_plane': len(planes),
            'plane_maximum_shortest_word': max(distance[k] for k in planes),
            'plane_farthest_state': decode(farthest_plane, m), 'plane_farthest_witness': pw,
            'H_residues': len(h_all), 'H_residues_with_b_congruent_one': len(h_ready),
        }

    receipt['CRT_overcount'] = {'local_product': 337*1311, 'actual_77': 429807, 'incompatible_phase_pairs': 12000}
    phantom = (5, 1, 72)
    require(phase7(phantom) == 1 and phase11(phantom) == 2, 'phantom phase changed')
    require(results[77][0][encode(phantom, 77)] < 0, 'phantom became reachable')
    receipt['CRT_hostile'] = {'state_mod77': phantom, 'local_witnesses': {str(m): results[m][1](tuple(z % m for z in phantom)) for m in (7, 11)}}

    checked = 0
    false_ready_lengths = []
    for length in range(7):
        for letters in product('LK', repeat=length):
            w = ''.join(letters)
            rows, root = physical(w)
            validate_physical(rows, root)
            v = lift(w)
            require(physical_state(rows, root) == v, f'literal Hamilton count failed on {w}')
            require((v[1] == 1) == (w.endswith('L') and set(w[:-1]) <= {'K'}), 'unit-defect characterization failed')
            if v[1] % 7 == 1 and (v[0]+v[2]) % 7 == 0:
                false_ready_lengths.append(length)
            if length <= 3:
                h_perm = sum(all(rows[p[i]] >> p[i+1] & 1 for i in range(len(rows)-1)) for p in permutations(range(len(rows))))
                require(h_perm == v[0]+v[2], 'independent permutation count failed')
            checked += 1
    receipt['literal_graph_controls'] = {'all_words_through_births': 6, 'states': checked, 'permutation_counts_through_births': 3}
    require(min(false_ready_lengths) == 6, 'false-ready first length changed')

    w = 'KLLLKL'
    v = lift(w)
    child = lift(w+'LL')
    require(v == (38, 15, 53) and child[0]+child[2] == 303, 'false-ready hostile changed')
    receipt['false_ready_hostile'] = {'word': w, 'state': v, 'H': 91, 'mod7_defect': 1, 'actual_next_H': 303, 'desired_next_H': 275}
    pure_ready = {}
    for m, period in ((7, 3), (11, 10), (77, 30)):
        pure_ready[str(m)] = sorted({(pow(2, j, m)+1) % m for j in range(2, 2+period)})
    receipt['exact_unit_defect_H_residues'] = pure_ready

    loops = {}
    for m in (7, 11, 77):
        unseen, cycles = set(range(m)), []
        while unseen:
            seed = min(unseen)
            v, cycle = seed, []
            while v not in cycle:
                cycle.append(v)
                unseen.remove(v)
                v = 3*v % m
            require(v == seed, 'target scalar map is not a permutation')
            cycles.append(cycle)
        loops[str(m)] = {'cycles': cycles, 'length_histogram': dict(sorted(Counter(map(len, cycles)).items()))}
    receipt['r_target_loops'] = loops
    receipt['full_state_L_orders'] = {str(m): matrix_order(m) for m in (7, 11, 14, 77)}
    print(json.dumps(receipt, sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
