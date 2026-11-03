"""Independent two-live-port audit; no repository imports or target scan.

Uses component chains for transfer and literal permutations/cuts for graphs.
Run normally and under -O; writes the adjacent deterministic JSON receipt.
"""
from collections import Counter, deque
from functools import lru_cache
from itertools import permutations, product
from pathlib import Path
import hashlib
import json
import sympy as sp


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def birth(rows, root):
    n = len(rows)
    return tuple(r | ((1 << n) if i == root else 0) for i, r in enumerate(rows)) + (((1 << n)-1) ^ (1 << root),)


def literal_covers(rows, ports):
    n = len(rows)
    actual = set()
    rank = {v: i for i, v in enumerate(ports)}
    for p in permutations(range(n)):
        for mask in range(1 << (n-1)):
            pieces, i, ok = [], 0, True
            for j in range(n-1):
                if mask >> j & 1:
                    pieces.append(p[i:j+1])
                    i = j+1
                elif not rows[p[j]] >> p[j+1] & 1:
                    ok = False
                    break
            if ok:
                actual.add(tuple(sorted(pieces + [p[i:]])))
    answer = Counter()
    for cover in actual:
        sig = tuple(sorted((rank.get(p[0], -1), rank.get(p[-1], -1)) for p in cover))
        if sum(e == -1 for s, e in sig) <= 1:
            answer[sig] += 1
    return answer


def transfer(sig, root, drop):
    """Build directed chains on old component nodes plus the newborn node."""
    k = len(sig)
    pairs = list(sig) + [(2, 2)]
    output = Counter()
    # Fixed-slot convention: a newborn replaces the discarded old slot.
    rename = {0: 0, 1: 1}
    if drop != 2:
        del rename[drop]
        rename[2] = drop
    incoming = [None] + [i for i, (s, e) in enumerate(sig) if e == root]
    outgoing = [None] + [i for i, (s, e) in enumerate(sig) if s != root]
    for pred, succ in product(incoming, outgoing):
        if pred is not None and pred == succ:
            continue
        nxt, back = {}, {}
        if pred is not None:
            nxt[pred], back[k] = k, pred
        if succ is not None:
            nxt[k], back[succ] = succ, k
        parts = []
        for i in range(k+1):
            if i in back:
                continue
            j = i
            while j in nxt:
                j = nxt[j]
            parts.append((rename.get(pairs[i][0], -1), rename.get(pairs[j][1], -1)))
        if sum(e == -1 for s, e in parts) <= 1:
            output[tuple(sorted(parts))] += 1
    return output


def advance(vec, root, drop):
    out = Counter()
    for sig, a in vec.items():
        for target, b in transfer(sig, root, drop).items():
            out[target] += a*b
    return out


def width(roots):
    return max([0] + [len({r for r in roots[t:] if r < 3+t}) for t in range(len(roots))])


def literal_h(rows):
    return sum(all(rows[p[i]] >> p[i+1] & 1 for i in range(len(p)-1))
               for p in permutations(range(len(rows))))


def recursive_h(rows):
    @lru_cache(None)
    def finish(last, remaining):
        if not remaining:
            return 1
        total = 0
        choices = rows[last] & remaining
        while choices:
            bit = choices & -choices
            choices -= bit
            total += finish(bit.bit_length()-1, remaining ^ bit)
        return total
    return sum(finish(i, ((1 << len(rows))-1) ^ (1 << i)) for i in range(len(rows)))


def decode(program):
    rows, ports = (2, 4, 1), [0, 1]
    for letter in program:
        root, drop = divmod(letter, 3)
        n = len(rows)
        rows = birth(rows, ports[root])
        if drop < 2:
            ports[drop] = n
    return rows


def modular_rank(rows, p=1000003):
    basis = {}
    for raw in rows:
        row = [int(x) % p for x in raw]
        for i, b in sorted(basis.items()):
            c = row[i]
            if c:
                row = [(x-c*y) % p for x, y in zip(row, b)]
        pivot = next((i for i, x in enumerate(row) if x), None)
        if pivot is not None:
            inv = pow(row[pivot], -1, p)
            basis[pivot] = [(x*inv) % p for x in row]
    return len(basis)


def common_static_subseries(matrices, u, f):
    """Exact finite-field comparison of the proper fixed-port sublanguage."""
    prime = 1000003

    def dimensions(ms, initial, final):
        def close(v, mm):
            rows = [v]
            for r in rows:
                for m in mm:
                    new = [int(x) % prime for x in list(sp.Matrix(1, len(r), r)*m)]
                    if modular_rank(rows+[new]) > len(rows):
                        rows.append(new)
            return rows
        reach = close(initial, ms)
        obs = close(final, [m.T for m in ms])
        hankel = [[sum(x*y for x, y in zip(r, c)) % prime for c in obs] for r in reach]
        return {'reachable': len(reach), 'observable': len(obs), 'Hankel': modular_rank(hankel)}

    new = dimensions([matrices[i] for i in (2, 5)], list(u), list(f))
    # Rebuild the old two-static-feedback gap automaton from its stated
    # empty/good-last/bad-last gap semantics; no old module is imported.
    gap_states = list(product(range(3), range(3), range(2)))
    ids = {s: i for i, s in enumerate(gap_states)}
    old = [sp.zeros(36) for _ in range(4)]
    start, finish = sp.zeros(1, 36), sp.zeros(36, 1)
    for block, p in enumerate(permutations(range(2))):
        off = 18*block
        start[off] = 1
        for s in gap_states:
            i = off+ids[s]
            if 2 not in s and (s[1] or p[0] == 0):
                finish[i] = 1
            for letter in range(4):
                for gap in range(3):
                    if s[gap] == 0 and gap and not (letter >> p[gap-1] & 1):
                        continue
                    t = list(s)
                    t[gap] = int(gap != 2 and bool(letter >> p[gap] & 1))+1
                    old[letter][i, off+ids[tuple(t)]] += 1
    # New chronology reads the old static chain backwards and ends with its
    # original third C3 vertex, whose static incidence letter is2.
    previous = dimensions([old[i].T for i in (1, 2)], list(old[2]*finish), list(start))
    need(new == {'reachable': 13, 'observable': 15, 'Hankel': 13}, 'fixed-port finite-field dimensions')
    need(previous == {'reachable': 15, 'observable': 20, 'Hankel': 13}, 'same static subseries dimensions')
    graph_checks = 0
    for length in range(4):
        for roots in product(range(2), repeat=length):
            new_rows = decode([3*r+2 for r in roots])
            letters = [1 << r for r in reversed(roots)] + [2]
            n = 2+len(letters)
            old_rows = [2, 0]+[0]*len(letters)
            for i, letter in enumerate(letters, 2):
                old_rows[i] |= sum(1 << j for j in range(i+1, n))
                for mark in range(2):
                    if letter >> mark & 1:
                        old_rows[mark] |= 1 << i
                    else:
                        old_rows[i] |= 1 << mark
            rename = [0, 1]+list(reversed(range(2, n)))
            transported = tuple(sum(1 << j for j in range(n) if new_rows[rename[i]] >> rename[j] & 1)
                                for i in range(n))
            need(transported == tuple(old_rows), 'literal common-language graph map')
            graph_checks += 1
    return {'field_prime': prime, 'new_fixed_port_dimensions': new,
            'old_reversed_mixed_chain_dimensions': previous, 'literal_graph_map_controls': graph_checks,
            'scope': 'A genuine proper common graph sublanguage and its exact rank13 over F_1000003. '
                     'This does not identify the full old/new languages or explain both full rational ranks17.'}


def rank_audit(signatures, initial):
    """Rebuild all matrices; check actual Hankel and exact observable closure."""
    path = Path(__file__).with_name('live_port.json')
    source_bytes = path.read_bytes()
    source = json.loads(source_bytes)
    states = tuple(sorted(signatures))
    index = {s: i for i, s in enumerate(states)}
    matrices = []
    sparse = []
    for root, drop in product(range(2), range(3)):
        m = sp.zeros(len(states))
        terms = []
        for i, sig in enumerate(states):
            for target, c in transfer(sig, root, drop).items():
                j = index[target]
                m[i, j] += c
                terms.append((i, j, c))
        matrices.append(m)
        sparse.append(sorted(terms))
    need(states == tuple(tuple(map(tuple, s)) for s in source['states']), 'producer support match')
    need(sparse == [sorted(map(tuple, m)) for m in source['matrices']], 'all six independently rebuilt matrices')
    u = sp.Matrix(1, len(states), [initial[s] for s in states])
    f = sp.Matrix([int(len(s) == 1) for s in states])
    need(list(u) == source['seed'] and list(f) == source['observer'], 'actual seed and observer')
    prefixes, suffixes = source['rank']['prefixes'], source['rank']['suffixes']
    reachable = []
    for word in prefixes:
        r = u
        for letter in word:
            r = r*matrices[letter]
        reachable.append(r)
    observations = []
    for word in suffixes:
        c = f
        for letter in reversed(word):
            c = matrices[letter]*c
        observations.append(c)
    R = sp.Matrix.vstack(*reachable)
    C = sp.Matrix.hstack(*observations)
    H = R*C
    need(H.tolist() == source['rank']['hankel'], 'frozen Hankel entries')
    # A different counting direction from the producer's forward subset-last DP.
    graph_counts = {}
    max_order = 0
    for i, prefix in enumerate(prefixes):
        for j, suffix in enumerate(suffixes):
            rows = decode(prefix+suffix)
            max_order = max(max_order, len(rows))
            if rows not in graph_counts:
                graph_counts[rows] = recursive_h(rows)
            need(graph_counts[rows] == H[i, j], 'actual graph Hankel entry')
    rr, cr, hr = (modular_rank(a.tolist()) for a in (R, C, H))
    need((rr, cr, hr) == (24, 17, 17), 'independent modular rank lower bounds')
    # Exact upper bound: all17 actual continuation columns span an invariant
    # rational column space containing f. Verify every entry, not a mod-p fit.
    row_pivots = list(C.T.rref()[1])
    B = C[row_pivots, :]
    need(len(row_pivots) == 17 and B.det() != 0, 'exact observer basis')
    inverse = B.inv()
    for m in matrices:
        image = m*C
        q = inverse*image[row_pivots, :]
        need(image == C*q, 'exact rational observable closure')
    need(C[:, 0] == f and suffixes[0] == [], 'observer included')
    hankel_rows = list(H.T.rref()[1])
    minor = H[hankel_rows, :]
    return {
        'producer_receipt_sha256': hashlib.sha256(source_bytes).hexdigest(),
        'six_matrices_rebuilt_independently': True,
        'actual_Hankel_entries_recounted': H.rows*H.cols,
        'distinct_actual_graphs_recounted': len(graph_counts), 'maximum_graph_order': max_order,
        'prime_for_lower_bounds': 1000003,
        'reachable_rank_mod_prime': rr, 'observable_rank_mod_prime': cr, 'Hankel_rank_mod_prime': hr,
        'exact_observable_closure_dimension': 17,
        'actual_Hankel_minor_row_indices': hankel_rows, 'actual_Hankel_minor_determinant': str(minor.det()),
        'common_static_subseries': common_static_subseries(matrices, u, f),
        'consequence': 'Exact rational reachable rank24, observable rank17, and scalar-series Hankel rank17. '
                       'No support saturation or all-target claim is inferred.',
    }


def main():
    seed, ports = (2, 4, 1), (0, 1)
    vec = literal_covers(seed, ports)
    seen, todo = set(vec), deque(vec)
    while todo:
        sig = todo.popleft()
        for root, drop in product(range(2), range(3)):
            for target in transfer(sig, root, drop):
                if target not in seen:
                    seen.add(target)
                    todo.append(target)
    need(len(seen) == 24, 'complete reachable signature closure')

    checks = 0
    layers = [(seed, ports, vec, (), ())]
    level_counts = []
    for depth in range(3):
        level_counts.append(len(layers))
        next_layer = []
        for rows, ports, v, roots, program in layers:
            need(literal_covers(rows, ports) == v, 'full literal cover vector')
            need(width(roots) <= 2, 'actual future-live width')
            need(sum(n for s, n in v.items() if len(s) == 1) == literal_h(rows), 'literal Hamilton count')
            checks += 1
            if depth == 2:
                continue
            for root, drop in product(range(2), range(3)):
                new = birth(rows, ports[root])
                newports = list(ports)
                if drop < 2:
                    newports[drop] = len(rows)
                next_layer.append((new, tuple(newports), advance(v, root, drop),
                                   roots + (ports[root],), program + (3*root+drop,)))
        layers = next_layer
    need(checks == 43 and level_counts == [1, 6, 36], 'complete short program universe')

    oriented_pairs = [(0, 1), (1, 2), (2, 0)]
    rotations = [literal_covers(seed, p) for p in oriented_pairs]
    need(all(v == rotations[0] for v in rotations), 'oriented C3 seed symmetry')

    controls = []
    for roots in [(0, 1, 2), (0, 0, 3, 3, 0)]:
        rows, degrees = seed, []
        for root in roots:
            degrees.append(len(rows)-1-rows[root].bit_count())
            rows = birth(rows, root)
        controls.append({'roots': roots, 'future_width': width(roots),
                         'root_indegrees': degrees, 'H_literal': literal_h(rows), 'rows': rows})
    need([(r['future_width'], r['H_literal']) for r in controls] == [(3, 37), (2, 127)],
         'actual width/degree separating controls')

    # Retain ports0 and3 after the first birth, use root3 repeatedly, then0.
    unbounded = []
    for k in [1, 2, 5, 12]:
        roots = (0,) + (3,)*k + (0,)
        rows = seed
        degrees = []
        for root in roots:
            degrees.append(len(rows)-1-rows[root].bit_count())
            rows = birth(rows, root)
        need(width(roots) == 2 and degrees[-1] == k+1, 'unbounded root degree in width2 family')
        unbounded.append({'k': k, 'future_width': width(roots), 'last_root_indegree': degrees[-1]})

    result = {
        'status': 'PASS',
        'scope': 'Independent component-chain transfer and literal graph covers. All six-action '
                 'words of length<=2; complete signature support closure; two selected literal '
                 'Hamilton controls; actual-Hankel recount and exact rank closure. No target saturation through999.',
        'letter_convention': '3*root_slot+drop; drop0/1 replaces that slot; drop2 discards newborn',
        'reachable_signatures': sorted(seen), 'signature_count': len(seen),
        'full_literal_vector_checks': checks, 'program_counts_by_length': level_counts,
        'oriented_seed_pairs': oriented_pairs, 'oriented_seed_vectors_equal': True,
        'width_degree_controls': controls, 'unbounded_degree_family_controls': unbounded,
        'rank_audit': rank_audit(seen, vec),
    }
    path = Path(__file__).with_name('live_port_independent_audit.json')
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': 'PASS', 'signature_count': len(seen),
                      'literal_full_vectors': checks, 'H_controls': [r['H_literal'] for r in controls]}))


if __name__ == '__main__':
    main()
