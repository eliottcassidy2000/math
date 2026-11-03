"""Independent twelve-command/common-observer audit; no producer imports.

No target-support scan. Source JSON supplies proposed certificates only.
"""
from collections import Counter, deque
from functools import lru_cache
from itertools import combinations_with_replacement, permutations, product
from pathlib import Path
import hashlib
import json
import sympy as sp


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def literal_covers(rows, ports):
    n = len(rows)
    actual = set()
    for word in permutations(range(n)):
        for cut in range(1 << (n-1)):
            parts, at = [], 0
            for i in range(n-1):
                if cut >> i & 1:
                    parts.append(word[at:i+1])
                    at = i+1
                elif not rows[word[i]] >> word[i+1] & 1:
                    break
            else:
                actual.add(tuple(sorted(parts+[word[at:]])))
    slot = {v: i for i, v in enumerate(ports)}
    out = Counter()
    for cover in actual:
        sig = tuple(sorted((slot.get(p[0], -1), slot.get(p[-1], -1)) for p in cover))
        if sum(e == -1 for s, e in sig) <= 1:
            out[sig] += 1
    return out


def step(sig, mask, drop):
    # Join a graph of component nodes, then read its chain endpoints.
    k = len(sig)
    ends = list(sig)+[(2, 2)]
    rename = {0: 0, 1: 1}
    if drop != 2:
        del rename[drop]
        rename[2] = drop
    ins = [None]+[i for i, (s, e) in enumerate(sig) if e >= 0 and mask >> e & 1]
    outs = [None]+[i for i, (s, e) in enumerate(sig) if s == -1 or not mask >> s & 1]
    out = Counter()
    for predecessor, successor in product(ins, outs):
        if predecessor is not None and predecessor == successor:
            continue
        nxt, prev = {}, {}
        if predecessor is not None:
            nxt[predecessor], prev[k] = k, predecessor
        if successor is not None:
            nxt[k], prev[successor] = successor, k
        parts = []
        for i in range(k+1):
            if i in prev:
                continue
            j = i
            while j in nxt:
                j = nxt[j]
            parts.append((rename.get(ends[i][0], -1), rename.get(ends[j][1], -1)))
        if sum(e == -1 for s, e in parts) <= 1:
            out[tuple(sorted(parts))] += 1
    return out


def advance(vec, mask, drop):
    out = Counter()
    for sig, coefficient in vec.items():
        for target, multiplicity in step(sig, mask, drop).items():
            out[target] += coefficient*multiplicity
    return out


def birth(rows, ports, mask, drop):
    n = len(rows)
    incoming = {ports[i] for i in range(2) if mask >> i & 1}
    child = tuple(r | ((1 << n) if i in incoming else 0) for i, r in enumerate(rows))
    child += (sum(1 << i for i in range(n) if i not in incoming),)
    ports = list(ports)
    if drop < 2:
        ports[drop] = n
    return child, tuple(ports)


def decode(word):
    rows, ports = (2, 0), (0, 1)
    for a in word:
        rows, ports = birth(rows, ports, *divmod(a, 3))
    return rows, ports


def H(rows):
    @lru_cache(None)
    def completion(last, remaining):
        if not remaining:
            return 1
        total, choices = 0, rows[last] & remaining
        while choices:
            bit = choices & -choices
            choices -= bit
            total += completion(bit.bit_length()-1, remaining ^ bit)
        return total
    return sum(completion(i, ((1 << len(rows))-1) ^ (1 << i)) for i in range(len(rows)))


def modular_rank(rows, p=1000003):
    basis = {}
    for raw in rows:
        v = [int(x) % p for x in raw]
        for pivot, b in sorted(basis.items()):
            c = v[pivot]
            if c:
                v = [(x-c*y) % p for x, y in zip(v, b)]
        pivot = next((i for i, x in enumerate(v) if x), None)
        if pivot is not None:
            z = pow(v[pivot], -1, p)
            basis[pivot] = [(x*z) % p for x in v]
    return len(basis)


def abstract_signatures():
    pairs = list(product((-1, 0, 1), repeat=2))
    allowed = set()
    for length in range(4):
        for sig in combinations_with_replacement(pairs, length):
            if sum(e == -1 for s, e in sig) > 1:
                continue
            valid = True
            for port in range(2):
                occurrences = sum(s == port for s, e in sig)+sum(e == port for s, e in sig)
                if occurrences > 1 and not (occurrences == 2 and (port, port) in sig):
                    valid = False
            if valid:
                allowed.add(sig)
    need(len(allowed) == 28, 'complete abstract two-port signature universe including degenerate covers')
    return allowed-{(), ((0, 0),), ((1, 1),)}


def main():
    here = Path(__file__).resolve().parent
    raw = (here/'two_port_envelope.json').read_bytes()
    source = json.loads(raw)
    initial = literal_covers((2, 0), (0, 1))
    def close(letters):
        found, todo = set(initial), deque(initial)
        while todo:
            s = todo.popleft()
            for mask, drop in letters:
                for t in step(s, mask, drop):
                    if t not in found:
                        found.add(t)
                        todo.append(t)
        return found
    signatures = close(list(product(range(4), range(3))))
    static = close([(a, 2) for a in range(4)])
    need(signatures == static == abstract_signatures(), 'all25 abstract signatures reached statically')
    states = tuple(sorted(signatures))
    ids = {s: i for i, s in enumerate(states)}
    matrices, sparse = [], []
    for mask, drop in product(range(4), range(3)):
        m, terms = sp.zeros(25), []
        for i, s in enumerate(states):
            for t, c in step(s, mask, drop).items():
                m[i, ids[t]] += c
                terms.append((i, ids[t], c))
        matrices.append(m)
        sparse.append(sorted(terms))
    need(states == tuple(tuple(map(tuple, s)) for s in source['states']), 'producer support')
    need(sparse == [sorted(map(tuple, m)) for m in source['matrices']], 'all12 independently rebuilt matrices')
    u = sp.Matrix(1, 25, [initial[s] for s in states])
    f = sp.Matrix([int(len(s) == 1) for s in states])
    need(list(u) == source['seed'] and list(f) == source['observer'], 'actual seed and acceptance')

    short, full_vectors = 0, 0
    todo = [((2, 0), (0, 1), initial, ())]
    while todo:
        rows, ports, vector, word = todo.pop()
        literal = sum(all(rows[p[i]] >> p[i+1] & 1 for i in range(len(p)-1))
                      for p in permutations(range(len(rows))))
        need(sum(c for s, c in vector.items() if len(s) == 1) == literal, 'literal Hamilton word count')
        short += 1
        if len(word) <= 2:
            need(vector == literal_covers(rows, ports), 'full literal cover distribution')
            full_vectors += 1
        if len(word) < 3:
            for a in range(12):
                mask, drop = divmod(a, 3)
                child, names = birth(rows, ports, mask, drop)
                todo.append((child, names, advance(vector, mask, drop), word+(a,)))
    need((short, full_vectors) == (1885, 157), 'complete short-word audit scopes')

    suffixes = source['rank']['suffixes']
    columns = []
    for word in suffixes:
        c = f
        for a in reversed(word):
            c = matrices[3*a+2]*c
        columns.append(c)
    C = sp.Matrix.hstack(*columns)
    need(C.cols == modular_rank(C.tolist()) == 17 and C[:, 0] == f, '17 independent static observation columns')
    for i, encoded in enumerate(source['rank']['query_matrices']):
        A = sp.Matrix([[sp.Rational(x) for x in row] for row in encoded])
        need(matrices[i]*C == C*A, 'all25 coefficients of common static-query update')
    need(list(u*C) == source['rank']['query_seed'], 'static-query seed')

    # Full rank of the STATIC reachable space: exact words chosen by independent
    # finite-field elimination give25 independent actual rows over Q as well.
    static_rows, static_words = [list(u)], [()]
    for r, word in zip(static_rows, static_words):
        for mask in range(4):
            v = list(sp.Matrix(1, 25, r)*matrices[3*mask+2])
            if modular_rank(static_rows+[v]) > len(static_rows):
                static_rows.append(v)
                static_words.append(word+(3*mask+2,))
    need(len(static_rows) == 25, 'static reachable rank25')
    static_hankel = sp.Matrix(static_rows)*C
    need(modular_rank(static_hankel.tolist()) == 17, 'static Hankel lower bound17')

    # Recount the producer's full-prefix/static-suffix425 entries on graphs.
    count_cache, entries, max_order = {}, 0, 0
    for i, prefix in enumerate(source['rank']['prefixes']):
        for j, suffix in enumerate(suffixes):
            rows, _ = decode(prefix+[3*a+2 for a in suffix])
            if rows not in count_cache:
                count_cache[rows] = H(rows)
            need(count_cache[rows] == source['rank']['hankel'][i][j], 'actual common Hankel entry')
            entries += 1
            max_order = max(max_order, len(rows))
    need(entries == 425, 'all common Hankel entries')

    old_graph_checks = 0
    for length in range(4):
        for word in product(range(4), repeat=length):
            rows, ports = decode([3*a+2 for a in word])
            old_word = tuple(reversed(word))
            n, old_rows = 2+length, [2, 0]+[0]*length
            for i, a in enumerate(old_word, 2):
                old_rows[i] |= sum(1 << j for j in range(i+1, n))
                for mark in range(2):
                    if a >> mark & 1:
                        old_rows[mark] |= 1 << i
                    else:
                        old_rows[i] |= 1 << mark
            rename = [0, 1]+list(reversed(range(2, n)))
            mapped = tuple(sum(1 << j for j in range(n) if rows[rename[i]] >> rename[j] & 1) for i in range(n))
            need(mapped == tuple(old_rows) and ports == (0, 1), 'old four-letter graph map')
            old_graph_checks += 1
    need(decode((8,)) == ((2, 4, 1), (0, 1)), 'literal C3 adaptive entrance')

    constant_controls = []
    for n in range(2, 7):
        vector = initial
        for _ in range(n):
            vector = advance(vector, 0, 2)
        coefficient = vector[((-1, -1), (-1, 0), (1, 1))]
        need(coefficient == 2**n-2 and sum(c for s,c in vector.items() if len(s) == 1) == 1,
             'constantH but unbounded cover vector hostile')
        constant_controls.append({'births': n, 'H': 1, 'selected_coefficient': coefficient})
    result = {
        'status': 'PASS', 'source_receipt_sha256': hashlib.sha256(raw).hexdigest(),
        'scope': 'Independent twelve-letter carrier and common static-query observer; no spectrum scan.',
        'entire_abstract_nondegenerate_signature_universe': len(signatures),
        'all_twelve_matrices_rebuilt': True, 'all_short_program_literal_H': short,
        'full_literal_cover_vectors': full_vectors, 'static_reachable_rank_lower_bound': 25,
        'field_prime': 1000003, 'common_static_query_rank': 17,
        'all_twelve_updates_verified_on_full25_coordinates': True,
        'actual_Hankel_entries': entries, 'distinct_Hankel_graphs': len(count_cache),
        'maximum_Hankel_graph_order': max_order, 'old_four_letter_graph_map_controls': old_graph_checks,
        'static_basis_prefix_words': static_words, 'static_queries': suffixes,
        'constantH_unbounded_weight_controls': constant_controls,
        'consequence': 'A common exact17-dimensional rational observer for all twelve operations, '
                       'expressed by static continuation counts, valid for every actual source with two '
                       'marked vertices. It identifies future count response, not graph languages or spectra.',
    }
    output = here/'two_port_envelope_independent_audit.json'
    output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'status': 'PASS', 'signatures': len(signatures), 'static_queries': C.cols,
                      'literal_programs': short, 'actual_Hankel_entries': entries, 'old_graph_maps': old_graph_checks}))


if __name__ == '__main__':
    main()
