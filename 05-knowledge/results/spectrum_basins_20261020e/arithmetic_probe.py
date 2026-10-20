"""Exact golden SCC word / parity-path matching dictionary. Standard library only."""
from collections import Counter
from functools import lru_cache
from hashlib import sha256
import json
from math import gcd
from pathlib import Path


def check(ok, detail):
    if not ok:
        raise RuntimeError(detail)


@lru_cache(None)
def words(n):
    if n == 0:
        return ((),)
    return tuple((a,) + w for a in (1, 3, 4) if a <= n for w in words(n-a))


def encode(w):
    offset = 0
    result = []
    for a in w:
        if a in (3, 4):
            result.append((offset, offset+2))
        if a == 4:
            result.append((offset+1, offset+3))
        offset += a
    return tuple(result)


def decode(n, edges):
    edges = set(edges)
    result, i = [], 0
    while i < n:
        if (i, i+2) not in edges:
            result.append(1)
            i += 1
        elif (i+1, i+3) in edges:
            result.append(4)
            i += 4
        else:
            result.append(3)
            i += 3
    check(sum(result) == n and set(encode(result)) == edges, ('decode', n, edges))
    return tuple(result)


def matchings(n):
    """Independent all-subsets recursion on edges i--i+2."""
    edges = tuple((i, i+2) for i in range(n-2))
    result = []
    def rec(j, used, chosen):
        if j == len(edges):
            result.append(tuple(chosen))
            return
        rec(j+1, used, chosen)
        a, b = edges[j]
        bits = (1 << a) | (1 << b)
        if used & bits == 0:
            rec(j+1, used | bits, chosen + [edges[j]])
    rec(0, 0, [])
    return result


def convolution(a, b):
    c = [0] * (len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i+j] += x*y
    return c


def main():
    fib = [0, 1]
    for _ in range(202):
        fib.append(fib[-1]+fib[-2])
    a = [1]
    for n in range(1, 201):
        a.append(sum(a[n-d] for d in (1, 3, 4) if n >= d))
        check(a[n] == fib[n//2+1]*fib[(n+1)//2+1], ('formula', n))
    matching_poly = [[1], [1]]
    for n in range(2, 19):
        x, y = matching_poly[-1], [0] + matching_poly[-2]
        matching_poly.append([(x[i] if i < len(x) else 0) + (y[i] if i < len(y) else 0)
                              for i in range(max(len(x), len(y)))])
    rows = []
    for n in range(19):
        ws, ms = words(n), matchings(n)
        encoded = {encode(w) for w in ws}
        check(len(encoded) == len(ws) and encoded == set(ms), ('bijection', n))
        check(all(decode(n, m) in ws for m in ms), ('inverse', n))
        potential = Counter(w.count(3)+2*w.count(4) for w in ws)
        check(potential == Counter(map(len, ms)), ('potential', n))
        poly = convolution(matching_poly[n//2], matching_poly[(n+1)//2])
        check([potential[j] for j in range(len(poly))] == poly, ('poly', n))
        check(len(ws) == a[n], ('count', n))
        fixed = 0
        for w in ws:
            reflected = tuple(sorted((n-1-b, n-1-a) for a,b in encode(w)))
            check(reflected == encode(w[::-1]), ('reflection equivariance', n, w))
            fixed += w == w[::-1]
        def symmetric_path_count(v):
            return fib[v//2+2] if v % 2 == 0 else fib[v//2+1]
        predicted_fixed = fib[n//2+1] if n % 2 == 0 else symmetric_path_count(n//2)*symmetric_path_count((n+1)//2)
        check(fixed == predicted_fixed, ('reflection fixed', n))
        rows.append({'n': n, 'count': len(ws), 'triangle_potential_coefficients': poly,
                     'self_converse_words': fixed, 'converse_merged_words': (len(ws)+fixed)//2,
                     'R4_count_distribution': dict(sorted(Counter(w.count(4) for w in ws).items())),
                     'feedback_inventory_distribution': dict(sorted(Counter(w.count(3)+w.count(4) for w in ws).items()))})
    farey = []
    for k in range(1, 101):
        p, q, r = fib[k+1:k+4]
        check(gcd(p, q) == 1 and abs(p*r-q*q) == 1, ('farey', k))
        check((q+p, 2*q+p) == (fib[k+3], fib[k+4]), ('two-step modular action', k))
        raw = [q*q-p*p, 2*p*q, q*q+p*p]
        content = gcd(gcd(raw[0], raw[1]), raw[2])
        check(content in (1, 2) and sum(x*x for x in raw[:2]) == raw[2]**2, ('euclid', k))
        if k <= 10:
            farey.append({'k': k, 'n': 2*k+1, 'p': p, 'q': q, 'pq': p*q,
                          'cassini_determinant': p*r-q*q,
                          'euclid_raw': raw, 'content': content,
                          'euclid_primitive': [v//content for v in raw]})
    for k in range(101):
        check(fib[k+1]*fib[k+2] == sum(fib[i]**2 for i in range(1, k+2)), ('squares', k))
        # Weighted sum of squares carries the number of deleted terminal edges.
        if k <= 8:
            target = convolution(matching_poly[k], matching_poly[k+1])
            rhs = [0]*len(target)
            for i in range(k+1):
                for j, value in enumerate(convolution(matching_poly[i], matching_poly[i])):
                    rhs[k-i+j] += value
            check(target == rhs, ('weighted squares', k))
    for key in ('triangle_potential_coefficients', 'R4_count_distribution', 'feedback_inventory_distribution'):
        values = rows[13][key]
        values = list(values.values()) if isinstance(values, dict) else values
        sums = {0}
        for value in values:
            sums |= {x+value for x in list(sums)}
        check(63 not in sums and 210 not in sums, ('unexpected inventory partition', key))
    # The missing coordinate is a positional edge mark, not a coarse inventory.
    n = 13
    marked, reflected = (3, 5), (7, 9)
    buckets = {(False,False):[], (False,True):[], (True,False):[], (True,True):[]}
    for w in words(n):
        edges = set(encode(w))
        flags = (marked in edges, reflected in edges)
        buckets[flags].append(w)
    check({k:len(v) for k,v in buckets.items()} == {(False,False):168,(False,True):42,(True,False):42,(True,True):21},
          'marked partition')
    union = [w for flags,ws in buckets.items() if any(flags) for w in ws]
    fixed_union = [w for w in union if w == w[::-1]]
    merged_union = {min(w,w[::-1]) for w in union}
    active_union = [w for w in union if 4 in w or w != w[::-1]]
    check(len(union) == 105 and len(fixed_union) == 3 and len(merged_union) == 54, 'union quotient')
    check(len(active_union) == 103 and len({min(w,w[::-1]) for w in active_union}) == 52,
          'active union quotient')
    expected_fixed = {(1,1,1,3,1,3,1,1,1), (3,3,1,3,3), (1,1,4,1,4,1,1)}
    check(set(fixed_union) == expected_fixed, 'fixed preimages')
    examples = []
    for edge_set in ((), (marked,), (reflected,), (marked,reflected), ((2,4),marked)):
        w = decode(n, edge_set)
        examples.append({'edges':edge_set,'word':w,'selected':marked in edge_set,'reflected':reflected in edge_set})
    marked_partition = {'n':13,'shorter_parity_path':[1,3,5,7,9,11],
                        'selected_edge_zero_based':marked,'reflected_edge_zero_based':reflected,
                        'selected':63,'avoid_selected':210,'both':21,'exactly_one':84,'union':105,'neither':168,
                        'selected_only':42,'reflected_only':42,'union_fixed_words':sorted(fixed_union),
                        'union_converse_merged':54,'active_union_ordinary':103,'active_union_converse_merged':52,
                        'examples':examples,
                        'sidecar':'Selected second edge of the shorter path; its reflected pair for the invariant union.',
                        'no_identification_claim':'The21+84split is not identified with Fano/junction orbits; no intertwining map has been given.'}
    payload = {'status': 'FINITE-EXACT', 'universe': 'all golden words and all parity-path matchings through vertex order 18',
               'source_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
               'formula_test_through_n': 200, 'odd_sequence': [a[n] for n in range(1, 20, 2)],
               'half_plane_action': {'matrix_on_q_p':[[2,1],[1,1]],'determinant':1,
                                     'mobius':'z -> (2z+1)/(z+1)', 'vertex_order_increment':4,
                                     'fixed_point_polynomial':[1,-1,-1],
                                     'scope':'Order-changing action on marked slopes, not a fixed-order black metagraph move'},
               'rows': rows, 'farey_controls': farey,
               '273': {'product': [13, 21], 'squares': [fib[i]**2 for i in range(1, 8)],
                       'triangle_potential_coefficients': rows[13]['triangle_potential_coefficients'],
                       'split_63_210': marked_partition},
               'claims_not_made': ['map of individual SCC words to distinct knots', 'black-move conjugacy with Farey moves',
                                   'primitive isosceles limit along Fibonacci parameters', 'finite Banach-Tarski paradox']}
    output = Path(__file__).with_suffix('.json')
    output.write_text(json.dumps(payload, indent=2, sort_keys=True)+'\n', encoding='utf-8')
    print(json.dumps({'all_checks': 'PASS', 'word_matching_pairs': sum(x['count'] for x in rows),
                      'odd_sequence': payload['odd_sequence'], '273': payload['273'],
                      'receipt_sha256': sha256(output.read_bytes()).hexdigest()}, sort_keys=True))


if __name__ == '__main__':
    main()
