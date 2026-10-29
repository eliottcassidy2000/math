"""Exact controls for fixed-port nilpotence, K4/Fano, reflection, and descent.

No removable assert gates. Same-agent independent constructions, not a referee.
Run from any directory; optional --output selects the deterministic JSON receipt.
"""
from collections import Counter
from itertools import combinations, permutations, product
from math import comb, gcd, prod
from pathlib import Path
import argparse
import importlib.util
import json
import sys

from flint import fmpz_mat
from sympy import Matrix, ZZ, expand, factorint, symbols, summation
from sympy.matrices.normalforms import smith_normal_form

HERE = Path(__file__).resolve().parent
GATES = 0


def need(ok, why):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(why)


def grammar(word):
    A = [[0] * 4 for _ in range(4)]
    for u, v in [(0, 1), (1, 2), (2, 3), (3, 0), (0, 2), (1, 3)]:
        A[u][v] = 1
    cycle = list(range(4))
    for letter in word:
        n = len(A)
        maxima = [i for i, row in enumerate(A) if sum(row) == n - 2]
        pred = lambda v: cycle[(cycle.index(v) - 1) % n]
        if len(maxima) == 1:
            choices = {'a': pred(maxima[0]), 'b': maxima[0]}
        else:
            need(len(maxima) == 2, 'classical type A/B')
            a, b = maxima
            if A[b][a]:
                a, b = b, a
            choices = {'c': pred(a), 'd': a, 'e': b}
        need(letter in choices, 'legal grammar port')
        p = choices[letter]
        for u in range(n):
            A[u].append(int(u == p))
        A.append([int(v != p) for v in range(n)] + [0])
        cycle.insert(cycle.index(p) + 1, n)
    return A, cycle


def observer(A, k):
    """Literal ordered path covers: a permutation and k-1 cuts, each exactly once."""
    n = len(A)
    tuples = list(permutations(range(n), k))
    index = {a: i for i, a in enumerate(tuples)}
    counts = [[0] * len(tuples) for _ in tuples]
    cuts = list(combinations(range(1, n), k - 1))
    for p in permutations(range(n)):
        bad = {i for i in range(1, n) if not A[p[i-1]][p[i]]}
        for cut in cuts:
            if not bad.issubset(cut):
                continue
            starts = tuple(p[i] for i in (0,) + cut)
            ends = tuple(p[i-1] for i in cut + (n,))
            counts[index[starts]][index[ends]] += 1
    return tuples, index, counts


def tournament_checks():
    records = []
    words = ['', 'c', 'd', 'e', 'ca', 'cb', 'dc', 'dd', 'de', 'ec', 'ed', 'ded']
    for word in words:
        A, cycle = grammar(word)
        n = len(A)
        succ = {cycle[i]: cycle[(i+1) % n] for i in range(n)}
        pred = {v: u for u, v in succ.items()}
        D = {(succ[u], v) for u in range(n) for v in range(n)
             if A[u][v] and v != succ[u]}
        order = []
        while len(order) < n:
            ready = [v for v in range(n) if v not in order and
                     all(u in order for u, w in D if w == v)]
            need(bool(ready), 'avoiding word has an auxiliary DAG')
            order.append(min(ready))
        f = {v: i for i, v in enumerate(order)}
        for k in ([1, 2, 3, 4] if n == 4 else [1, 2, 3] if n <= 6 else [2]):
            tuples, index, M = observer(A, k)
            phi = []
            for a in tuples:
                ends = []
                for u in a:
                    v = succ[u]
                    while v not in a:
                        v = succ[v]
                    ends.append(pred[v])
                phi.append(index[tuple(ends)])
            need(len(set(phi)) == len(tuples), 'canonical port matching is bijective')
            rows = [[M[i][phi[j]] for j in range(len(tuples))] for i in range(len(tuples))]
            weights = [sum(f[v] for v in a) for a in tuples]
            for i, row in enumerate(rows):
                need(row[i] == 1, 'one canonical equality cover')
                need(all(not v or j == i or weights[j] > weights[i]
                         for j, v in enumerate(row)), 'every other entry raises potential')
            B = fmpz_mat(rows)
            size = len(tuples)
            I = fmpz_mat([[int(i == j) for j in range(size)] for i in range(size)])
            zero = fmpz_mat(size, size)
            N = B - I
            bound = k * (n-k)
            powers = [I]
            nil_index = None
            for j in range(1, bound + 2):
                powers.append(powers[-1] * N)
                if powers[-1] == zero:
                    nil_index = j
                    break
            need(nil_index is not None, 'N^(k(n-k)+1)=0')
            # Independent path-length computation in the positive off-diagonal support.
            longest = [0] * size
            for i in sorted(range(size), key=lambda i: weights[i]):
                for j, value in enumerate(rows[i]):
                    if i != j and value:
                        longest[j] = max(longest[j], longest[i] + 1)
            need(max(longest) + 1 == nil_index, 'nilpotence agrees with support longest path')
            inverse = sum(((-1)**j * P for j, P in enumerate(powers)), zero)
            need(B * inverse == I and inverse * B == I, 'integral inverse finite series')
            for depth in [0, 1, 2, 5, 11]:
                polynomial = sum((comb(depth, j) * P for j, P in enumerate(powers)
                                  if j <= depth), zero)
                need(B**depth == polynomial, 'exact binomial composition recurrence')
            records.append({'word': word, 'n': n, 'k': k, 'dimension': size,
                            'degree_bound': bound, 'actual_nilpotence_index': nil_index})
            if not word and k == 2:
                candidates = []
                for vertex_swap in permutations(range(4)):
                    pi = [index[tuple(vertex_swap[v] for v in a)] for a in tuples]
                    remarked = fmpz_mat([[rows[pi[i]][pi[j]] for j in range(size)] for i in range(size)])
                    pair = B * remarked
                    trace = sum(int(pair[i, i]) for i in range(size))
                    candidates.append((trace, vertex_swap))
                trace, vertex_swap = max(candidates)
                need(trace > size, 'native vertex re-marking destroys common nilpotent filtration')
                native_hostile = {'word': '', 'k': 2, 'vertex_swap': vertex_swap,
                                  'dimension': size, 'product_trace': trace,
                                  'all_24_trace_histogram': dict(sorted(Counter(t for t, p in candidates).items()))}
    upper = Matrix([[1, 1], [0, 1]])
    swap = Matrix([[0, 1], [1, 0]])
    mixed = upper * swap * upper * swap
    need(mixed == Matrix([[2, 1], [1, 1]]), 're-marked shear hostile')
    need((mixed - Matrix.eye(2)).det() == -1, 'mixed product minus I not nilpotent')
    A, cycle = grammar('ee')
    tuples, index, M = observer(A, 2)
    wedges = []
    for a, s, r in permutations(range(6), 3):
        if all(row[index[a, r]]-row[index[r, a]]-row[index[s, r]]+row[index[r, s]] == 0
               for row in M):
            wedges.append([a, s, r])
    need(bool(wedges), 'forbidden word ee has native bridge-star wedge kernel')
    return {'native_observers': records, 'remarking_hostile': mixed.tolist(),
            'native_remarking_hostile': native_hostile,
            'ee_wedge_witness': wedges[0], 'ee_wedge_witness_count': len(wedges)}


def symmetry_checks():
    # E = even-weight F2^4; basis e0+e3,e1+e3,e2+e3.
    lift = lambda v: (v & 7) | ((v.bit_count() % 2) << 3)
    bits = [lift(v) for v in range(8)]
    need(len(set(bits)) == 8 and all(v.bit_count() % 2 == 0 for v in bits), 'even space')
    def linear(columns, v):
        out = 0
        for i in range(3):
            if (v >> i) & 1:
                out ^= columns[i]
        return out
    gl = [tuple(linear(cols, v) for v in range(8))
          for cols in product(range(8), repeat=3)
          if len({linear(cols, v) for v in range(8)}) == 8]
    stabilizer = {g for g in gl if g[7] == 7}
    sheet = set()
    for p in permutations(range(4)):
        sheet.add(tuple(sum(((lift(v) >> i) & 1) << p[i] for i in range(4)) & 7
                        for v in range(8)))
    need(len(gl) == 168 and len(stabilizer) == 24 and stabilizer == sheet,
         'K4 sheet action exactly Fano point stabilizer')
    edges = list(combinations(range(4), 2))
    edge_bits = [(1 << i) | (1 << j) for i, j in edges]
    need(set(edge_bits) == set(bits) - {0, 15}, 'six edges plus distinguished seventh vector')
    pairs = {tuple(sorted((v, v ^ 15))) for v in edge_bits}
    need(len(pairs) == 3, 'opposite edges are Fano lines through p')
    # Oriented boundary, NOT the unsigned S4 edge module in THM-3045.
    incidence = Matrix([[int(v == j) - int(v == i) for i, j in edges] for v in range(4)])
    need(incidence.rank() == 3 and len(incidence.nullspace()) == 3, 'wedge boundary dimensions')
    treecols = [edges.index((0, v)) for v in [1, 2, 3]]
    cyclecols = []
    for i, j in combinations([1, 2, 3], 2):
        vec = [0] * 6
        for e, a in [((0, i), 1), ((i, j), 1), ((0, j), -1)]:
            vec[edges.index(e)] = a
        cyclecols.append(Matrix(vec))
    Z = Matrix.hstack(*cyclecols)
    cut = incidence.T[:, 1:]
    need(incidence * Z == Matrix.zeros(4, 3), 'cycle basis exact')
    whole = cut.row_join(Z)
    smith = smith_normal_form(whole, domain=ZZ)
    diag = [abs(int(smith[i, i])) for i in range(6)]
    need(diag == [1, 1, 1, 1, 4, 4], 'oriented cut-cycle lattice index 16')
    states = ['A', 'B0', 'B1', 'FA', 'FB', 'X']
    rules = {'A': {'a': 'A', 'b': 'B1'}, 'B0': {'c': 'A', 'd': 'B0', 'e': 'B1'},
             'B1': {'c': 'A', 'd': 'B0', 'e': 'FB'}, 'FA': {'a': 'FA', 'b': 'FB'},
             'FB': {'c': 'FA', 'd': 'FB', 'e': 'FB'}, 'X': {}}
    delta = [[states.index(rules[s].get(a, 'X')) for a in 'abcde'] for s in states]
    output = [1, 1, 1, 0, 0, -1]
    automorphisms = [p for p in permutations(range(6))
                     if all(output[i] == output[p[i]] and
                            all(p[delta[i][a]] == delta[p[i]][a] for a in range(5))
                            for i in range(6))]
    need(len(automorphisms) == 1, 'named-port six-state machine has trivial automorphism group')
    need(len({row[0] for row in delta}) == 3, 'a operation not an edge permutation')
    return {'GL3F2_order': len(gl), 'distinguished_point_stabilizer': len(stabilizer),
            'opposition_pairs': sorted(pairs), 'oriented_cut_cycle_smith': diag,
            'six_state_port_automorphism_order': len(automorphisms),
            'letter_image_sizes': dict(zip('abcde', [len({r[a] for r in delta}) for a in range(5)]))}


def field_clock_checks():
    # F169 = F13[t]/(t^2-2), since 2 is nonsquare mod13.
    need(2 not in {x*x % 13 for x in range(13)}, 'irreducible quadratic')
    def mul(x, y):
        a, b = x; c, d = y
        return ((a*c+2*b*d) % 13, (a*d+b*c) % 13)
    def power(x, n):
        result = (1, 0)
        while n:
            if n & 1:
                result = mul(result, x)
            x = mul(x, x)
            n //= 2
        return result
    nonzero = [x for x in product(range(13), repeat=2) if x != (0, 0)]
    fibres = {}
    for x in nonzero:
        need(power(x, 168) == (1, 0), 'F169 multiplicative clock')
        norm, ratio = power(x, 14), power(x, 12)
        need(norm[1] == 0 and power(ratio, 14) == (1, 0), 'split/norm-one targets')
        need(power(norm, 6) == power(ratio, 7), 'CRT shared-sign compatibility')
        fibres.setdefault((norm, ratio), []).append(x)
    need(len(fibres) == 84 and set(map(len, fibres.values())) == {2}, 'one missing sign bit')
    need(all(set(xs) == {xs[0], tuple(-v % 13 for v in xs[0])}
             for xs in fibres.values()), 'all fibres exactly +/- x')
    primitive = next(x for x in nonzero if all(power(x, 168//p) != (1, 0) for p in [2, 3, 7]))
    return {'field_size': 169, 'primitive_element': primitive, 'nonzero_clock': 168,
            'norm_values': len({key[0] for key in fibres}), 'ratio_values': len({key[1] for key in fibres}),
            'compatible_pairs': len(fibres), 'each_fibre_size': 2}


def menu_checks():
    """Exhaust all pairs of binary off-diagonal supports on three states."""
    edges = list(permutations(range(3), 2))
    matrices = []
    for bits in range(64):
        A = [[int(i == j) for j in range(3)] for i in range(3)]
        for t, (i, j) in enumerate(edges):
            A[i][j] = (bits >> t) & 1
        matrices.append(A)
    def mul(A, B):
        return [[sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    identity = [[int(i == j) for j in range(3)] for i in range(3)]
    acyclic = 0
    cyclic = 0
    for A, B in product(matrices, repeat=2):
        union = [[int(i != j and (A[i][j] or B[i][j])) for j in range(3)] for i in range(3)]
        N3 = mul(mul(union, union), union)
        cycle = None
        for length in [2, 3]:
            for seq in permutations(range(3), length):
                if all(union[seq[i]][seq[(i+1) % length]] for i in range(length)):
                    cycle = seq
                    break
            if cycle is not None:
                break
        need((cycle is None) == (not any(map(any, N3))), 'cycle search equals union nilpotence')
        if cycle is None:
            acyclic += 1
            for word in product([A, B], repeat=3):
                value = identity
                for factor in word:
                    value = mul(value, factor)
                need(sum(value[i][i] for i in range(3)) == 3, 'common-rank product control')
        else:
            cyclic += 1
            value = identity
            for i, u in enumerate(cycle):
                v = cycle[(i+1) % len(cycle)]
                value = mul(value, A if A[u][v] else B)
            need(sum(value[i][i] for i in range(3)) > 3, 'explicit short cyclic-word trace witness')
    return {'universe': 'all 4096 ordered pairs of 3x3 binary off-diagonal supports',
            'acyclic_menus': acyclic, 'cyclic_menus': cyclic, 'maximum_witness_length': 3}


def arithmetic_checks():
    c, k, j = symbols('c k j', integer=True)
    signed = lambda power, step: c**power + summation((c-step*j)**power-(c+step*j)**power, (j, 1, k))
    need(expand(signed(2, 1)-c*signed(1, 2)) == 0, 'symbolic reflection identity')
    s = k*(k+1)/2
    need(expand(signed(3, 1)-(c*c*(c-6*s)-2*s*s)) == 0, 'symbolic cubic obstruction')
    total = summation((2*k*(k+1)+j)**2, (j, 1, k))
    need(expand(total-k*(k+1)*(2*k+1)*(12*k*(k+1)+1)/6) == 0, 'all-k square total')
    rows = []
    for k in range(1, 301):
        c = 2*k*(k+1)
        linear = sum(c+2*j for j in range(-k, 1)) - sum(c+2*j for j in range(1, k+1))
        square = sum((c+j)**2 for j in range(-k, 1)) - sum((c+j)**2 for j in range(1, k+1))
        need(linear == square == 0, 'linear-square reflection family')
        s = k*(k+1)//2
        cubic = lambda c: c**3-6*s*c*c-2*s*s
        need(cubic(6*s) < 0 < cubic(6*s+1), 'cubic root between consecutive integers')
        if k <= 6:
            rows.append({'k': k, 'center': c, 'linear_start': 2*k*k,
                         'square_base_start': k*(2*k+1), 'square_total': sum((c+j)**2 for j in range(1, k+1))})
    offsets_cases = list(product(range(1, 6), repeat=3))
    for offsets in offsets_cases:
        c = 4*sum(offsets)
        need(c*c + sum((c-d)**2-(c+d)**2 for d in offsets) == 0, 'arbitrary-offset polarization')
    def rad(n):
        return prod(factorint(n))
    examples = []
    for a in range(1, 61):
        for b in range(a+1, 61):
            if gcd(a, b) != 1:
                continue
            c = a+b
            g = gcd((b-a)**2, 4*a*b)
            A, B, C = (b-a)**2//g, 4*a*b//g, c*c//g
            need(g in [1, 4] and A+B == C and gcd(A, B) == 1, 'primitive ABC square lift')
            odd = b-a
            while odd % 2 == 0:
                odd //= 2
            need(rad(A*B*C) == rad(a*b*c)*rad(odd), 'exact new-prime radical budget')
            need(C >= c and (C == c) == ((a, b) == (1, 3)), 'square lift grows except fixed 1+3=4')
            if (a, b) in [(1, 2), (1, 3), (1, 8)]:
                examples.append({'source': [a, b, c], 'lift': [A, B, C], 'radicals': [rad(a*b*c), rad(A*B*C)]})
    return {'reflection_checked_k': [1, 300], 'rows': rows,
            'arbitrary_offset_cases': len(offsets_cases), 'abc_universe': '1<=a<b<=60, gcd(a,b)=1',
            'abc_examples': examples}


def descent_checks():
    def odd_map(n, m=3):
        y = m*n+1
        while y % 2 == 0:
            y //= 2
        return y
    for length in range(1, 81):
        n = 2**(length+1)-1
        x = n
        for j in range(length+1):
            need(x == 3**j*2**(length+1-j)-1, 'actual arbitrary long expanding prefix')
            if j < length:
                need((3*x+1) % 4 == 2, 'valuation exactly one')
                x = odd_map(x)
    need([odd_map(x, 5) for x in [13, 33, 83]] == [33, 83, 13], '5x+1 cycle rejects ranks')
    # Finite certificates only: first descent, retaining each exact valuation.
    receipts = []
    for n in range(3, 2001, 2):
        x = n
        exponents = []
        for _ in range(10000):
            y = 3*x+1
            e = 0
            while y % 2 == 0:
                y //= 2; e += 1
            exponents.append(e); x = y
            if x < n:
                break
        need(x < n, 'finite first-descent source certificate')
        receipts.append({'n': n, 'target': x, 'exponents': exponents})
    return {'expanding_prefix_lengths': [1, 80], 'five_x_plus_one_cycle': [13, 33, 83],
            'scope': 'FINITE-EXACT for odd sources 3..1999; no all-source coverage theorem',
            'max_first_descent_steps': max(len(r['exponents']) for r in receipts), 'receipts': receipts}


def inherited_hostile():
    path = HERE.parent / 'four_vertex_glue_fermat_20261029/verify.py'
    spec = importlib.util.spec_from_file_location('inherited_81_27', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    result = mod.hostile_control()
    need(result['automorphisms'] == [81, 27], '81/27 remains distinct')
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=HERE/'certificate.json')
    args = parser.parse_args()
    out = {'tournaments': tournament_checks(), 'menus': menu_checks(), 'symmetry': symmetry_checks(),
           'field_clocks': field_clock_checks(), 'arithmetic': arithmetic_checks(),
           'descent': descent_checks(), 'hostile_81_27': inherited_hostile()}
    out['gates'] = GATES
    args.output.write_text(json.dumps(out, indent=2, default=int)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({'gates': GATES, 'observers': len(out['tournaments']['native_observers']),
                      'descent_sources': len(out['descent']['receipts']), 'output': str(args.output)}))


if __name__ == '__main__':
    main()
