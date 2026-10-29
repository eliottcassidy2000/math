"""Replay integer minors and exact null vectors, with independent universe checks."""
from collections import Counter
from itertools import combinations, permutations, product
from math import comb, factorial, gcd
import json
import sys

import numpy as np
from sympy import Matrix
from sympy.polys.matrices import DomainMatrix

from atom_precision import automorphisms, hostile
from complement_kernel import (HERE, complement, exact_rank_certificate, native,
                               need, paths, rank_mod, rational_kernel, unit_certificate)
from paley_pencil import diagonal, local_layers


def strong(rows):
    full = (1 << len(rows))-1
    for start in range(len(rows)):
        seen = 1 << start
        while True:
            more = seen
            for i, row in enumerate(rows):
                if seen >> i & 1:
                    more |= row
            if more == seen:
                break
            seen = more
        if seen != full:
            return False
    return True


def prime_atom(rows):
    n = len(rows)
    for size in range(2, n):
        for subset in combinations(range(n), size):
            if all(len({rows[v] >> u & 1 for u in subset}) == 1
                   for v in range(n) if v not in subset):
                return False
    return True


def canonical(rows):
    groups = [[v for v in range(len(rows)) if rows[v].bit_count() == degree]
              for degree in sorted({r.bit_count() for r in rows})]
    best = None
    for parts in product(*(permutations(group) for group in groups)):
        order = sum((list(part) for part in parts), [])
        code = tuple(sum((rows[u] >> v & 1) << j for j, v in enumerate(order)) for u in order)
        best = code if best is None else min(best, code)
    return best


def verify_rank(M, cert):
    p, r, K = cert['prime'], cert['rank'], cert['kernel']
    need(rank_mod(M, p) == r, 'modular lower rank certificate')
    need(len(K) == len(M[0])-r, 'null vector upper bound has matching size')
    if K:
        need(rank_mod(K, p) == len(K), 'independent integer null vectors')
        need(not np.any(np.array(M, dtype=object) @ np.array(K, dtype=object).T),
             'exact integer annihilation, no machine-integer overflow')


def verify_minors(M, minors):
    divisor = 0
    for item in minors:
        indices = item['rows']
        need(len(indices) == len(M[0]) and len(set(indices)) == len(indices), 'square row minor')
        need(all(0 <= i < len(M) for i in indices), 'minor rows within actual observer')
        value = int(DomainMatrix.from_Matrix(Matrix([M[i] for i in indices])).det())
        need(value == item['determinant'], 'exact determinant replay')
        divisor = gcd(divisor, value)
    need(divisor == 1, 'Bezout certificate for integral recovery at every modulus')


def cycle_histogram(rows):
    ports = list(permutations(range(len(rows)), 3))
    index = {p: i for i, p in enumerate(ports)}
    hist = Counter()
    for g in automorphisms(rows):
        sigma = [index[tuple(g[v] for v in p)] for p in ports]
        unseen, count = set(range(len(ports))), 0
        while unseen:
            start = next(iter(unseen))
            u = start
            while u in unseen:
                unseen.remove(u)
                u = sigma[u]
            count += 1
        hist[count] += 1
    return dict(sorted(hist.items()))


def details():
    cases = []
    for label, rows in [('T5', [6, 12, 24, 17, 3]), ('B5', [16, 5, 9, 3, 14])]:
        _, ports, M = native(rows)
        records = []
        for path in paths(rows):
            target = complement(rows, path)
            N = native(target)[2]
            snf = diagonal(M+N)
            records.append(dict(path=path, target=target, target_strong=strong(target), smith=snf))
        histogram = Counter(tuple(x['smith']) for x in records)
        if label == 'T5':
            need(histogram == Counter({(1,)*60: 10, (1,)*59+(2,): 5}), 'T5 complete integral path classification')
        else:
            need(sum(v for k,v in histogram.items() if k == (1,)*60) == 9, 'nine unit B5 paths')
            need(histogram[(1,)*59+(3,)] == 3, 'three ternary B5 paths')
            need(sum(v for k,v in histogram.items() if k.count(0) == 2) == 3, 'three rank58 B5 paths')
        cases.append(dict(label=label, rows=rows, paths=records, port_cycles=cycle_histogram(rows)))
    rows, P = [16, 5, 9, 3, 14], (1, 2, 3, 0, 4)
    _, ports, M = native(rows)
    N = native(complement(rows, P))[2]
    K = rational_kernel(M+N)
    need(len(K) == 2 and all(sum(v) == 0 for v in K), 'two equal-mass hidden directions')
    index = {p:i for i,p in enumerate(ports)}
    need(all(v[index[tuple(g[u] for u in port)]] == v[i]
             for g in automorphisms(rows) for v in K for i,port in enumerate(ports)),
         'hidden vectors already invariant under every vertex automorphism')
    # Run action on this kernel: its character is the standard S3 character.
    B = Matrix(K).T
    pivot_rows = list(B.T.rref()[1])
    inverse = B[pivot_rows, :].inv()
    traces = []
    for g in permutations(range(3)):
        moved = Matrix([[v[index[tuple(port[g[j]] for j in range(3))]] for v in K] for port in ports])
        action = inverse*moved[pivot_rows, :]
        need(B*action == moved, 'entire run action preserves the hidden plane')
        traces.append(int(action.trace()))
    need(Counter(traces) == Counter({2:1, 0:3, -1:2}), 'standard two-dimensional run representation')
    v = K[0]
    plus, minus = [max(x,0) for x in v], [max(-x,0) for x in v]
    need(sum(plus) == sum(minus) == 6, 'nonnegative six-unit weight laws')
    need(all(sum(a*x for a,x in zip(row,v)) == 0 for row in M+N), 'both original and bad complement are blind')
    good = next(x for x in cases[1]['paths'] if x['smith'] == [1]*60)
    G = native(good['target'])[2]
    witness_row = next(i for i,row in enumerate(G) if sum(a*x for a,x in zip(row,v)))
    witness = dict(bad_path=P, kernel=K, run_character=traces, plus=plus, minus=minus,
                   good_path=good['path'], start=ports[witness_row],
                   before=[sum(a*x for a,x in zip(M[witness_row],w)) for w in (plus,minus)],
                   after=[sum(a*x for a,x in zip(G[witness_row],w)) for w in (plus,minus)])
    paley = [sum(1 << ((v+d)%7) for d in (1,2,4)) for v in range(7)]
    need(cycle_histogram(paley) == {30:6,70:14,210:1}, 'Paley affine/Frobenius port cycles')
    # Explicit F8* chart with alpha^3=alpha+1. No normal-basis assumption.
    def mul(a,b):
        out=0
        while b:
            if b & 1: out ^= a
            b >>= 1
            a <<= 1
            if a & 8: a ^= 11
        return out
    powers=[1]
    for _ in range(6): powers.append(mul(powers[-1],2))
    need(len(set(powers)) == 7 and all(mul(powers[i],powers[i]) == powers[2*i%7] for i in range(7)),
         'Frobenius squaring is exponent doubling')
    affine = {tuple((b+pow(2,j,7)*i)%7 for i in range(7)) for b in range(7) for j in range(3)}
    need(affine == {tuple(g) for g in automorphisms(paley)}, 'actual 21-element field action')
    hostile_cycles=[cycle_histogram(hostile(i)) for i in (False,True)]
    need(hostile_cycles == [{148:36,440:18,444:8,520:12,776:6,1320:1},
                            {444:8,520:12,776:6,1320:1}], '81/27 native three-port cycle control')
    need({tuple(g) for g in automorphisms(hostile(True))} <=
         {tuple(g) for g in automorphisms(hostile(False))}, 'literal labelled subgroup control')
    for hist in [case['port_cycles'] for case in cases]+[cycle_histogram(paley)]+hostile_cycles:
        order=sum(hist.values())
        for modulus in range(1,21):
            need(sum(mult*modulus**degree for degree,mult in hist.items()) % order == 0,
                 'Burnside integrality including nonunit group-order characteristics')
        for p in (2,3,5):
            poly=[1]
            for degree in hist:
                root=p**degree
                new=[0]*(len(poly)+1)
                for i,c in enumerate(poly): new[i]-=root*c;new[i+1]+=c
                poly=new
            values=[sum(mult*p**(degree*k) for degree,mult in hist.items())//order for k in range(10)]
            need(all(sum(c*values[k+i] for i,c in enumerate(poly)) == 0
                     for k in range(len(values)-len(poly)+1)), 'all-precision recurrence check')
    return dict(five_vertex=cases, positive_witness=witness, paley_rows=paley,
                paley_port_cycles=cycle_histogram(paley), field8_powers=powers,
                hostile_port_cycles=hostile_cycles)


def verify_scout():
    def blowup(blocks):
        offsets, size = [], 0
        for block in blocks:
            offsets.append(size)
            size += len(block)
        rows=[0]*size
        for b,block in enumerate(blocks):
            for v,row in enumerate(block):
                rows[offsets[b]+v]=row << offsets[b]
                c=(b+1)%3
                rows[offsets[b]+v] |= ((1 << len(blocks[c]))-1) << offsets[c]
        return rows
    c3=[2,4,1]
    expected=[blowup([c3,c3,[2,0]]),blowup([[2,4,9,19,39,15],[0],[0]]),
              blowup([[6,12,24,17,3],[2,0],[0]]),blowup([c3,c3,c3])]
    records=json.loads((HERE/'complement_scout.json').read_text())
    need(len(records) == 4,'four specified larger examples, not a census')
    for record,rows in zip(records,expected):
        need(record['rows'] == rows and strong(rows),'literal strong substitution source')
        cert=record['certificate']
        M=native(rows)[2]
        N=native(complement(rows,cert['path']))[2]
        verify_minors(M+N,cert['minors'])
        if len(rows)==9:
            need(len(automorphisms(rows)) == 81,'recursive nine-vertex source group')
            need([rank_mod(M,p) for p in (2,3,1009)] == [432,255,504], 'nine-vertex source precision controls')
            B=native(complement(rows,tuple(range(9))))[2]
            need([rank_mod(M+B,p) for p in (2,3,1009)] == [502,504,504], 'default path leaves two binary directions')


def verify():
    ranks = json.loads((HERE/'complement_kernel.json').read_text())
    units = json.loads((HERE/'complement_units.json').read_text())
    need(len(ranks) == len(units) == 396, 'complete declared source bank')
    counts = Counter(x['n'] for x in ranks)
    need(counts == Counter({3:1,4:1,5:6,6:35,7:353}), 'strong class counts')
    seen, labelled, depth0 = {}, Counter(), Counter()
    integral_obstructions = Counter()
    for number,(record,unit) in enumerate(zip(ranks, units), 1):
        n, rows = record['n'], record['source']
        need(rows == unit['source'] and n == unit['n'], 'shared source across receipts')
        need(len(rows) == n and all(0 <= r < 1 << n and not(r >> i & 1) for i,r in enumerate(rows)),
             'native vertex set and loop-free rows')
        need(all(((rows[i] >> j) ^ (rows[j] >> i)) & 1 for i in range(n) for j in range(i+1,n)),
             'exactly one directed arc per unordered pair')
        need(strong(rows) and prime_atom(rows) == record['prime_atom'], 'strong versus modular-prime predicates')
        key = canonical(rows)
        need(key not in seen.setdefault(n,set()), 'pairwise nonisomorphic representatives')
        seen[n].add(key)
        aut = len(automorphisms(rows))
        need(aut == record['aut'], 'actual source automorphism order')
        labelled[n] += factorial(n)//aut
        _, ports, M = native(rows)
        need(len(ports) == record['size'], 'port universe')
        verify_rank(M, record['source_certificate'])
        selected_paths = paths(rows)
        need(len(selected_paths) == record['path_count'], 'Hamilton source mass')
        if record['source_certificate']['rank'] < len(ports) and n <= 6:
            need([list(p) for p in selected_paths] == [x['path'] for x in record['complements']],
                 'all legal complements of every deficient source through order six')
        for c in record.get('complements', []):
            target = complement(rows, c['path'])
            need(target == c['target'] and complement(target,c['path']) == rows, 'literal involution with path retained')
            verify_rank(M+native(target)[2],c['certificate'])
        if unit['path'] is None:
            observer = M
            depth0[n] += 1
        else:
            need(tuple(unit['path']) in selected_paths,'unit witness is a legal path')
            observer = M+native(complement(rows,unit['path']))[2]
            bad_prime = next((p for p in (2,3,5,7,11,13,17,19,23,29,31,37) if rank_mod(M,p)<len(ports)),None)
            need(bad_prime is not None, 'source alone cannot recover at all moduli')
            integral_obstructions[(n,bad_prime)] += 1
        verify_minors(observer,unit['minors'])
        if number % 50 == 0: print('verified',number,flush=True)
    # Independent exhaustiveness: unique first strong component of a tournament.
    strong_labelled={1:1}
    for n in range(2,8):
        strong_labelled[n]=2**comb(n,2)-sum(comb(n,k)*strong_labelled[k]*2**comb(n-k,2) for k in range(1,n))
    need(all(labelled[n] == strong_labelled[n] for n in counts), 'orbit masses exhaust every labelled strong tournament')
    need(depth0 == Counter({3:1,4:1,5:3,6:7,7:16}), 'sharp integral depth-zero census')
    saved=json.loads((HERE/'complement_details.json').read_text())
    current=json.loads(json.dumps(details()))
    need(current == saved,'five-vertex, run-action, field and hostile details replay')
    verify_scout()
    print('PASS',dict(counts),'integral depth0',dict(depth0),'labelled',dict(labelled))


if __name__ == '__main__':
    if '--write-details' in sys.argv:
        (HERE/'complement_details.json').write_text(json.dumps(details(),separators=(',',':'))+'\n')
    else:
        verify()
