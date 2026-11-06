#!/usr/bin/env python3
"""Independent module-gap binomial census and literal path/cycle DFS."""
from itertools import permutations
from math import comb
from collections import Counter
import json

CASES = [
    ('A', [2,24,3,37,13,23], 15, [2,8,67,101,15,87,19], 79, (1,4)),
    ('B', [32,5,1,39,47,70,31], 5, [0,5,129,39,47,199,159,27], 5, (0,5)),
    ('C', [2,4,41,3,15,91,31], 27, [0,5,41,131,15,219,159,23], 155, (0,1)),
]
EXPECTED = {
    'A': [[25,59,41,8],[51,119,82,16],[2,3,1],[3,5,2]],
    'B': [[25,45,20],[51,91,40],[1,1],[0]],
    'C': [[25,53,27],[51,107,54],[1,1],[0]],
}

def need(b, s):
    if not b:
        raise RuntimeError(s)

def valid(a):
    n = len(a)
    for i, row in enumerate(a):
        need(0 <= row < 1 << n and not row & (1 << i), 'bad row')
        for j in range(i):
            need(((row >> j) & 1) + ((a[j] >> i) & 1) == 1, 'bad pair')

def restrict(a, order):
    return [sum(1 << j for j, v in enumerate(order) if a[u] & (1 << v)) for u in order]

def inflate(a, mask, k):
    n = len(a)
    result = a[:] + [0] * k
    for v in range(n, n + k):
        result[v] = mask | sum(1 << w for w in range(v + 1, n + k))
    for v in range(n):
        if not mask & (1 << v):
            result[v] |= sum(1 << w for w in range(n, n + k))
    valid(result)
    return result

def gap_census(a, mask, cycle=False):
    n = len(a)
    coeff = [0] * (n + 2)
    hist = Counter()
    count = 0
    words = ((0,) + p for p in permutations(range(1,n))) if cycle else permutations(range(n))
    for p in words:
        count += 1
        allowed = 0 if cycle else int(bool(mask & (1 << p[0]))) + int(not mask & (1 << p[-1]))
        mandatory = 0
        pairs = zip(p, p[1:] + p[:1]) if cycle else zip(p, p[1:])
        for left, right in pairs:
            fillable = not mask & (1 << left) and bool(mask & (1 << right))
            allowed += int(fillable)
            if not a[left] & (1 << right):
                if not fillable:
                    break
                mandatory += 1
        else:
            hist[allowed, mandatory] += 1
    # Select r occupied gaps, including all mandatory ones.
    for (allowed, mandatory), multiplicity in hist.items():
        for r in range(mandatory, allowed + 1):
            coeff[r] += multiplicity * comb(allowed - mandatory, r - mandatory)
    while len(coeff) > 1 and coeff[-1] == 0:
        coeff.pop()
    exp = [0] * len(coeff)
    for r, value in enumerate(coeff):
        for j in range(r + 1):
            exp[r - j] += value * (-1) ** j * comb(r,j)
    return {'coefficients': coeff, 'power_coefficients': exp,
            'orders_examined': count,
            'allowed_mandatory_histogram': {str(x): y for x,y in sorted(hist.items())}}

def literal(a):
    full = (1 << len(a)) - 1
    calls = 0
    def paths(last, used, root):
        nonlocal calls
        calls += 1
        if used == full:
            return (1, int(root == 0 and bool(a[last] & 1)))
        h = c = 0
        choices = a[last] & ~used
        while choices:
            bit = choices & -choices
            choices -= bit
            nh, nc = paths(bit.bit_length()-1, used|bit, root)
            h += nh
            c += nc
        return h, c
    h = c = 0
    for root in range(len(a)):
        nh, nc = paths(root, 1 << root, root)
        h += nh
        c += nc
    return h, c, calls

def evaluate(power, k):
    return sum(c * r ** k for r,c in enumerate(power))

def main():
    results = []
    total_graphs = total_calls = 0
    for name, parent, pm, child, tm, flip in CASES:
        valid(parent); valid(child)
        n = len(parent)
        need(len(child) == n + 1 and (tm & ((1 << n)-1)) == pm, 'core/mask sizes')
        changed = [(i,j) for i in range(n) for j in range(i+1,n)
                   if bool(parent[i] & (1 << j)) != bool(child[i] & (1 << j))]
        need(changed == [flip], 'wrong core edit')
        need(parent[flip[0]] & (1 << flip[1]), 'wrong directed reversal')
        stats = [gap_census(parent,pm), gap_census(child,tm),
                 gap_census(parent,pm,True), gap_census(child,tm,True)]
        need([s['coefficients'] for s in stats] == EXPECTED[name], 'coefficient discrepancy')
        if name != 'A':
            need(child[0] == 0 and tm & 1, 'not a universal sink')
            reduced = restrict(child, list(range(1,n+1)))
            reduced_stats = [gap_census(reduced,tm>>1), gap_census(reduced,tm>>1,True)]
            need(reduced_stats[0]['coefficients'] == stats[1]['coefficients'], 'sink path count')
            need(reduced_stats[1]['coefficients'] == [3,3], 'replacement cycles')
        else:
            reduced_stats = None
        samples = []
        for k in range(6):
            p = inflate(parent,pm,k)
            t = inflate(child,tm,k)
            # Identify retained old labels in child despite moving the new vertex ahead of the module.
            old_order = list(range(n)) + list(range(n+1,n+k+1))
            old_child = restrict(t,old_order)
            physical_changes = [(i,j) for i in range(len(p)) for j in range(i+1,len(p))
                                if bool(p[i] & (1 << j)) != bool(old_child[i] & (1 << j))]
            need(physical_changes == [flip], 'inflation changed extra old arcs')
            ph,pc,workp = literal(p)
            th,tc,workt = literal(t)
            total_graphs += 2; total_calls += workp + workt
            need((ph,pc) == (evaluate(stats[0]['power_coefficients'],k), evaluate(stats[2]['power_coefficients'],k)), 'parent literal control')
            need((th,tc) == (evaluate(stats[1]['power_coefficients'],k), evaluate(stats[3]['power_coefficients'],k)), 'child literal control')
            need(th == 2*ph+1, 'doubling control')
            sample = {'k':k,'parent_H':ph,'parent_C':pc,'child_H':th,'child_C':tc}
            if name != 'A':
                reduced = restrict(t,list(range(1,len(t))))
                rh,rc,workr = literal(reduced)
                total_graphs += 1; total_calls += workr
                need((rh,rc) == (th,3), 'replacement literal control')
                need(len(reduced) == len(p), 'replacement order')
                sample.update({'replacement_H':rh,'replacement_C':rc})
            samples.append(sample)
        results.append({'name':name,'parent_core':parent,'child_core':child,
                        'parent_mask':pm,'child_mask':tm,'flip':flip,
                        'path_parent':stats[0],'path_child':stats[1],
                        'cycle_parent':stats[2],'cycle_child':stats[3],
                        'reduced_stats':reduced_stats,'literal_samples':samples})
    print(json.dumps({'status':'PASS','families':results,
                      'literal_graphs':total_graphs,'literal_DFS_calls':total_calls},indent=2,sort_keys=True))

if __name__ == '__main__':
    main()
