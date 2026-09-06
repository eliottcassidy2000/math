#!/usr/bin/env python3
"""Bounded complete K8 signed-cycle profile scout, all 2^21 switching classes."""
from collections import Counter
from itertools import combinations
import json
import sys
import numpy as np
import sympy as sp
from even_graph_d5_six_vertex_synthesis_sep05 import cycles, transform, single_edge_masks


def check(condition, label):
    if not condition:
        raise RuntimeError(label)


def main():
    n = 8
    index = {edge: j for j, edge in enumerate(combinations(range(1, n), 2))}
    size = 1 << len(index)
    counts = []
    for k in range(3, 9):
        masks = cycles(n, k, index)
        f = np.zeros(size, dtype=np.int64)
        f[masks] = 1
        counts.append((len(masks)-transform(f))//2)
    total = sum(counts)
    minimum = int(total[1:].min())
    minimizers = set(map(int, np.flatnonzero(total == minimum)))
    check(minimizers == single_edge_masks(n, index), 'D7 K8 equality classes')
    records = {'n': n, 'classes': size, 'D7_minimum': minimum,
               'D7_minimizers': len(minimizers)}
    c3,c4,c5,c6,c7,c8 = counts
    for k,c in ((4,c4),(6,c6),(8,c8)):
        m=int(c[1:-1].min())
        args=set(map(int,np.flatnonzero(c==m)))
        edges=single_edge_masks(n,index)
        both=edges|{(size-1)^h for h in edges}
        records[f'even_layer_{k}']={'minimum':m,'ties':len(args),
                                   'exact_signed_edge_classes':args==both}
    check(int(c3[0]) == int(c8[0]) == 0, 'balanced control')
    check(int(c4[-1]) == int(c6[-1]) == int(c8[-1]) == 0, 'antibalanced even control')
    for label, values in (
            ('3c8-c6',3*c8-c6), ('c8-c6',c8-c6),
            ('3c8-2c6',3*c8-2*c6), ('3c8-4c6',3*c8-4*c6),
            ('c8-3c4-c6',c8-3*c4-c6),
            ('c7+3c8-5c4-3c5-c6',c7+3*c8-5*c4-3*c5-c6)):
        m = int(values.min())
        arg = int(values.argmin())
        records[label] = {'minimum': m, 'mask': arg,
                          'profile': [int(c[arg]) for c in counts]}
    check(int((3*c8-4*c6).min()) == 0, 'exact local ratio comparison')
    profiles = Counter(tuple(map(int, row)) for row in np.stack(counts, axis=1))
    records['profile_count'] = len(profiles)
    positive = c6>0
    indices = np.flatnonzero(positive)
    ratios = c8[positive]/c6[positive]
    arg = int(indices[int(ratios.argmin())])
    records['min_c8_c6'] = {'numerator': int(c8[arg]), 'denominator':int(c6[arg]),
                           'profile':[int(c[arg]) for c in counts]}
    x,t = sp.symbols('n t')
    A = sum(sp.prod(x-j for j in range(2,k)) for k in range(3,9))
    for q in (sp.Rational(13,2),sp.Rational(20,3),sp.Rational(27,4),7):
        gap = sp.expand(x*A.subs(x,x-1)-(x-q)*A-(q-3)*x*(x-1)*(x-2)/6)
        records[f'induction_center_{q}'] = {
            'gap':str(sp.factor(gap)),
            'at9':str(gap.subs(x,9)),
            'shift9':str(sp.Poly(sp.expand(gap.subs(x,9+t)),t).as_expr())}
    print(json.dumps(records,indent=2,sort_keys=True))
    print('RESULT=PASS')


if __name__ == '__main__':
    main()
