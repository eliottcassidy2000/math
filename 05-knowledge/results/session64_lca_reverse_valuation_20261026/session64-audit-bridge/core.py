#!/usr/bin/env python3
"""Independent finite audit: raw forward-pair Gaussian equations, no imports.

The only producer data read are found.json's twenty literal certificates.
Every source and every old word is reconstructed independently from the
published source definition and displayed base ledgers in LOW_ORDERS_REVIEW.
"""
from collections import Counter
from hashlib import sha256
from itertools import combinations, product, permutations
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
INPUT = HERE / 'found.json'

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def pair(a, b):
    return (min(a, b), max(a, b))

def path(q):
    return {pair(a, b) for a, b in zip(q, q[1:])}

class Linear:
    def __init__(self, n):
        self.n = n
        self.pairs = list(combinations(range(n), 2))
        self.index = {p: i for i, p in enumerate(self.pairs)}
        self.m = len(self.pairs)
        self.one = 1 << self.m
        self.mask = self.one - 1
        self.basis = {}
        self.equations = []
        self.ok = True

    def copy(self):
        c = Linear(self.n)
        c.basis = self.basis.copy()
        c.equations = self.equations[:]
        c.ok = self.ok
        return c

    def arc(self, a, b, flips=frozenset()):
        p = pair(a, b)
        return (1 << self.index[p]) ^ (self.one if (a > b) ^ (p in flips) else 0)

    def reduce(self, row):
        for pivot in sorted(self.basis, reverse=True):
            if row & (1 << pivot):
                row ^= self.basis[pivot]
        return row

    def add(self, row):
        self.equations.append(row)
        row = self.reduce(row)
        if row & self.mask:
            self.basis[(row & self.mask).bit_length() - 1] = row
        elif row:
            self.ok = False

    def force(self, a, b, value):
        self.add(self.arc(a, b) ^ (self.one if value else 0))

    def anti(self, q, flips=frozenset()):
        r = dict(zip(q, reversed(q)))
        for a, b in self.pairs:
            self.add(self.arc(a, b, flips) ^ self.arc(r[b], r[a], flips))

    def dimension(self):
        need(self.ok, 'dimension of inconsistent family')
        return self.m - len(self.basis)

    def models(self):
        need(self.ok, 'models of inconsistent family')
        free = [v for v in range(self.m) if v not in self.basis]
        for choice in range(1 << len(free)):
            value = sum(((choice >> j) & 1) << v for j, v in enumerate(free))
            for pivot, row in sorted(self.basis.items()):
                parity = ((row & self.mask & value).bit_count() & 1) ^ bool(row & self.one)
                if parity:
                    value |= 1 << pivot
            need(all(((row & self.mask & value).bit_count() & 1) == bool(row & self.one)
                     for row in self.equations), 'full raw equation solution')
            rows = [0] * self.n
            for i, (a, b) in enumerate(self.pairs):
                u, v = (a, b) if value >> i & 1 else (b, a)
                rows[u] |= 1 << v
            yield tuple(rows)

def source(n, incoming, opened, chord):
    s = Linear(n)
    s.anti(tuple(range(n)))
    for i in range(n - 1):
        s.force(i, i + 1, 1)
    s.force(n - 1, 0, 1)
    for i in range(1, n - 1):
        s.force(0, i, i not in incoming)
    s.force(1, n - 2, opened)
    s.force(2, n - 2, chord)
    need(s.ok, ('source inconsistent', n, incoming, opened, chord))
    return s

def old_words(n, I, opened, chord):
    a, b, x, y, c, z = 0, 1, 2, n-3, n-2, n-1
    M = list(range(2, c))
    if opened:
        if y not in I:
            return 'open_release', [[b]+M+[c,z,a], [z]+M+[c,b,a], [a,z]+M+[c,b]]
        if c not in I:
            return 'open_penultimate', [[a,b,c,z]+M, [b,c,a,z]+M, [b]+M+[c,a,z]]
        h = next(iter(I-{y,c}))
        cut = y-int(h == y-1)
        L, R = list(range(2,cut)), list(range(cut,c))
        return 'open_endrun', [[a,b,c,z]+L+R, [b,c,z]+L+[a]+R, [b]+L+[a]+R+[c,z]]
    if c in I:
        need(y not in I, 'closed-endpair family is separate')
        cut = min({z-v for v in I}-{1})
        L, R = list(range(1,cut)), list(range(cut,c))
        return 'closed_last_movable', [[a]+L+[z]+R+[c], L+[z]+R+[a,c], L+R+[a,c,z]]
    if not chord:
        return 'closed_zero', [[c]+M+[b,z,a], [z,c]+M+[b,a], [a,z,b]+M+[c]]
    if x not in I:
        return 'closed_one', [[c,b,z,a]+M, [z,c,b,a]+M, [a,z,b]+M+[c]]
    if y not in I:
        cut = min(z-v for v in I if z-v > 2)
        name = 'closed_second_movable'
    elif 3 not in I:
        cut, name = y, 'closed_second_fixed'
    else:
        cut, name = n-4, 'closed_second_tail4'
    L, R = list(range(2,cut)), list(range(cut,c))
    return name, [[a,b]+L+[z]+R+[c], [a]+L+[z]+R+[b,c], [a]+L+R+[b,c,z]]

def symbolic_ledger(s, words):
    P = list(range(s.n))
    flips = set()
    blue_families = []
    for q in words:
        need(sorted(q) == P, 'symbolic word is a permutation')
        for a, b in zip(q, q[1:]):
            need(s.reduce(s.arc(a,b,flips)) == s.one, ('unforced join',q,a,b))
        t = s.copy()
        t.anti(q, flips)
        blue_families.append(t)
        flips ^= set(s.pairs) ^ path(q)
    need(flips == set(s.pairs) ^ path(P), 'complete symbolic labelled K target')
    need(path(words[0]) ^ path(words[1]) ^ path(words[2]) == path(P), 'edge XOR')
    return blue_families

def endpair_subcells(s, h, chord):
    """Literal complete complementary splits from the inherited h-family proofs.

    Return (name, guards, words). Guards branch only on genuine source arcs;
    each leaf's cardinality is checked below against the unsplit source cell.
    """
    n = s.n
    a,b,x,y,c,z = 0,1,2,n-3,n-2,n-1
    D = list(range(3,n-3))
    U = list(range(4,n-3))
    def forward(cut):
        L,R = list(range(3,cut)),list(range(cut,n))
        return [[a,b]+L+R+[x],[b]+L+[a]+R+[x],[b,x]+L+[a]+R]
    def backward(start):
        V = list(range(start,n))
        tail = list(range(3,start))
        return [[a]+V+[b,x]+tail,[x,a]+V+[b]+tail,[x,a,b]+tail+V]
    if h == 2:
        if chord:
            return [('h2_xc',[],[[a,b,x]+D+[z,y,c],
                                  [a,x]+D+[z,y,b,c],[a,x]+D+[y,b,c,z]])]
        return [('h2_cx_3c',[(3,c,1)],[[z,a,b,x]+D+[y,c],
                    [z,a,x]+D+[b,y,c],[a,x]+D+[b,y,c,z]]),
                ('h2_cx_c3',[(3,c,0)],[[b,x,a,3]+U+[y,c,z],
                    [b,x,a,3,c,z]+U+[y],[a,b,x,3,c,z]+U+[y]])]
    if h == 3:
        W = list(range(4,n-3))
        R = [y,c,z]
        return [('h3_forward',[(1,3,1)],forward(y)),
                ('h3_backward_4b',[(1,3,0),(1,4,0)],backward(4)),
                ('h3_backward_b4',[(1,3,0),(1,4,1)],
                 [[a,b]+W+R+[x,3],[b]+W+[a]+R+[x,3],[b,x,3]+W+[a]+R])]
    return [(f'h{h}_forward',[(1,3,1)],forward(h)),
            (f'h{h}_backward',[(1,3,0)],backward(5 if h == 4 else 4))]

def arc(rows, a, b):
    return (rows[a] >> b) & 1

def valid(rows):
    n = len(rows)
    need(all(0 <= row < 1 << n for row in rows), 'row bit range')
    need(all(not arc(rows,a,a) for a in range(n)), 'no loops')
    need(all(arc(rows,a,b) ^ arc(rows,b,a) == 1 for a,b in combinations(range(n),2)), 'tournament')

def mismatch(rows, q):
    r = dict(zip(q, reversed(q)))
    for a, b in combinations(range(len(rows)),2):
        if arc(rows,a,b) != arc(rows,r[b],r[a]):
            return {'ordered_pair':[a,b], 'reflected_ordered_pair':[r[b],r[a]],
                    'arc_values':[arc(rows,a,b),arc(rows,r[b],r[a])]}
    return None

def move(rows, q):
    E = path(q)
    n = len(rows)
    return tuple(sum((arc(rows,a,b) ^ (pair(a,b) not in E)) << b
                     for b in range(n) if a != b) for a in range(n))

def twins(rows):
    n = len(rows)
    return [[a,b] for a,b in combinations(range(n),2)
            if all(arc(rows,a,v) == arc(rows,b,v) for v in range(n) if v not in (a,b))]

def literal_source(rows, I, opened, chord):
    n = len(rows)
    valid(rows)
    need(all(arc(rows,v,v+1) for v in range(n-1)), 'source path arcs')
    need(arc(rows,n-1,0) and mismatch(rows,list(range(n))) is None, 'source closed blue')
    need({v for v in range(1,n-1) if arc(rows,v,0)} == set(I), 'exact source inports')
    need({v for v in range(n-1) if arc(rows,n-1,v)} == {0}|{n-1-v for v in I}, 'exact reflected outports')
    need(rows[0].bit_count() == n-5 and rows[-1].bit_count() == 4, 'endpoint scores')
    need(arc(rows,1,n-2) == opened and arc(rows,2,n-2) == chord, 'literal branch bits')

def literal_ledger(rows, words, all_black):
    n = len(rows)
    target = move(rows,list(range(n)))
    trace, flags, xor = [], [], set()
    for q in words:
        need(sorted(q) == list(range(n)), 'literal word permutation')
        need(all(arc(rows,a,b) for a,b in zip(q,q[1:])), 'literal directed joins')
        witness = mismatch(rows,q)
        if all_black:
            need(witness is not None, 'actual word must be black')
        flags.append(witness is None)
        trace.append({'word':q,'source_rows':rows,'black_witness':witness,
                      'actual_first_inports':[v for v in range(n) if v != q[0] and arc(rows,v,q[0])],
                      'actual_last_outports':[v for v in range(n) if v != q[-1] and arc(rows,q[-1],v)],
                      'literal_twin_pairs':twins(rows)})
        after = move(rows,q)
        valid(after)
        need((mismatch(after,q) is None) == (witness is None), 'blue invariance across own K edge')
        rows = after
        xor ^= path(q)
    need(xor == path(list(range(n))) and rows == target, 'literal complete target and XOR')
    return flags, trace, target

def solver_controls():
    """Exhaustively cross-check this elimination/decoder on small raw systems."""
    state = 104729
    consistent = inconsistent = 0
    for case in range(96):
        s = Linear(4)
        rows = []
        for _ in range(case % 13):
            state = (1103515245 * state + 12345) & 0x7fffffff
            row = (state >> 8) & ((1 << 7)-1)
            rows.append(row)
            s.add(row)
        brute = {value for value in range(64)
                 if all(((row & 63 & value).bit_count() & 1) == bool(row & 64) for row in rows)}
        need(s.ok == bool(brute), 'brute Gaussian consistency control')
        if s.ok:
            decoded = {sum(arc(A,a,b) << i for i,(a,b) in enumerate(s.pairs)) for A in s.models()}
            need(decoded == brute and len(decoded) == 1 << s.dimension(), 'brute exhaustive model control')
            consistent += 1
        else:
            inconsistent += 1
    need(consistent and inconsistent, 'both solver control polarities')
    return {'raw_systems':96,'brute_assignments_per_system':64,
            'consistent_systems':consistent,'inconsistent_systems':inconsistent}

