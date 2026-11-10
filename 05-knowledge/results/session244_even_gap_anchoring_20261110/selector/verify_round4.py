#!/usr/bin/env python3
"""Two-index anchors, all-tail normal form, finite fibers and graph rebasing.

All checks use explicit exceptions and remain active under Python -O.
"""
from collections import deque
from fractions import Fraction
import importlib.util
import itertools
import json
import math
from pathlib import Path

p = Path(__file__).with_name('verify_round2.py')
sp = importlib.util.spec_from_file_location('round2', p)
r = importlib.util.module_from_spec(sp)
sp.loader.exec_module(r)
v = r.v
require = v.require
LOW = 'OJJJJOZ'
HIGH = 'OJOZJOZ'


def run(s, word):
    for c in word:
        s = r.jstep(s) if c == 'J' else v.step(s, c)
    return s


def graph(a, root, word):
    for c in word:
        a, root = r.graph_step(a, root, c)
    return a, root


def tail(m):
    # Chronological least-significant-first reading of m-1.
    return ''.join('Z' if b == '1' else 'O' for b in bin(m-1)[:1:-1]) if m > 1 else ''


def lower(m, n):
    return (18*m+3*v.g(m)-4)*v.g(n) + (44*m+6*v.g(m)-10)*n


def block_data(word):
    e = [int(c == 'Z') for c in word]
    a = sum(x * 2**i for i, x in enumerate(e))
    c = sum(e[i]*e[j]*2**(i-j-1) for i in range(len(e)) for j in range(i))
    c += sum(e[i]*e[i-1] for i in range(1, len(e)))
    return a, c


def general_tail(s, m):
    d0, e, u, d = s
    a = (m-1) % 2
    return d0 + (m-1)*u - (m-v.g(m)-a)*d - a*e


def fiber(cap=27783):
    start = (1, 0, 2, 1)
    best = {start: 0}
    q = deque([start])
    while q:
        s = q.popleft()
        for c in 'OZJ':
            cost = best[s] + (c == 'J')
            t = run(s, c)
            require(t[0]+2*t[2] > s[0]+2*s[2], 'strict balanced readout')
            if cost > 3 or t[0]+2*t[2] > cap:
                continue
            if cost < best.get(t, 4):
                best[t] = cost
                (q.append if c == 'J' else q.appendleft)(t)
    f = sorted((s, j) for s, j in best.items() if s[0]+2*s[2] == cap)
    return len(best), f


def modular_orbit(mod):
    start = (2 % mod, 1)
    seen = {start}
    todo = [start]
    while todo:
        n, g = todo.pop()
        for t in [((2*n) % mod, g), ((2*n-1) % mod, (g+n-n % 2) % mod)]:
            if t not in seen:
                seen.add(t)
                todo.append(t)
    return seen


def main():
    basis = [tuple(int(i == j) for i in range(4)) for j in range(4)]
    lm = tuple(zip(*(run(s, LOW) for s in basis)))
    hm = tuple(zip(*(run(s, HIGH) for s in basis)))
    require(lm == ((17,0,40,0),(3,0,6,0),(21,0,50,0),(3,0,6,0)), 'lower matrix')
    require(hm == ((17,0,42,0),(3,0,6,0),(21,0,50,0),(3,0,6,0)), 'upper matrix')
    blocks = 0
    seeds = [(1,0,2,1),(15,1,59,9),(13,2,18,3),(97,15,121,15)]
    for k in range(12):
        for w in itertools.product('OZ', repeat=k):
            w = ''.join(w)
            a, c = block_data(w)
            m = a+1
            require(c+a % 2 == m-v.g(m), ('quadratic collapse', w))
            for s in seeds:
                require(run(s,w)[0] == general_tail(s,m), ('general tail',s,w))
            blocks += 1
    for n in range(1,129):
        s = (v.g(n), n % 2, n, 1)
        lo, hi = run(s, LOW), run(s, HIGH)
        require(tuple(y-x for x,y in zip(lo,hi)) == (2*n,0,0,0), ('pure D',n))
        for m in range(1,129):
            w = tail(m)
            require(run(lo,w)[0] == lower(m,n), ('lower',m,n))
            require(run(hi,w)[0] == lower(m,n)+2*n, ('upper',m,n))
            require(lower(m,n) >= 40*m*n+17, ('inverse bound',m,n))
            qn = math.gcd(n,3*v.g(n))
            require((lower(m,n)+v.g(n)) % (2*qn) == 0, ('fixed-input residue',m,n))
            common = math.gcd(n,v.g(n))
            if common > 1:
                require(lower(m,n) % common == 0 and lower(m,n) > common, ('prime filter',m,n))
    literal = 0
    for n in range(1,5):
        a = (0,) if n == 1 else v.graph_word(v.word(n))
        for m in range(1,5):
            for word, inc in [(LOW,0),(HIGH,2*n)]:
                b, root = graph(a,0,word+tail(m))
                prof = v.profile(b,root)
                require(prof[0] == lower(m,n)+inc, ('literal anchor',m,n,word))
                literal += 1
    mod_sizes = {}
    for mod in range(2,65,2):
        orbit = modular_orbit(mod)
        expected = {(n,g) for n in range(mod) for g in range(1,mod,2)}
        require(orbit == expected, ('inherited joint modular theorem',mod))
        image = {(17*g+40*n) % mod for n,g in orbit}
        require(image == set(range(1,mod,2)), ('new family modular consequence',mod))
        mod_sizes[mod] = len(orbit)
        for n in range(1,33):
            gn = v.g(n)
            a,b,c = 18*gn+44*n,3*gn+6*n,-4*gn-10*n
            image_n = {(a*m+b*gm+c) % mod for m,gm in orbit}
            congruence = math.gcd(2*math.gcd(n,3*gn),mod)
            require(image_n == {x for x in range(mod) if (x+gn) % congruence == 0}, ('exact fixed-input modular image',n,mod))
    # All M,N are included for each target under this cap by the proved MN bound.
    cap = 2001
    image = set()
    params = []
    for n in range(1,(cap-17)//40+1):
        for m in range(1,(cap-17)//(40*n)+1):
            h = lower(m,n)
            for e in [0,1]:
                x = h+2*e*n
                if x <= cap:
                    image.add(x)
                    params.append((m,n,e,x))
    require(57 in image and 59 in image, 'small boundary controls')
    require(61 not in image and 173 not in image, 'all-order family hostiles')
    # A rational certificate for H_32768 < 11, then the bank density budget.
    h32 = sum((Fraction(1,j) for j in range(1,33)), Fraction(0))
    log2_upper = sum((Fraction(2,(2*j+1)*3**(2*j+1)) for j in range(3)), Fraction(0)) + Fraction(1,6804)
    require(log2_upper < Fraction(347,500), 'rational logarithm bound')
    hbound = h32+10*Fraction(347,500)
    require(hbound < 11, 'harmonic bank bound')
    pair_coefficient = 2*(Fraction(1,94)+Fraction(1,88))
    require(pair_coefficient == Fraction(91,2068), 'density tail bound')
    density = 11*pair_coefficient
    omitted = 1-2*density
    require(density == Fraction(1001,2068) and omitted == Fraction(33,1034), 'bank omission')
    state_count, f = fiber()
    us = sorted({s[2] for s,j in f})
    require(state_count == 216030 and len(f) == 33 and len(us) == 20, 'complete fiber size')
    require(not any(u+1 in us for u in us), 'fiber has no consecutive responses')
    # Actual rooted isomorphisms, not equality of scalar recurrences.
    isomorphisms = 0
    for length in range(4):
        for w in itertools.product('OZJ', repeat=length):
            a, root = graph((2,0),0,''.join(w))
            old = len(a)
            for k in range(9):
                b, mark = graph(a,root,'J'*k)
                order = list(range(len(b)))
                for j in range(0,k-1,2):
                    order[old+j],order[old+j+1] = order[old+j+1],order[old+j]
                pos = {vertex:i for i,vertex in enumerate(order)}
                require(all(pos[i] == i for i in range(old)), 'old graph fixed')
                expected = set()
                if k:
                    expected.add(frozenset((root,old)))
                    for j in range(1,k-1,2):
                        expected.add(frozenset((old+j,old+j+1)))
                reversals = {frozenset((i,j)) for i in range(len(b)) for j in range(i+1,len(b))
                             if j >= old and bool(b[j]>>i&1) == (pos[i]<pos[j])}
                require(reversals == expected, ('matching reversal isomorphism',w,k))
                degrees = {}
                for edge in reversals:
                    for vertex in edge:degrees[vertex] = degrees.get(vertex,0)+1
                require(max(degrees.values(),default=0) <= 1, 'matching degrees')
                expected_mark = root if not k else old+k-1-(k%2 == 0)
                require(pos[mark] == expected_mark, 'root transported, not reset')
                isomorphisms += 1
    print(json.dumps({
        'status':'PROVED two-index form, complete inverse, bounded-bank density obstruction; finite controls',
        'binary_blocks':blocks,'general_seed_profiles':seeds,
        'two_index_arithmetic_checks':128**2,'literal_anchor_profiles':literal,
        'lower_matrix':lm,'upper_matrix':hm,
        'modular_source':'INHERITED Session199/congruence-intervals/PROOF.md Theorem 1',
        'modular_orbit_sizes':mod_sizes,'fixed_input_modular_image_checks':32**2,
        'complete_anchor_image_cap':cap,'complete_anchor_image_size':len(image),
        'complete_anchor_parameter_records':len(params),'all_order_family_hostiles':[61,173],
        'density_bank_max_N':32768,'density_upper':str(density),
        'odd_omitted_fraction_lower':str(omitted),'harmonic_upper_certificate':str(hbound),
        'balanced_fiber_target':27783,'source_states_under_target':state_count,
        'balanced_fiber_profiles':f,'balanced_fiber_responses':us,
        'balanced_fiber_has_consecutive_responses':False,
        'literal_rooted_matching_isomorphisms':isomorphisms,
        'unresolved':'Unbounded N family coverage/density, four-block all-target lemma, and Chan-Pak remain OPEN.'
    },indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
