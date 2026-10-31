#!/usr/bin/env python3
"""Bounded Session126 addendum: A5 subgroups, F/U/S, and exact seed coalescence.

All gates survive python -O. Group/arithmetic/program-prefix checks are stdlib,
from definitions. The all-time seed2 balance conclusion is an explicit imported
proof dependency, not inferred from the prefix replay in this file.
"""
from collections import Counter, deque
from itertools import combinations, permutations
from pathlib import Path
import hashlib
import json

gates = 0
def need(ok, label):
    global gates
    gates += 1
    if not ok:
        raise RuntimeError(label)

def factor(n):
    out = {}
    d = 2
    while d*d <= n:
        while n % d == 0:
            out[d] = out.get(d,0)+1
            n //= d
        d += 1
    if n > 1:
        out[n] = out.get(n,0)+1
    return out

def prime(n):
    return n > 1 and factor(n) == {n:1}

def fus(n):
    need(n >= 2, 'FUS domain')
    divs = [d for d in range(2,n) if n % d == 0]
    ps = [d for d in divs if prime(d)]
    sq = [d for d in divs if all(e == 1 for e in factor(d).values())]
    return {'N':n,'factors':factor(n),'proper_divisors':divs,
            'prime_divisors':ps,'squarefree_divisors':sq,
            'F':len(divs),'U':len(ps),'S':len(sq),
            'defect':len(divs)-len(ps)-len(sq)}

def constructible(n):
    return all(p == 2 or (e == 1 and (p-1)&(p-2) == 0)
               for p,e in factor(n).items())

P = [p for p in permutations(range(5))
     if sum(p[i] > p[j] for i,j in combinations(range(5),2)) % 2 == 0]
idx = {p:i for i,p in enumerate(P)}
identity = idx[tuple(range(5))]
def compose(a,b):
    return tuple(a[b[i]] for i in range(5))
mul = [[idx[compose(a,b)] for b in P] for a in P]
inv = [next(j for j in range(60) if mul[i][j] == identity) for i in range(60)]
need(len(P) == 60, 'A5 order')

def closure(generators):
    seen = {identity}
    todo = [identity]
    while todo:
        x = todo.pop()
        for g in generators:
            y = mul[x][g]
            if y not in seen:
                seen.add(y)
                todo.append(y)
    return frozenset(seen)

# Every subgroup can be reached by repeatedly adjoining one generator.
initial = frozenset([identity])
subgroups = {initial:()}
queue = deque([initial])
while queue:
    H = queue.popleft()
    gens = subgroups[H]
    for g in range(60):
        if g not in H:
            K = closure(gens+(g,))
            if K not in subgroups:
                subgroups[K] = gens+(g,)
                queue.append(K)
hist = dict(sorted(Counter(map(len,subgroups)).items()))
need(hist == {1:1,2:15,3:10,4:5,5:6,6:10,10:6,12:5,60:1},
     'complete A5 subgroup census')
for H in subgroups:
    need(60 % len(H) == 0, 'subgroup Lagrange check')
    need(all(mul[h][k] in H for h in H for k in H), 'subgroup multiplication closure')
    need(all(inv[h] in H for h in H), 'subgroup inverse closure')

def cycle(*labels):
    out = list(range(5))
    for a,b in zip(labels,labels[1:]+labels[:1]):
        out[a] = b
    return tuple(out)

a = idx[cycle(0,1,2)]
s = idx[compose(cycle(0,1),cycle(3,4))]
S3 = closure((a,s))
need(len(S3) == 6 and mul[s][mul[a][s]] == inv[a], 'S3 embedding')
a5cycle = idx[cycle(0,1,2,3,4)]
dreflection = idx[compose(cycle(1,4),cycle(2,3))]
need(len(closure((a5cycle,dreflection))) == 10 and
     mul[dreflection][mul[a5cycle][dreflection]] == inv[a5cycle], 'D10 embedding')
C3 = closure((a,))
def conjugate(g,h):
    return mul[g][mul[h][inv[g]]]
normalizer = [g for g in range(60)
              if frozenset(conjugate(g,h) for h in C3) == C3]
need(frozenset(normalizer) == S3, 'normalizer of C3 is explicit S3')
need(len({frozenset(conjugate(g,h) for h in C3) for g in range(60)}) == 10,
     'ten Sylow3 subgroups')
for H in subgroups:
    if len(H) in (3,5):
        normalizer_size = sum(frozenset(conjugate(g,h) for h in H) == H
                              for g in range(60))
        need(normalizer_size == 2*len(H), 'C3/C5 normalizers')

conjugacy = {frozenset(conjugate(g,h) for g in range(60)) for h in range(60)}
class_sizes = sorted(map(len,conjugacy))
need(class_sizes == [1,12,12,15,20], 'A5 conjugacy sizes')
normal_candidates = sorted({1+sum(chosen) for r in range(5)
                           for chosen in combinations(class_sizes[1:],r)
                           if 60 % (1+sum(chosen)) == 0})
need(normal_candidates == [1,60], 'simplicity by normal class unions and Lagrange')

rotation = [fus(n) for n in (12,24,60)]
full = [fus(n) for n in (24,48,120)]
need([r['defect'] for r in rotation] == [-1,1,0], 'rotation balance table')
need([r['defect'] for r in full] == [1,3,4], 'full symmetry balance table')
f60 = fus(60)
realized = sorted(k for k in hist if 1 < k < 60)
missing = sorted(set(f60['proper_divisors'])-set(realized))
need(realized == [2,3,4,5,6,10,12] and missing == [15,20,30], 'seven/three orders')
need(f60['squarefree_divisors'] == [2,3,5,6,10,15,30], 'squarefree seven')
need(f60['squarefree_divisors'] != realized, 'two sevens are different sets')
need(all(not(sorted(factor(n).values()) == [1,1,2]) for n in range(2,60)),
     'no smaller mixed balanced integer')

fractions = [(17,91),(78,85),(19,51),(23,38),(29,33),(77,29),(95,23),
             (77,19),(1,17),(11,13),(13,11),(15,2),(1,7),(55,1)]
def step(n):
    for a,b in fractions:
        if n % b == 0:
            return n//b*a
    raise RuntimeError('unreachable: final denominator1')
def trace(n,steps):
    states = [n]
    for _ in range(steps):
        states.append(step(states[-1]))
    return states
standard = trace(2,145)
need(standard[1] == 15, '15 is standard seed next state')
need(trace(23,11)[-1] == 2, '23 joins seed2')
need(trace(16,14)[-1] == standard[145] == 97540625, '16 joins standard C(5,4)')
need(trace(17,4) == [17,1,55,65,55], '17 eventually two-cyclic')
need(trace(25,3) == [25,1375,1625,1375], '25 eventually two-cyclic')
need(standard[69:72] == [8,60,450], '60 balance exit witness')
need(fus(450)['defect'] != 0, 'balance is not preserved from60')

matrix = []
for n in (15,16,17,23,25):
    ar = fus(n)
    matrix.append({'N':n,'factors':ar['factors'],'FUS':[ar['F'],ar['U'],ar['S']],
                   'constructible_polygon':constructible(n),
                   'joins_standard_primegame_orbit':n in (15,16,23),
                   'balanced':ar['defect'] == 0})
need([r['N'] for r in matrix if r['constructible_polygon']] == [15,16,17],
     'constructible triple')
need([r['N'] for r in matrix if r['balanced']] == [17,23], 'balanced pair')

ROOT = Path(__file__).parent
pins = []
for rel in ['dependencies/balance_dynamics/README.md','dependencies/balance_dynamics/seed_balance_window.json',
            'dependencies/fractran/FRACTRAN_PRIME_EXPONENT_REPORT.md','dependencies/fractran/input_seed_results.json']:
    path = ROOT/rel
    need(path.exists(), 'explicit imported dependency '+rel)
    pins.append({'path':rel,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
balance = json.loads((ROOT/'dependencies/balance_dynamics/seed_balance_window.json').read_text())
need(balance['balanced_visit_count'] == 16, 'imported complete balance count')
need(balance['balanced_visits'][-1]['n'] == 60 and
     balance['balanced_visits'][-1]['step'] == 70, 'imported last balance event')

out = {'status':'PASS','checks':gates,'A5_subgroup_order_histogram':hist,
       'A5_conjugacy_class_sizes':class_sizes,'possible_normal_subgroup_orders':normal_candidates,
       'proper_nontrivial_subgroup_order_values':realized,'missing_proper_divisor_orders':missing,
       'C3_normalizer_order':len(normalizer),'Sylow3_subgroup_count':10,
       'S3_generators_zero_based':[P[a],P[s]],
       'rotational_group_balance':rotation,'full_symmetry_group_balance':full,
       'D60_squarefree_divisors':f60['squarefree_divisors'],'five_input_features':matrix,
       'seed_coalescence':{'15':{'seed2_update':1},'23':{'updates_to_seed2':11},
                          '16':{'updates_to_common_state':14,'seed2_update':145,'common_state':97540625}},
       'seed2_last_balanced_state':{'N':60,'update':70,'status':'IMPORTED PROOF; exact prefix checked'},
       'source_pins':pins}
print(json.dumps(out,indent=2,sort_keys=True))
