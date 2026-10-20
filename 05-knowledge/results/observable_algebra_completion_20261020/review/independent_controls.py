#!/usr/bin/env python3
"""Tiny exact finite-state algebra controls, independent of the source verifier.
No repository imports or writes; emits JSON to stdout.
"""
from fractions import Fraction as F
from itertools import product
from collections import Counter
import json

counts = Counter()
def require(name, truth):
    if not truth:
        raise RuntimeError('Failed independent control: ' + name)
    counts[name] += 1

def rank(rows):
    if not rows:
        return 0
    a = [[F(v) for v in row] for row in rows]
    r = 0
    for c in range(len(a[0])):
        pivot = next((i for i in range(r,len(a)) if a[i][c]),None)
        if pivot is None:
            continue
        a[r],a[pivot] = a[pivot],a[r]
        d = a[r][c]
        a[r] = [v/d for v in a[r]]
        for i in range(len(a)):
            if i != r:
                d = a[i][c]
                a[i] = [u-d*v for u,v in zip(a[i],a[r])]
        r += 1
        if r == len(a):
            break
    return r

def basis(rows):
    out = []
    for row in rows:
        if rank(out+[row]) > len(out):
            out.append(row)
    return out

def algebra(n, observables, transitions):
    out = basis([[1]*n]+observables)
    while True:
        candidates = out + [[f[t[x]] for x in range(n)] for f in out for t in transitions]
        candidates += [[u*v for u,v in zip(a,b)] for i,a in enumerate(out) for b in out[i:]]
        nxt = basis(candidates)
        if len(nxt) == len(out):
            return out
        out = nxt

def quotient(n, observables, transitions):
    signatures = [tuple(f[x] for f in observables) for x in range(n)]
    def labels(sig):
        keys = sorted(set(sig))
        return [keys.index(x) for x in sig]
    old = labels(signatures)
    while True:
        new = labels([(old[x],)+tuple(old[t[x]] for t in transitions) for x in range(n)])
        if all((old[x]==old[y]) == (new[x]==new[y]) for x in range(n) for y in range(n)):
            return new
        old = new

# Finite nonlinear example: four observational classes among five states.
t = [1,2,2,4,0]
obs = [[0,1,1,0,1]]
alg = algebra(5,obs,[t]); blocks = quotient(5,obs,[t])
require('finite_algebra_quotient', len(alg) == len(set(blocks)) == 4)
require('finite_algebra_quotient', blocks[1] == blocks[2] and len({blocks[i] for i in (0,1,3,4)}) == 4)
for block in set(blocks):
    indicator = [int(v==block) for v in blocks]
    require('quotient_indicators_in_algebra', rank(alg+[indicator]) == len(alg))

# Three-bit cyclic register: state-observable, but timewise means have rank four.
X = list(product((0,1),repeat=3))
shift = [X.index((x[1],x[2],x[0])) for x in X]
first = [(-1)**x[0] for x in X]
f = first[:]; time_rows = [[1]*8]
for _ in range(6):
    time_rows.append(f)
    f = [f[shift[i]] for i in range(8)]
alg = algebra(8,[first],[shift]); blocks = quotient(8,[first],[shift])
require('cyclic_register', len(alg) == len(set(blocks)) == 8)
require('cyclic_register', rank(time_rows) == 4)
parity = [(-1)**sum(x) for x in X]
even = [F(1+p,8) for p in parity]; odd = [F(1-p,8) for p in parity]
require('same_timewise_means', sum(even) == sum(odd) == 1 and all(a*b==0 for a,b in zip(even,odd)))
require('same_timewise_means', all(sum((a-b)*v for a,b,v in zip(even,odd,row)) == 0 for row in time_rows))
require('joint_observation_distinguishes', all(time_rows[1][i]*time_rows[2][i]*time_rows[3][i] == parity[i] for i in range(8)))
require('joint_observation_distinguishes', sum(p*w for p,w in zip(parity,even)) == 1 and sum(p*w for p,w in zip(parity,odd)) == -1)

# Binary affine system with one invertible and one singular transition.
matrices = [((0,1,0),(0,0,1),(1,0,0)), ((1,0,0),(0,0,0),(0,0,1))]
translations = [(1,0,1),(0,1,0)]
def dot(a,b): return sum(x*y for x,y in zip(a,b)) % 2
def transpose_label(A,l): return tuple(sum(A[j][i]*l[j] for j in range(3))%2 for i in range(3))
def affine(A,b,x): return tuple((dot(row,x)+v)%2 for row,v in zip(A,b))
T = [[X.index(affine(A,b,x)) for x in X] for A,b in zip(matrices,translations)]
labels = {(0,0,0),(1,0,0)}
while True:
    enlarged = labels | {transpose_label(A,l) for A in matrices for l in labels} | {tuple((u+v)%2 for u,v in zip(a,b)) for a in labels for b in labels}
    if enlarged == labels:
        break
    labels = enlarged
chars = {l:[(-1)**dot(l,x) for x in X] for l in labels}
require('affine_character_span', len(labels) == 8 and rank(list(chars.values())) == 8)
require('affine_character_span', len(algebra(8,[first],T)) == len(labels))
for A,b,tmap in zip(matrices,translations,T):
    for l,f in chars.items():
        image = chars[transpose_label(A,l)]
        require('affine_signed_pullback', [f[tmap[i]] for i in range(8)] == [(-1)**dot(l,b)*v for v in image])

# Strict lower-rank case: observe first coordinate, no mixing of labels.
A = ((1,0,0),(0,0,0),(0,0,1)); b = (1,1,0)
tmap = [X.index(affine(A,b,x)) for x in X]
require('annihilator_cosets', len(algebra(8,[first],[tmap])) == 2)
blocks = quotient(8,[first],[tmap])
require('annihilator_cosets', all((blocks[i]==blocks[j]) == (X[i][0]==X[j][0]) for i in range(8) for j in range(8)))

# Scope guard: arbitrary Fourier sums need not expose each support character.
Y = list(product((0,1),repeat=2))
combined = [(-1)**x[0]+(-1)**x[1] for x in Y]
require('individual_character_scope_guard', len(algebra(4,[combined],[list(range(4))])) == 3)
require('individual_character_scope_guard', len(algebra(4,[[(-1)**x[0] for x in Y],[(-1)**x[1] for x in Y]],[list(range(4))])) == 4)

print(json.dumps({'status':'PASS','exact_control_total':sum(counts.values()),'checks':dict(sorted(counts.items())), 'scope':'Tiny exact finite-state controls, independent rational linear algebra and independent partition refinement. No broad census or source verifier imports.'},indent=2,sort_keys=True))
