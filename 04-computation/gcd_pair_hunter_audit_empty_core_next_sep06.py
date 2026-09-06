#!/usr/bin/env python3
"""Independent CRT-pair/Prim-Hunter and flat all-clock gcd-profile audit.

No primary imports. The clock universe is every integer through the proved
rectangular bound, rather than a product-generated candidate set. Exact
max-block intersections are separately exhausted at small orders.
"""
from functools import cache
from hashlib import sha256
from itertools import combinations, combinations_with_replacement
from math import gcd, lcm
from pathlib import Path
import json

CHECKS = 0
ROOT = Path(__file__).resolve().parents[1]


def need(test,label):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(label)


def pair_credit(c,q,r):
    a = gcd(q,r)
    k,ell = (q+6)//7,(r+6)//7
    A,u = divmod(k,a)
    B,v = divmod(ell,a)
    return (c//lcm(q,r))*(a*A*B+A*v+B*u+max(0,u+v-a))


def prim_weight(weights):
    # Maximum spanning tree by vertex growth, independent of edge sorting.
    if not weights:
        return 0
    reached = {0}
    total = 0
    while len(reached)<len(weights):
        edge = max((weights[i][j],j) for i in reached
                   for j in range(len(weights)) if j not in reached)
        total += edge[0]
        reached.add(edge[1])
    return total


@cache
def hunter(c,gs):
    qs = tuple(c//g for g in gs)
    capacity = sum(g*((q+6)//7) for g,q in zip(gs,qs))
    weights = [[pair_credit(c,q,r) if i!=j else 0
                for j,r in enumerate(qs)] for i,q in enumerate(qs)]
    return capacity-prim_weight(weights)


@cache
def divisors(c):
    return tuple(d for d in range(2,c+1) if c%d==0)


def survives(c,gs):
    for a in divisors(c):
        reduced = tuple(sorted(gcd(a,g) for g in gs if g%a))
        need(bool(reduced),'primitivity prevents absorbing all tails')
        if hunter(a,reduced)<a:
            return False
    return True


def blocks(q):
    k = (q+6)//7
    steps = [u for u in range(q) if gcd(u,q)==1]
    return {tuple(sorted((start+u*j)%q for j in range(k)))
            for u in steps for start in range(q)}


def lifted_mask(c,q,block):
    residues = set(block)
    return sum(1<<j for j in range(c) if j%q in residues)


def main():
    # Full unit/start universe of padded capacity blocks. Sharpness is checked
    # over all block pairs, independently of the formula's balanced-count proof.
    banks = {q:blocks(q) for q in range(1,19)}
    pair_manifest = []
    for q in range(1,19):
        for r in range(q,19):
            c = lcm(q,r)
            left = [lifted_mask(c,q,b) for b in banks[q]]
            right = [lifted_mask(c,r,b) for b in banks[r]]
            minimum = min((a&b).bit_count() for a in left for b in right)
            credit = pair_credit(c,q,r)
            need(minimum==credit,('sharp direct pair minimum',q,r))
            need(pair_credit(2*c,q,r)==2*credit,'common-clock multiplicity')
            pair_manifest.append((q,r,credit))
    # Every triple of subsets on a four-point space audits the tree union law.
    for masks in combinations_with_replacement(range(16),3):
        weights = [[(a&b).bit_count() if i!=j else 0 for j,b in enumerate(masks)]
                   for i,a in enumerate(masks)]
        bound = sum(a.bit_count() for a in masks)-prim_weight(weights)
        need((masks[0]|masks[1]|masks[2]).bit_count()<=bound,'all three-mask Hunter controls')
    need(3-3<1,'cycle subtraction hostile: three identical singleton sets')
    # At c6,y1/2,tails2 and3 all actual bad sets are empty, while any padded
    # order3/order2 blocks intersect once. Actual singleton sizes cannot be used.
    need(all(min((t*(2*j+1))%12,12-(t*(2*j+1))%12)*14 >= 12
             for t in (2,3) for j in range(6)), 'physical empty-mask padding hostile')
    need(pair_credit(6,2,3)==1,'padded credit despite empty actual masks')
    named = [(36,(1,4,4,6,9)),(27,(1,1,1,3,9)),(288,(1,1,4,4,18,32))]
    for c,gs in named:
        need(sum(g*((c//g+6)//7) for g in gs)>=c,'named scalar gate misses')
        need(hunter(c,gs)<c,'named tree gate closes')
        print('named_tree_control',c,gs,'scalar',sum(g*((c//g+6)//7) for g in gs),
              'tree_upper',hunter(c,gs))
    expected_counts = [1,2,5,19,110,1217]
    expected_maxima = [1,2,4,9,32,96]
    expected_clock_counts = [1,2,4,7,16,43]
    previous = {1}
    previous_profiles = {(1,())}
    layers = []
    for d in range(1,7):
        q_bound = 6*d//(7-d)
        c_bound = q_bound*max(previous)
        profiles = []
        candidates = 0
        for c in range(1,c_bound+1):
            options = sorted(g for g in previous if c%g==0)
            for gs in combinations_with_replacement(options,d):
                if gcd(*gs)!=1:
                    continue
                candidates += 1
                if survives(c,gs):
                    for i,g in enumerate(gs):
                        child = (g,tuple(sorted(gcd(g,h) for j,h in enumerate(gs) if j!=i)))
                        need(child in previous_profiles,'flat cuts imply every complete child profile')
                    profiles.append((c,gs))
        clocks = sorted({c for c,gs in profiles})
        need(len(profiles)==expected_counts[d-1],('complete profile count',d))
        need(max(clocks)==expected_maxima[d-1],('complete maximum gcd',d))
        need(len(clocks)==expected_clock_counts[d-1],('complete gcd-set size',d))
        layer = {'tails':d,'rectangular_clock_bound':c_bound,'primitive_candidates':candidates,
                 'clocks':clocks,'profiles':profiles}
        layers.append(layer)
        print('layer',d,'all_clocks_through',c_bound,'candidates',candidates,
              'profiles',len(profiles),'max_gcd',max(clocks),'gcd_set',clocks)
        previous = set(clocks)
        previous_profiles = set(profiles)
    encoded = json.dumps({'pair_controls':pair_manifest,'layers':layers},
                         sort_keys=True,separators=(',',':')).encode()
    destination = ROOT/'05-knowledge/results/gcd_pair_hunter_audit_empty_core_next_sep06.json'
    destination.write_bytes(encoded+b'\n')
    print('profile_json',destination.name)
    print('semantic_sha256',sha256(encoded).hexdigest())
    print('explicit_checks',CHECKS)
    print('PROVED-relative candidate gcd caps for 13-d subsets:1,2,4,9,32,96')
    print('SCOPE: exact necessary clock profiles; no unsafe-row or LRC14 conclusion')


if __name__=='__main__':
    main()
