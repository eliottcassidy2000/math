#!/usr/bin/env python3
"""Exact finite controls; no producer imports, dependencies, or removable asserts."""
import itertools
import json
from collections import Counter


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


PAIRS = list(itertools.combinations(range(4), 2))


def algebra_checks():
    checked = triples = gauges = 0
    for parity in itertools.product((0, 1), repeat=4):
        for pattern in range(64):
            chi = [[0]*4 for _ in range(4)]
            for k, (i, j) in enumerate(PAIRS):
                chi[i][j] = 1 if pattern >> k & 1 else -1
                chi[j][i] = chi[i][j] * (-1)**(parity[i]*parity[j])
            members = [[i for i in range(4) if mask >> i & 1] for mask in range(16)]
            gauge = []
            for A in range(16):
                g = 1
                for i, j in itertools.combinations(members[A], 2):
                    g *= chi[i][j]
                gauge.append(g)
            table = []
            for A in range(16):
                row = []
                for B in range(16):
                    if A & B:
                        row.append((0, 0))
                        continue
                    s = 1
                    for i in members[A]:
                        for j in members[B]:
                            s *= chi[i][j]
                    row.append((s, A | B))
                    inversions = sum(parity[i]*parity[j]
                                     for i in members[A] for j in members[B] if i > j)
                    need(gauge[A]*gauge[B]*s*gauge[A | B] == (-1)**inversions,
                         'ordered basis removes all but the graded sign')
                    gauges += 1
                table.append(row)
            for A, B, C in itertools.product(range(16), repeat=3):
                s, AB = table[A][B]
                t, ABC = table[AB][C]
                u, BC = table[B][C]
                v, ABC2 = table[A][BC]
                need(s*t == u*v and (s*t == 0 or ABC == ABC2), 'associativity')
                triples += 1
            checked += 1
    return {'systems': checked, 'associativity_triples': triples, 'nonzero_gauge_pairs': gauges}


def tournament_checks():
    canonical = Counter()
    names = {(0,1,2,3): 'TT4', (1,1,1,3): 'source_triangle',
             (0,2,2,2): 'triangle_sink', (1,1,2,2): 'strong'}
    functions = {'TT4': (0,1), 'source_triangle': (0,0),
                 'triangle_sink': (1,1), 'strong': (1,0)}
    counts = Counter()
    witnesses = {}
    for word in range(64):
        adj = [[0]*4 for _ in range(4)]
        for k, (i,j) in enumerate(PAIRS):
            adj[i][j] = (word >> k) & 1
            adj[j][i] = 1-adj[i][j]
        orbit = [tuple(adj[p[i]][p[j]] for i,j in PAIRS)
                 for p in itertools.permutations(range(4))]
        canonical[min(orbit)] += 1
        name = names[tuple(sorted(map(sum, adj)))]
        reverse = names[tuple(sorted(3-sum(row) for row in adj))]
        f = functions[name]
        need(functions[reverse] == (1-f[1], 1-f[0]), 'converse equals flip conjugation')
        counts[name] += 1
        witnesses.setdefault(name, word)
    need(len(canonical) == 4 and sorted(canonical.values()) == [8,8,24,24], 'full orbit census')
    for a,b in itertools.product((0,1), repeat=2):
        f = tuple(((1+a+b)*x+a)%2 for x in (0,1))
        swapped = tuple(((1+a+b)*x+b)%2 for x in (0,1))
        need(swapped == (1-f[1],1-f[0]), 'parity swap equivariance')
    return {'classes':len(canonical), 'labelled_counts':dict(counts), 'least_mask_witnesses':witnesses}


def affine_checks():
    total = 0
    histogram = Counter()
    maps = tuple(itertools.product((0,1), repeat=2))  # (multiplier, constant)
    for length in range(1,7):
        for word in itertools.product(maps, repeat=length):
            M,C = 1,0
            values = [0,1]
            for a,b in word:
                M,C = a*M, (a*C+b)%2
                values = [(a*x+b)%2 for x in values]
            need(values == [C,(M+C)%2], 'literal versus affine composition')
            fixed = sum(values[x] == x for x in (0,1))
            need(fixed == (1 if M == 0 else 2*(1-C)), 'reset / identity / flip fixed points')
            need((M == 0) == any(a == 0 for a,b in word), 'reset persistence')
            total += 1
            histogram[f'{M},{C}'] += 1
    for bits in itertools.product((0,1),repeat=6):
        parity=sum(bits)%2
        solutions=0
        for x in (0,1):
            y=x
            for b in bits:y ^= b
            solutions += y == x
        need(solutions==2*(1-parity),'cycle holonomy versus sections')
    f=(0,2,1)
    need([i for i in range(3) if f[i]==i]==[0] and f[f[1]]==1 and f[1]!=1,
         'unique fixed point does not imply convergence')
    return {'words':total,'return_histogram':dict(histogram),'cycle_words':64,'fixed_point_hostile':f}


def trace(initial, rule):
    row=list(initial)
    out=[row[len(row)//2]]
    while len(row)>1:
        row=[(row[i-1] ^ (row[i] | row[i+1])) if rule==30
             else row[i-1]^row[i] for i in range(1,len(row)-1)]
        out.append(row[len(row)//2])
    return tuple(out)


def trace_checks():
    rows=[]
    total=0
    for rule in (30,60):
        for n in range(1,7):
            seen={}
            histogram=Counter()
            for integer in range(1 << (2*n+1)):
                initial=tuple((integer>>i)&1 for i in range(2*n+1))
                y=trace(initial,rule)
                key=(initial[n:],y)
                need(key not in seen,'same right boundary and trace have unique left input')
                seen[key]=initial[:n]
                histogram[y]+=1
            need(len(histogram)==1 << (n+1),'all trace words appear')
            need(set(histogram.values())=={1 << n},'exact uniform cylinder counts')
            total+=len(seen)
            rows.append({'rule':rule,'depth':n,'initial_blocks':len(seen),
                         'trace_words':len(histogram),'fibre_size':1 << n})
    seed=(0,)*32+(1,)+(0,)*32
    sixty=trace(seed,60);thirty=trace(seed,30)
    need(sixty==(1,)*33 and len(set(thirty))==2,'fixed-seed hostile differs despite free law')
    return {'initial_blocks':total,'cases':rows,'single_seed_rule30_prefix':thirty,
            'single_seed_rule60_prefix':sixty}


if __name__=='__main__':
    print(json.dumps({'status':'FINITE-EXACT PASS','algebras':algebra_checks(),
        'tournaments':tournament_checks(),'affine_returns':affine_checks(),
        'trace_cones':trace_checks(),
        'scope':'finite controls for the stated elementary proofs; no fixed-seed asymptotic conclusion'},
        indent=2,sort_keys=True))
