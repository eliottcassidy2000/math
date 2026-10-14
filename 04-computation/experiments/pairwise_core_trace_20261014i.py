"""Signed injective graph sums via a guarded three-matrix trace.

Inspired by Lemma6 of attached arXiv2610.15426v1. The implementation uses
classical exact-object matrix multiplication, not a fast-omega backend.
"""
from itertools import combinations, permutations
from math import isqrt, prod
from pathlib import Path
import json
import numpy as np


def exact_int(x, minimum=0):
    if type(x) is not int or x < minimum:
        raise ValueError('exact integer out of range')


def prime(q):
    return q >= 2 and all(q % d for d in range(2, isqrt(q) + 1))


def validate(q, n, edges):
    exact_int(q, 3); exact_int(n, 2)
    if q % 2 == 0 or not prime(q): raise ValueError('odd prime required')
    if type(edges) is not tuple: raise ValueError('tuple edges required')
    for edge in edges:
        if type(edge) is not tuple or len(edge) != 2: raise ValueError('edge pair required')
        a, b = edge
        exact_int(a); exact_int(b)
        if a >= n or b >= n or a == b: raise ValueError('distinct labelled vertices required')


def character(q):
    squares = {x*x % q for x in range(1, q)}
    return [0] + [1 if x in squares else -1 for x in range(1, q)]


def trace_sum(q, n, edges, max_entries=2_000_000, cross_guard=True):
    """Returns full injective signed sum. cross_guard=False is hostile-only."""
    validate(q, n, edges); exact_int(max_entries, 1)
    if type(cross_guard) is not bool: raise ValueError('Boolean guard flag required')
    if n > q: return 0, {'empty': True}
    if len(edges) % 2: return 0, {'odd_edge_scaling_cancellation': True}
    chi = character(q)
    # Affine normalization fixes two marked, distinct vertices, not just an edge.
    sizes = [(n-2)//3 + (i < (n-2)%3) for i in range(3)]
    blocks, first = [], 2
    for size in sizes:
        blocks.append(tuple(range(first,first+size))); first += size
    count = lambda size: prod(range(q-2-size+1,q-1))
    counts = [count(size) for size in sizes]
    entries = sum(counts[i]*counts[(i+1)%3] for i in range(3))
    if entries > max_entries: raise ValueError('explicit matrix-entry cap exceeded')
    states = [[dict(zip(block,state)) for state in permutations(range(2,q),len(block))]
              for block in blocks]
    owner = {v:i for i,block in enumerate(blocks) for v in block}
    fixed = {0:0,1:1}
    constant, internal, crossing = 1, [[],[],[]], [[],[],[]]
    cycle_pairs = [(0,1),(1,2),(2,0)]
    for a,b in edges:
        deps = {owner[v] for v in (a,b) if v >= 2}
        if not deps: constant *= chi[(fixed[b]-fixed[a])%q]
        elif len(deps)==1: internal[next(iter(deps))].append((a,b))
        else:
            pair = next(i for i,p in enumerate(cycle_pairs) if set(p)==deps)
            crossing[pair].append((a,b))
    unary = []
    for i in range(3):
        vals=[]
        for state in states[i]:
            values=fixed|state
            vals.append(prod(chi[(values[b]-values[a])%q] for a,b in internal[i]))
        unary.append(vals)
    matrices=[]
    for pair,(i,j) in enumerate(cycle_pairs):
        matrix=np.empty((len(states[i]),len(states[j])),dtype=object)
        for ai,sa in enumerate(states[i]):
            for bj,sb in enumerate(states[j]):
                if cross_guard and set(sa.values()) & set(sb.values()):
                    matrix[ai,bj]=0; continue
                values=sa|sb
                matrix[ai,bj]=unary[i][ai]*prod(chi[(values[b]-values[a])%q]
                                               for a,b in crossing[pair])
        matrices.append(matrix)
    ab=matrices[0]@matrices[1]
    normalized=int(np.sum(ab*matrices[2].T))
    return constant*q*(q-1)*normalized, {
        'block_sizes':sizes,'state_counts':counts,'matrix_entries':entries,
        'normalized_sum':constant*normalized,'cross_distinctness':cross_guard}


def brute(q,n,edges):
    """Independent direct affine-normalized sum; even edges only."""
    validate(q,n,edges)
    if n > q: return 0
    if len(edges)%2: raise ValueError('brute affine control requires even edges')
    chi=character(q); total=0
    for tail in permutations(range(2,q),n-2):
        x=(0,1)+tail
        total += prod(chi[(x[b]-x[a])%q] for a,b in edges)
    return total*q*(q-1)


WHEEL=tuple([(0,i) for i in range(1,5)]+[(1,2),(2,3),(3,4),(4,1)])
NEXT_CORE=tuple(e for e in combinations(range(6),2)
                if e not in {(0,1),(1,2),(2,3),(3,4),(4,5)})


def main():
    checks=0
    def check(condition):
        nonlocal checks
        checks+=1
        if not condition: raise RuntimeError(f'failed check {checks}')
    rows=[]
    for name,n,edges in [('wheel',5,WHEEL),('K6-minus-P6',6,NEXT_CORE),
                         ('K4',4,tuple(combinations(range(4),2)))]:
        for q in (3,5,7,11,19,31):
            value,receipt=trace_sum(q,n,edges)
            check(value==brute(q,n,edges))
            reverse=(tuple(reversed(edges[0])),)+edges[1:]
            reverse_value,_=trace_sum(q,n,reverse)
            check(reverse_value==(-1 if q%4==3 else 1)*value)
            total=prod(range(q-n+1,q+1)) if n<=q else 0
            check((total+value)%2==0 and abs(value)<=total)
            rows.append({'name':name,'q':q,'sum':value,'positive':(total+value)//2,
                         'negative':(total-value)//2,'receipt':receipt})
    for q,value in ((7,336),(11,-12320),(19,-76608)):
        check(next(r['sum'] for r in rows if r['name']=='K6-minus-P6' and r['q']==q)==value)
    # Omitting cross-block inequality is a real quotient, not a harmless shortcut.
    bad,_=trace_sum(7,5,WHEEL,cross_guard=False)
    good,_=trace_sum(7,5,WHEEL)
    check(bad!=good)
    check(next(r for r in rows if r['name']=='K4' and r['q']==7)['positive']==420)
    # Keep semantic/API guards, including exact prime field and resource bounds.
    for action in [lambda:trace_sum(9,5,WHEEL),lambda:trace_sum(True,5,WHEEL),
                   lambda:trace_sum(7,5,WHEEL,max_entries=1),
                   lambda:trace_sum(7,5,((0,0),)),
                   lambda:trace_sum(7,5,WHEEL,cross_guard=1)]:
        try: action()
        except ValueError: check(True)
        else: check(False)
    # Independent audit: full injections, without affine normalization or trace.
    full_controls=0
    for q in (3,5,7):
        chi=character(q)
        for n in (2,3,4):
            available=tuple(combinations(range(n),2))
            for mask in range(1<<len(available)):
                edges=tuple(e for i,e in enumerate(available) if mask>>i&1)
                full=sum(prod(chi[(x[b]-x[a])%q] for a,b in edges)
                         for x in permutations(range(q),n))
                check(trace_sum(q,n,edges)[0]==full)
                full_controls+=1
    # Two parallel edges square to one on every injective assignment.
    check(trace_sum(7,3,((0,2),(0,2)))[0]==7*6*5)
    check(trace_sum(7,2,((0,1),(1,0)))[0]==-7*6)
    out={'checks':checks,'status':'PROVED trace identity; FINITE-EXACT stated prime graph universe',
         'arithmetic':'Python integers/object arrays; classical matrix backend',
         'independent_full_injection_controls':full_controls,
         'rows':rows,'hostile_missing_cross_guard':{'q':7,'guarded':good,'unguarded':bad}}
    path=Path(__file__).resolve().parents[2]/'05-knowledge/results/pairwise_core_trace_20261014i.json'
    path.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(out,indent=2))


if __name__=='__main__': main()
