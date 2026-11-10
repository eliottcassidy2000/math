#!/usr/bin/env python3
"""Independent small controls for primitive named deletion cubes.

Universe: every labelled tournament of orders 2..5, roots named 0 and 1.
Counts use literal permutations; marked classes quotient only the other names.
"""
from collections import Counter, defaultdict
from itertools import permutations
from math import comb, gcd
from functools import reduce
from hashlib import sha256
from pathlib import Path
import json

def require(value, message):
    if not value:
        raise RuntimeError(message)

def graph(n, bits):
    a = [[0]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1,n):
            a[i][j] = (bits >> k) & 1
            a[j][i] = 1-a[i][j]
            k += 1
    return tuple(map(tuple,a))

def h(a, vertices):
    return sum(all(a[p[i]][p[i+1]] for i in range(len(p)-1))
               for p in permutations(vertices))

def components(a):
    n = len(a)
    reach = []
    for v in range(n):
        seen = {v}; todo = [v]
        while todo:
            x = todo.pop()
            for y in range(n):
                if a[x][y] and y not in seen:
                    seen.add(y); todo.append(y)
        reach.append(seen)
    return {frozenset(w for w in range(n) if w in reach[v] and v in reach[w])
            for v in range(n)}

def canonical(a, vertices):
    # The two named roots are fixed, not quotiented against each other.
    rest = tuple(v for v in vertices if v not in (0,1))
    return min(tuple(a[p[i]][p[j]] for i in range(len(p)) for j in range(i+1,len(p)))
               for q in permutations(rest) for p in [(0,1)+q])

def rows(a):
    return [sum(x<<j for j,x in enumerate(row)) for row in a]

def main():
    classes = defaultdict(dict); tests = 0; primitive = 0
    same_component_zero = None
    for n in range(2,6):
        for bits in range(1<<(n*(n-1)//2)):
            a = graph(n,bits); V = tuple(range(n))
            c = tuple(h(a,tuple(v for v in V if not (mask >> v & 1)))
                      for mask in range(4))
            parts = components(a)
            touched = [p for p in parts if p & {0,1}]
            outside = [p for p in parts if not p & {0,1}]
            factor = 1
            for p in outside: factor *= h(a,tuple(p))
            d = reduce(gcd,c)
            require(d % factor == 0, 'unmarked factor divides whole cube')
            determinant = c[0]*c[3]-c[1]*c[2]
            if len(touched)==2:
                require(determinant==0,'separate roots give a rank-one tensor')
            if len(parts)==1 and determinant==0 and same_component_zero is None:
                same_component_zero = {'rows':rows(a),'deletion_cube':c,'order':n}
            if d==1:
                primitive += 1
                require(all(len(p)==1 for p in outside),'primitive forces unmarked singletons')
                keep = tuple(sorted(set().union(*touched)))
                k = len(touched)
                require(5**(len(keep)-k) <= c[0]**3, 'exact integer Busch core bound')
                if determinant:
                    require(k==1,'nonzero determinant forces one marked SCC')
                key = canonical(a,V)
                value = (c,canonical(a,keep),len(keep),k)
                require(key not in classes[n] or classes[n][key]==value,'canonical consistency')
                classes[n][key] = value
            tests += 1
    actual = defaultdict(Counter); kernels = defaultdict(set)
    for n,bank in classes.items():
        for c,key,b,k in bank.values():
            actual[c][n] += 1
            kernels[c].add((key,b,k))
    for c,polynomial in actual.items():
        predicted = Counter()
        for key,b,k in kernels[c]:
            for n in range(b,6):
                predicted[n] += comb(n-b+k,k)
        require(predicted==polynomial,('joint rational series',c))
    require(same_component_zero is not None,'rank-one converse hostile exists')
    # Explicit padding checks beyond the census; keep both named ports in one C3.
    base = graph(3,5) # 0->1->2->0.
    require(h(base,(0,1,2))==3,'cyclic positive control')
    padded_checks = 0
    for left in range(4):
        for right in range(4):
            # Reindex only after building: roots still vertices0 and1 in base.
            n = 3+left+right
            a = [[0]*n for _ in range(n)]
            order = list(range(3,3+left))+[0,1,2]+list(range(3+left,n))
            pos = {v:i for i,v in enumerate(order)}
            for i in range(n):
                for j in range(n):
                    a[i][j] = base[i][j] if i<3 and j<3 else int(pos[i]<pos[j])
            c = tuple(h(a,tuple(v for v in range(n) if not(mask >> v & 1))) for mask in range(4))
            require(c==(3,1,1,1),'joint padding')
            padded_checks += 1
    print(json.dumps({
        'status':'PASS', 'universe':'all labelled tournaments n=2..5, two distinct fixed roots0,1',
        'labelled_graphs':tests,'primitive_labelled_cubes':primitive,
        'primitive_marked_classes_by_order':{n:len(bank) for n,bank in classes.items()},
        'joint_rational_series_checked':len(actual),'padding_controls':padded_checks,
        'strong_rank_one_hostile':same_component_zero,
        'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest()
    },indent=2,sort_keys=True))

if __name__=='__main__':main()
