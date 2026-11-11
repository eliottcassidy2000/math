#!/usr/bin/env python3
"""Independent product, explicit-vector-subspace, and literal-path audit."""
import ast
import hashlib
import itertools
import json
from collections import deque
from functools import lru_cache, reduce
from math import gcd
from pathlib import Path

HERE = Path(__file__).resolve().parent
TARGET = HERE.parent/'golden'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def main():
    source = TARGET/'check.py'
    assignments = {}
    for node in ast.parse(source.read_text()).body:
        if isinstance(node, ast.Assign) and len(node.targets)==1 and isinstance(node.targets[0],ast.Name):
            name=node.targets[0].id
            if name in ('M','WORDS'):
                assignments[name]=ast.literal_eval(node.value)
    matrices, words = assignments['M'],assignments['WORDS']
    receipt = json.loads((TARGET/'CHECK.json').read_text())
    ident=((1,0,0),(0,1,0),(0,0,1))
    def multiply(a,b):
        return tuple(tuple(sum(x*y for x,y in zip(row,col)) for col in zip(*b)) for row in a)
    def vector(v,a):
        return tuple(sum(x*y for x,y in zip(v,col)) for col in zip(*a))
    def determinant(a):
        x,y,z=a
        return x[0]*(y[1]*z[2]-y[2]*z[1])-x[1]*(y[0]*z[2]-y[2]*z[0])+x[2]*(y[0]*z[1]-y[1]*z[0])
    def snf(a):
        d1=reduce(gcd,[x for row in a for x in row])
        minors=[a[i][j]*a[k][l]-a[i][l]*a[k][j] for i,k in itertools.combinations(range(3),2)
                for j,l in itertools.combinations(range(3),2)]
        d2=reduce(gcd,minors)
        return [d1,d2//d1,abs(determinant(a))//d2]
    require({k:determinant(v) for k,v in matrices.items()}==receipt['determinants'],'determinants')
    normalized={}
    for p,word in words.items():
        product=reduce(multiply,[matrices[letter] for letter in word],ident)
        data=receipt['primitive_contractions'][str(p)]
        require(product==tuple(map(tuple,data['matrix'])),'word product')
        require(snf(product)==data['smith'],'Smith data')
        require(all(x%p==0 for row in product for x in row),'entry divisibility')
        u=tuple(tuple(x//p for x in row) for row in product)
        normalized[p]=u
        require(snf(u)==data['normalized_smith'],'normalized Smith data')
        require(gcd(determinant(u),p)==1,'primitive p-adic unit')
        require(determinant(u)==data['normalized_determinant'],'normalized determinant')
    require(determinant(normalized[5])==1,'integer unimodular five normalization')
    for v in ((0,1,1),(34,21,61),(-5,0,15),(25,-10,0),(9,18,27)):
        out=vector(v,tuple(tuple(5*x for x in row) for row in normalized[5]))
        require(reduce(gcd,out)==5*reduce(gcd,v),'exact integer gcd sample')
    # Separate representation from the target's RREF algorithm: retain EVERY vector.
    graphs={}
    for p in (3,5,7,11):
        full=frozenset(itertools.product(range(p),repeat=3))
        queue=deque([(full,0)])
        seen={full}
        shortest=None
        while queue:
            space,depth=queue.popleft()
            if len(space)==1 and shortest is None:
                shortest=depth
            for matrix in matrices.values():
                child=frozenset(tuple(x%p for x in vector(v,matrix)) for v in space)
                if child not in seen:
                    seen.add(child);queue.append((child,depth+1))
        by_dimension=[sum(len(s)==p**d for s in seen) for d in range(4)]
        expected=receipt['complete_subspace_graphs'][str(p)]
        require(len(seen)==expected['states'] and by_dimension==expected['ranks'],'explicit subspace universe')
        want=expected['shortest_letter_word_to_zero']
        require(shortest==(None if want is None else len(want)),'shortest zero word')
        graphs[p]={'states':len(seen),'dimension_histogram':by_dimension,'shortest_zero_length':shortest}
    @lru_cache(None)
    def path_count(rows,subset):
        return sum(all(rows[p[i]]>>p[i+1]&1 for i in range(len(p)-1))
                   for p in itertools.permutations(subset))
    def literal_state(rows):
        n=len(rows);port=n-1
        paths=[p for p in itertools.permutations(range(n))
               if all(rows[p[i]]>>p[i+1]&1 for i in range(n-1))]
        e=sum(p[-1]==port for p in paths)
        c=0
        for mask in range(1<<n):
            if not mask>>port&1:continue
            rem=tuple(i for i in range(n) if not mask>>i&1)
            other=[i for i in range(n-1) if mask>>i&1]
            start_count=sum(all(rows[p[i]]>>p[i+1]&1 for i in range(len(p)-1))
                            for tail in itertools.permutations(other) for p in [(port,)+tail])
            c+=start_count*path_count(rows,rem)
        return (len(paths)-e,e,c)
    twins=[]
    for data in receipt['H55_twins']:
        rows=tuple(data['rows']);v=literal_state(rows)
        require(v==tuple(data['state']),'literal H55 state')
        phase=tuple(x%7 for x in vector(v,normalized[7]))
        require(phase==tuple(data['normalized_phase_mod7']),'normalized H55 phase')
        twins.append({'literal_state':v,'phase_mod7':phase})
    basis=((0,1,1),(1,0,2),(2,1,3))
    require(determinant(basis)==2,'actual prefix basis index2')
    output={'status':'PASS','independent_methods':['explicit finite-field vector-set BFS, no RREF',
            'integer product and 2x2-minor Smith arithmetic','literal permutation-and-subset cover counts for H55 twins'],
            'graphs':graphs,'H55_twins':twins,
            'source_sha256':{str(p.relative_to(HERE.parent)):hashlib.sha256(p.read_bytes()).hexdigest()
                             for p in (source,TARGET/'NOTE.md',TARGET/'CHECK.json')},
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    print(json.dumps(output,indent=2,sort_keys=True))


if __name__=='__main__':
    main()
