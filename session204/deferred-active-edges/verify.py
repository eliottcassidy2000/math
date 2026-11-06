#!/usr/bin/env python3
"""Standalone exact checks; no inherited matrices, bank, or third-party modules."""
from collections import Counter
from itertools import combinations,product,permutations
import json,random
from deferred import patterns,update,complete,read
from literal_forest import forest

def need(ok,what):
    if not ok:raise RuntimeError(what)

def partials(n,width):
    edges=[(a,b) for a,b in combinations(range(n),2) if b>=width]
    for bits in range(1<<len(edges)):
        rows=[0]*n
        for bit,(a,b) in enumerate(edges):
            if bits>>bit&1:rows[a]|=1<<b
            else:rows[b]|=1<<a
        yield rows

def physical(rows,ports,drop,mask,permutation=None):
    n=len(rows);rows=rows[:]+[0];active=list(ports)+[n]
    for v in range(n):
        if v not in ports:rows[v]|=1<<n
    forgotten=active[drop];remaining=[i for i in range(len(active)) if i!=drop]
    for bit,label in enumerate(remaining):
        other=active[label]
        if mask>>bit&1:rows[forgotten]|=1<<other
        else:rows[other]|=1<<forgotten
    if permutation is not None:remaining=[remaining[i] for i in permutation]
    return rows,[active[i] for i in remaining]

def physical_complete(rows,ports,bits):
    rows=rows[:]
    for bit,(a,b) in enumerate(combinations(ports,2)):
        if bits>>bit&1:rows[a]|=1<<b
        else:rows[b]|=1<<a
    return rows

def main():
    dimensions={};columns=checks=closures=permutation_checks=0
    # Regenerate every sparse matrix column, including all port permutations.
    for width,p in product((2,3),(1,2)):
        pats=patterns(width,p);known=set(pats);dimensions[f'W{width}_p{p}']=len(pats)
        for drop,mask,perm in product(range(width+1),range(1<<width),permutations(range(width))):
            for ps in pats:
                col=update({ps:1},width,p,drop,mask,perm)
                need(set(col)<=known,'state closure')
                need(all(type(c) is int and c>=0 for c in col.values()),'nonnegative exact column')
                columns+=1
        # All partial digraphs through order four, every incidence/forget choice.
        for n in range(width,5):
            for rows in partials(n,width):
                ports=list(range(width));state,_=forest(rows,ports,p)
                for drop,mask in product(range(width+1),range(1<<width)):
                    got=update(state,width,p,drop,mask)
                    child,remaining=physical(rows,ports,drop,mask)
                    expected,_=forest(child,remaining,p)
                    need(got==expected,('literal transition',width,p,n,drop,mask))
                    checks+=1
                # Every final orientation; p=2 also checks designated-root readout.
                for bits in range(1<<(width*(width-1)//2)):
                    got=complete(state,width,bits)
                    full=physical_complete(rows,ports,bits)
                    expected,starts=forest(full,ports,p)
                    need(got==expected,('literal completion',width,p,n,bits))
                    need(read(got)[0]==sum(starts),'Hamilton readout')
                    need(read(got)[1]==starts[0],'root-start readout')
                    closures+=1
    rng=random.Random(204)
    for width,p in product((2,3),(1,2)):
        for _ in range(12):
            n=6;rows=[0]*n
            for a,b in combinations(range(n),2):
                if b<width:continue
                if rng.randrange(2):rows[a]|=1<<b
                else:rows[b]|=1<<a
            ports=list(range(width));state,_=forest(rows,ports,p)
            drop=rng.randrange(width+1);mask=rng.randrange(1<<width)
            perm=list(range(width));rng.shuffle(perm)
            got=update(state,width,p,drop,mask,perm)
            child,remaining=physical(rows,ports,drop,mask,perm)
            expected,_=forest(child,remaining,p)
            need(got==expected,'permuted literal transition');permutation_checks+=1
    print(json.dumps(dict(status='PASS',dimensions=dimensions,regenerated_sparse_columns=columns,
                         exhaustive_literal_transitions=checks,literal_completions=closures,
                         extra_permuted_transitions=permutation_checks),indent=2))
if __name__=='__main__':main()
