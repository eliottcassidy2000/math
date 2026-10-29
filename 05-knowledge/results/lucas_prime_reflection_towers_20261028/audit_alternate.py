#!/usr/bin/env python3
"""Same-agent alternate exact paths, no producer imports. Requires SymPy."""
import hashlib
import itertools
import json
import math
import random
from pathlib import Path
from sympy.combinatorics import Permutation, PermutationGroup


def compose(p,q):
    return tuple(p[x] for x in q)


def cycles(p):
    todo=set(range(len(p))); out=[]
    while todo:
        row=[]; x=min(todo)
        while x in todo:
            todo.remove(x); row.append(x); x=p[x]
        out.append(row)
    return out


def jordan_certificate(t,r):
    """Produce a short explicit prime-cycle word; never infer from random success alone."""
    size=len(t); rng=random.Random(20261028+size)
    g=tuple(range(size)); word=[]
    for step in range(1,10001):
        a=rng.randrange(size); word.append(a)
        g=compose(r,tuple((v+a)%size for v in g))
        cs=cycles(g)
        for row in cs:
            q=len(row)
            if not(2<=q<=size-3) or any(q%d==0 for d in range(2,math.isqrt(q)+1)):
                continue
            others=[len(c) for c in cs if c is not row]
            e=math.lcm(*others)
            if e%q:
                h=list(range(size))
                for c in cs:
                    for i,v in enumerate(c):h[v]=c[(i+e)%len(c)]
                lengths=sorted(map(len,cycles(h)))
                if lengths!=[1]*(size-q)+[q]:raise RuntimeError('cycle extraction')
                return dict(word_exponents=word,cycle_lengths=sorted(map(len,cs)),
                            isolating_power=str(e),prime_cycle=q,fixed_points=size-q)
    raise RuntimeError('no bounded Jordan certificate found')


def negation(k):
    return tuple(sum((-(v//3**i)%3)*3**i for i in range(k)) for v in range(3**k))


def automorphisms(k):
    if not k:
        return [tuple([0])]
    child=automorphisms(k-1); s=3**(k-1)
    return [tuple(((b+r)%3)*s+parts[b][j] for b in range(3) for j in range(s))
            for r in range(3) for parts in itertools.product(child,repeat=3)]


def centralizer(k):
    if not k:
        return [tuple([0])]
    t=negation(k-1); s=3**(k-1)
    return [tuple(c)+tuple(s+x for x in a)+tuple(2*s+x for x in compose(t,compose(a,t)))
            for c in centralizer(k-1) for a in automorphisms(k-1)]


def main():
    primary=json.loads(Path(__file__).with_name('receipt.json').read_text())
    hashes=[]
    for k in [1,2,3]:
        c=centralizer(k)
        digest=hashlib.sha256('\n'.join(','.join(map(str,g)) for g in sorted(c)).encode()).hexdigest()
        expected=primary['marking_holonomy'][k-1]
        if digest != expected['holonomy_sha256']:
            raise RuntimeError('independent recursive centralizer differs from loop subgroup')
        hashes.append(dict(depth=k,size=len(c),sha256=digest))
    groups=[]
    for n,p,root,edge in [(5,11,4,False),(7,29,24,False),(8,47,None,True),(11,199,62,False)]:
        # Generate independent subsets directly as tuples, not integer-mask filtering.
        w=[v for v in itertools.product([0,1],repeat=n) if all(not(v[i] and v[(i+1)%n]) for i in range(n))]
        if len(w)!=p: raise RuntimeError('word census')
        labels=({v:sum(v[i]*pow(root,i,p) for i in range(n))%p for v in w}
                if root is not None else {v:i for i,v in enumerate(w)})
        if set(labels.values())!=set(range(p)):raise RuntimeError('alternate label bijection')
        if any(p%d==0 for d in range(2,math.isqrt(p)+1)):raise RuntimeError('prime degree')
        r=[0]*p
        for v in w:
            vprime=tuple(v[((1 if edge else 0)-i)%n] for i in range(n))
            r[labels[v]]=labels[vprime]
        t=[(i+1)%p for i in range(p)]
        if compose(r,r)!=tuple(range(p)):raise RuntimeError('alternate mirror involution')
        odd=sum(len(c)-1 for c in cycles(r))%2
        if odd!=(n==11):raise RuntimeError('alternate reflection parity')
        target=660 if n==5 else math.factorial(p)//(1 if n==11 else 2)
        if p<=29:
            order=int(PermutationGroup(Permutation(t),Permutation(r)).order())
            if order!=target: raise RuntimeError(f'group order failed at n={n}')
            groups.append(dict(n=n,p=p,exact_order=str(order),agrees=True,method='Schreier-Sims'))
        else:
            cert=jordan_certificate(t,r)
            groups.append(dict(n=n,p=p,exact_order=str(target),agrees=True,
                               method='explicit Jordan prime-cycle certificate plus reflection sign',certificate=cert))
    print(json.dumps(dict(scope='same-agent alternate computation, not independent-agent audit',
                          recursive_centralizers=hashes,schreier_sims_groups=groups),indent=2,sort_keys=True))


if __name__=='__main__': main()
