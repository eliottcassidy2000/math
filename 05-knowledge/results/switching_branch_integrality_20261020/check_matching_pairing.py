#!/usr/bin/env python3
"""Bounded hostile for extending matching-flip pairing to n divisible by 4."""
from itertools import permutations
from collections import deque
from fractions import Fraction
import json

def setup(n):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    pos={ij:k for k,ij in enumerate(pairs)}
    cuts=[sum(1<<k for k,(i,j) in enumerate(pairs) if v in (i,j)) for v in range(n)]
    inside=[k for k,(i,j) in enumerate(pairs) if i]
    def normalize(x):
        cut=0
        for v in range(1,n):
            if (x>>pos[0,v])&1:
                cut^=cuts[v]
        return x^cut
    def full(q):
        return sum(((q>>j)&1)<<k for j,k in enumerate(inside))
    def compressed(x):
        return sum(((x>>k)&1)<<j for j,k in enumerate(inside))
    def transform(q,p):
        x=full(q); y=0
        for k,(i,j) in enumerate(pairs):
            a,b=p[i],p[j]
            bit=((x>>k)&1) ^ (a>b)
            if bit:
                y|=1<<pos[tuple(sorted((a,b))) ]
        return compressed(normalize(y))
    def flip(q,g):
        matching=sum(1<<pos[i,g[i]] for i in range(n) if i<g[i])
        return compressed(normalize(full(q)^matching))
    return transform,flip

def compose(p,q):
    return tuple(p[q[i]] for i in range(len(p)))

def inv(p):
    q=[0]*len(p)
    for i,j in enumerate(p):q[j]=i
    return tuple(q)

def order(p):
    from math import lcm
    seen=set(); answer=1
    for i in range(len(p)):
        if i in seen:continue
        j=i; k=0
        while j not in seen:
            seen.add(j);k+=1;j=p[j]
        answer=lcm(answer,k)
    return answer

def run(n):
    transform,flip=setup(n)
    perms=list(permutations(range(n)))
    identity=tuple(range(n))
    canonical={}; classes=[]
    for q in range(1<<((n-1)*(n-2)//2)):
        if q in canonical:continue
        orbit={transform(q,p) for p in perms}
        for y in orbit:canonical[y]=q
        aut=[p for p in perms if transform(q,p)==q]
        odd=sum(order(p)%2 for p in aut)
        involutions=[p for p in aut if p!=identity and compose(p,p)==identity]
        flips=[]
        for g in involutions:
            target=flip(q,g)
            witnesses=[p for p in perms if transform(q,p)==target]
            commuting=[p for p in witnesses if compose(p,g)==compose(g,p)]
            flips.append({'g':g,'target':target,'target_type':min(transform(target,p) for p in perms),
                          'isomorphic':bool(witnesses),'commuting_witness':commuting[0] if commuting else None,
                          'commuting_witness_order':order(commuting[0]) if commuting else None})
        classes.append({'q':q,'aut_order':len(aut),'rho':str(Fraction(odd,len(aut))),
                        'number_involutions':len(involutions),'flips':flips})
    return {'n':n,'classes':classes,'N_odd':str(sum((Fraction(x['rho']) for x in classes),Fraction())),
            'N_lev':str(sum((1-Fraction(x['rho']) for x in classes),Fraction()))}

if __name__=='__main__':
    result={'records':[run(4),run(6)]}
    n4,n6=result['records']
    if len(n4['classes'])!=2 or len(n6['classes'])!=6:
        raise AssertionError('class-count control failed')
    if (n4['N_odd'],n4['N_lev'])!=('1','1') or (n6['N_odd'],n6['N_lev'])!=('4','2'):
        raise AssertionError('branch-count control failed')
    if not all(f['isomorphic'] for c in n4['classes'] for f in c['flips']):
        raise AssertionError('n4 hostile failed')
    if any(f['isomorphic'] for c in n6['classes'] for f in c['flips']):
        raise AssertionError('n6 positive control failed')
    if sum(bool(c['flips']) for c in n6['classes'])!=4:
        raise AssertionError('n6 even-Aut carrier failed')
    inherited_n8_weights=[Fraction(1,2)]*19+[Fraction(3,4)]*3+[Fraction(1,4)]+[Fraction(7,8)]*3+[Fraction(3,8)]
    if len(inherited_n8_weights)!=27 or sum(inherited_n8_weights)!=15:
        raise AssertionError('inherited n8 arithmetic control failed')
    result['inherited_n8_arithmetic']={'even_aut_types':27,'N_lev':15,'carrier_cardinality_is_odd':True}
    result['status']='ALL CHECKS PASSED'
    with open('matching_pairing_probe.json','w') as f:json.dump(result,f,indent=2);f.write('\n')
    print(result['status'])
    for r in result['records']:
        print('n',r['n'],'types',len(r['classes']),'N_odd',r['N_odd'],'N_lev',r['N_lev'])
        for c in r['classes']:
            if c['flips']:
                print('  type',c['q'],'Aut',c['aut_order'],'rho',c['rho'],
                      'flip target types',sorted(set(f['target_type'] for f in c['flips'])))
