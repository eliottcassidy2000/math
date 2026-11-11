#!/usr/bin/env python3
"""New local-image theorem and exact short-word automaton controls.
No inherited census is rerun. Python standard library only.
"""
from collections import deque
from functools import lru_cache
from math import gcd
import json

@lru_cache(None)
def f(n):
    if n < 2: return 1
    q,b = divmod(n,2)
    return f(q) + (q+q%2 if b==0 else 0)

def image(h,e,r,L):
    return h*f(L-1)+r*L-e*(L%2)

def predicted(h,e,r,m):
    # Work modulo 2m, then project. The lift retains physical parity.
    M=2*m
    d=gcd(gcd(2*h,2*r),M)
    return {a%m for a in range(M) if (a-h)%d==0 or (a-h-r+e)%d==0}

def automaton(M):
    assert M%2==0
    start=(0,1%M)
    prev={start:None}; dist={start:0}; todo=deque([start])
    while todo:
        n,z=todo.popleft()
        for b in (0,1):
            nxt=((2*n+b)%M,(z+(n+n%2 if b==0 else 0))%M)
            if nxt not in prev:
                prev[nxt]=((n,z),b)
                dist[nxt]=dist[(n,z)]+1
                todo.append(nxt)
    assert len(prev)==M*M//2
    return prev,dist

def word(prev,state):
    bits=[]
    while prev[state] is not None:
        state,b=prev[state];bits.append(str(b))
    return ''.join(reversed(bits))

def run():
    roots=[(1,1,1),(3,1,3),(5,1,5),(9,1,15)]
    stats=[]; checks=0
    for M in range(2,130,2):
        prev,dist=automaton(M)
        stats.append({'M':M,'states':len(prev),'max_shortest_bits':max(dist.values())})
        for h,e,r in roots:
            actual={(h*z+r*((n+1)%M)-e*((n+1)%2))%M for n,z in prev}
            assert actual==predicted(h,e,r,M)
            checks+=1
        # Reconstruct one farthest state at every modulus and check integer recurrence.
        state=max(dist,key=dist.get); bits=word(prev,state); n=int(bits or '0',2)
        assert (n%M,f(n)%M)==state
    # Primitive C3 prefix: exact local gate 0 or 2 modulo 3.
    assert predicted(3,1,3,3)=={0,2}
    # Singleton prefix: no finite-modulus odd obstruction, but an exact gap at 13.
    assert all(image(1,1,1,L)!=13 for L in range(1,14))
    # Larger L cannot hit 13: G(L)>=1 and L-parity(L)>=14 for L>=14.
    gaps=[N for N in range(1,102,2) if not any(image(1,1,1,L)==N for L in range(1,N+2))]
    print(json.dumps({'local_image_checks':checks,'diameters':stats,'singleton_gaps_through_101':gaps,
      'scope':'FINITE-EXACT shortest-word diameters through modulus 128; local-image identity proved in README'},indent=2))
if __name__=='__main__':run()
