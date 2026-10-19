#!/usr/bin/env python3
"""Exact, dependency-free controls for prescribed-edge square-sum absorbers."""
import json
from math import gcd, isqrt

checks = 0

def check(ok, message):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(message)

def square(n):
    return n >= 0 and isqrt(n)**2 == n

def ceilroot(n):
    r = isqrt(n)
    return r + (r*r < n)

def cycle_check(vertices, expected):
    check(len(vertices) == len(set(vertices)), 'cycle repeats a vertex')
    check(set(vertices) == set(expected), 'cycle vertex set mismatch')
    for x,y in zip(vertices,vertices[1:]+vertices[:1]):
        check(x != y and square(x+y), 'invalid cycle edge')

def pair(x,y,A=5,B=5,C=7):
    check(x != y and square(x+y), 'invalid omitted pair')
    check(A*A+B*B == C*C+1, 'coefficient relation')
    u,v = A*A*(x+y)-x, B*B*(x+y)-y
    check(max(x,y) < min(u,v) and u != v, 'host distinctness')
    check(all(square(z) for z in (u+x,x+y,y+v,u+v)), 'insertion squares')
    return u,v

pell=[]
A,C=1,1
for _ in range(12):
    check(C*C-2*A*A == -1, 'Pell identity')
    if A > 1:
        u,v=pair(7,9,A,A,C)
        check(u-v == 2 and u+v == 16*C*C, 'gap-preserving Pell lift')
    pell.append([A,C])
    A,C=3*A+2*C,4*A+3*C

for x in range(1,100001):
    p=ceilroot(2*x)+2
    p += (p%2 == 0)
    q=(p*p-2*x-1)//2
    u,v=p*p-x,q*q-x
    check(2*q == p*p-2*x-1, 'singleton integrality')
    check(x < u < v < 196*x, 'singleton linear budget')
    check(u+x == p*p and v+x == q*q and u+v == (q+1)**2,
          'singleton squares')

windows=[]
for s in (3,5,17,57,59,71,101,201,301):
    S=s*s; h=(S-1)//2; q=14*s
    host=[]; omitted=[]
    for x in range(1,h+1):
        u,v=pair(x,S-x)
        check((u,v) == (25*S-x,24*S+x), 'host formula')
        host.extend((u,v)); omitted.extend((x,S-x))
    check(sorted(host)==list(range(24*S+1,25*S)), 'host interval partition')
    check(sorted(omitted)==list(range(1,S)), 'omitted interval partition')
    check(gcd(q-1,2*q)==1, 'rotation coprimality')
    for b in range(1,h-2*q+2,2*q):
        indices=[]; x=b
        for _ in range(2*q):
            indices.append(x)
            y=b+((x-b+q-1)%(2*q))
            k=-1 if y-x == q-1 else 1
            check(y-x == -q*k-k*k, 'connector offset')
            check((24*S+x)+(25*S-y) == (7*s+k)**2, 'connector square')
            x=y
        check(x==b and set(indices)==set(range(b,b+2*q)), 'rotation cycle')
        base=[]; full=[]; partial=[]; protected=[]
        for x in indices:
            u,v=25*S-x,24*S+x
            base += [u,v]
            full += [u,x,S-x,v]
            if (x-b)%3 == 0:
                partial += [u,x,S-x,v]
            else:
                partial += [u,v]
                protected.append((u,v))
        cycle_check(base, [z for x in range(b,b+2*q) for z in (25*S-x,24*S+x)])
        cycle_check(full, base+[z for x in indices for z in (x,S-x)])
        cycle_check(partial, base+[z for x in indices if (x-b)%3==0 for z in (x,S-x)])
        partial_edges={frozenset((a,c)) for a,c in zip(partial,partial[1:]+partial[:1])}
        check(all(frozenset(e) in partial_edges for e in protected), 'protected host edge lost')
        ks=[]
        for x,y in zip(indices,indices[1:]+indices[:1]):
            ks.append(-1 if y-x==q-1 else 1)
        check(q*sum(ks)+sum(k*k for k in ks)==0, 'closed-cycle energy')
        windows.append({'s':s,'first_index':b,'pairs':2*q,'absorbed_vertices':len(full)})

# Exact fixed-edge descent controls for the bounded-endpoint lift obstruction.
check(49*8+24 == 416, 'glue endpoint bound')
check(482+3362 == 62**2 and 62%7 != 0 and 482 > 416, 'fixed seed obstruction')
u,v=pair(482,3362)
check((u,v)==(95618,92738), 'seed absorber host')
check(u==49*1951+19 and v==49*1893-19, 'centered affine descent')
check(1951+1893==62**2 and min(1951,1893)>416, 'descended host obstruction')

print(json.dumps({'status':'passed','checks':checks,'singleton_range':[1,100000],
                  'pell_pairs':pell,'rotation_windows':windows,
                  'scope':'local identities and prescribed-host conditional absorption; no new unrestricted interval theorem'},
                 sort_keys=True,indent=2))
