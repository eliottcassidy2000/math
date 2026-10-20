"""Independent directed-set verification of the complete n>=11 partition.

No imports from companion code. The proof is structural; these are exact
small universes and deterministic controls, not an enumeration-based theorem.
"""
from collections import Counter
from itertools import combinations
from pathlib import Path
import hashlib
import json
import random


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def make_source(n, incoming, rng, mask=None):
    arcs, seen, index = set(), set(), 0
    for i, j in combinations(range(1, n-1), 2):
        if (i,j) in seen:
            continue
        orbit = {(i,j),(n-1-j,n-1-i)}
        seen.update(orbit)
        if j == i+1:
            forward = True
        else:
            forward = bool((mask >> index) & 1) if mask is not None else bool(rng.getrandbits(1))
            index += 1
        arcs.update((u,v) if forward else (v,u) for u,v in orbit)
    outgoing = {n-1-u for u in incoming}
    for v in range(1,n-1):
        arcs.add((v,0) if v in incoming else (0,v))
        arcs.add((n-1,v) if v in outgoing else (v,n-1))
    arcs.add((n-1,0))
    return arcs, index


def path(arcs, q):
    return len(set(q)) == len(q) and all(pair in arcs for pair in zip(q,q[1:]))


def blue(arcs, q):
    r = dict(zip(q,q[::-1]))
    return all((r[v],r[u]) in arcs for u,v in arcs)


def move(arcs, q):
    retained = {frozenset(e) for e in zip(q,q[1:])}
    return {(u,v) if frozenset((u,v)) in retained else (v,u) for u,v in arcs}


def ledger(n, incoming, arcs):
    a,b,x,y,c,z = 0,1,2,n-3,n-2,n-1
    d,w = tuple(range(3,n-2)),tuple(range(2,n-2))
    if (b,c) in arcs:
        if y not in incoming:
            return 'open_release', ((b,)+w+(c,z,a), (z,)+w+(c,b,a), (a,z)+w+(c,b))
        if c not in incoming:
            return 'open_penultimate', ((a,b,c,z,x)+d, (b,c,a,z,x)+d, (b,x)+d+(c,a,z))
        interior = tuple(range(3,n-3))
        return 'open_endpair', ((a,b,c,z,x)+interior+(y,), (b,c,z,x)+interior+(a,y), (b,x)+interior+(a,y,c,z))
    if c in incoming:
        if y not in incoming:
            h = next(u for u in incoming if u != c)
            k = n-1-h
            left,right = tuple(range(1,k)),tuple(range(k,n-2))
            return 'closed_movable_cut', ((a,)+left+(z,)+right+(c,), left+(z,)+right+(a,c), left+right+(a,c,z))
        if (x,c) in arcs:
            q = (a,x)+d+(c,z,b)
            s = (a,b,x,c)+d[::-1]+(z,)
            v = (b,z)+d+(c,x,a)
            return 'closed_endpair_two_block', (q,s[::-1],v)
        interior = tuple(range(3,n-3))
        q1 = (a,x)+interior+(y,c,z,b)
        if (y,x) in arcs:
            return 'closed_endpair_yx', (q1, interior+(a,x,y,c,z,b), interior+(a,b,x,y,c,z))
        return 'closed_endpair_xy', (q1, (z,)+interior+(y,x,b,a,c), (b,z)+interior+(y,x,a,c))
    if (c,x) in arcs:
        return 'closed_epsilon_zero', ((c,)+w+(b,z,a), (z,c)+w+(b,a), (a,z,b)+w+(c,))
    if x not in incoming:
        return 'closed_epsilon_one', ((c,b,z,a)+w, (z,c,b,a)+w, (a,z,b)+w+(c,))
    if 3 not in incoming:
        interior = tuple(range(3,n-3))
        return 'closed_second', ((a,b,x)+interior+(z,y,c), (a,x)+interior+(z,y,b,c), (a,x)+interior+(y,b,c,z))
    interior,e = tuple(range(3,n-4)),n-4
    return 'closed_second_pair', ((a,b,x)+interior+(z,e,y,c), (a,x)+interior+(z,e,y,b,c), (a,x)+interior+(e,y,b,c,z))


def replay(n, incoming, arcs, strict, counts, digest):
    p = tuple(range(n))
    need(len(arcs)==n*(n-1)//2 and path(arcs,p) and blue(arcs,p), 'valid source')
    need(sum(u==0 for u,v in arcs)==n-4, 'a score')
    need(sum(u==n-1 for u,v in arcs)==3, 'z score')
    name,words = ledger(n,incoming,arcs)
    now,flags = arcs,[]
    for q in words:
        need(set(q)==set(p) and path(now,q), ('directed',n,incoming,name,q,sorted(arcs)))
        flags.append(blue(now,q))
        now = move(now,q)
    need(now==move(arcs,p), ('labelled target',n,incoming,name))
    if strict:
        need(not any(flags), ('black theorem',n,incoming,name,flags,sorted(arcs)))
    counts[(n,name,'all_black' if not any(flags) else 'some_blue')] += 1
    digest.update(json.dumps([n,incoming,sorted(arcs),name,flags],separators=(',',':')).encode())


def main():
    rng,counts,digest = random.Random(2026102003),Counter(),hashlib.sha256()
    # Every source in these small universes; low-order templates need at least9.
    small = {}
    for n in [9]:
        _,bits = make_source(n,(2,3),rng)
        number = 0
        for mask in range(1<<bits):
            for incoming in combinations(range(2,n-1),2):
                arcs,_ = make_source(n,incoming,rng,mask)
                replay(n,incoming,arcs,False,counts,digest)
                number += 1
        small[n] = number
    larger = {}
    for n in list(range(11,23))+[25,30,41,60]:
        choices = list(combinations(range(2,n-1),2)) if n<=22 else sorted(set([(h,n-2) for h in range(2,n-2)]+[(h,n-3) for h in range(2,n-3)]+[(2,3),(2,4),(4,5)]))
        number = 0
        for incoming in choices:
            for _ in range(8):
                arcs,_ = make_source(n,incoming,rng)
                replay(n,incoming,arcs,True,counts,digest)
                number += 1
        larger[n] = number
    # R1 sharpness is an infinite construction, sampled for k=2..16.
    sharp = []
    for k in range(2,17):
        m,n = 2*k-1,4*k-2
        arcs = set()
        for i,j in combinations(range(m),2):
            u,v = (i,j) if 1<=(j-i)%m<=k-1 else (j,i)
            arcs.update({(2*u,2*v),(2*u+1,2*v+1)})
        for i in range(m):
            for j in range(m):
                arcs.add((2*j+1,2*i) if j==(i-1)%m else (2*i,2*j+1))
        p,q = tuple(range(n)),tuple(range(1,n))+(0,)
        need(path(arcs,p) and blue(arcs,p) and path(arcs,q) and blue(arcs,q), 'sharp cyclic cut')
        need(sum(u==0 for u,v in arcs)==n-1-k and sum(u==n-1 for u,v in arcs)==k,'sharp scores')
        sharp.append([k,n])
    receipt = {'status':'PASS','universe':'standard marked closed blue sources; no strong-partner, mixedness, or isomorphism filters','complete_small':small,'deterministic_larger':larger,'sharp_cyclic_cut_controls':sharp,'replays':sum(counts.values()),'sha256':digest.hexdigest(),'branch_counts':[{'n':n,'branch':b,'outcome':o,'count':c} for (n,b,o),c in sorted(counts.items())]}
    Path(__file__).with_suffix('.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k!='branch_counts'},indent=2))


if __name__=='__main__':
    main()
