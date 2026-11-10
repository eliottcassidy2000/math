#!/usr/bin/env python3
"""Self-contained exact regression for the all-positive-index fan chart."""
from collections import defaultdict
from itertools import combinations, permutations
from math import gcd, lcm
import json
from pathlib import Path


def require(ok, detail):
    if not ok:
        raise RuntimeError(detail)


def fan(n):
    if n == 1:
        return [2, 4, 1]
    k = 2*n
    size = k+2
    losses = {n, n+1, 2*n+1}
    rows = [0]*size
    for i in range(size):
        for j in range(i+1, size):
            u, v = (j, i) if i == 0 and j in losses else (i, j)
            rows[u] |= 1 << v
    return rows


def paths(rows):
    n = len(rows)
    layers = {(1 << v, v): 1 for v in range(n)}
    for _ in range(1, n):
        new = defaultdict(int)
        for (mask, end), count in layers.items():
            choices = rows[end] & ~mask
            while choices:
                bit = choices & -choices
                choices -= bit
                new[mask | bit, bit.bit_length()-1] += count
        layers = new
    return sum(layers.values())


def h(n):
    return 4**n-1


def phi(n):
    result, p = n, 2
    while p*p <= n:
        if n % p == 0:
            result -= result//p
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        result -= result//n
    return result


def embedding(small, large):
    k = len(small)
    for vertices in combinations(range(len(large)), k):
        induced_degrees = {v: sum((large[v] >> w)&1 for w in vertices) for v in vertices}
        target_degrees = [r.bit_count() for r in small]
        groups = [[v for v in vertices if induced_degrees[v] == d] for d in target_degrees]
        if any(not g for g in groups):
            continue
        for perm in permutations(vertices):
            if any(perm[i] not in groups[i] for i in range(k)):
                continue
            if all(((small[i]>>j)&1) == ((large[perm[i]]>>perm[j])&1) for i in range(k) for j in range(k)):
                return list(perm)
    return None


def main():
    actual = []
    for n in range(1,7):
        rows=fan(n)
        hp=paths(rows)
        require(hp==h(n), ('actual fan',n,hp))
        actual.append({'n':n,'order':len(rows),'rows':rows,'H':hp})
    for a in range(1,201):
        for b in range(1,201):
            require(gcd(h(a),h(b))==h(gcd(a,b)), ('gcd',a,b))
            require((h(b)%h(a)==0)==(b%a==0), ('divisibility',a,b))
            require(gcd(h(a)//3,h(b)//3)==h(gcd(a,b))//3, ('normalized gcd',a,b))
            require(h(a+b)==h(a)*h(b)+h(a)+h(b), ('addition',a,b))
    totient=[]
    for n in range(1,201):
        chart_count=sum(gcd(h(a),h(n))==3 for a in range(1,n+1))
        require(chart_count==phi(n), ('index totient',n))
        divisors=[a for a in range(1,n+1) if h(n)%h(a)==0]
        require(sum(phi(a) for a in divisors)==n, ('divisor totient',n))
        if n<=16:
            image_phi=phi(h(n))
            require((image_phi & (image_phi-1)==0)==(n in {1,2,4,8,16}), ('image totient',n))
            totient.append({'n':n,'H':h(n),'phi_index':phi(n),'phi_H':image_phi})
    primes=[p for p in range(3,102) if all(p%d for d in range(2,int(p**.5)+1))]
    order_checks=0
    for p in primes:
        t=1
        while pow(4,t,p)!=1:
            t+=1
        for n in range(1,201):
            require((h(n)%p==0)==(n%t==0), ('prime visibility',p,n,t))
            order_checks+=1
    witness=embedding(fan(2),fan(3))
    require(witness is not None, 'expected induced-containment hostile')
    out={'status':'PASS','universe':{'adjacency_fans':'n=1..6','arithmetic_pairs':40000,'totient_indices':'n=1..200','image_totients':'n=1..16','prime_order_checks':order_checks},'actual_fans':actual,'totient_rows':totient,'hostiles':{'induced_Q2_in_Q3':witness,'index_divisibility':False,'H_divisibility':False,'ordinary_lcm_H2_H3':lcm(h(2),h(3)),'chart_join_H6':h(6),'H2_times_H3':h(2)*h(3),'normalized_R3':h(3)//3,'R3_is_forbidden_H':True},'limits':['The normalized chart scalar R_n is not asserted to be realizable as H.','Finite tests supplement all-index proofs; no unrestricted Chan-Pak interval is claimed.']}
    text=json.dumps(out,indent=2,sort_keys=True)+'\n'
    (Path(__file__).parent/'exact.json').write_text(text)
    print(text,end='')

if __name__=='__main__':
    main()
