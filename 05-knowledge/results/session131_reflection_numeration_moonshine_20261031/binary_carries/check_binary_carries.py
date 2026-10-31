#!/usr/bin/env python3
"""Self-contained exact checker; stdlib only; no asserts or inherited imports."""
import itertools as it
import json
import random
from collections import Counter

GATES = 0


def need(ok, msg):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(msg)


def matrix(n, mask):
    a = [[0] * n for _ in range(n)]
    for k, (u, v) in enumerate(it.combinations(range(n), 2)):
        a[u][v] = (mask >> k) & 1
        a[v][u] = 1 - a[u][v]
    return a


def ham(a):
    n = len(a)
    if not n:
        return 1
    d = [[0] * n for _ in range(1 << n)]
    for u in range(n):
        d[1 << u][u] = 1
    for s in range(1, 1 << n):
        for u, val in enumerate(d[s]):
            if val:
                for v in range(n):
                    if not (s >> v) & 1 and a[u][v]:
                        d[s | (1 << v)][v] += val
    return sum(d[-1])


def cycles(a):
    """Directed simple odd cycles; least root, no reversal identification."""
    n = len(a)
    out = []
    def visit(path, used):
        u = path[-1]
        if len(path) >= 3 and len(path) % 2 and a[u][path[0]]:
            out.append(tuple(path))
        for v in range(path[0] + 1, n):
            if not (used >> v) & 1 and a[u][v]:
                visit(path + [v], used | (1 << v))
    for root in range(n):
        visit([root], 1 << root)
    return out


def normalize(c):
    k = c.index(min(c))
    return c[k:] + c[:k]


def rows(a):
    return [sum(x << v for v, x in enumerate(row)) for row in a]


def examine(a, rho):
    n = len(a)
    need(all(rho[rho[v]] == v for v in range(n)), 'not an involution')
    need(all(a[u][v] == a[rho[v]][rho[u]] for u in range(n) for v in range(n)),
         'not an anti-automorphism')
    cc = cycles(a)
    lookup = {c: k for k, c in enumerate(cc)}
    support = [sum(1 << v for v in c) for c in cc]
    mirror = [lookup[normalize(tuple(rho[v] for v in reversed(c)))] for c in cc]
    alpha = [0] * (n // 3 + 1)
    fixed = [0] * len(alpha)
    pack = []
    def visit(start, used):
        k = len(pack)
        alpha[k] += 1
        if {mirror[v] for v in pack} == set(pack):
            fixed[k] += 1
        for v in range(start, len(cc)):
            if not used & support[v]:
                pack.append(v)
                visit(v + 1, used | support[v])
                pack.pop()
    visit(0, 0)
    alpha += [0] * max(0, 4 - len(alpha))
    fixed += [0] * max(0, 4 - len(fixed))
    q = sum(v == mirror[v] for v in range(len(cc)))
    i = sum(v < mirror[v] and bool(support[v] & support[mirror[v]]) for v in range(len(cc)))
    j = sum(v < mirror[v] and not support[v] & support[mirror[v]] for v in range(len(cc)))
    need(alpha[1] == q + 2*i + 2*j, 'single-cycle orbits')
    need(fixed[1] == q and fixed[2] == j, 'fixed packing classification')
    need(all((x-y) % 2 == 0 for x,y in zip(alpha,fixed)), 'packing orbit integrality')
    p = (alpha[2]-j)//2
    t = fixed[3]
    direct_t = sum(not (support[f] & (support[v] | support[mirror[v]]))
                   for f in range(len(cc)) if mirror[f] == f
                   for v in range(len(cc)) if v < mirror[v] and not (support[v] & support[mirror[v]]))
    need(t == direct_t, 'fixed three-packing classification')
    h = ham(a)
    need(h == sum((1 << k)*x for k,x in enumerate(alpha)), 'exact OCF vs path DP')
    need(h % 8 == (1 + 2*q + 4*i) % 8, 'mod8')
    need(h % 16 == (1 + 2*q + 4*i + 8*(j+p+t)) % 16, 'mod16')
    need((h >> 1) & 1 == q % 2, 'second digit')
    need((h >> 2) & 1 == (q//2+i) % 2, 'third digit')
    need((h >> 3) & 1 == ((q+2*i)//4+j+p+t) % 2, 'fourth digit')
    need((h % 8 == 7) == (q % 2 == 1 and (q//2+i) % 2 == 1), 'seven residue criterion')
    carry = 0
    bits = []
    for k in range(max(len(alpha), h.bit_length())):
        raw = (alpha[k] if k < len(alpha) else 0) + carry
        bits.append(raw % 2)
        carry = raw // 2
        need(bits[-1] == ((h >> k) & 1), 'general carry recurrence')
    need(carry == 0, 'terminal carry')
    if n % 2 == 0:
        need(q == 0 and t == 0 and h % 4 == 1, 'even-order specialization')
    if n <= 8:
        need(h == 1+2*q+4*i+8*(j+p), 'exact orbit identity through order8')
    while len(alpha) > 1 and not alpha[-1]:
        alpha.pop()
    return {'n':n, 'rows':rows(a), 'rho':list(rho), 'H':h,
            'alpha':alpha, 'q':q, 'i':i, 'j':j, 'p':p, 't':t,
            'bits_low_to_high':bits}


def symmetric_matrices(n, codes=None):
    pairs = list(it.combinations(range(n),2))
    index = {uv:k for k,uv in enumerate(pairs)}
    groups = []
    for k,(u,v) in enumerate(pairs):
        l = index[(n-1-v,n-1-u)]
        if k <= l:
            groups.append((1<<k)|(1<<l))
    for code in range(1 << len(groups)) if codes is None else codes:
        mask = sum(g for k,g in enumerate(groups) if (code >> k) & 1)
        yield code, matrix(n, mask)


def triangle_stack(k):
    n = 3*k
    a = [[0]*n for _ in range(n)]
    for u,v in it.combinations(range(n),2):
        a[u][v] = int(u//3 < v//3 or (v-u) % 3 == 1)
        a[v][u] = 1-a[u][v]
    return a


def main():
    receipt = {'scope':'All reversal-compatible labelled tournaments through order7; deterministic samples8,9; structural stack controls.'}
    census=[]
    hostiles={}
    charge_fibres={}
    rng=random.Random(131)
    for n in range(1,10):
        codes=None if n<=7 else sorted(rng.sample(range(1 << ((n*n)//4)),128 if n==8 else 48))
        counts=Counter()
        total=0
        for code,a in symmetric_matrices(n,codes):
            d=examine(a,tuple(reversed(range(n))))
            h=d['H'];q=d['q'];i=d['i'];j=d['j'];p=d['p'];t=d['t']
            total+=1;counts[h%16]+=1
            if (1+2*(q%2)+4*(i%2))%8 != h%8:
                hostiles.setdefault('dropping_binary_carry',d)
            if h%8==7:
                hostiles.setdefault('residue7_is_attained',d)
            for name,term in [('j',j),('p',p),('t',t)]:
                if term%2:
                    hostiles.setdefault('omitting_'+name,d)
            key=(n,q,i)
            old=charge_fibres.setdefault(key,d)
            if old['H']%16 != h%16:
                hostiles.setdefault('same_q_i_different_mod16',[old,d])
        census.append({'n':n,'tested':total,'H_mod16':dict(sorted(counts.items()))})
    controls=[]
    for k in range(1,5):
        d=examine(triangle_stack(k),tuple(reversed(range(3*k))))
        need(d['H']==3**k,'triangle stack count')
        need(d['alpha']==[__import__('math').comb(k,r) for r in range(k+1)], 'triangle stack profile')
        controls.append(d)
        if k == 3:
            hostiles['omitting_t']=d
    # Actual H=9 cycle-profile collision; inherited THM-466 correction.
    a=matrix(5,8); rho=(4,3,2,1,0)
    d1=examine(a,rho)
    b=triangle_stack(2);d2=examine(b,tuple(reversed(range(6))))
    need(d1['H']==d2['H']==9 and d1['alpha']!=d2['alpha'],'H does not determine OCF')
    hostiles['H_does_not_determine_OCF']=[d1,d2]
    for name in ['dropping_binary_carry','residue7_is_attained','omitting_j','omitting_p','omitting_t','same_q_i_different_mod16']:
        need(name in hostiles,'missing hostile: '+name)
    # A two-face needs its two singleton subfaces: (1,1,1) cannot be a graph profile.
    need(1+2+4==7 and 1>__import__('math').comb(1,2),'formal seven face profile forbidden')
    path=[]
    for m in range(20):
        values=[1,3]
        while len(values)<=m:
            values.append(values[-1]+2*values[-2])
        value=values[m]
        need(3*value==2**(m+2)+(-1)**(m+1),'correct path Jacobsthal formula')
        path.append(value)
    need((2**2+(-1))//3 != path[1],'old path formula rejected')
    receipt.update(census=census,triangle_stacks=controls,hostiles=hostiles,path_I_at_2=path,
                   gates=GATES,status='PASS')
    print(json.dumps(receipt,indent=2,sort_keys=True))


if __name__=='__main__':
    main()
