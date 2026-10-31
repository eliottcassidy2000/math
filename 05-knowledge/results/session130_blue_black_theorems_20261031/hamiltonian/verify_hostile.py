#!/usr/bin/env python3
"""Independent direct-word verification of the equivariant OCF failure."""
from itertools import permutations, combinations
import json


def need(ok, msg):
    if not ok:
        raise RuntimeError(msg)


def matrix(n, code):
    edges = list(combinations(range(n), 2))
    a = [[0]*n for _ in range(n)]
    for k,(u,v) in enumerate(edges):
        a[u][v] = (code//(2**k)) % 2
        a[v][u] = 1-a[u][v]
    return a


def word_code(a,p):
    return sum(a[p[u]][p[v]]*2**k for k,(u,v) in
               enumerate(combinations(range(len(a)),2)))


def hpaths(a):
    return [p for p in permutations(range(len(a)))
            if all(a[p[k]][p[k+1]] for k in range(len(a)-1))]


def cyc_normal(c):
    return min(c[k:]+c[:k] for k in range(len(c)))


def cycles(a):
    n=len(a)
    return sorted({cyc_normal(c) for k in range(3,n+1,2)
                   for c in permutations(range(n),k)
                   if all(a[c[j]][c[(j+1)%k]] for j in range(k))})


def witness(code):
    a=matrix(5,code)
    ps=list(permutations(range(5)))
    aut=[p for p in ps if word_code(a,p)==code]
    anti=[p for p in ps if all(a[p[u]][p[v]]==a[v][u]
                              for u in range(5) for v in range(5) if u!=v)]
    need(len(aut)==len(anti)==1, 'rigidity/unique anti')
    rho=anti[0]
    pp=hpaths(a)
    cc=cycles(a)
    image=[cc.index(cyc_normal(tuple(rho[v] for v in reversed(c)))) for c in cc]
    fixed=[p for p in pp if tuple(rho[v] for v in p[::-1])==p]
    need(len(pp)==9 and len(cc)==4,'path/cycle count')
    need(all(set(c)&set(d) for c,d in combinations(cc,2)), 'conflict K4')
    need(sorted(len(cc[k]) for k,j in enumerate(image) if k==j)==[3,5],
         'fixed length signature')
    need(sorted(len(cc[k]) for k,j in enumerate(image) if k!=j)==[3,3],
         'swapped length signature')
    return dict(code=code,adjacency=a,rho=rho,cycles=[dict(cycle=c,image=cc[image[k]])
                for k,c in enumerate(cc)],fixed_paths=fixed,all_paths=pp)


def smaller_orders():
    rows=[]
    for n in range(1,5):
        pp=list(permutations(range(n)))
        seen=set()
        sc=[]
        for code in range(2**(n*(n-1)//2)):
            a=matrix(n,code)
            c=min(word_code(a,p) for p in pp)
            if c in seen:
                continue
            seen.add(c)
            aa=matrix(n,c)
            conv=[[aa[v][u] for v in range(n)] for u in range(n)]
            if min(word_code(conv,p) for p in pp)==c:
                sc.append(dict(code=c,H=len(hpaths(aa))))
        need(len({x['H'] for x in sc})==len(sc),'smaller selfconverse H collision')
        rows.append(dict(n=n,labelled_universe=2**(n*(n-1)//2),
                         classes=len(seen),selfconverse=sc))
    return rows


def main():
    cases=[witness(8),witness(10)]
    need([len(c['fixed_paths']) for c in cases]==[1,3],'blue hostile')
    # Explicit length-preserving, reflection-intertwining conflict-graph map.
    cc0=[tuple(x['cycle']) for x in cases[0]['cycles']]
    cc1=[tuple(x['cycle']) for x in cases[1]['cycles']]
    perm=next(p for p in permutations(range(4))
              if all(len(cc0[k])==len(cc1[p[k]]) and
                     p[cc0.index(tuple(cases[0]['cycles'][k]['image']))]==
                     cc1.index(tuple(cases[1]['cycles'][p[k]]['image']))
                     for k in range(4)))
    print(json.dumps(dict(status='PASS',cases=cases,equivariant_cycle_bijection=perm,
                          smaller_order_minimality=smaller_orders()),indent=2,sort_keys=True))


if __name__=='__main__':
    main()
