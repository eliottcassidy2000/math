"""Exact replay of the closed-middle endpoint-two repairs; no assertions."""
from collections import Counter
import json
from pathlib import Path
from random import Random
from theory_probe import G,endpoint_two

HERE=Path(__file__).resolve().parent


def reflected_product(p,q):
    rp={u:v for u,v in zip(p,p[::-1])};rq={u:v for u,v in zip(q,q[::-1])}
    g={u:rq[rp[u]] for u in p};remaining=set(p);cycles=[]
    while remaining:
        c=[];u=min(remaining)
        while u in remaining:c.append(u);remaining.remove(u);u=g[u]
        cycles.append(c)
    return cycles


def ledger(n,u,epsilon):
    a,b,x,y,c,z=0,1,2,n-3,n-2,n-1
    W=tuple(range(2,n-2));U=tuple(range(3,n-3))
    if not epsilon:
        return 'epsilon0',[(c,)+W+(b,z,a),(z,c)+W+(b,a),(a,z,b)+W+(c,)]
    if u==x:
        return 'second_port',[(a,b,x)+U+(z,y,c),(a,x)+U+(z,y,b,c),(a,x)+U+(y,b,c,z)]
    return 'epsilon1',[(c,b,z,a)+W,(z,c,b,a)+W,(a,z,b)+W+(c,)]


def replay(a,words):
    original=a;n=len(a)
    for p in words:
        G.need(G.is_path(a,p),'ledger directedness')
        G.need(not G.blue(a,p),'ledger blackness')
        a=G.complement(a,p)
    G.need(a==G.complement(original,tuple(range(n))),'exact labelled target')


def random_blue_middle(r,rng):
    out=[0]*r;seen=set()
    for u in range(r-1):out[u]|=1<<(u+1)
    for u in range(r):
        for v in range(u+2,r):
            if (u,v) in seen:continue
            orb={(u,v),(r-1-v,r-1-u)};seen.update(orb)
            forward=bool(rng.getrandbits(1))
            if (0,r-1) in orb:forward=False
            for x,y in orb:
                if not forward:x,y=y,x
                out[x]|=1<<y
    return tuple(out)


def main():
    permutation_checks=0
    for n in range(6,201):
        P=tuple(range(n));q=(n-2,1,n-1,0)+tuple(range(2,n-2))
        sizes=sorted(map(len,reflected_product(P,q)))
        expected=([n] if n%2 else ([n//2,n//2] if n%4==0 else
                  [(n-2)//4,(n-2)//4,(n+2)//4,(n+2)//4]))
        G.need(sizes==sorted(expected),'first-word reflection product')
        q=(0,n-1)+tuple(range(1,n-1))
        sizes=sorted(map(len,reflected_product(P,q)))
        G.need(sizes==([n] if n%2 else [n//2,n//2]),'final-word reflection product')
        permutation_checks+=2
    results=[]
    for n in range(8,11):
        counts=Counter()
        for m in G.blue_tilings(n-2):
            if not m[-1]&1:continue
            epsilon=bool(m[0]>>(n-4)&1)
            for u in range(2,n-2):
                a=endpoint_two(m,u)
                G.need(G.blue(a,tuple(range(n))),'source blue')
                G.need(a[0].bit_count()==n-3 and a[-1].bit_count()==2,'endpoint layer')
                branch,words=ledger(n,u,epsilon);replay(a,words);counts[branch]+=1
        results.append(dict(n=n,counts=dict(counts)))
    small_second=0
    for m in G.blue_tilings(5):
        if m[-1]&1 and m[0]>>3&1:
            a=endpoint_two(m,2);replay(a,ledger(7,2,True)[1]);small_second+=1
    rng=Random(20261020);large=Counter()
    for n in range(11,32):
        for _ in range(32):
            m=random_blue_middle(n-2,rng);u=rng.randrange(2,n-2)
            epsilon=bool(m[0]>>(n-4)&1)
            a=endpoint_two(m,u);branch,words=ledger(n,u,epsilon)
            G.need(G.blue(a,tuple(range(n))) and a[-1]&1,'large blue closed source')
            replay(a,words);large[branch]+=1
    out=dict(status='FINITE-EXACT bounded controls plus deterministic large samples',
             complete_universe='closed standard-blue middles, total n8..10, every extra incoming port u except first/last middle',
             complete=results,second_port_n7=small_second,
             permutation_identity_checks=permutation_checks,
             large_sample_orders='11..31; 32 per order; seed20261020',large_samples=dict(large))
    (HERE/'theory_verify.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
