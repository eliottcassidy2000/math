"""Fresh local-gap, named-normalization and full-star observer-loss controls."""
from itertools import product,permutations,combinations
from pathlib import Path
from math import gcd
import json,hashlib
from check import need,H,full_state,induced

HERE=Path(__file__).resolve().parent


def build(gaps):
    components=[];n=0
    for i,g in enumerate(gaps):
        for _ in range(g):components.append([n]);n+=1
        if i<len(gaps)-1:components.append([n,n+1,n+2]);n+=3
    T=[0]*n
    for i,c in enumerate(components):
        for v in c:
            for d in components[i+1:]:
                for u in d:T[v]|=1<<u
        if len(c)==3:
            for j in range(3):T[c[j]]|=1<<c[(j+1)%3]
    return tuple(T),components


def flip(T,a,b):
    R=list(T);R[a]^=1<<b;R[b]^=1<<a
    return tuple(R)


def sink(T):
    n=len(T)
    return tuple(row|(1<<n) for row in T)+(0,)


def star(T,mask):
    n=len(T)
    return tuple(row|((1<<n) if not(mask>>i&1) else 0) for i,row in enumerate(T))+(mask,)


def normalize(gaps):
    T,components=build(gaps);crossings=0
    while True:
        # Deliberately use the rightmost inversion, unlike the producer.
        slots=[i for i in range(len(components)-1) if len(components[i])==3 and len(components[i+1])==1]
        if not slots:break
        at=slots[-1];tri=components[at];v=components[at+1][0]
        for u in reversed(tri):T=flip(T,u,v)
        components[at],components[at+1]=components[at+1],components[at]
        crossings+=1
    order=[v for c in components for v in c]
    expected,_=build([sum(gaps)]+[0]*(len(gaps)-1))
    need(induced(T,order)==expected,'same named normal form under alternate crossing schedule')
    need(crossings==sum(i*g for i,g in enumerate(gaps)),'independent inversion count')
    return T


def main():
    local,_=build((0,1));local_traces=[]
    for order in permutations((0,1,2)):
        T=local;trace=[H(T)]
        for v in order:T=flip(T,v,3);trace.append(H(T))
        need(trace==[3,5,5,3],'every local three-flip order')
        low=flip(local,2,0);lowtrace=[H(local),H(low)]
        for v in(2,1,0):low=flip(low,v,3);lowtrace.append(H(low))
        low=flip(low,2,0);lowtrace.append(H(low))
        need(low==T and lowtrace==[3,1,1,1,1,3],'same named low route endpoint')
        local_traces.append({'order':order,'upper':trace,'lower':lowtrace})
    normalizations=[]
    for r in range(1,4):
        for gaps in product(range(3),repeat=r+1):
            if sum(gaps)>3:continue
            T=normalize(gaps)
            need(H(T)==3**r,'independent normal-form full count')
            normalizations.append(gaps)
    need(len(normalizations)==56,'same complete finite gap universe')
    old,_=build((0,1,0));new=normalize((0,1,0))
    neighbors=[{'edge':(i,j),'H':H(flip(old,i,j))} for i,j in combinations(range(7),2)]
    need(min(rec['H'] for rec in neighbors if rec['H']>=9)==15,'independent exact floor barrier')
    rooted=[sink(old),sink(new)];contents=[]
    for T in rooted:
        need(full_state(T,7)==[9,9,0,9,288,9,288],'independent seven-state collision')
        values=[H(star(T,mask)) for mask in range(1<<len(T))]
        content=0
        for h in values:content=gcd(content,h)
        contents.append({'C1':content,'all256_extension_values':values})
    need([r['C1'] for r in contents]==[3,1],'complete literal star-image gcds')
    for r in(2,3):
        separated,_=build([0]+[1]*(r-1)+[0]);collapsed=normalize([0]+[1]*(r-1)+[0])
        expected=[3**r,3**r,0,3**r,12**r*2**(r-1),3**r,12**r*2**(r-1)]
        need(full_state(sink(separated),len(separated))==full_state(sink(collapsed),len(collapsed))==expected,
             'larger observer family')
    paths=[HERE.parent/'fibres'/name for name in('GAP_TRANSPORT.md','ROOT_OBSERVER_LOSS.md','gap_transport.py')]
    result={'status':'PASS','scope':'All6 local flip orders;56 alternate-schedule named normalizations;all21 minimal-source neighbors;all512 literal one-birth extensions of the two eight-vertex rooted graphs;seven-state larger controls r2,r3.',
            'local_traces':local_traces,'gap_vectors':normalizations,'neighbors':neighbors,
            'star_contents':contents,'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (HERE/'gap_audit.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print('PASS: independent named gap macros, sharp floor barrier, full-seven-state equality, and complete star gcd3 versus1')


if __name__=='__main__':main()
