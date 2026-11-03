#!/usr/bin/env python3
"""Independent exact wrapper and repeated-operation controls."""
from extended import *

def legalwords(queries,runs):
    names=['Y']*runs+list('ab'+'cd'[:queries]);allowed={('Y','a'),('a','b'),('b','Y')}
    old=['Y','a','b']
    for vertex,target in [('c','a'),('d','c')][:queries]:
        for v in old:allowed.add((vertex,v) if v==target else (v,vertex))
        old.append(vertex)
    return sorted(''.join(p) for p in set(permutations(names)) if all((u,v) in allowed for u,v in zip(p,p[1:])))

def main():
    expected=[(['Yab','abY','bYa'],['YabY']),(['Yabc','Ycab','abYc','bYca','cabY'],['YabYc','YcabY']),(['Yabdc','Ycabd','Ydcab','abYdc','bYadc','bYcad','bYdca','cabYd','dcabY'],['YabYdc','YcabYd','YdcabY'])]
    for q in range(3):
        for runs in [1,2]:need(legalwords(q,runs)==sorted(expected[q][runs-1]),('run words',q,runs,legalwords(q,runs)))
    print('All six collapsed-run word lists independently verified',flush=True)
    count=0
    for n in range(1,5):
        for code in range(1<<comb(n,2)):
            Y=tour(n,code);P=perm_h(Y);R=2*P+sum(map(sum,direct_covers(Y)))
            need(obs(wrap(Y),n)==(P+R,R//2,P+R//2),('universal wrapper finite',n,code));count+=1
    print('Exact wrapper observer: all',count,'labeled nonempty Y through order four PASS',flush=True)
    gs=gates();seed=tour(2,1);source=wrap(seed)
    need(obs(source,2)==(5,2,3),'initial source')
    for t0,t1 in [(0,0),(0,1),(1,0),(1,1)]:
        Y=ordinal([seed,gs[0,t0],((False,),),gs[0,t1],((False,),)])
        z=wrap(Y);k=t0+t1;r=2-k;P=9**r*11**k;R=4*168**r*200**k
        v=obs(z,len(Y));need(v==(P+R,R//2,P+R//2),('two repeated ops',t0,t1,v))
        need(v3(gcd(*v))==r,'closed valuation');need(len(sccs(z))==1,'closed strong')
        need(sub(z,[0,1,len(Y),len(Y)+1])==source,'induced initial source')
        print('Two literal closed operations:',(t0,t1),'observer',v,'PASS',flush=True)
    print('PASS: exact wrapper compression and closed operations',flush=True)

if __name__=='__main__':main()
