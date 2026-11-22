"""Literal perfect-center pairs and their small birth incidences."""
import sys
from collections import Counter
from itertools import permutations, product
from pathlib import Path

from check_core_grammar import ham, pairs, permutation_action, transform, cuts

sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'chanpak_gap_star_ports_20261122'))
from replay import profile, count_from_profile, star_rows, strong, flipped, p

EXACT = {
 3: ((30,52,56,50,32,1),(1,5)),
 5: ((160,145,67,7,13,30,59,124),(0,2)),
 7: ((988,433,74,34,236,645,938,780,60,286),(8,9)),
}


def center(j): return (1<<(j-1))*((1<<j)-1)


def tobits(A):
    return sum(((A[v]>>u)&1)<<j for j,(u,v) in enumerate(pairs(len(A))))


def exact():
    for j,(A,arc) in EXACT.items():
        B=flipped(A,*arc)
        assert strong(A) and strong(B)
        assert [ham(len(A),tobits(D)) for D in (A,B)]==[center(j)-1,center(j)+1]
        print("p =",j,"literal order",len(A),center(j)-1,center(j)+1)
    A=(0,33,27,211,163,141,55,71); tag=55; edge=(3,7); y=23955
    Bs=(A,flipped(A,*edge))
    P=[profile(B,tag) for B in Bs]
    assert P==[(111,2856,(633,744),(184,245),309),
               (111,2856,(633,744),(186,246),311)]
    assert p(y)==5080
    hosts=[star_rows(B,tag,0,y) for B in Bs]
    assert all(len(B)==24 and strong(B) for B in hosts)
    assert flipped(hosts[0],*edge)==hosts[1]
    assert [count_from_profile(row,0,y) for row in P]==[center(13)-1,center(13)+1]
    print("p = 13 exact profile, literal order 24:",center(13)-1,center(13)+1)


def census(path):
    seen=Counter(); maxH=Counter(); byH=Counter(); joint=Counter(); lowers=Counter()
    for line in open(path):
        ni,hi,rows=line.split()
        n,h=int(ni),int(hi)
        seen[n]+=1
        maxH[n]=max(maxH[n],h)
        byH[n,h]+=1
        if (n,h) not in {(6,27),(8,495)}:
            continue
        A=tuple(map(int,rows.split(','))); x=tobits(A)
        trace=[k for k in range(len(pairs(n))) if ham(n,x^(1<<k))==h+2]
        if not trace:continue
        assert not tuple(cuts(n,x))
        assert all(not tuple(cuts(n,x^(1<<k))) for k in trace)
        scores=[a.bit_count() for a in A]
        cells=[[i for i in range(n) if scores[i]==v] for v in sorted(set(scores))]
        aut=0
        for labels in product(*(permutations(c) for c in cells)):
            q=list(range(n))
            for c,js in zip(cells,labels):
                for i,j in zip(c,js):q[i]=j
            aut+=transform(x,permutation_action(q))==x
        assert len(trace)%aut==0
        joint[n,h]+=len(trace)//aut
        lowers[n,h]+=1
    assert [seen[i] for i in range(1,10)]==[1,1,2,4,12,56,456,6880,191536]
    assert [maxH[i] for i in [5,7,9]]==[15,189,3357]
    for n,h,lo,hi,s,a in [(6,27,2,4,2,2),(8,495,22,33,18,14)]:
        assert (byH[n,h],byH[n,h+2],joint[n,h],lowers[n,h])==(lo,hi,s,a)
        print("joint census",n,h,lo,hi,s,a)


if __name__=="__main__":
    exact()
    if len(sys.argv)>1:census(sys.argv[1])
