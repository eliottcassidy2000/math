"""Exact finite controls for the Singer module, triple route rule, and count pulses."""
from collections import Counter
from itertools import combinations, product
from pathlib import Path
from runpy import run_path

def orbit_cases():
    def diffs(S,p):return Counter((a-b)%p for a in S for b in S if a!=b)
    for p,S in [(7,(0,1,3)),(13,(0,1,3,9)),(21,(0,1,4,14,16)),
                (57,(1,6,7,9,19,38,42,49))]:
        assert diffs(S,p)==Counter(range(1,p))
    p=13
    D=(0,1,3,9)
    norm={(0,*S) for S in combinations(range(1,13),3)
          if len(diffs((0,*S),p))==12}
    All={tuple(sorted((a*x+b)%p for x in D))
         for a in range(1,p) for b in range(p)}
    assert len(All)==52 and len(norm)==16
    assert norm=={S for S in All if S[0]==0}
    gap=lambda j:min(j%13,-j%13)
    g=lambda S:tuple(sorted(gap(a-b) for a,b in combinations(S,2)))
    bank={}
    for a,b in combinations(range(1,13),2):
        S=(0,a,b);G=g(S)
        if len(set(G))==3:bank.setdefault(G,[]).append(S)
    assert sorted(bank)==[(1,2,3),(1,3,4),(1,4,5),(1,5,6),
                          (2,3,5),(2,4,6),(2,5,6),(3,4,6)]
    all6=set(range(1,7))
    assert [G for G in sorted(bank) if tuple(sorted(all6-set(G))) in bank]==[(1,3,4),(2,5,6)]
    assert [G for G in bank if not {1,2}&set(G)]==[(3,4,6)]
    assert {tuple(sorted(gap(s*i) for i in (3,4,6))) for s in all6}==set(bank)-{(1,3,4),(2,5,6)}
    assert all(len(all6-{gap(x) for x in (s,2*s,t,2*t,s+t,s-t)})<=2
               for s,t in combinations(range(1,7),2))
    def add(A,B):return {(a+b)%13 for a in A for b in B}
    for a,c in product(bank[1,3,4],repeat=2):
        for b,d in product(bank[2,5,6],repeat=2):
            assert add(add(a,c),add(b,d))==set(range(13))
    print('52 Singer contacts; alternating triple four-cycles all excluded')

def pulses():
    hp=run_path(str(Path(__file__).with_name('thirteen-alphabet-check.py')))['hp']
    H=[]
    for mask in range(64):
        D=[(j if (mask>>(j-1))&1 else -j) for j in range(1,7)]
        H.append(hp([sum(1<<((i+d)%13) for d in D) for i in range(13)]))
    assert sorted(Counter(H).items())==[(3669497,12),(3683797,12),(3703011,4),
                                       (3704857,12),(3707483,12),(3711175,12)]
    a,b=([H[base]]+[H[base^(1<<i)] for i in range(6)] for base in (0,7))
    assert len(set(zip(a,b)))==7
    def v(h,p):
        k=0
        while h%p==0:h//=p;k+=1
        return k
    assert [v(h,13) for h in a]==[1,2,1,2,1,1,1]
    assert [v(h,5) for h in a]==[2,0,0,0,0,0,0]
    print('marked pulses',a,b)

def rank_packet():
    assert len({pow(13,j,125) for j in range(100)})==100
    assert pow(13,100,125)==1 and 25%13!=0
    out=[101+12*h for h in range(101,118)]
    assert (out[0],out[-1])==(1313,1505)
    print('full crossing rank:',out)

if __name__=='__main__':
    orbit_cases()
    rank_packet()
    pulses()
