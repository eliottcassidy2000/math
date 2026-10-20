#!/usr/bin/env python3
"""Exhaust all circulant presentations at odd orders 3..15.
Find a non-automorphic Hamiltonian-cycle successor and a bad cut reversal.
No isomorphism filter or randomized sampling. DFS stops at the first certificate.
Control counts use independent subset dynamic programming.
"""
from functools import lru_cache
from itertools import permutations
from math import factorial
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parent

def circulant(n, mask):
    S={d if mask>>(d-1)&1 else n-d for d in range(1,(n+1)//2)}
    return S, [[int(i != j and (j-i)%n in S) for j in range(n)] for i in range(n)]

def cycle_count(A):
    n=len(A)
    @lru_cache(None)
    def f(mask,last):
        if mask==(1<<n)-1: return A[last][0]
        return sum(f(mask|1<<v,v) for v in range(1,n) if not(mask>>v&1) and A[last][v])
    return f(1,0)

def bad_pair(A,p,anti=False):
    n=len(A)
    return next(((i,j) for i in range(n) for j in range(i+1,n)
                 if A[i][j] != (A[p[j]][p[i]] if anti else A[p[i]][p[j]])),None)

def first_bad_cycle(A):
    n=len(A); c=[0]; visited=0
    def dfs(mask):
        nonlocal visited
        if len(c)==n:
            if not A[c[-1]][0]: return None
            visited += 1
            succ=[0]*n
            for j in range(n): succ[c[j]]=c[(j+1)%n]
            pair=bad_pair(A,succ)
            if pair is not None:
                for cut in range(n):
                    p=[c[(cut+j)%n] for j in range(n)]
                    rev=[0]*n
                    for j in range(n): rev[p[j]]=p[-1-j]
                    rp=bad_pair(A,rev,anti=True)
                    if rp is not None:
                        return dict(cycle=c.copy(), successor_bad_pair=pair,
                                    bad_cut_path=p, reversal_bad_pair=rp,
                                    cycles_checked=visited)
                raise AssertionError('odd cycle proof failure')
            return None
        for v in range(1,n):
            if not mask>>v&1 and A[c[-1]][v]:
                c.append(v); ans=dfs(mask|1<<v); c.pop()
                if ans is not None:return ans
        return None
    ans=dfs(1)
    return ans,visited

def verify_certificate(n,S,w):
    # Separate presentation-level verifier, not using the matrix helper.
    arc=lambda a,b:(b-a)%n in S
    C=w['cycle']; P=w['bad_cut_path']
    assert sorted(C)==list(range(n))==sorted(P)
    assert all(arc(C[i],C[(i+1)%n]) for i in range(n))
    assert all(arc(P[i],P[i+1]) for i in range(n-1))
    succ=dict(zip(C,C[1:]+C[:1])); rev=dict(zip(P,P[::-1]))
    u,v=w['successor_bad_pair']; assert arc(u,v)!=arc(succ[u],succ[v])
    u,v=w['reversal_bad_pair']; assert arc(u,v)!=arc(rev[v],rev[u])

def aut_count(A):
    # Tiny controls only; direct complete permutation test.
    return sum(bad_pair(A,p) is None for p in permutations(range(len(A))))

def main():
    rows=[]; summaries=[]
    for n in range(3,16,2):
        bad=0; max_checked=0; control_hcs=[]
        for mask in range(1<<((n-1)//2)):
            S,A=circulant(n,mask); w,checked=first_bad_cycle(A)
            max_checked=max(max_checked,checked)
            if n>=7:
                assert w is not None
                verify_certificate(n,S,w);bad+=1
            else:
                assert w is None
                assert checked==cycle_count(A)==(1 if n==3 else 2)
            rows.append(dict(n=n,mask=mask,S=sorted(S),certificate=w))
        summaries.append(dict(n=n,presentations=1<<((n-1)//2),excluded=bad,
                              max_cycles_examined=max_checked))
    carousel=[]
    for n in range(3,16,2):
        _,A=circulant(n,(1<<((n-1)//2))-1)
        carousel.append(dict(n=n,HC=cycle_count(A)))
    assert [r['HC'] for r in carousel[:6]]==[1,2,17,222,5109,166562]
    # C3 lexicographic C3: blocks 0,1,2 in cyclic order; each block cyclic.
    A=[[int(i!=j and (((j//3-i//3)%3==1) if i//3!=j//3 else (j-i)%3==1))
        for j in range(9)] for i in range(9)]
    lex=dict(HC=cycle_count(A),Aut=aut_count(A))
    assert lex==dict(HC=207,Aut=81)
    # Both named order-five tournaments checked without using their aliases.
    _,reg=circulant(5,3)
    E={(0,1),(1,2),(2,3),(3,0),(0,4),(2,4),(4,1),(4,3),(0,2),(1,3)}
    exception=[[int((i,j) in E) for j in range(5)]for i in range(5)]
    fives=[]
    for label,A in [('regular_C5',reg),('Darbinyan_T5_exception',exception)]:
        hps=[p for p in permutations(range(5)) if all(A[p[i]][p[i+1]] for i in range(4))]
        fives.append(dict(label=label,scores=list(map(sum,A)),HP=len(hps),
                          Aut=aut_count(A),HC=cycle_count(A),
                          bypass_paths=[list(p) for p in hps if A[p[0]][p[-1]]]))
    assert fives[0]['HP']==15 and fives[0]['Aut']==5
    assert fives[1]['scores']==[3,2,2,1,2] and not fives[1]['bypass_paths']
    # Exact bound comparison and induction base.
    base_L=4**22*factorial(22);base_U=3**71
    assert base_L>base_U
    comparison=[dict(n=n,k=(n-7)//6,lower=4**((n-7)//6)*factorial((n-7)//6),
                     upper=3**((n-1)//2)-1) for n in range(7,144,2)]
    exceptional_exclusions=[r['n'] for r in comparison if r['n']<139 and r['lower']>r['upper']]
    assert exceptional_exclusions==[127,133,135]
    assert 2*4**21*factorial(21)>3**68
    assert [n for n in range(7,133,2) if 2*4**((n-7)//6)*factorial((n-7)//6)>3**((n-1)//2)-1]==[121,127,129]
    for m in range(3,101):
        n=2*m+1; S=set(range(1,m+1))
        C=list(range(m-1))+list(range(m,2*m-1))+[2*m,m-1,2*m-1]
        assert sorted(C)==list(range(n))
        assert all((C[(i+1)%n]-C[i])%n in S for i in range(n))
        succ=dict(zip(C,C[1:]+C[:1]));assert (m-1) in S and (succ[m-1]-succ[0])%n not in S
    payload=dict(universe=summaries,carousel_HC=carousel,lex_C3_C3=lex,order_five=fives,
                 bound_base=dict(k=22,lower=base_L,upper_at_n143=base_U),
                 below139_counting_exclusions=exceptional_exclusions,certificates=rows)
    (ROOT/'cycle_certificates.json').write_text(json.dumps(payload,indent=2)+'\n')
    for summary in summaries:print(summary)
    print('carousel_HC',carousel); print('lex_C3_C3',lex);print('order_five',fives)
    print('k22 exact bound',base_L,'>',base_U)
    print('extra counting exclusions',exceptional_exclusions)
    print('ALL CHECKS PASSED')

if __name__=='__main__':main()
