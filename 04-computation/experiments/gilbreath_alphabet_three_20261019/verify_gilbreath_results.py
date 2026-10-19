#!/usr/bin/env python3
"""Exact certificates for the alphabet-three shield/first-zero-run research.
No assert statements; all checks remain live under python -O.
"""
from collections import deque
from fractions import Fraction
from itertools import product
import json, hashlib
from pathlib import Path

def require(v,msg):
    if not v:raise RuntimeError(msg)
def D(w):return tuple(abs(x-y)for x,y in zip(w,w[1:]))
def edge(w):
    ans=[]
    while w:ans.append(w[0]);w=D(w)
    return ans
def good(w):return max(edge(w),default=0)<=1
def pascal_front(w):
    # Diagonal entries of the pure binary triangle, right to left.
    a=w;out=[]
    while a:out.append(a[-1]);a=D(a)
    return tuple(out)
def scan_step(s,b):
    out=[]
    for x in s:b=abs(b-x);out.append(b)
    return tuple(out)
def dfa(age):
    states=[(3,)*(age+1)];ids={states[0]:0};E=[]
    for s in states:
        row=[]
        for b in (0,1):
            v=scan_step(s,b)
            if v not in ids:ids[v]=len(states);states.append(v)
            row.append(ids[v])
        E.append(tuple(row))
    A={i for i,s in enumerate(states)if max(s[:-1],default=0)<=1 and s[-1]>=2}
    P=[[]for _ in states]
    for i,e in enumerate(E):
        for j in e:P[j].append(i)
    C=set(A);q=deque(C)
    while q:
        for j in P[q.popleft()]:
            if j not in C:C.add(j);q.append(j)
    return states,E,A,C

def main():
    # Independent direct shrinking triangle vs rotated scan, all binary seas.
    mapped=0
    for F in range(1,12):
        fronts=set()
        for w in product((0,1),repeat=F-1):
            a=pascal_front(w);fronts.add(a)
            s=(3,)*7
            for b in a:s=scan_step(s,b)
            ee=edge(w+(3,)+(0,)*6)
            require(s==tuple(ee[F-1:F+6]),f'rotation mismatch {w}')
            mapped+=1
        require(len(fronts)==2**(F-1),'front transform not bijective')
    # Independent all-language product check for age2 regex DFA.
    states,E,A,C=dfa(2)
    def regex_step(q,b):
        if q==7:return 7
        if q<4:return (q+1)%4 if b==0 else (4 if q==0 else 7)
        if q==4:return 5 if b else 7
        if q==5:return 6 if b else 7
        return 6 if b==0 else 7
    seen={(0,0)};q=deque(seen)
    while q:
        i,j=q.popleft()
        require((i in A)==(j in (5,6)),'language equivalence failed')
        for b in (0,1):
            pair=(E[i][b],regex_step(j,b))
            if pair not in seen:seen.add(pair);q.append(pair)
    require(len(states)==23 and len(C)==7,'unexpected age2 automaton')
    # All-F threshold counts for fixed ages, using Cayley-Hamilton certificate:
    # actual generating function has denominator degree <=number coaccessible
    # states; proposed function C*x^h/((1-x)*(1-x^(2h))) degree2h+1.
    # Checking m+2h+1 initial coefficients forces equality identically.
    constants={2:(2,1),4:(4,3),6:(4,1),8:(8,7),10:(8,2),12:(8,3),14:(8,1)}
    rows=[];slope=Fraction(0)
    for r,(h,c) in constants.items():
        S,T,Acc,Co=dfa(r);m=len(Co);terms=m+2*h+1
        counts={0:1} if 0 in Co else {}
        for L in range(terms):
            val=sum(n for i,n in counts.items()if i in Acc)
            expected=c*((L+h)//(2*h))
            require(val==expected,f'GF coefficient failure age{r} length{L}')
            nxt={}
            for i,n in counts.items():
                for j in T[i]:
                    if j in Co:nxt[j]=nxt.get(j,0)+n
            counts=nxt
        rows.append({'age':r,'reachable_states':len(S),'coaccessible_states':m,'coefficients_checked':terms,'gf_numerator_coefficient':c,'gf_numerator_degree':h,'gf_denominator_cycle':2*h})
        slope+=Fraction(c, h*(2**r))
    odd=[]
    for r in range(1,15,2):
        S,T,A,C=dfa(r);require(not A,f'odd threshold present atage{r}');odd.append({'age':r,'reachable_states':len(S),'accepting':0})
    # Minimal alphabet-four hostile to the shield.
    require(good((1,2,3,4)) and edge((1,2,3,4,1))==[1,1,0,0,2],'hostile control')
    receipt={'rotated_direct_seeds':mapped,'age2_dfa_states':len(states),'age2_product_states':len(seen),'age2_language':'(0000)*11 OR (0000)*1110*','age2_table':[{'id':i,'state':s,'next':E[i],'accepting':i in A}for i,s in enumerate(states)],'all_F_count_certificates':rows,'odd_age_controls':odd,'c6_linear_liminf_lower_bound':str(slope),'probability_ratio_liminf_lower_bound':str(1+slope/2),'ratio_lower_decimal':float(1+slope/2),'scope':'All-F exact threshold counts only for ages2,4,6,8,10,12,14. No all-age asymptotic claimed.'}
    # Rebuild age2 accepting flags after loop variable reuse.
    _,_,AA,_=dfa(2)
    for row in receipt['age2_table']:row['accepting']=row['id']in AA
    print(json.dumps(receipt,indent=2,sort_keys=True))
if __name__=='__main__':main()
