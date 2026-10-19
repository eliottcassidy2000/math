#!/usr/bin/env python3
"""Exact controls for arbitrary terminal-upset/initial-ideal Gold contexts."""
from itertools import product
from fractions import Fraction
from pathlib import Path
import json
from gold_terminal_chain import count,extensions,gold,deficits,require
from gold_four_collar import CORE,PORTS,flags,corner,table_pred
FIRST=(5,7);SECOND=(4,7)
EXPECTED=[[[325,0],[201,138]],[[112,0],[84,84]],[[80,0],[48,32]],[[28,0],[21,21]]]

def subsets(p,up):
    n=len(p)
    for s in range(1<<n):
        if up:
            ok=all(not(s>>i&1) or all(not(p[j]>>i&1) or s>>j&1 for j in range(n)) for i in range(n))
        else:ok=all(not(s>>i&1) or p[i]&s==p[i] for i in range(n))
        if ok:yield s

def insert(A,U,C,V):
    na,nb=len(A),len(CORE)
    p=list(A)+[(x<<na)|((1<<na)-1) for x in CORE]
    p += [(x<<(na+nb))|((1<<(na+nb))-1) for x in C]
    p[na+PORTS[0]]&=~U
    for z in range(len(C)):
        if V>>z&1:p[na+nb+z]&=~(1<<(na+PORTS[1]))
    return tuple(p)

def tail(w,S,rev):
    k=0
    for x in reversed(w) if rev else w:
        if not(S>>x&1):break
        k+=1
    return k

def hand():
    F=[[[0,0],[0,0]] for _ in range(4)]
    for h in range(3,7):
        for q in range(2,h):
            for t in range(1,q):
                i=0 if h==6 else 1
                j=0 if h>=5 else 1
                K=36-t*(t+1)//2
                F[0][i][j]+=K
                F[2][i][j]+=8
                if t==1:
                    F[1][i][j]+=28
                    F[3][i][j]+=7
    return F

def main():
    F=flags(FIRST,SECOND)
    require(F==EXPECTED==hand(),'basis tables')
    require(all(gold(t) for t in F),'coefficient cone')
    contexts=[(),(0,),(0,1),(0,0),(0,0,1),(0,1,3)]
    tests=0;maxrs=[Fraction(0),Fraction(0)];minmargin=None
    for A,C in product(contexts,repeat=2):
        ea,ec=count(A),count(C)
        aw=list(extensions(A));cw=list(extensions(C))
        for U,V in product(subsets(A,True),subsets(C,False)):
            r=Fraction(sum(tail(w,U,True) for w in aw),ea)
            s=Fraction(sum(tail(w,V,False) for w in cw),ec)
            p=insert(A,U,C,V)
            got=table_pred(p,tuple(x+len(A) for x in FIRST),tuple(x+len(A) for x in SECOND))
            expected=[[ea*ec*x for x in row] for row in corner(F,r,s)]
            require(got==expected,'general collar table')
            require(min(sum(deficits(got),[]))>0,'strict Gold')
            maxrs=[max(maxrs[0],r),max(maxrs[1],s)]
            margin=Fraction(min(sum(deficits(got),[])),sum(sum(got,[])))
            minmargin=margin if minmargin is None else min(minmargin,margin)
            tests+=1
    chains=[]
    for r,s in product((0,1,2,5,10,30),repeat=2):
        A=tuple((1<<i)-1 for i in range(r));C=tuple((1<<i)-1 for i in range(s))
        p=insert(A,(1<<r)-1,C,(1<<s)-1)
        got=table_pred(p,tuple(x+r for x in FIRST),tuple(x+r for x in SECOND))
        require(got==corner(F,r,s),'long chain collar')
        require(min(sum(deficits(got),[]))>0,'long collar strictness')
        chains.append(dict(left=r,right=s,table=got,minimum_deficit=min(sum(deficits(got),[]))))
    out={'status':'proved nonnegative flag-cone Gold certificate; exact independent controls passed','first':FIRST,'second':SECOND,'flag_tables':F,'basis_deficits':[deficits(t) for t in F],'general_context_count':tests,'maximum_observed_parameters':[str(x) for x in maxrs],'minimum_observed_normalized_slack':str(minmargin),'chain_contexts':chains}
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
