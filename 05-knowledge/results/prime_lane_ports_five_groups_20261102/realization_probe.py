#!/usr/bin/env python3
"""Exact source-preserving component swaps and lane-candidate families.
No repository code is imported. Labels and apex cut sets are literal.
"""
from collections import Counter
from itertools import product, combinations
from functools import lru_cache
import json


def check(ok, message):
    if not ok: raise RuntimeError(message)


def mm(A, B):
    return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(4)) for j in range(4)) for i in range(4))


def mv(A,x):
    return tuple(sum(a*b for a,b in zip(row,x)) for row in A)


def sub(A,B): return tuple(tuple(a-b for a,b in zip(x,y)) for x,y in zip(A,B))
def scale(a,A): return tuple(tuple(a*x for x in row) for row in A)
I=tuple(tuple(int(i==j) for j in range(4)) for i in range(4))
ZERO=scale(0,I)


def singleton(e):
    return ((1,1-e,1-e,0),(0,2,0,e),(0,0,0,e),(0,0,0,1))


def triangle(a):
    return ((3,9-3*a,3-a,a*(3-a)),(0,12,0,3*a),(0,0,0,a),(0,0,0,3))


def kernel(u,d):
    return ((0,d,-u,d+u),(0,0,0,d),(0,0,0,-u),(0,0,0,0))


ALPHABET=tuple((1,e) for e in (0,1))+tuple((3,a) for a in range(4))

def mat(letter):
    return singleton(letter[1]) if letter[0]==1 else triangle(letter[1])


def wordstate(word):
    state=(1,0,0,1)
    for letter in word: state=mv(mat(letter),state)
    return state


def literal(word):
    rows=[]; S=0
    for n,a in word:
        old=len(rows)
        rows=[r|(((1<<n)-1)<<old) for r in rows]+[0]*n
        if n==3:
            for i in range(3): rows[old+i]|=1<<(old+(i+1)%3)
        else: check(n==1,'bad block')
        S|=((1<<a)-1)<<old
    return tuple(rows),S


@lru_cache(None)
def counts(rows):
    n=len(rows); full=(1<<n)-1
    dp=[[0]*n for _ in range(full+1)]
    H=[0]*(full+1); H[0]=1
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1,full+1):
        for v,c in enumerate(dp[mask]):
            if c:
                rest=rows[v]&~mask
                while rest:
                    b=rest&-rest; rest-=b; w=b.bit_length()-1
                    dp[mask|b][w]+=c
        H[mask]=sum(dp[mask])
    return H,dp


def actualstate(word):
    rows,S=literal(word); n=len(rows); full=(1<<n)-1
    H,ends=counts(rows)
    opposite=tuple(full^r^(1<<i) for i,r in enumerate(rows))
    _,starts=counts(opposite)
    q=sum(H[full^Q]*sum(starts[Q][v] for v in range(n) if S>>v&1)
          for Q in range(1,full+1))
    g=sum(ends[full][v] for v in range(n) if S>>v&1)
    child=tuple(r|((1<<n) if not S>>i&1 else 0) for i,r in enumerate(rows))+(S,)
    F=counts(child)[0][-1]
    return F,q,g,H[-1]


def actual_block(rows,S):
    n=len(rows); full=(1<<n)-1
    Hs,ends=counts(rows)
    _,starts=counts(tuple(full^row^(1<<i) for i,row in enumerate(rows)))
    H=Hs[-1]
    E=sum(Hs[L]*Hs[full^L] for L in range(full+1))
    P=sum(sum(ends[L][i] for i in range(n) if not S>>i&1)*Hs[full^L] for L in range(1,full+1))
    Q=sum(sum(starts[L][i] for i in range(n) if S>>i&1)*Hs[full^L] for L in range(1,full+1))
    st=sum(starts[full][i] for i in range(n) if S>>i&1)
    t=sum(ends[full][i] for i in range(n) if not S>>i&1)
    child=tuple(row|((1<<n) if not S>>i&1 else 0) for i,row in enumerate(rows))+(S,)
    f=counts(child)[0][-1]
    B=((H,P,H-st,f-H),(0,E,0,Q),(0,0,0,H-t),(0,0,0,H))
    return B,(H,E,P,Q,st,t,f)


def strong(rows):
    if not rows: return False
    reach=[r|(1<<i) for i,r in enumerate(rows)]
    for k in range(len(rows)):
        for i in range(len(rows)):
            if reach[i]>>k&1: reach[i]|=reach[k]
    return all(r==(1<<len(rows))-1 for r in reach)


def isprime(n):
    if n<2: return False
    d=2
    while d*d<=n:
        if n%d==0: return False
        d+=1
    return True


def power_two_r(n):
    if n<=1 or n&(n-1): return None
    return n.bit_length()-1


def family_membership(h):
    recipes=[]
    for a,cA,cB,offset in ((1,15,18,7),(2,9,15,1)):
        for name,c in (('U',cA),('V',cB)):
            if (h+offset)%c==0:
                r=power_two_r((h+offset)//c)
                if r is not None: recipes.append((name,a,r))
    return recipes


def main():
    checks=Counter()
    for u,d,n in ((1,3,3),(3,21,5),(25,259,7),(27,261,7)):
        K=kernel(u,d)
        square=tuple(tuple((u*u+d*d) if i==0 and j==3 else 0 for j in range(4)) for i in range(4))
        check(mm(K,K)==square and mm(mm(K,K),K)==ZERO,'nilpotent kernel')
        for a in range(n+1):
            # f cancels, so test several arbitrary f coefficients as an algebraic control.
            for f in (0,1,17):
                B=((n*u,(n-a)*d,(n-a)*u,f),(0,n*(u+d),0,a*d),(0,0,0,a*u),(0,0,0,n*u))
                for e in (0,1):
                    check(sub(mm(singleton(e),B),mm(B,singleton(e)))==scale(a-n*e,K),'general commutator')
                    checks['general_matrix_checks']+=1
    # Every actual labelled block/cut through order four; no symmetry filter.
    for n in range(1,5):
        es=tuple(combinations(range(n),2))
        for code in range(1<<len(es)):
            rows=[0]*n
            for bit,(i,j) in enumerate(es):
                u,v=(i,j) if code>>bit&1 else (j,i)
                rows[u]|=1<<v
            for S in range(1<<n):
                B,(H,E,P,Q,st,t,f)=actual_block(tuple(rows),S)
                D=E-H-P
                C0=((0,D,-st,H+Q-t),(0,0,0,Q),(0,0,0,t-H),(0,0,0,0))
                C1=((0,P,st-H,H+P-st),(0,0,0,E-H-Q),(0,0,0,-t),(0,0,0,0))
                check(sub(mm(singleton(0),B),mm(B,singleton(0)))==C0,'arbitrary block S0')
                check(sub(mm(B,singleton(1)),mm(singleton(1),B))==C1,'arbitrary block S1')
                check(D>=0 and E-H-Q>=0,'actual cover guards')
                if n>=2:
                    check(Q-st>=S.bit_count() and P-t>=n-S.bit_count(),'proper-cover strictness')
                checks['arbitrary_block_cut_profiles']+=1
    for length in range(5):
        for word in product(ALPHABET,repeat=length):
            if sum(n for n,a in word)>8: continue
            check(actualstate(word)==wordstate(word),('actual state',word))
            checks['direct_full_state_words']+=1
    prefixes=[w for k in range(3) for w in product(ALPHABET,repeat=k)]
    suffixes=[()]+[(a,) for a in ALPHABET]
    for pre in prefixes:
        state=wordstate(pre)
        for suffix in suffixes:
            W=I
            for letter in suffix: W=mm(mat(letter),W)
            for a in range(4):
                for e in (0,1):
                    U=pre+((1,e),(3,a))+suffix
                    V=pre+((3,a),(1,e))+suffix
                    diff=wordstate(V)[0]-wordstate(U)[0]
                    expected=mv(W,mv(scale(a-3*e,kernel(1,3)),state))[0]
                    check(diff==expected,('transported commutator',pre,suffix,a,e))
                    check((diff>0)-(diff<0)==(a-3*e>0)-(a-3*e<0),'strict sign')
                    checks['actual_word_swap_controls']+=1
    records=[]
    for r in range(1,31):
        for a,cA,cB,off in ((1,15,18,7),(2,9,15,1)):
            U=((1,1),)*r+((1,0),(3,a))
            V=((1,1),)*r+((3,a),(1,0))
            vals=[wordstate(U)[0],wordstate(V)[0]]
            check(vals==[cA*2**r-off,cB*2**r-off],'family formula')
            for name,h in zip(('U','V'),vals): check((name,a,r) in family_membership(h),'exact family inverse')
            if r<=5:
                for word,h in zip((U,V),vals):
                    check(actualstate(word)[0]==h,'family direct graph')
                    rows,S=literal(word); n=len(rows)
                    child=tuple(x|((1<<n) if not S>>i&1 else 0) for i,x in enumerate(rows))+(S,)
                    check(strong(child),'strong family')
                records.append({'r':r,'a':a,'U':vals[0],'V':vals[1],'vertices':r+5})
            if a==1:
                check(all(h%7 for h in vals),'a1 avoids7')
                check((vals[1]%5==0)==(r%4==2),'a1V mod5')
            else:
                check((vals[0]%7==0)==(r%3==2),'a2U mod7')
                check((vals[1]%7==0)==(r%3==0),'a2V mod7')
                check((vals[0]%5==0)==(r%4==2),'a2U mod5')
            checks['family_formula_and_residue_controls']+=1
    check(not family_membership(13),'family completeness hostile13')
    check(all(isprime(p) for p in (17,23,29,5,59,71,41,479)),'prime examples')
    check(9*2**2-1==7*5 and 15*2**2-1==59,'lane35-59')
    check(9*2**3-1==71 and 15*2**3-1==7*17,'lane71-119')
    check(9*2**5-1==7*41 and 15*2**5-1==479,'lane287-479')
    return {'status':'PROVED formulas; FINITE-EXACT implementation controls, no all-prime coverage',
            'checks':dict(checks),'small_families':records,
            'H13_family_membership':family_membership(13),
            'same_port_prime_lanes':[[23,29],[35,59],[71,119],[287,479]],
            'port_convention':'apex beats first a triangle vertices, all prefix vertices; loses to zero singleton'}


if __name__=='__main__': print(json.dumps(main(),indent=2,sort_keys=True))
