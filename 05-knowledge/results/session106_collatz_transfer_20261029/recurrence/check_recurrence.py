#!/usr/bin/env python3
"""Standard-library, explicit checks; normal and -O have identical semantics."""
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
import json

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def v2(n):
    need(isinstance(n, int) and n > 0, 'valuation domain')
    return (n & -n).bit_length()-1

def native(b, n):
    x=b*n+1
    a=v2(x)
    return x//(1<<a), a

@dataclass(frozen=True)
class State:
    r:int=1
    M:int=2
    A:int=1
    t:int=1
    lo:int=0
    hi:object=None

    def child(self,b,a):
        need(b>=3 and b%2 and a>=1, 'branch domain')
        P=1<<a
        s=((P//2-(b*self.t+1)//2)*pow(b*self.A,-1,P))%P
        lo=max(0, -((s-self.lo)//P))
        hi=None if self.hi is None else (self.hi-s)//P
        if hi is not None and lo>hi:
            return None
        numerator=b*self.t+1+2*b*self.A*s
        need(numerator%P==0 and (numerator//P)%2==1, 'exact division guard')
        return State(self.r+self.M*s,self.M*P,b*self.A,numerator//P,lo,hi)

    def survive(self):
        # endpoint - source = (t-r)+(2A-M)q >= 0
        d=2*self.A-self.M
        h=self.r-self.t
        lo,hi=self.lo,self.hi
        if d>0:
            lo=max(lo,-((-h)//d))
        elif d<0:
            bound=h//d
            hi=bound if hi is None else min(hi,bound)
        elif h>0:
            return None
        if hi is not None and lo>hi:
            return None
        return State(self.r,self.M,self.A,self.t,lo,hi)

    def values(self):
        need(self.hi is not None,'finite materialization only')
        return [self.r+self.M*q for q in range(self.lo,self.hi+1)]

    def mass(self,z):
        start=z**(self.r+self.M*self.lo)
        if self.hi is None:
            return start/(1-z**self.M)
        return start*(1-z**(self.M*(self.hi-self.lo+1)))/(1-z**self.M)

def run(b,N,depth):
    states=[State(lo=1,hi=(N-1)//2)] # positive odd sources >1
    expected={n:n for n in range(3,N+1,2)}
    rows=[]
    for k in range(1,depth+1):
        nxt=[]
        for st in states:
            max_y=st.t+2*st.A*st.hi
            for a in range(1,(b*max_y+1).bit_length()):
                child=st.child(b,a)
                if child is not None:
                    child=child.survive()
                    if child is not None:
                        nxt.append(child)
        expected={n:native(b,y)[0] for n,y in expected.items()}
        expected={n:y for n,y in expected.items() if y>=n}
        actual={}
        for st in nxt:
            for q in range(st.lo,st.hi+1):
                n=st.r+st.M*q
                need(n not in actual,'source partition disjointness')
                actual[n]=st.t+2*st.A*q
        need(actual==expected,'bounded-source sieve vs native replay')
        rows.append({'depth':k,'survivors':len(actual),'cells':len(nxt)})
        states=nxt
    return rows

def main():
    count=0
    for b in (3,5,7):
        for k in range(1,5):
            for w in product(range(1,5),repeat=k):
                st=State()
                A=B=1; C=0
                for a in w:
                    st=st.child(b,a)
                    A,B,C=b*A,B*(1<<a),b*C+B
                need(st.M==2*B and st.A==A,'clock agreement')
                need(st.r==((B-C)*pow(A,-1,2*B))%(2*B),'independent affine cylinder')
                need(st.t==(A*st.r+C)//B,'endpoint carry')
                for q in (0,1,2,7):
                    n=st.r+st.M*q; y=n; got=[]
                    for _ in w:
                        y,a=native(b,y); got.append(a)
                    need(tuple(got)==w and y==st.t+2*A*q,'native exact word')
                # Infinite interval survives iff each sampled q satisfies current endpoint test.
                sv=st.survive()
                for q in range(80):
                    member=sv is not None and q>=sv.lo and (sv.hi is None or q<=sv.hi)
                    need(member==((st.t-st.r)+(2*A-st.M)*q>=0),'interval linear cut')
                count+=1
    # Independently test inherited hostile append/reclose rather than temporal orbit.
    def fixed(b,w):
        x=Fraction(0); slope=Fraction(1)
        for a in w:
            x=(b*x+1)/(1<<a); slope=b*slope/(1<<a)
        return x/(1-slope)
    need(fixed(3,(1,1,5))==Fraction(19,101),'101 hostile 1')
    need(fixed(3,(1,3,3))==Fraction(31,101),'101 hostile 2')
    need(fixed(3,(1,1,5,1))==Fraction(37,35),'35 hostile')
    need(fixed(3,(1,3,3,1))==Fraction(221,175),'175 hostile')
    need([native(5,n)[0] for n in (13,33,83)]==[33,83,13],'5-cycle')
    st=State(lo=1)
    for a in (1,1,5)*4:
        st=st.child(5,a)
        need(st is not None,'cycle branch exists')
        st=st.survive()
        need(st is not None,'cycle remains unpaid')
        need((13-st.r)%st.M==0,'cycle source cylinder')
        q=(13-st.r)//st.M
        need(q>=st.lo and (st.hi is None or q<=st.hi),'cycle retained')
    # Finite geometric generating functions independently sum literal source atoms.
    st=State(lo=1,hi=35); z=Fraction(1,2)
    total=Fraction(0)
    for a in range(1,(3*(st.t+2*st.A*st.hi)+1).bit_length()):
        ch=st.child(3,a)
        if ch is not None:
            total+=ch.mass(z)
    need(total==st.mass(z),'source-weight extension conservation')
    need(st.mass(z)==sum((z**n for n in st.values()),Fraction(0)),'geometric mass formula')
    result={'status':'PASS','word_states':count,'bounded_sieves':{},
            'hostiles':{'five_cycle':[13,33,83],'same_denominator_then_reclose':[101,35,175]}}
    for b in (3,5,7):
        result['bounded_sieves'][str(b)]=run(b,4095,18)
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__':
    main()
