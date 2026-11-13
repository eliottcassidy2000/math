#!/usr/bin/env python3
from check import H,join,tt,tournaments,endpoint
from itertools import permutations
from math import comb
from functools import lru_cache
import json

def require(x,msg):
    if not x:raise AssertionError(msg)

def exterior(A,B,E):
    na=len(A);u=list(join(A,B))
    for a,b in E:u[a]^=1<<(na+b);u[na+b]^=1<<a
    return tuple(u)
def coefficients(A,B,E):
    na=len(A);n=na+len(B);out=[0]*(min(na,len(B))+2);u=exterior(A,B,E)
    for p in permutations(range(n)):
        s=int(p[0]>=na)+int(p[-1]<na);r=0;ok=True
        for a,b in zip(p,p[1:]):
            if a<na<=b:
                s+=1;r+=int((a,b-na) in E)
            elif not u[a]>>b&1:ok=False;break
        if not ok:continue
        require(s<len(out),'slot bound')
        for k in range(r,s+1):out[k]+=comb(s-r,k-r)
    return out

def covers(M):
    n=len(M)
    @lru_cache(None)
    def f(mask,k):
        if k==0:return int(mask==0)
        if mask.bit_count()<k:return 0
        S=mask;total=0
        while S:
            total+=H(M,S)*f(mask^S,k-1);S=(S-1)&mask
        return total
    return [f((1<<n)-1,k) for k in range(n+1)]
def graph(A,M,B,E):
    u=list(join(A,M,B));off=len(A)+len(M)
    for a,b in E:u[a]^=1<<(off+b);u[off+b]^=1<<a
    return tuple(u)

def main():
    middle=[()]+[M for n in range(1,4) for M in tournaments(n)]
    cases=single=0
    for na in range(1,3):
      for nb in range(1,3):
       for A in tournaments(na):
        for B in tournaments(nb):
         pairs=[(a,b) for a in range(na) for b in range(nb)]
         for mask in range(1<<len(pairs)):
          E={e for i,e in enumerate(pairs) if mask>>i&1};d=coefficients(A,B,E)
          if E:
           require(d[2]>=len(E)*H(A)*H(B),'positive second cover')
          for M in middle:
           C=covers(M);expected=sum(x*y for x,y in zip(d,C));actual=H(graph(A,M,B,E))
           require(actual==expected,(A,B,M,E,d,C,actual,expected))
           require(actual>=H(graph(A,tt(len(M)),B,E)),'middle lower bound')
           cases+=1
           if E and M:
            W=sum(H(A,((1<<na)-1)^(1<<a))*H(B,((1<<nb)-1)^(1<<b)) for a,b in E)
            full=(1<<len(M))-1
            odd=sum(H(M,S)*H(M,full^S) for S in range(1,full+1) if S.bit_count()%2)
            require(actual>=H(A)*H(B)*H(M)+2*W*odd,'OCF feedback lower bound')
            require(actual>=H(A)*H(B)+(1<<len(M))*W,'sharp size lower bound')
           if len(E)==1 and M:
            a,b=next(iter(E));lr=endpoint(A,a,True)*endpoint(B,b,False);c2=C[2] if len(C)>2 else 0
            require(actual==H(A)*H(B)*H(M)+lr*(2*H(M)+c2),'single formula')
            require(actual>=H(A)*H(B)+(1<<len(M))*lr,'single bound')
            require((actual==H(A)*H(B)+(1<<len(M))*lr)==(H(M)==1),'equality')
            single+=1
    out={'arbitrary_middle_exact_cases':cases,'single_backarc_cases':single,'status':'PASS'}
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
