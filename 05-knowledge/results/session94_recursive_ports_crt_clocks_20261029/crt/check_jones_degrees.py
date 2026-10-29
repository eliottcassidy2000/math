"""Exact degree exclusion + generator signs; Jones classification is external."""
from math import isqrt
import json

def require(x,m):
    if not x:raise RuntimeError(m)

def fib(n):
    a,b=0,1
    for _ in range(n):a,b=b,a+b
    return a

out=[]
for n in (25,55):
    N=fib(n-1)+fib(n+1);fixed=fib((n+3)//2)
    require(N%11==0 and N>11 and N not in (11,23),'nonprime/nonexceptional degree')
    q=N-1
    require(q%2==0 and q&(q-1)!=0,'d2: even non-power-of-two')
    matches=[]
    for q in range(2,isqrt(N)+1):
        total=1+q+q*q;d=3
        while total<=N:
            if total==N:matches.append([q,d])
            total=total*q+1;d+=1
    require(not matches,'no projective repunit degree')
    transpositions=(N-fixed)//2
    require(N%2==1 and (N-fixed)%2==0 and transpositions%2==0,'even generators')
    out.append({'period':n,'degree':N,'reflection_fixed':fixed,'reflection_transpositions':transpositions,'integer_bases_checked':[2,isqrt(N)],'projective_matches':matches,'consequence_given_sparse_primitivity_and_Jones':'alternating group of stated degree'})
print(json.dumps(out,indent=2))
