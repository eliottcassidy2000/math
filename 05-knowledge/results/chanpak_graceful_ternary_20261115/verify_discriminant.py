from math import isqrt
from random import Random
from verify_scaled import h

def inverse(S,N,ori):
    j=2*N%3
    B=432*S+(129 if ori==0 else -129)
    C=(129 if ori==0 else 258)*S+135+2*j
    m=S//3
    def base(t):return -432*t*t+B*t+C
    if N<base(1) or m<1:return None
    if N>=base(m):t=m
    else:
        disc=B*B-1728*(N-C)
        z=isqrt(disc)
        t=(B-z-(z*z<disc))//864
    x,y=(S-t,t) if ori==0 else (t,S-t)
    return (x,y,j) if h(x,y,j)==N else None

def need(x):
    if not x:raise RuntimeError("inverse discrepancy")

if __name__=='__main__':
    rng=Random(432)
    for _ in range(3000):
        S=rng.randrange(50,1<<260)
        t=rng.randrange(1,S//3+1)
        j=rng.randrange(3)
        for ori,(x,y) in enumerate(((S-t,t),(t,S-t))):
            N=h(x,y,j)
            need(inverse(S,N,ori)==(x,y,j))
            need(inverse(S,N+6,ori) is None)
    print("6,000 direct orientation inverses and near misses")
