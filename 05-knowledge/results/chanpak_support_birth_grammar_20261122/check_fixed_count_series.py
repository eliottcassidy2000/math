#!/usr/bin/env python3
"""Check exact transients and parity tails using the inherited finite inventory."""
from functools import cache
from collections import Counter
from math import comb

# Global strong inventories for these counts are inherited from Session 24;
# tuples list orders with repetitions for distinct cores.
A={3:(3,),5:(4,),9:(5,5),11:(5,),15:(5,5,6,6),25:(6,6,7),
   27:(6,6,7,7,7),33:(6,6,6,6,6,7,7,7)}
C={1:(0,1),3:(3,),5:(4,),9:(5,5),15:(5,5),25:(7,),27:(7,),33:(6,7)}

def need(v,msg):
    if not v:
        raise ArithmeticError(msg)

@cache
def actual(h,n):
    if n<0:
        return 0
    if n==0:
        return int(h==1)
    return actual(h,n-1)+sum(actual(h//d,n-m) for d,orders in A.items()
                              if h%d==0 for m in orders)

@cache
def core(h):
    out=Counter({(0,0):1}) if h==1 else Counter()
    for d,orders in A.items():
        if h%d==0:
            for m in orders:
                for (k,s),b in core(h//d).items():
                    out[k+1,m+s]+=b
    return out

def bpoly(q,k):
    if q>=0: return comb(q,k) if q>=k else 0
    return (-1)**k*comb(k-q-1,k)

def eventual(h,n):
    return sum(b*bpoly(n-m+k,k) for (k,m),b in core(h).items())

def error(h,n):
    return sum((-1)**(k+1)*b*comb(m-n-1,k)
               for (k,m),b in core(h).items() if n<m-k)

def selfconverse(h,n):
    return sum(actual(j,(n-m)//2) for j in range(1,h+1) if h%(j*j)==0
               for m in C.get(h//(j*j),()) if (n-m)%2==0)

def main():
    for h in (1,3,5,7,9,11,15,25,27,33):
        for n in range(30):
            need(actual(h,n)==eventual(h,n)+error(h,n),"transient replay")
    need([n for n in range(20) if error(15,n)]==[0,1,2],"exact H=15 onset")
    need(next(n for n in range(20) if actual(15,n))==5,"exact H=15 birth")
    for n in range(10,31):
        expected=(n-6)//2 if n%2==0 else (n-5)//2
        need(selfconverse(25,n)==expected,"partial negative pole")
        need(selfconverse(33,n)==1,"complete negative pole cancellation")
    print("exact onset; H=15 birth; partial and complete negative-pole controls")
if __name__=="__main__":
    main()
