#!/usr/bin/env python3
"""Exact prime-target prime-power fiber checks, stdlib only."""
from math import isqrt

def need(ok,why):
    if not ok: raise ValueError(why)
def primes(N):
    a=bytearray(b'\1')*(N+1); a[:2]=b'\0\0'
    for p in range(2,isqrt(N)+1):
        if a[p]: a[p*p:N+1:p]=b'\0'*((N-p*p)//p+1)
    return [i for i in range(2,N+1) if a[i]]
def rep(p,a): return (p**a-1)//(p-1)
def inverse(q):
    out=[]
    for a in primes((q+1).bit_length()):
        if (q-1)%a: continue
        lo,hi=2,q
        while lo<=hi:
            p=(lo+hi)//2; r=rep(p,a)
            if r<q:lo=p+1
            elif r>q:hi=p-1
            else:
                if all(p%d for d in range(2,isqrt(p)+1)):out.append((p,a,p**a))
                break
    return sorted(out)
def main():
    N=100000; ps=primes(N); pset=set(ps); direct={q:[] for q in ps}
    for p in ps:
        a=2; r=1+p
        while r<=N:
            if r in pset:
                direct[r].append((p,a,p**a))
                need(a in pset,'exponent prime'); need((r-1)%a==0,'order divides')
                need(pow(p,a,r)==1 and p%r!=1,'exact prime order')
            a+=1;r=r*p+1
    for q in ps:need(inverse(q)==sorted(direct[q]),f'inverse {q}')
    collisions={q:v for q,v in direct.items() if len(v)>1}
    need(collisions=={31:[(2,5,32),(5,3,125)]},'finite collisions')
    need(inverse(3)==[(2,2,4)],'3'); need(inverse(7)==[(2,3,8)],'7')
    need(inverse(13)==[(3,3,27)],'13'); need(inverse(11)==[],'11')
    need(inverse(8191)==[(2,13,8192)],'8191 prime-base guard')
    need(rep(90,3)==8191,'composite base repunit')
    # sigma(90^3)=sigma(2^3*3^6*5^3), not a geometric sum in composite base90.
    sigma90cube=(2**4-1)*(3**7-1)//2*(5**4-1)//4
    need(sigma90cube-90**3!=8191,'composite base hostile')
    print('PASS prime rows through',N,'prime rows',len(ps),'collisions',collisions)
if __name__=='__main__':main()
