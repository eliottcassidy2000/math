#!/usr/bin/env python3
from fractions import Fraction
from math import isqrt,prod

def need(b,m):
    if not b:raise RuntimeError(m)
def primes(N):
    a=bytearray(b'\1')*(N+1);a[:2]=b'\0\0'
    for p in range(2,isqrt(N)+1):
        if a[p]:a[p*p:N+1:p]=b'\0'*((N-p*p)//p+1)
    return [p for p in range(2,N+1) if a[p]]
def bracket(p):
    j=(isqrt(p-1)+1)//2
    return (2*j-1)**2,(2*j+1)**2
def main():
    ps=primes(100000)
    def E(x):return [p for p in ps if bracket(p)[0]<x*p<=bracket(p)[1]]
    for x in [Fraction(2),Fraction(25,13)+Fraction(1,100000),Fraction(25,11)]:
        need(E(x)==[2,3,11],f'plateau {x}')
    need(E(Fraction(25,13))==[2,3,11,13],'lower boundary')
    need(E(Fraction(25,11)+Fraction(1,100000))==[2,3],'upper boundary')
    pairs=[(k,p) for k in range(2,10) for p in ps if bracket(p)[0]<k*p<=bracket(p)[1]]
    need(pairs==[(2,2),(2,3),(2,11),(3,2),(3,3),(4,2)],'integer pairs')
    need(prod(E(Fraction(2)))==66,'product')
    def sigma(n):return sum(d for d in range(1,n+1) if n%d==0)
    need(sigma(22)==36 and sigma(66)==144,'square lift')
    need([sigma(n) for n in [132,198,726]]==[336,468,1596],'promotion totals')
    need(sigma(364)==784,'compatible elsewhere')
    for n,expected in [(66,(6,3,6)),(132,(10,3,7)),(198,(10,3,7)),(726,(10,3,7))]:
        ds=[d for d in range(2,n) if n%d==0]
        U=sum(d in ps for d in ds)
        S=sum(not any(d%(p*p)==0 for p in ps if p*p<=d) for d in ds)
        need((len(ds),U,S)==expected,'FUS literal count')
    print('PASS bracket plateau and integer pairs; 66 square lift and three non-square promotions')
if __name__=='__main__':main()
