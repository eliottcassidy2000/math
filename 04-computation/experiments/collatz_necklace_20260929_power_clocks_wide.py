"""Hostile probe: perfect-power clocks |2^K - q^X| = m^r (r>=2) with X>=2, over ALL odd q <= Qmax, 2<=X<=K<=Kmax.
Anything outside {Catalan, Pythagorean family q=2^(K-2)+1 with X=2, the FC solutions 2^5+7^2=3^4, 7^3+13^2=2^9, 2^7+17^3=71^2 (in all prime-power readings)} would be a new Fermat-Catalan-type identity."""
import gmpy2, sys
from gmpy2 import mpz
Qmax=int(sys.argv[1]); Kmax=int(sys.argv[2])
found=[]
two=[mpz(1)<<K for K in range(Kmax+1)]
for q in range(3,Qmax+1,2):
    qp=mpz(q)
    P=[mpz(1)]
    for X in range(1,Kmax+1): P.append(P[-1]*qp)
    for X in range(2,Kmax+1):
        for K in range(X,Kmax+1):
            g=two[K]-P[X]; a=abs(g)
            if a==1 or gmpy2.is_power(a):
                r=None
                for rr in range(64,1,-1):
                    m,ex=gmpy2.iroot(a,rr)
                    if ex: r=(int(m),rr); break
                found.append((q,K,X,int(g),r))
print("odd q<=%d, 2<=X<=K<=%d: %d perfect-power clocks"%(Qmax,Kmax,len(found)))
fam=0; other=[]
for (q,K,X,g,r) in found:
    if X==2 and q==(1<<(K-2))+1: fam+=1; continue
    other.append((q,K,X,g,r))
print("Pythagorean-family members:",fam)
print("all others:")
for o in other: print(" ",o)
