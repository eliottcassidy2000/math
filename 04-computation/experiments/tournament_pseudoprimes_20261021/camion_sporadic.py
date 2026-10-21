# Sporadic zero classes: for odd primes 7 <= p <= 400, count r in {4,...,p-1} with p | D(r); compare with Poisson.
from math import exp, factorial
def smod(M,p):
    T=[pow(2,j*(j-1)//2,p) for j in range(M+1)]; s=[0]*(M+1); row=[1]
    for k in range(1,M+1):
        row=[1]+[(row[i-1]+row[i])%p for i in range(1,k)]+[1]
        acc=0
        for i in range(1,k): acc=(acc+row[i]*s[i]*T[k-i])%p
        s[k]=(T[k]-acc)%p
    return s
primes=[p for p in range(7,401) if all(p%d for d in range(2,int(p**.5)+1))]
counts=[]; expct=[]
for p in primes:
    s=smod(p,p)
    z=[r for r in range(4,p) if (s[r]-(pow(2,(r-1)//2,p) if r%2 else 0))%p==0]
    counts.append(len(z)); expct.append((p-4)/p)
lam=sum(expct)/len(expct)
from collections import Counter
C=Counter(counts); N=len(primes)
print('primes',N,'mean sporadic classes',sum(counts)/N,'model mean',lam)
for k in range(0,6): print(k, C[k], round(N*exp(-lam)*lam**k/factorial(k),1))
