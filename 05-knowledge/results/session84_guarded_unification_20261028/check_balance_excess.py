#!/usr/bin/env python3
from math import isqrt

def need(b,m):
    if not b:raise RuntimeError(m)
def main():
    N=100000;sig=[0]*(N+1);tau=[0]*(N+1);U=[0]*(N+1);S=[0]*(N+1)
    prime=bytearray(b'\1')*(N+1);prime[:2]=b'\0\0'
    sqfree=bytearray(b'\1')*(N+1)
    for p in range(2,isqrt(N)+1):
        if prime[p]:
            prime[p*p:N+1:p]=b'\0'*((N-p*p)//p+1)
            sqfree[p*p:N+1:p*p]=b'\0'*((N-p*p)//(p*p)+1)
    ps=[p for p in range(2,N+1) if prime[p]]
    for d in range(1,N+1):
        for n in range(d,N+1,d):sig[n]+=d;tau[n]+=1
        if d>1:
            for n in range(2*d,N+1,d):U[n]+=int(prime[d]);S[n]+=int(sqfree[d])
    balanced=[n for n in range(2,N+1) if tau[n]-2==U[n]+S[n]]
    E={350,490,550,572,650,748,836};families=[]
    for c,lo in [(12,5),(20,7),(28,11),(18,5)]:families.append({c*p for p in ps if p>=lo and c*p<=N})
    families.append({6*p*p for p in ps if p>=5 and 6*p*p<=N})
    expected=set(E)
    for f in families:need(not expected&f,'disjoint families');expected|=f
    actual={n for n in balanced if sig[n]>2*n}
    need(actual==expected,'complete balanced abundance list')
    need([n for n in balanced if sig[n]==2*n]==[],'no perfect balanced')
    need([n for n in balanced if sig[n]==2*n+2]==[650],'unique 650')
    pi=[0]*(N+1)
    for n in range(1,N+1):pi[n]=pi[n-1]+int(prime[n])
    count=0
    for X in range(N+1):
        count+=X in actual
        predicted=sum(max(pi[X//c]-k,0) for c,k in[(12,2),(20,3),(28,4),(18,2)])+max(pi[isqrt(X//6)]-2,0)+sum(e<=X for e in E)
        need(count==predicted,'counting formula')
    need(sig[650]==1302 and tau[650]==12 and U[650]==3 and S[650]==7,'named counts')
    print('PASS balanced abundance through',N,'sources',len(actual),'balanced count',len(balanced),'unique excess2 650; all cutoff counts agree')
if __name__=='__main__':main()
