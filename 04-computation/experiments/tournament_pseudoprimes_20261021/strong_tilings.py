# sigma(n): strong tilings = tournaments on [n] containing the path 1->2->...->n that are strong.
# Composition law: F(x) = 1 + sum_{n>=1} 2^C(n-1,2) x^n = 1/(1 - Sigma(x)).
N=40
b=[1]+[2**((n-1)*(n-2)//2) for n in range(1,N+1)]
sig=[0]*(N+1)
for n in range(1,N+1):
    sig[n]=b[n]-sum(sig[k]*b[n-k] for k in range(1,n))
print('sigma(1..12)=',sig[1:13])
# check against hs(n)/n! for n<=5: 1,0,1,5,50
for p in [3,5,7,11,13,17,19,23]:
    print(p, 'sigma(p) mod p =', sig[p]%p, ' sigma(p-1) mod p=', sig[p-1]%p, ' sigma(p+1) mod p=', sig[p+1]%p if p+1<=N else None)
