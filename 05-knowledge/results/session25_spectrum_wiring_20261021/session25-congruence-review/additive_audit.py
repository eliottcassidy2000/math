"""Independent recurrence checks of additive prime-power quasi-periodicity."""
def require(ok,msg):
    if not ok: raise RuntimeError(msg)
def strong_counts(N,mod):
    T=[pow(2,n*(n-1)//2,mod) for n in range(N+1)]
    s=[0]*(N+1); row=[1]
    for n in range(1,N+1):
        row=[1]+[(row[k-1]+row[k])%mod for k in range(1,n)]+[1]
        s[n]=(T[n]-sum(row[k]*s[k]*T[n-k] for k in range(1,n)))%mod
    return s
checks=0
for p,maxr in [(3,6),(5,4),(7,3),(11,3),(13,2),(17,2)]:
    for r in range(1,maxr+1):
        mod=p**r; P=p**(r-1)*(p-1); threshold=2*r-1
        eps=1 if pow(2,(p-1)//2,p)==1 else -1
        vals=strong_counts(P+threshold+49,mod)
        for m in range(threshold,threshold+50):
            require((vals[m+P]-eps*vals[m])%mod==0,f'p={p},r={r},m={m}')
            checks+=1
        failures=[m for m in range(1,threshold) if (vals[m+P]-eps*vals[m])%mod]
        print('p,r,modulus,shift,proved-threshold,earlier-failures:',p,r,mod,P,threshold,failures)
print('PASS:',checks,'additive shift checks; earlier indices reported, not used to assert universal sharpness')
# Arbitrary unit bases: an independent inverse recurrence, including a negative base.
def inverse_component_counts(N,mod,base):
    T=[pow(base,n*(n-1)//2,mod) for n in range(N+1)]
    s=[0]*(N+1); row=[1]
    for n in range(1,N+1):
        row=[1]+[(row[k-1]+row[k])%mod for k in range(1,n)]+[1]
        s[n]=(T[n]-sum(row[k]*s[k]*T[n-k] for k in range(1,n)))%mod
    return s
base_checks=0
for base in [-2,3,5]:
    for p in [3,5,7]:
        if base%p==0:continue
        for r in [1,2,3]:
            P=p**(r-1)*(p-1); mod=p**r; threshold=2*r-1
            vals=inverse_component_counts(P+threshold+19,mod,base)
            eps=1 if pow(base,(p-1)//2,p)==1 else -1
            for n in range(threshold,threshold+20):
                require((vals[n+P]-eps*vals[n])%mod==0,f'base={base},p={p},r={r},n={n}')
                base_checks+=1
print('PASS:',base_checks,'arbitrary-unit-base checks')
from math import comb
connected_checks=0
for p in [3,5,7]:
    for r in [1,2,3]:
        P=p**(r-1)*(p-1);mod=p**r;threshold=2*r;N=P+threshold+19
        T=[pow(2,n*(n-1)//2,mod) for n in range(N+1)];b=[0]*(N+1)
        for n in range(1,N+1):
            b[n]=(T[n]-sum(comb(n-1,k-1)*b[k]*T[n-k] for k in range(1,n)))%mod
        eps=1 if pow(2,(p-1)//2,p)==1 else -1
        for n in range(threshold,threshold+20):
            require((b[n+P]-eps*b[n])%mod==0,f'connected p={p},r={r},n={n}')
            connected_checks+=1
print('PASS:',connected_checks,'connected-graph prime-power tail checks')
# Full CRT criterion against independent direct recurrence at all composites <=1000.
from math import lcm
limit=1000; modulus=lcm(*range(1,limit+1)); vals=strong_counts(limit,modulus)
def discrepancy(n):return vals[n]-(pow(2,(n-1)//2,modulus) if n%2 else 0)
def factor(n):
    factors=[];p=2
    while p*p<=n:
        if n%p==0:
            r=0
            while n%p==0:n//=p;r+=1
            factors.append((p,r))
        p+=1
    if n>1:factors.append((n,1))
    return factors
crt_checks=0
for n in range(4,limit+1):
    fs=factor(n)
    if len(fs)==1 and fs[0]==(n,1):continue
    conditions=[]
    for p,r in fs:
        if p==2:continue
        power=p**r;period=p**(r-1)*(p-1);threshold=2*r-1
        t=threshold+(n-threshold)%period
        require(threshold<=t<=power<=n,f'CRT representative n={n},p={p},r={r}')
        conditions.append(discrepancy(t)%power==0)
    require((discrepancy(n)%n==0)==all(conditions),f'full CRT n={n}')
    crt_checks+=1
print('PASS:',crt_checks,'full CRT checks on every composite n<=1000, including unrestricted odd prime powers')
