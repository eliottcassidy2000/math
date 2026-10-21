# Positive control: strong semicomplete digraphs (3^C(m,2), linear condensation) -> quasi-period (3/p).
# Hostile control: connected labeled graphs (2^C(m,2), unordered components, exp/log) -> expect no such law.
from math import comb
def strong_mod(base,M,p):
    T=[pow(base,j*(j-1)//2,p) for j in range(M+1)]; s=[0]*(M+1); row=[1]
    for k in range(1,M+1):
        row=[1]+[(row[i-1]+row[i])%p for i in range(1,k)]+[1]
        s[k]=(T[k]-sum(row[i]*s[i]*T[k-i] for i in range(1,k)))%p
    return s
def connected_mod(M,p):
    # g(m) = sum_{k=1}^m C(m-1,k-1) c(k) g(m-k)  (component of vertex 1), g(m)=2^C(m,2)
    g=[pow(2,j*(j-1)//2,p) for j in range(M+1)]; c=[0]*(M+1); row=[1]
    for m in range(1,M+1):
        # row = C(m-1, .)
        if m>1: row=[1]+[(row[i-1]+row[i])%p for i in range(1,m-1)]+[1]
        c[m]=(g[m]-sum(row[k-1]*c[k]*g[m-k] for k in range(1,m)))%p
    return c
def leg(a,p): return 1 if pow(a,(p-1)//2,p)==1 else -1
for p in [5,7,11,13,17,19,23]:
    M=4*p; s3=strong_mod(3,M,p); e=leg(3,p)
    ok3=all((s3[m+p-1]-e*s3[m])%p==0 for m in range(1,M-p+2))
    cg=connected_mod(M,p)
    okc=[L for L in range(1,2*p) if all((cg[m+L]-cg[m])%p==0 for m in range(1,M-L+1)) or all((cg[m+L]+cg[m])%p==0 for m in range(1,M-L+1))]
    print(f'p={p}: semicomplete strong quasi-period (3/p)={e:+d}: {ok3};  connected graphs: (anti)periods L<2p: {okc}')
