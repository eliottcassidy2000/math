# spacer_identity.py (wave 22). Usage: python3 spacer_identity.py pairs.txt (from fcollide.c)
# Verify the spacer factorization identity on every collision pair of F(u)=908u+9g(u) up to 2^24:
#   u = 2^m V + (2^q - 1) 2^p + c,  u' = 2^m V + 2^(q+p) + c',  m = q+p+1, (c_{p-1},c'_{p-1}) != (1,0)
#   (2V+1) * (2^(q+1) D + 3) = 1834 t + 2(Psi(c)-Psi(c')) - 2 rho + 2(c0-c0') - 2 v0 + kappa
#   D = rev_p(c') - rev_p(c) - 1 (never 0),  t = (2^p + c' - c)/9,  rho = rev_p(c),
#   kappa = 2 c_{p-1} - 1 if q >= 1, else 1 - 2 c'_{p-1}.
import sys, collections
def Psi(x):
    b=[r for r in range(x.bit_length()) if (x>>r)&1]; s=0
    for i in range(len(b)):
        for k in range(i+1,len(b)):
            d=b[k]-b[i]; s+= 2 if d==1 else 1<<(d-1)
    return s
def rev(x,n): return int(format(x,'0%db'%n)[::-1],2) if n>0 else 0
pairs=[tuple(map(int,l.split())) for l in open(sys.argv[1] if len(sys.argv)>1 else 'coll24.txt')]
ok=0; Dzero=0; qdist=collections.Counter()
for u,up in pairs:
    x=u^up; m=x.bit_length(); V=u>>m; k=m-1
    W=u&((1<<k)-1); Wp=up&((1<<k)-1); q=0
    while q<k and (W>>(k-1-q))&1 and not (Wp>>(k-1-q))&1: q+=1
    qdist[q]+=1
    p=k-q; c=W&((1<<p)-1); cp=Wp&((1<<p)-1)
    assert u==(V<<m)+(((1<<q)-1)<<p)+c and up==(V<<m)+(1<<(q+p))+cp
    rho=rev(c,p); D=rev(cp,p)-rho-1; Dzero+=(D==0)
    t,rem=divmod((1<<p)+cp-c,9); assert rem==0
    kap=(2*((c>>(p-1))&1)-1) if q>=1 else (1-2*((cp>>(p-1))&1))
    lhs=(2*V+1)*((1<<(q+1))*D+3)
    rhs=1834*t+2*(Psi(c)-Psi(cp))-2*rho+2*((c&1)-(cp&1))-2*(V&1)+kap
    ok+=(lhs==rhs)
print(f"identity holds on {ok} of {len(pairs)} collision pairs; D=0 occurs {Dzero} times")
print("spacer length q distribution:",sorted(qdist.items()))
