# (a) even-n Korselt check: n | s(n) predicted iff for every odd p|n: s(n mod (p-1) [in 1..p-1]) = 0 mod p
#     (n squarefree odd part; 2-part automatic).  (b) D(m) mod p^2 universality tests.
from math import comb
import sys
def smod(M,m):
    T=[pow(2,j*(j-1)//2,m) for j in range(M+1)]; s=[0]*(M+1); row=[1]
    for k in range(1,M+1):
        row=[1]+[(row[i-1]+row[i])%m for i in range(1,k)]+[1]
        acc=0
        for i in range(1,k): acc=(acc+row[i]*s[i]*T[k-i])%m
        s[k]=(T[k]-acc)%m
    return s
def factor(n):
    f={};d=2;x=n
    while d*d<=x:
        while x%d==0: f[d]=f.get(d,0)+1; x//=d
        d+=1
    if x>1: f[x]=f.get(x,0)+1
    return f
circ=lambda m,mod: (pow(2,(m-1)//2,mod) if m%2 else 0)
# zero classes for odd primes p<=100 (all residues r in 1..p-1)
Zc={}
for p in [q for q in range(3,101) if all(q%d for d in range(2,q))]:
    s=smod(p,p); Zc[p]={r for r in range(1,p) if (s[r]-circ(r,p))%p==0}
agree=dis=0
for n in range(4,401,2):
    f=factor(n); odd=[p for p in f if p>2]
    if any(f[p]>1 for p in odd) or any(p>100 for p in odd): continue
    pred=all(((n-1)%(p-1))+1 in Zc[p] for p in odd)
    s=smod(n,n); direct=(s[n]%n==0)
    if pred==direct: agree+=1
    else: dis+=1; print('DIS even',n,pred,direct)
print('even n<=400 (odd part squarefree): agree',agree,'disagree',dis)
print('universal zero classes include 1,2,3:', all({1,2,3}<=Zc[p] or p<5 for p in Zc))
for q in [9,25,27,49,121,169,289,361]:
    M=4*q if q<100 else 300
    s=smod(M,q)
    bad=[m for m in range(1,M+1,2) if (s[m]-circ(m,q))%q]
    print(f'mod {q}: odd m<={M} with D(m) != 0: {len(bad)}; first few {bad[:12]}')
# Lemma B of the independent audit (MISTAKE-651 corrected the pre-audit claim "D is not periodic mod p^2",
# which had tested periods from m = 1): s(m + p(p-1)) = (2/p) s(m) (mod p^2) for all m >= 3, failing only at m = 1, 2.
for p in [3,5,7,11,13]:
    q=p*p; L=p*(p-1); M=3*L+10; s=smod(M,q); eps=1 if p%8 in (1,7) else -1
    fails=[m for m in range(1,M-L+1) if (s[m+L]-eps*s[m])%q]
    D=[None]+[(s[m]-circ(m,q))%q for m in range(1,M+1)]
    failsD=[m for m in range(1,M-L+1) if (D[m+L]-eps*D[m])%q]
    print(f'Lemma B check p={p}: period L={L}, eps={eps:+d}: s fails at m={fails}, D fails at m={failsD}')
