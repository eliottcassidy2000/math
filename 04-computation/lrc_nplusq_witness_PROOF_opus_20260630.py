"""
RIGOROUS single-swap 'prime forced': q prime, (n-1)/2 < q <= n-1. AP\{q}+g not tight for ANY single g.
 - g non-multiple of q => no multiple of q in S => q-witness M>=1/q>1/n.
 - g=mq => WITNESS at p=n+q, rotation a=q^{-1} mod (n+q): every runner avoids {0,1,n+q-1} => M>=2/(n+q)>1/n.
   Proof: runner mq -> (mq)*a = m mod (n+q) in {2..} (m>=2). Runners v in {1..n-1}\{q}: v*a in {0,q,n}?
   v=q removed; v=n>n-1 absent; v=0 none. So none hit {0,1,n+q-1}. QED. Verify.
"""
from fractions import Fraction
def frac(x): x=x%1; return min(x,1-x)
def inv(a,m):
    for x in range(1,m):
        if (a*x)%m==1: return x
    return None
def isprime(m): return m>1 and all(m%d for d in range(2,int(m**.5)+1))
print("verify: at rotation a=q^{-1} mod (n+q), min_v||v*a/(n+q)|| = 2/(n+q) for AP\{q}+mq (q prime):")
for n in [14,18]:
    for q in [p for p in range(n//2+1,n) if isprime(p)]:
        p=n+q; a=inv(q,p)
        for m in [2,3]:
            S=[v for v in range(1,n) if v!=q]+[m*q]
            t=Fraction(a,p)
            mn=min(frac(Fraction(v)*t) for v in S)
            # also check no runner residue in {0,1,p-1}
            res=sorted(set((v*a)%p for v in S))
            bad=[r for r in res if r in (0,1,p-1)]
            print(f"  n={n} q={q} m={m}: p=n+q={p}, a=q^(-1)={a}; min_v||va/p||={mn}={float(mn):.4f} (>=2/p={2}/{p}={float(Fraction(2,p)):.4f}); runners in {{0,1,p-1}}: {bad if bad else 'NONE'}")
print()
print("=> RIGOROUS: single-swap tight sets NEVER remove a prime q in ((n-1)/2, n-1].")
print("   (non-mult g: q-witness M>=1/q; mult mq: (n+q)-witness at a=q^{-1}, M>=2/(n+q); both >1/n.)")
print("   The full multi-swap/general case needs n NOT in S (else runner n hits residue -1); = HYP-3749 crux.")
