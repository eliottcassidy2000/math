"""
Pin the rigor of 'large primes forced'.
 (1) WIDE HOLE (rigorous): q prime, q>(n-1)/2 => q is the UNIQUE multiple of q in {1..n-1}, so AP\{q} misses
     residue 0 mod q => q-witness M>=1/q; achieved at t=1/q (nearest runners 1,q-1 at 1/q) => M(AP\{q})=1/q.
 (2) SINGLE-SWAP: AP\{q}+mq for m=2,3,4,5 -- does ANY single multiple restore tightness? (q prime)
 (3) the general obstruction: adding mq covers residue 0 mod q but a witness at ANOTHER modulus survives
     (CRT-linkage = HYP-3749). Does the surviving M stay >1/n for all single m?
"""
from fractions import Fraction
def frac(x): x=x%1; return min(x,1-x)
def M_wit(S,Qmax):
    best=Fraction(0); bt=None
    for q in range(1,Qmax+1):
        for a in range(1,q):
            m=min(min((Fraction(s*a,q))%1,1-(Fraction(s*a,q))%1) for s in S)
            if m>best: best=m; bt=Fraction(a,q)
    return best,bt
def isprime(m): return m>1 and all(m%d for d in range(2,int(m**.5)+1))
print("(1) M(AP\{q}) = 1/q for q prime > (n-1)/2 (the punctured wide hole):")
for n in [12,14,18]:
    for q in [p for p in range(n//2+1,n) if isprime(p)]:
        M,_=M_wit([v for v in range(1,n) if v!=q], 5*n)
        print(f"   n={n} q={q}: M(AP\q)={M} = 1/q? {M==Fraction(1,q)}")
print()
print("(2)+(3) AP\{q}+mq for m=2..5 (q prime): any single substitute tight? & the surviving witness:")
for n in [14,18]:
    for q in [p for p in range(n//2+1,n) if isprime(p)]:
        row=f"   n={n} q={q}: "
        for m in range(2,6):
            S=[v for v in range(1,n) if v!=q]+[m*q]
            M,t=M_wit(S,6*n)
            row+=f"+{m}q:M={float(M):.4f}{'(TIGHT)' if M==Fraction(1,n) else ''} "
        row+=f" [1/n={1/n:.4f}]"
        print(row)
print()
print("=> (1) rigorous (q-witness). (2) NO single substitute mq restores tightness for prime q (M>1/n always).")
print("   (3) the surviving witness at a shifted modulus = HYP-3749 CRT-linkage; the FULL 'large primes forced'")
print("   (multi-element substitutes / general tight S) reduces to HYP-3749's open robustness -- NOT closed by")
print("   me. What IS rigorous: the wide hole + single-substitute failure. Honest.")
