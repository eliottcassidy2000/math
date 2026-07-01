"""
Search 'AP with one element replaced' tight sets: S = ({1..n-1} minus k) plus g, M(S)=1/n?
Scan k in {1..n-1}, g in [n, 4n]. Directly find the sporadic family + the (k,g) pattern.
Confirm tightness with large Qmax.
"""
from fractions import Fraction
def fr(x): x=x%1.0; return min(x,1-x)
def Mfloat(S,Qmax):
    best=0.0
    for q in range(1,Qmax+1):
        for a in range(1,q):
            m=min(fr(s*a/q) for s in S)
            if m>best: best=m
    return best
def Mexact(S,Qmax):
    best=Fraction(0)
    for q in range(1,Qmax+1):
        for a in range(1,q):
            m=min(min((Fraction(s*a,q))%1,1-(Fraction(s*a,q))%1) for s in S)
            if m>best: best=m
    return best
print("AP-with-one-element-replaced tight sets  S=({1..n-1}\{k}) U {g}, M=1/n:")
for n in range(5,13):
    base=list(range(1,n)); found=[]
    for k in base:
        for g in range(n, 4*n+1):
            if g in base: continue
            S=[x for x in base if x!=k]+[g]
            if len(set(x%1 for x in S))<n-1: pass
            if Mfloat(S,5*n) <= 1.0/n + 1e-9:
                if Mexact(S,10*n)==Fraction(1,n):
                    found.append((k,g))
    # summarize the (k,g) pairs
    print(f"  n={n}: (remove k, add g) tight -> {found}")
print()
print("interpret: for each n, which (k,g)? Is g a MULTIPLE of k, or g = k + n*m (a lift mod n)?")
for n in range(5,13):
    base=list(range(1,n)); found=[]
    for k in base:
        for g in range(n, 4*n+1):
            if g in base: continue
            S=[x for x in base if x!=k]+[g]
            if Mfloat(S,5*n) <= 1.0/n + 1e-9 and Mexact(S,10*n)==Fraction(1,n):
                rel=[]
                if g%k==0: rel.append(f"g={g//k}*k")
                if (g-k)%n==0: rel.append(f"g=k+{(g-k)//n}*n")
                found.append((k,g,rel))
    if found: print(f"  n={n}: {found}")
