from fractions import Fraction as F
def nrm(x): f=x-(x.numerator//x.denominator); return min(f,1-f)
def M_exact(S):
    C=set()
    for i in range(len(S)):
        for j in range(len(S)):
            for d in (S[i]+S[j],abs(S[i]-S[j])):
                if d:
                    for m in range(d+1): C.add(F(m,d))
        for m in range(2*S[i]+1): C.add(F(m,2*S[i]))
    best=F(0);arg=F(0)
    for t in C:
        if 0<t<1:
            v=min(nrm(s*t) for s in S)
            if v>best: best,arg=v,t
    return best,arg
def is_cov(S,n): return all(any(s%q==0 for s in S) for q in range(2,n+1))
print("FIRST-PYRAMID construction {1,..,n-2,(n-1)n} vs M_min formula n/(n^2-n+1):")
print(f"{'n':>3} {'killer':>7} {'M':>12} {'n/(n^2-n+1)':>13} {'match':>6} {'witness':>10}")
for n in range(4,15):
    S=list(range(1,n-1))+[(n-1)*n]
    M,t=M_exact(S); formula=F(n,n*n-n+1)
    print(f"{n:>3} {(n-1)*n:>7} {str(M):>12} {str(formula):>13} {str(M==formula):>6} {str(t):>10}")
print(f"\nSECOND pyramid (squares): row k center 2k(k+1). k=13 -> 2*182=364=2*killer. start k(2k+1).")
print("first pyramid (linear) row k split = pronic k(k+1); second (squares) center = 2k(k+1).")
# the two pyramids: linear split = killer; squares center = 2*killer
for k in [13]:
    print(f"  k={k}=n-1: 1st-pyramid split (pronic) = {k*(k+1)} = killer; 2nd-pyramid center = {2*k*(k+1)}")
# does the second pyramid (squares identity) give a DUAL covering construction? squares -> different killer?
print("\ncheck: is {1..n-2, (n-1)n} the global min? compare to nearby covering constructions at n=14:")
ap=list(range(1,15))
for desc,S in [("drop unit 13, killer 182 (1st-pyramid)",list(range(1,13))+[182]),
               ("drop 13, killer 364=2*182",list(range(1,13))+[364]),
               ("drop 11&13, killers 154,182",list(range(1,11))+[12,154,182]),
               ("drop unit 1, killer 14",list(range(2,13))+[13,14])]:
    if len(set(S))==13 and is_cov(S,14):
        M,t=M_exact(S); print(f"   {desc:>40}: M={M}={float(M):.6f}")
