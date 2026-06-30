from fractions import Fraction
def Mc(S,c,Qmax):
    best=Fraction(0)
    for q in range(1,Qmax+1):
        for a in range(q):
            v=(Fraction(s*a,q)-c for s in S)
            m=min(min((x-x.__floor__()),1-(x-x.__floor__())) for x in [Fraction(s*a,q)-c for s in S])
            if m>best: best=m
    return best
# fast: L via envelope is analytic; just confirm (L-1/4)*n rises toward 1/2 with Qmax, for n=10
n=10; AP=list(range(1,n)); G=20
for Q in [2*n,4*n,8*n]:
    L=sum(float(Mc(AP,Fraction(i,2*G),Q)) for i in range(0,2*G))/(2*G)
    print(f"n={n} Qmax={Q}: (L-1/4)*n = {(L-0.25)*n:.4f}  (3/8=0.375, 1/2=0.500)")
print("envelope-integral L = 1/4 + 1/(2n): (L-1/4)*n = 1/2 EXACTLY (analytic).")
