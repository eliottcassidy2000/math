"""
CRITICAL CHECK: is M_c(AP_n) = 1/n + c*(n-2)/n EXACTLY (the linear envelope through the spectrum points)?
If so: L = 1/4 + 1/(2n) (coeff 1/2, NOT 3/8); spread regime = ALL of [0,1/2); binding c*=1/2 (antipode).
The earlier 3/8 and c*=1/3 would be Qmax-under-estimation ARTIFACTS. Test with growing Qmax.
"""
from fractions import Fraction
def Mc(S,c,Qmax):
    best=Fraction(0); bw=None
    for q in range(1,Qmax+1):
        for a in range(q):
            m=min(min((Fraction(s*a,q)-c)-((Fraction(s*a,q)-c).__floor__()), 1-((Fraction(s*a,q)-c)-((Fraction(s*a,q)-c).__floor__()))) for s in S)
            if m>best: best=m; bw=(a,q)
    return best,bw
def env(n,c): return Fraction(1,n)+c*Fraction(n-2,n)
n=14; AP=list(range(1,n))
print(f"n={n}: M_c vs envelope 1/n+c(n-2)/n, with GROWING Qmax (does M_c rise to the envelope?):")
print(f"   {'c':>6} {'env':>8} {'M(Q=2n)':>9} {'M(Q=5n)':>9} {'M(Q=12n)':>10} {'=env?':>6}")
for cnum,cden in [(1,28),(1,14),(1,7),(3,28),(1,4),(2,7),(3,7),(9,28)]:
    c=Fraction(cnum,cden); e=env(n,c)
    m2,_=Mc(AP,c,2*n); m5,_=Mc(AP,c,5*n); m12,b=Mc(AP,c,12*n)
    print(f"   {str(c):>6} {str(e):>8} {str(m2):>9} {str(m5):>9} {str(m12):>10} {str(m12==e):>6}")
print()
# the 'trivial' c=0.4-region: does big Qmax show M_c > c (spread), refuting the trivial regime?
print("the supposed 'trivial regime' c>1/3: does large Qmax show M_c = envelope > c (spread, not trivial)?")
for cnum,cden in [(2,5),(3,7),(11,28),(0,1)]:
    if cden==1: continue
    c=Fraction(cnum,cden); e=env(n,c)
    m3,_=Mc(AP,c,3*n); m12,_=Mc(AP,c,12*n)
    print(f"   c={str(c):>5}: env={str(e)}={float(e):.4f}; M(Q=3n)={float(m3):.4f}; M(Q=12n)={float(m12):.4f}; >c? {m12>c}")
print()
print("If M(Q=12n)=env at all c: the LINEAR envelope is EXACT, L=1/4+1/(2n) (coeff 1/2), spread=[0,1/2).")
print("=> the earlier 3/8 (Qmax=2n) and c*=1/3 (Qmax=3n) were UNDER-ESTIMATION ARTIFACTS. Honest correction.")
