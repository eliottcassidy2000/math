"""
Explore the UPPER bound min_v ||vt-c|| <= env for ALL t. Generate/test hypotheses.
 H1: the bound holds for all t (sup = env).
 H2: at the maximizing t, the min is achieved at an EXTREME runner v in {1, n-1}.
 H3: runners {vt} are symmetric about nt/2 (v <-> n-v); the big gap containing c reflects to a big gap at nt-c.
 H4: the # of 'big gaps' (>= 2*env) is constrained; c-gap + its reflection forces structure.
"""
from fractions import Fraction
def frac(x): x=x%1; return min(x,1-x)
n=9; AP=list(range(1,n))
print(f"n={n}. For each c, scan t=a/q (q<=Qmax), find max_t min_v||vt-c|| and the argmin runner:")
Qmax=6*n
for cn,cd in [(0,1),(1,18),(1,9),(1,6),(2,9),(1,4),(7,18),(4,9)]:
    c=Fraction(cn,cd)
    best=Fraction(-1); bt=None; bv=None
    for q in range(1,Qmax+1):
        for a in range(q):
            t=Fraction(a,q)
            dv=[(frac(Fraction(v)*t-c),v) for v in AP]
            m,vmin=min(dv)
            if m>best: best=m; bt=t; bv=vmin
    env=Fraction(1,n)+c*Fraction(n-2,n) if c<=Fraction(1,2) else Fraction(1,n)+(1-c)*Fraction(n-2,n)
    print(f"  c={str(c):>5}: M_c={str(best):>7}={float(best):.4f} env={float(env):.4f} {'OK' if best<=env else 'VIOLATION!!'}; argmin runner v={bv} (extreme={bv in (1,n-1)}); t*={bt}")
print()
print("H2 test: at the MAX, is argmin always an extreme runner v in {1,n-1}? (above)")
print()
# H3: symmetry about nt/2
print("H3: runners {vt} symmetric about nt/2? check for a random t:")
t=Fraction(3,7); pts=sorted(frac(Fraction(v)*t)*0+ (Fraction(v)*t)%1 for v in AP)
ctr=(Fraction(n)*t/2)%1
refl=sorted((2*ctr - p)%1 for p in pts)
print(f"   n={n}, t={t}: center nt/2={ctr}; runners symmetric about it: {pts==refl}")
print()
# H1 broader: random irrational-ish t, verify <= env
import random; random.seed(0); viol=0
for _ in range(3000):
    q=random.randint(2,200); a=random.randint(1,q-1); t=Fraction(a,q)
    c=Fraction(random.randint(0,n),2*n)
    m=min(frac(Fraction(v)*t-c) for v in AP)
    env=Fraction(1,n)+c*Fraction(n-2,n)
    if m>env+Fraction(1,10**9): viol+=1
print(f"H1: random (t,c) pairs with min_v||vt-c|| > env: {viol}/3000 (0 = bound holds)")
