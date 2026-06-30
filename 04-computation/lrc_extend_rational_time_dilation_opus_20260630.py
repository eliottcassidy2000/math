from fractions import Fraction
from math import gcd
import random
def frac(x): x=x%1; return min(x,1-x)
def M_exact(S, Qmax):
    best=Fraction(0); bt=None
    for q in range(1,Qmax+1):
        for a in range(q):
            m=min(frac(Fraction(s*a,q)) for s in S)
            if m>best: best=m; bt=Fraction(a,q)
    return best,bt
n=14
print("(A) RATIONAL-TIME LEMMA: if n does NOT divide any s_i, t=1/n gives min_i||s_i/n|| >= 1/n:")
random.seed(1)
for trial in range(4):
    S=sorted(random.sample([x for x in range(1,40) if x%n!=0], n-1))
    m=min(frac(Fraction(s,n)) for s in S)
    print(f"   S={S[:6]}...: min_i||s_i/14|| = {m} >= 1/14: {m>=Fraction(1,n)}")
print("   => PROVEN: M(S)>=1/n whenever n does not divide s_i for all i. LRC reduces to 'some s_i divisible by n'.")
print()
print("(B) DILATION M_c(dS)=M_c(S): homogeneous AP result holds for all d*{1..n-1}:")
for d in [2,3,5]:
    S=[d*k for k in range(1,n)]; m,t=M_exact(S, 3*n)
    print(f"   d={d}: M(d*AP)={m}=1/n? {m==Fraction(1,n)} (t={t})")
print()
print("(C) GENERAL AP speed sets {a,a+e,...} (n=6): M=1/n still?")
nn=6
for (a,e) in [(1,1),(2,3),(1,2),(3,5),(2,5)]:
    S=[a+k*e for k in range(nn-1)]
    if len(set(S))<nn-1: continue
    m,t=M_exact(S, 6*nn)
    print(f"   a={a},e={e}: S={S}: M={m}={float(m):.4f}, 1/n={1/nn:.4f}, M>=1/n:{m>=Fraction(1,nn)}, =1/n:{m==Fraction(1,nn)}")
print()
print("(D) HARD CASE (some s_i divisible by n), n=6: t=1/6 fails; is there a working t with M>=1/6?")
nn=6
for S in [[6,1,2,3,4],[6,2,3,4,5],[6,1,2,4,5],[12,1,2,3,5]]:
    m,t=M_exact(S, 12*nn)
    print(f"   S={S}: M={m}={float(m):.4f} >= 1/6={1/nn:.4f}? {m>=Fraction(1,nn)} (best t={t})")
print("   => divisible case still has M>=1/n but needs a non-1/n time = the hard core (the open LRC).")
