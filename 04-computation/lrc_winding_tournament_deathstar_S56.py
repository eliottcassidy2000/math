from fractions import Fraction as F
from math import gcd, comb
def M_arg(fam):
    Q=2*max(fam)+2; best=F(0); arg=None
    for q in range(2,Q+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            if F(m,q)>best: best=F(m,q); arg=(a,q)
    return best,arg
def winding_tournament(fam, a, q):
    # observer speed 0 included: 14 vertices. i->j iff frac((s_i - s_j) a/q) in (0,1/2)
    S=[0]+list(fam); n=len(S)
    # positions numerator mod q: (s * a) % q ; frac in (0,1/2) <=> residue in (0, q/2)
    scores=[0]*n
    c3check=0
    for i in range(n):
        for j in range(n):
            if i==j: continue
            r=((S[i]-S[j])*a)%q
            # frac in (0,1/2): 0 < r < q/2  (handle r=0 tie: skip/orient by index)
            if 0<2*r<q:
                scores[i]+=1
    return scores,n
def analyze(name,fam):
    fam=sorted(fam); M,(a,q)=M_arg(fam)
    scores,n=winding_tournament(fam,a,q)
    c3=comb(n,3)-sum(comb(s,2) for s in scores)
    sc=sorted(scores)
    regular = (max(scores)-min(scores)<=1)   # near-regular
    print(f"{name}: M={float(M):.4f}={M} @t={a}/{q}  n={n}")
    print(f"   scores(sorted)={sc}")
    print(f"   c3(3-cycles)={c3}  (regular C(n,3)*.. ) near-regular={regular}  score-range={max(scores)-min(scores)}")
analyze("AP {1..13}", list(range(1,14)))
analyze("GW", [1,2,3,4,5,6,7,8,9,10,11,13,24])
analyze("deep well {1..12,182}", list(range(1,13))+[182])
analyze("covering min", [2,5,6,8,11,13,14,16,19,20,24,25,27])
analyze("covering {2..14}", list(range(2,15)))
