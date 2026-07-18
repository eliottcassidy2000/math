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
def comparison_tournament(R,q):
    # i->j iff (r_i - r_j) mod q in (0,q/2). scores = out-degrees. return score-sequence + c3
    n=len(R); scores=[0]*n
    for i in range(n):
        for j in range(n):
            if i==j: continue
            d=(R[i]-R[j])%q
            if 0<2*d<q: scores[i]+=1
    c3=comb(n,3)-sum(comb(s,2) for s in scores)
    return sorted(scores),c3
print("M<1/13 families: is a=val (t=M fixed point)? tournament of residues regular (R_12+1)?")
fams=[("deep well",list(range(1,13))+[182]),
      ("ladder m=2",list(range(1,13))+[364]),
      ("ladder m=3",list(range(1,13))+[546]),
      ("GW-ish core {1..11,13,182}",list(range(1,12))+[13,182])]  # non-{1..12} core test
for name,fam in fams:
    (M,(a,q))=M_arg(fam); val=int(round(float(M)*q))
    R=sorted((v*a)%q for v in fam)
    sc,c3=comparison_tournament(R,q)
    ap12 = (val>1 and sorted(r for r in R if r%val==0)==[val*i for i in range(1,13)])
    print(f"  {name}: M={float(M):.4f} a={a} val={val}  a=val? {a==val}  1 in V? {1 in fam}")
    print(f"     residue-tournament scores={sc}  c3={c3}  (regular R_12+1 => 12 scores equal + 1 outlier)")
    print(f"     core=val*{{1..12}}? {ap12}")

# Is a=val NECESSARY for M<1/13? check near-1/13 families (M slightly above) -- do they have a!=val?
print("\nNear-1/13 (M>=1/13): is a=val? (if a=val characterizes M<1/13, that's the hook)")
for name,fam in [("{1..11,13,84} M>1/13",list(range(1,12))+[13,84]),
                 ("{2,4..24,25}",[2,4,6,8,10,12,14,16,18,20,22,24,25]),
                 ("AP {1..13} tight",list(range(1,14)))]:
    (M,(a,q))=M_arg(fam); val=int(round(float(M)*q))
    print(f"  {name}: M={float(M):.4f} a={a} val={val} a=val? {a==val}")
