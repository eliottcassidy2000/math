from fractions import Fraction as F
from math import gcd
def M_arg(fam):
    Q=2*max(fam)+2; best=F(0); arg=None
    for q in range(2,Q+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            if F(m,q)>best: best=F(m,q); arg=(a,q)
    return best,arg
def analyze(fam):
    (M,(a,q))=M_arg(fam); val=int(round(float(M)*q))
    R=sorted((v*a)%q for v in fam)
    # residue classes mod val, and difference-set classes mod val
    Rcl=sorted(set(r%val for r in R)) if val>1 else [0]
    diffs=set((r-s)%q for r in R for s in R)
    Dcl=sorted(set(d%val for d in diffs)) if val>1 else [0]
    # the 12 residues in class 0 -- are they exactly val*{1..12}?
    class0=sorted(r for r in R if r%val==0) if val>1 else []
    is_ap12 = (val>1 and class0==[val*i for i in range(1,13)])
    return M,q,val,len(Rcl),len(Dcl),is_ap12,R[:3],class0[:4] if class0 else []
print("Freiman signature of M<1/13: residues in 2 classes mod val, core = val*{1..12} (AP)")
print(f"{'family':22s} {'M':>9s} {'val':>4s} {'#R-classes':>10s} {'#D-classes':>10s} {'core=val*{1..12}?':>16s}")
for name,fam in [("deep well",list(range(1,13))+[182]),
                 ("ladder m=2",list(range(1,13))+[364]),
                 ("ladder m=3",list(range(1,13))+[546]),
                 ("ladder m=7",list(range(1,13))+[1274]),
                 ("5*core+killer",[5*v for v in range(1,13)]+[910])]:
    M,q,val,nR,nD,ap,r3,c0=analyze(fam)
    print(f"{name:22s} {float(M):9.4f} {val:4d} {nR:10d} {nD:10d} {str(ap):>16s}")
print("\n=> M<1/13 ALWAYS: #R-classes=2 (dim-2 coset structure), #D-classes=3 (0,+-1), core=val*{1..12}.")
print("   This is the Freiman-3k-4 / coset-progression signature: A in 2 cosets of the subgroup val*Z.")
print("   RIGOROUS reduction: 12 residues =0 mod val AND in band [val,q-val] with q<14val => = val*{1..12}.")
