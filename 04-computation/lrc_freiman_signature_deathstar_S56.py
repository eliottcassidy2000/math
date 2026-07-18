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
def residue_structure(fam):
    (M,(a,q))=M_arg(fam); val=M.numerator*(q//M.denominator) if False else int(M*q)
    # val = M*q (integer). residues r_i = v_i*a mod q, folded to [0,q); they lie in [val,q-val]
    val=M.numerator; qq=M.denominator  # M=val/q reduced; but need actual q from arg
    a,q=arg=(a,q); val=int(round(float(M)*q))
    R=sorted((v*a)%q for v in fam)
    # classes mod val
    classes={}
    for r in R: classes.setdefault(r%val,[]).append(r)
    return M,a,q,val,R,{c:len(v) for c,v in sorted(classes.items())}

fams={
 "deep well {1..12,182}": list(range(1,13))+[182],
 "ladder {1..12,364}": list(range(1,13))+[364],
 "ladder {1..12,546}": list(range(1,13))+[546],
 "3*deepwell {3,6..,546}": [3*v for v in range(1,13)]+[546],
}
print("M<1/13 families: residue classes mod val (expect 12 in one class + 1 in another = AP+perturbation):")
for name,fam in fams.items():
    M,a,q,val,R,cl=residue_structure(fam)
    print(f"  {name}: M={M}={float(M):.4f} q={q} val={val}  #classes mod val={len(cl)}  class-sizes={cl}")

print("\nNEAR-MISS (M slightly > 1/13): more classes mod val? (breaking the rigidity):")
near=[("{1..12,84}",list(range(1,13))+[84]),
      ("{1..11,13,84}",list(range(1,12))+[13,84]),
      ("{1..12,26}",list(range(1,13))+[26]),
      ("{2,4..24,25}",[2,4,6,8,10,12,14,16,18,20,22,24,25])]
for name,fam in near:
    M,a,q,val,R,cl=residue_structure(fam)
    print(f"  {name}: M={float(M):.4f} (>1/13? {M>F(1,13)}) val={val} #classes mod val={len(cl)} sizes={cl}")
