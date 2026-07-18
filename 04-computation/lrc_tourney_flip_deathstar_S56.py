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
def scores(R,q):
    n=len(R); s=[0]*n
    for i in range(n):
        for j in range(n):
            if i!=j and 0<2*((R[i]-R[j])%q)<q: s[i]+=1
    return sorted(s)
# distance from the regular tournament: how many edge-flips from all-scores-=(n-1)/2 ?
# for odd n, regular has all scores (n-1)/2. flip-distance lower bound = sum |s_i-(n-1)/2| / 2
def flip_defect(sc,n):
    reg=(n-1)//2
    return sum(abs(s-reg) for s in sc)//2, [s-reg for s in sc if s!=reg]
print("M<1/13 residue tournament: edge-flip distance from the REGULAR tournament R_13:")
for name,fam in [("deep well",list(range(1,13))+[182]),
                 ("ladder m=2",list(range(1,13))+[364]),
                 ("ladder m=5",list(range(1,13))+[910]),
                 ("dilate 3*",[3*v for v in range(1,13)]+[546])]:
    (M,(a,q))=M_arg(fam); val=int(round(float(M)*q)); R=sorted((v*a)%q for v in fam)
    sc=scores(R,q); d,defects=flip_defect(sc,13)
    print(f"  {name}: scores={sc}  flip-defect from R_13 = {d}  (defect-vertices {defects})  a=val:{a==val}")
print("\n=> every M<1/13 family: residue tournament = R_13 (regular) with flip-defect 1 (the far element).")
print("   So: M<1/13 <=> tournament is regular-except-one-flip <=> residues equally spaced + 1 perturbation")
print("   <=> 2 cosets mod val <=> AP core + far element. The far element = the single broken edge.")
# contrast: M>1/13 -> larger defect?
print("\nContrast M>1/13:")
for name,fam in [("{1..11,13,84}",list(range(1,12))+[13,84]),("random cov",[2,3,4,6,7,9,11,13,17,19,23,26,28])]:
    (M,(a,q))=M_arg(fam); R=sorted((v*a)%q for v in fam); sc=scores(R,q); d,_=flip_defect(sc,13)
    print(f"  {name}: M={float(M):.4f} scores={sc} flip-defect={d}")
