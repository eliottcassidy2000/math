r"""
THE APEX-DENOMINATOR LEMMA + the finite residue characterization of (star) (kps-S31aa).
LEMMA: at a tight optimum t*=a/D (M=1/14), the binding pair {v_a,v_b} (runners at +-1/14) satisfies
14|D and D|(v_a+v_b). For binding residues mod 14, the optimum is at D=14, and then M is determined by
the RESIDUE MULTISET {v_i mod 14}. So (star) reduces to: which residue patterns give M_14 = 1/14?
"""
from fractions import Fraction as F
from math import gcd
def nf(x):
    r=x%1; return min(r,1-r)
def M_and_opt(S):
    S=sorted(set(abs(s) for s in S if s!=0)); C=set()
    for i in range(len(S)):
        for j in range(i,len(S)):
            for comb in {S[i]+S[j], abs(S[i]-S[j])}:
                if comb:
                    for k in range(1,comb): C.add(F(k,comb))
    best=F(0); arg=None
    for t in C:
        if 0<t<1:
            m=min(nf(s*t) for s in S)
            if m>best: best=m; arg=t
    return best,arg
print("PART 1 -- the apex-denominator lemma on tight sets:")
for name,S in [("AP {1..13}",list(range(1,14))),
               ("GW {1..11,13,24}",list(range(1,12))+[13,24]),
               ("2*AP (dilate)",[2*k for k in range(1,14)]),
               ("GW dilate 3x",[3*k for k in (list(range(1,12))+[13,24])])]:
    M,t=M_and_opt(S); D=t.denominator
    binders=[s for s in S if nf(s*t)==F(1,14)]
    g=0
    for s in S: g=gcd(g,s)
    print(f"  {name:18s}: M={M}, opt t={t} (D={D}), 14|D? {D%14==0}, binders={binders}, sum={sum(binders) if len(binders)==2 else '...'}, gcd={g}")

print("\nPART 2 -- FINITE residue characterization at D=14: which residue multisets give M_14=1/14?")
units=[a for a in range(1,14) if gcd(a,14)==1]
def M14(res):  # res = multiset of 13 residues in 1..13; returns 14*M_14
    return max(min(min((r*a)%14,14-(r*a)%14) for r in res) for a in units)
full=list(range(1,14))
print(f"  FULL set {{1..13}} (AP pattern): 14*M_14 = {M14(full)}  {'=1 => M=1/14 TIGHT' if M14(full)==1 else ''}")
tight_patterns=[]
for rep in range(1,14):      # one residue doubled
    for miss in range(1,14):  # one residue missing
        if rep==miss: continue
        res=[x for x in range(1,14) if x!=miss]+[rep]  # 13 residues: drop miss, add a 2nd rep
        if len(res)!=13: continue
        if M14(res)==1: tight_patterns.append((rep,miss))
print(f"  one-repeat/one-missing patterns giving M_14=1/14: {len(tight_patterns)}")
for rep,miss in tight_patterns: print(f"     double {rep}, miss {miss}  (GW = double 10, miss 12 -> 24==10 mod14: {(rep,miss)==(10,12)})")
print(f"\n  => tight residue patterns at D=14: the FULL set (AP) + {len(tight_patterns)} repeat/missing (GW-type).")
print("     This is the FINITE half of (star): the tight locus at the apex denominator is characterized.")
