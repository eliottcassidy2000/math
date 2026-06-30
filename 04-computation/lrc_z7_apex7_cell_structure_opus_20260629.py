"""
THE Z_7 / apex-7 SECONDARY obstruction for LRC(14).
Key: at the 7-scale t=a/7, ||s_i a/7|| < 1/14  <=>  7 | s_i a  <=> (a!=0) 7|s_i.
So the 6 points {a/7: a=1..6} (the Z_7-orbit) are BLOCKED exactly by the mult-of-7 runners
(covering sets HAVE them: 7 and 14 are in {2..14}). The lonely set lives in the 7 open cells
[a/7,(a+1)/7); the Z_7 rotation t->t+1/7 permutes the cells (shifting runner i by s_i/7).
Q: is the per-cell lonely structure Z_7-forced? Compare to the Z_2 mirror (which gave even, no force).
"""
from fractions import Fraction as F
def nrm(x): f=x-(x.numerator//x.denominator); return min(f,1-f)
def lonely_intervals(S,n=14,lo=F(0),hi=F(1)):
    B=set([lo,hi])
    for s in S:
        k=0
        while F(n*k-1,n*s) < hi+1 and k < n*int(hi)+s+5:
            for t in (F(n*k+1,n*s),F(n*k-1,n*s)):
                if lo<t<hi: B.add(t)
            k+=1
    B=sorted(B); ivs=[]
    for i in range(len(B)-1):
        mid=(B[i]+B[i+1])/2
        if min(nrm(s*mid) for s in S)>=F(1,n): ivs.append((B[i],B[i+1]))
    return ivs
def is_cov(S): return all(any(s%q==0 for s in S) for q in range(2,15))
def per_cell(S):
    iv=lonely_intervals(S); cells=[0]*7
    for (a,b) in iv:
        mid=(a+b)/2; c=int(mid*7); 
        if c>6: c=6
        cells[c]+=1
    return cells, len(iv)
import random; random.seed(11)
sets={"{2..14}":list(range(2,15)),
      "{1..11,13,84}":[1,2,3,4,5,6,7,8,9,10,11,13,84],
      "{1..12,182}":[1,2,3,4,5,6,7,8,9,10,11,12,182]}
cc=0
while len(sets)<10 and cc<50000:
    cc+=1; S=tuple(sorted(random.sample(range(1,30),13)))
    if is_cov(S): sets["r"+str(S[:2])]=list(S)
print("lonely intervals per 1/7-cell [a/7,(a+1)/7), a=0..6 (the Z_7 structure):")
print(f"{'set':28s} {'cell0':>5}{'cell1':>6}{'cell2':>6}{'cell3':>6}{'cell4':>6}{'cell5':>6}{'cell6':>6}  total")
for nm,S in sets.items():
    cells,tot=per_cell(S)
    mult7=[s for s in S if s%7==0]
    print(f"{nm:28s} " + "".join(f"{c:>6}" for c in cells) + f"   {tot}  (mult-of-7 in set: {mult7})")
print()
print("Z_2 mirror: cell a <-> cell 6-a (t->1-t maps [a/7,(a+1)/7) to (6-a)/7 cell). Check symmetry:")
for nm,S in list(sets.items())[:4]:
    cells,_=per_cell(S)
    sym = all(cells[a]==cells[6-a] for a in range(7))
    print(f"  {nm:28s}: cells={cells}  mirror-symmetric(a<->6-a): {sym}")
print("\nb1^- (metagraph R-odd Betti, the SECONDARY obstruction) = 1, 7, 119 for n=4,5,6; factor:")
for n,b in [(4,1),(5,7),(6,119)]:
    fac=[]; x=b; d=2
    while x>1:
        while x%d==0: fac.append(d); x//=d
        d+=1
    print(f"   n={n}: b1^-={b} = {'*'.join(map(str,fac)) if fac else '1'}   apex-7 divides: {b%7==0}")
