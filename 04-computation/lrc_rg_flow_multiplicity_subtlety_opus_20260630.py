"""Check the multiplicity subtlety: random A's g=0.061 is from a MULTISET odd core (repeated residues mod 7),
outside the THM-590 distinct-subset gap list {0,0.198,0.308,1,2}. Clean AP cores are distinct-residue."""
import cmath, math, random
def gmult(res,p=7):
    w=cmath.exp(2j*cmath.pi/p)
    return min(abs(sum(w**((k*x)%p) for x in res))**2 for k in range(1,p))
random.seed(1)
A=sorted(random.sample(range(1,50),13))
O0=[x for x in A if x%2==1]
res=[x%7 for x in O0]
from collections import Counter
print(f"random A = {A}")
print(f"  level-0 odd core = {O0}, residues mod7 = {res}")
print(f"  residue multiset = {dict(Counter(res))}  (distinct? {len(set(res))==len(res)})")
print(f"  g(multiset) = {gmult(res):.4f}  -- below 0.198 BECAUSE of multiplicity (repeated residues)")
print(f"  g(distinct set) = {gmult(sorted(set(res))):.4f}  -- the THM-590 value for the underlying SET")
print()
print("HONEST SUBTLETY: THM-590's floor 4cos^2(3pi/7)=0.198 is for DISTINCT-residue cores (subsets of Z_7).")
print("  When several speeds share a residue mod 7 (a MULTISET core), |sum omega^kx|^2 can dip below 0.198.")
print("  The AP and the covering constructions have DISTINCT-residue cores at every level (clean flow),")
print("  so the doublet floor applies there. General multiplicity is a separate (looser-bound) regime.")
print()
# the AP cores ARE distinct at every level (verify for p=7)
ch=[];cur=list(range(1,14))
while cur:
    O=[x for x in cur if x%2==1]
    if O: ch.append(O)
    cur=[x//2 for x in cur if x%2==0]
print("AP {1..13} cores distinct-residue mod 7 at every level:")
for j,O in enumerate(ch):
    r=[x%7 for x in O]; print(f"   level {j}: {O} -> mod7 {r} distinct={len(set(r))==len(r)}")
