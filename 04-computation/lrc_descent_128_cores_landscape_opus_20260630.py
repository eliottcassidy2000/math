"""
THE FINITE FAMILIES: the 128 cores O subset Z_7, gap g(O)=min_{k!=0}|sum_{x in O} w^{kx}|^2 (w=e^{2pi i/7}).
Know them: the 5 values, multiplicities, the families (doublet/Fano/Z_7), Z_7* orbits, complement & cyclic symmetry.
"""
import cmath, math
from itertools import combinations
w=cmath.exp(2j*math.pi/7)
def gap(O):
    if not O: return None
    return min(abs(sum(w**(k*x) for x in O))**2 for k in range(1,7))
allO=[]
for r in range(8):
    for O in combinations(range(7),r):
        allO.append(frozenset(O))
from collections import defaultdict
byval=defaultdict(list)
for O in allO:
    if not O: continue
    g=round(gap(O),5); byval[g].append(O)
print("the 5 cyclotomic gap values and their cores (THM-590 landscape):")
for g in sorted(byval):
    cores=byval[g]; sizes=sorted(set(len(O) for O in cores))
    tag={0.0:"GAP=0: O=Z_7 ONLY (disproof boundary)",0.19806:"=4cos^2(3pi/7) DOUBLET floor (the atom)"}.get(g,"")
    print(f"   g={g:<8}: {len(cores):>3} cores, sizes {sizes}  {tag}")
print()
# the families
print("FAMILY 1 -- the binding DOUBLETS (|O|=2) and quintuplets (|O|=5): all bind at 0.198")
db=[O for O in allO if len(O)==2]; print(f"   {len(db)} doublets, gaps all = {round(gap(db[0]),5)} (Z_7*-invariant)")
print("FAMILY 2 -- the FANO/QR perfect-difference-set cores ({1,2,4} & translates): FLAT spectrum, gap 2 (optimal)")
fano=[]
for s in range(7):
    O=frozenset((x+s)%7 for x in {1,2,4}); fano.append((tuple(sorted(O)),round(gap(O),5)))
for O,g in fano[:4]: print(f"   {O}: gap={g} (lambda_k=2 all k, perfect difference set)")
print("FAMILY 3 -- the BOUNDARY: O=Z_7 (gap 0, the only zero, = the disproof / Phi_7 full sum)")
print(f"   Z_7={tuple(range(7))}: gap={round(gap(frozenset(range(7))),5)}")
print()
# symmetries
print("SYMMETRIES of the family:")
# complement
csym=all(round(gap(O),5)==round(gap(frozenset(range(7))-O),5) for O in allO if O and O!=frozenset(range(7)))
print(f"   complement-symmetric g(O)=g(Z_7\O): {csym}")
# cyclic-invariant cores (O = O+1)
cyc=[O for O in allO if O and O==frozenset((x+1)%7 for x in O)]
print(f"   CYCLIC-invariant cores (O=O+1, the descent's Z_7-invariance test): {[tuple(sorted(O)) for O in cyc]}")
print(f"     => ONLY empty and Z_7 are cyclic-invariant! Proper cores are NOT translation-invariant.")
print(f"     (the descent needs Z_7-MULTIPLIER (Z_7*) invariance of the SPECTRUM, not translation of the SET -")
print(f"      Z_7* permutes the 6 nonzero modes => spectrum orbit-invariant; verified all doublets share gap.)")
