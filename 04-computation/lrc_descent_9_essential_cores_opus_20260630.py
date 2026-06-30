"""
The ESSENTIAL finite families: cores up to the affine group AGL(1,7)=Z_7 |x| Z_7* (order 42).
gap is AGL-invariant. Find the orbits, their gaps, sizes. Connect the binding DOUBLET to my binding pairs.
"""
import cmath, math
from itertools import combinations
w=cmath.exp(2j*math.pi/7)
def gap(O):
    if not O: return None
    return round(min(abs(sum(w**(k*x) for x in O))**2 for k in range(1,7)),5)
def affine_orbit(O):
    orb=set()
    for u in [1,2,3,4,5,6]:
        for s in range(7):
            orb.add(frozenset((u*x+s)%7 for x in O))
    return orb
allO=[frozenset(O) for r in range(8) for O in combinations(range(7),r)]
seen=set(); orbits=[]
for O in allO:
    if O in seen or not O: continue
    orb=affine_orbit(O); seen|=orb
    orbits.append((len(next(iter([O]))), gap(O), len(orb), tuple(sorted(O))))
print("AGL(1,7)-orbits of cores (the ESSENTIAL finite families), by gap:")
print(f"{'|O|':>4} {'gap g(O)':>10} {'orbit size':>11} {'representative':>16}  role")
for sz,g,osz,rep in sorted(orbits, key=lambda t:(t[1],t[0])):
    role=""
    if g==0.0: role="BOUNDARY (disproof: O=Z_7)"
    elif g==0.19806: role="<<< THE FLOOR (doublet/quintuplet = binding pair)"
    elif g==2.0: role="OPTIMAL (Fano/perfect difference set = octonion)"
    print(f"{sz:>4} {g:>10} {osz:>11} {str(rep):>16}  {role}")
print(f"\n  => only {len(orbits)} essential cores (AGL-orbits). The whole descent is a {len(orbits)}-row table.")
print("  THE FLOOR is the doublet orbit (gap 0.198=4cos^2(3pi/7)); the BOUNDARY is the single orbit O=Z_7 (gap 0).")
print()
# my binding pairs mod 14 descend to doublets mod 7
print("my binding pairs mod 14 -> doublets mod 7 (all bind at the floor 0.198):")
for pair in [(1,13),(5,9),(3,11)]:
    O=frozenset(x%7 for x in pair); print(f"   {{{pair[0]},{pair[1]}}} mod14 -> {tuple(sorted(O))} mod7, gap={gap(O)}, diff={ (max(O)-min(O))}")
print("\n  => the descent's binding core IS my binding pair = THM-578's doublet. The floor binds at the same object.")
print("  THE DESCENT WALL (precise): show the descended core O_j is never the FULL Z_7 (the only gap-0 orbit).")
print("  Everything else (127 of 128 cores) already clears 0.198 by THM-590 (finite, PROVED).")
