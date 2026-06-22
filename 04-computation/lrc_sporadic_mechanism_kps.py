from fractions import Fraction as F
# The witness-point structure: at the tight t, the 13 points {s t mod 1} must avoid (-1/14,1/14)
# around 0, with one point AT 1/14 (tightness). AP = equal-spacing; GW = ?
def points(S, t):
    return sorted(set((s*t) % 1 for s in S)), [ (s*t)%1 for s in S]
def show(name, S, t):
    distinct, allp = points(S, t)
    res = sorted((s, (F(s)*t*14)) for s in S)   # residue *14
    mult = {}
    for s in S:
        r = (s*t) % 1
        mult[r] = mult.get(r, 0) + 1
    occupied = sorted(int(r*14) for r in mult)   # which k/14 are hit
    missing = [k for k in range(1,14) if k not in occupied]
    collisions = sorted(int(r*14) for r,c in mult.items() if c>1)
    print(f"{name}: t={t}")
    print(f"  points (x14): {sorted(int((s*t%1)*14) for s in S)}")
    print(f"  occupied k/14: {occupied}")
    print(f"  MISSING k/14 : {missing}   COLLISIONS (k/14 with >1 speed): {collisions}")
    print(f"  min ||s t|| = {min(min(r,1-r) for r in allp)} = 1/14? {min(min(r,1-r) for r in allp)==F(1,14)}")
show("AP {1..13}", list(range(1,14)), F(3,14))
print()
show("GW {1..11,13,24}", list(range(1,12))+[13,24], F(5,14))
print("\n=> AP: all 13 residues 1..13 occupied, equally spaced, NO gap/collision (difference-closed).")
print("   GW: a GAP (missing residue) compensated by a COLLISION (two speeds same residue).")
print("   The sporadic mechanism trades equal-spacing for a balanced gap+collision that still")
print("   isolates 0 at exactly 1/14. Finiteness of such patterns = OPEN-Q-108 (sporadic residue).")
