from fractions import Fraction as F
LO=F(1,14); HI=F(3,41)
print("M = D/s lies in (1/14, 3/41)  <=>  41D/3 < s < 14D")
print()
print("  D    integer range for s        admissible s -> M")
cands=[]
for D in range(1,16):
    lo=F(41*D,3); hi=14*D
    ss=[s for s in range(int(lo)+1, hi) if LO<F(D,s)<HI]
    cands.extend((D,s) for s in ss)
    rng=f"({float(lo):7.2f}, {hi:3d})"
    tag=", ".join(f"{s} -> {F(D,s)}" for s in ss) if ss else "NONE"
    print(f"  {D:2d}   {rng}   {tag}")
print()
print(f"  D = 1,2,3 give NO integer s  =>  the interval FORCES D >= 4.")
print(f"  (boxeph-S123's stratification gave D >= 3 for the wider Farey interval;")
print(f"   the attained interval (1/14,3/41) pushes it one step further.)")
print()
print(f"  total candidate (D,s) pairs with D <= 15: {len(cands)}")
print(f"  first few: {cands[:8]}")
print()
print("  Note the list is INFINITE (D unbounded), since s ~ 14D always admits")
print("  an integer once D >= 4.  So emptiness cannot be settled by finite check")
print("  alone -- it needs an argument that bounds D for a realising family.")
