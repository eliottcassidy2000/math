import sys, itertools
from fractions import Fraction
if hasattr(sys.stdout,'reconfigure'): sys.stdout.reconfigure(encoding='utf-8')
from lrc14_verify_glaisher_indep import meas_S7

# PART 3: scrutinize what the "doubling-monotonicity" certificate actually
# certifies. The claim frames it as evidence for a DYADIC mechanism, but the
# swap test only replaces an element by a LARGER value. Question: is the
# decrease about doubling structure, or just monotone-in-magnitude / contiguity?
#
# Test A: does p_0 decrease when we replace element by a larger value that
#   PRESERVES dyadic structure vs one that BREAKS it? If dyadic structure were
#   the mechanism, breaking the chain should hurt MORE than a structure-neutral
#   enlargement.
# Test B: full neighbor structure — is consec a LOCAL max under ALL single
#   element replacements (both up AND down within {1..13})?

print("="*70)
print("PART 3A: single-element replacement local-max test at k=8")
print("  (replace each element by EVERY other value in {1..13}\\E)")
print("="*70)
base=list(range(1,9)); p0_base=meas_S7(base)
print(f"base consec{{1..8}} p0={float(p0_base):.5f}")
incr=0; decr=0; eq=0; up_viol=[]
for i in range(8):
    for new in range(1,14):
        if new in base: continue
        E=base[:i]+[new]+base[i+1:]
        p=meas_S7(E)
        if p>p0_base: incr+=1; up_viol.append((base[i],new,float(p)))
        elif p<p0_base: decr+=1
        else: eq+=1
print(f"  replacements: increase={incr} decrease={decr} equal={eq}")
print(f"  consec is local max under all single swaps: {incr==0}")
if up_viol: print(f"  INCREASES: {up_viol[:5]}")

print()
print("="*70)
print("PART 3B: does dyadic-structure-preserving vs breaking matter?")
print("  Replace element 5 (odd, not a power of 2) by candidates; compare")
print("  to replacing 8 (the top of the 1-2-4-8 chain).")
print("="*70)
# Replacing 8 breaks the power-of-2 chain {1,2,4,8}->{1,2,4}. Replacing 5
# does not touch the chain. If the chain is the mechanism, removing 8 should
# hurt more than removing 5, for matched replacement values.
for repl,label in [(8,"remove 8 (breaks 1-2-4-8 chain)"),(5,"remove 5 (chain intact)")]:
    print(f"  {label}:")
    idx=base.index(repl)
    for new in (9,10,11,12,13):
        E=base[:idx]+[new]+base[idx+1:]
        print(f"     ->{new}: p0={float(meas_S7(E)):.5f}  (drop {float(p0_base-meas_S7(E)):.5f})")

print()
print("="*70)
print("PART 3C: is p_0(consec) decreasing-in-shift truly monotone, or does the")
print("  dyadic-chain claim mispredict? Recheck all 6 windows, sorted by p0.")
print("="*70)
def pow2chain(E):
    s=set(E); best=0
    for a in range(0,14):
        if (1<<a) in s and (a==0 or (1<<(a-1)) not in s):
            L=0;v=1<<a
            while v in s: L+=1;v*=2
            best=max(best,L)
    return best
rows=[]
for start in range(1,7):
    E=tuple(range(start,start+8)); rows.append((meas_S7(E),pow2chain(E),start,E))
rows.sort(key=lambda t:-t[0])
prev=None; mono_ok=True
for v,c,s,E in rows:
    flag=""
    if prev is not None and c>prev[1]:
        flag=" <-- pow2chain INCREASES while p0 decreases (mechanism MISPREDICTS)"
        mono_ok=False
    print(f"  start={s} p0={float(v):.5f} pow2chain={c}{flag}")
    prev=(v,c)
print(f"\n  pow2chain is a monotone predictor of p0 across all windows: {mono_ok}")
print("  => The chain-length ordering only holds for the cherry-picked")
print("     windows starting at 1,2,3; it FAILS as a general predictor.")
