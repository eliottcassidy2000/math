import sys, itertools
from fractions import Fraction
if hasattr(sys.stdout,'reconfigure'): sys.stdout.reconfigure(encoding='utf-8')
from lrc14_verify_glaisher_indep import meas_S7, odd_part, maxrun, richC2

# PART 2: scrutinize the window chain-length interpretation and the b=1 claim,
# plus stress-test the "consec = global max" claim at k=9,10 (the OTHER binding
# cluster sizes) and check whether 'start at 1' really beats all windows.

def pow2_chain_in(E):
    """longest run of consecutive powers of two {2^a, 2^(a+1), ...} all in E."""
    s = set(E); best = 0
    for a0 in range(0, 14):
        if (1<<a0) not in s: continue
        # only count as a chain-start if 2^(a0-1) not in s
        if (a0>0) and ((1<<(a0-1)) in s): continue
        L=0; v=1<<a0
        while v in s and v<=13:
            L+=1; v*=2
        best=max(best,L)
    return best

print("="*70)
print("PART 2A: window chain length interpretation (windows of size 8)")
print("="*70)
for start in range(1,7):
    E=tuple(range(start,start+8))
    p2=[e for e in E if (e&(e-1))==0]  # powers of 2 in E
    print(f"  {E}: powers-of-2={p2} pow2chain={pow2_chain_in(E)} p0={float(meas_S7(E)):.5f}")

print()
print("="*70)
print("PART 2B: ALL contiguous windows of size 8 ordered by p0")
print("="*70)
wins=[]
for start in range(1,7):
    E=tuple(range(start,start+8))
    wins.append((meas_S7(E),pow2_chain_in(E),E))
wins.sort(key=lambda t:-t[0])
for v,c,E in wins:
    print(f"  p0={float(v):.5f}  pow2chain={c}  {E}")
mono = all(wins[i][1]>=wins[i+1][1] for i in range(len(wins)-1))
print(f"  p0-order matches pow2chain-order (weakly decreasing): {mono}")

print()
print("="*70)
print("PART 2C: consec global-max at k=9 and k=10 (exhaustive)")
print("="*70)
for k in (9,10):
    allk=list(itertools.combinations(range(1,14),k))
    vals=[(meas_S7(E),E) for E in allk]
    vals.sort(key=lambda t:-t[0])
    consec=tuple(range(1,k+1))
    pc=meas_S7(consec)
    rank=[i for i,(v,E) in enumerate(vals) if E==consec][0]+1
    n_at_max=sum(1 for v,E in vals if v==vals[0][0])
    print(f"  k={k}: #sets={len(vals)} consec p0={float(pc):.5f} rank={rank}/{len(vals)} #at-max={n_at_max} top={vals[0][1]}")

print()
print("="*70)
print("PART 2D: counterexample hunt — does ANY non-consec set beat consec at k=8,9,10?")
print("         (i.e. is consec strictly the unique max at every binding k?)")
print("="*70)
for k in (8,9,10):
    allk=list(itertools.combinations(range(1,14),k))
    consec=tuple(range(1,k+1))
    pc=meas_S7(consec)
    beats=[E for E in allk if E!=consec and meas_S7(E)>=pc]
    print(f"  k={k}: sets with p0 >= consec (excl consec): {len(beats)}  -> {'CONSEC UNIQUE MAX' if not beats else beats[:3]}")
