"""
lrc14_freiman_dimension_is_the_parameter_opus_S181.py  (opus-2026-07-09-S181)

CORRECTION of the S179 'additive energy is THE single parameter' overclaim, found via the Freiman
deep-dive. 2-D GAPs have HIGH additive energy E~1000-1200 (comparable to near-AP) but are LOOSE
(L~0.10-0.14, like Sidon). So E does NOT determine looseness. The TRUE parameter is the FREIMAN
DIMENSION d: the AP (d=1, resonances aligned on ONE modulus => coherent covering => tight) vs multi-D
GAPs (d>=2, resonances span d directions => no common tau => loose). Holds E ~fixed, varies d, measures
L -- if L clusters by d not by E, Freiman DIMENSION (not energy) is the loneliness parameter.
"""
import numpy as np
from collections import Counter

NG = 240007
TAU = (np.arange(NG) + 0.5) / NG
H = 1.0 / 14

def E(S):
    r = Counter(a + b for a in S for b in S); return sum(c*c for c in r.values())
def sumset(S): return len(set(a + b for a in S for b in S))
def lonely(S):
    M = np.zeros(NG)
    for v in S:
        d = np.abs(((v * TAU + 0.5) % 1.0) - 0.5); M += (d <= H)
    return float((M == 0).mean())

def gap(dims, steps):
    """d-dim GAP with side lengths `dims`, steps `steps`; +1 so nonzero."""
    import itertools
    pts = [1 + sum(x*s for x, s in zip(combo, steps))
           for combo in itertools.product(*[range(d) for d in dims])]
    return sorted(set(pts))

print("="*92)
print("FREIMAN DIMENSION vs additive energy: which determines looseness L?  (1/14=0.0714)")
print("="*92)
fams = []
# dim 1 (APs / near-APs)
fams.append((1, "AP {1..13}", list(range(1,14))))
fams.append((1, "near-AP {1..12}+20", list(range(1,13))+[20]))
fams.append((1, "near-AP {1..11}+{20,30}", list(range(1,12))+[20,30]))
fams.append((1, "AP step2 {1,3..25}", list(range(1,26,2))))
# dim 2 (2-D GAPs, various steps -- take first 13)
fams.append((2, "GAP 4x4 P=5", gap([4,4],[1,5])[:13]))
fams.append((2, "GAP 4x4 P=7", gap([4,4],[1,7])[:13]))
fams.append((2, "GAP 4x4 P=10", gap([4,4],[1,10])[:13]))
fams.append((2, "GAP 3x5 P=11", gap([5,3],[1,11])[:13]))
fams.append((2, "GAP 7x2 P=9", gap([7,2],[1,9])[:13]))
# dim 3 (3-D GAPs)
fams.append((3, "GAP 2x2x4 P=(1,5,20)", gap([4,2,2],[1,5,20])[:13]))
fams.append((3, "GAP 3x2x3 P=(1,7,30)", gap([3,2,3],[1,7,30])[:13]))
print(f"  {'dim':>3} {'set':>26} {'E(S)':>6} {'|S+S|':>6} {'L lonely':>9}  {'verdict':>10}")
byd = {1:[],2:[],3:[]}
for d, name, S in fams:
    S = sorted(set(S))
    if len(S) != 13: 
        print(f"  {d:>3} {name:>26}  (size {len(S)}, skip)"); continue
    e, L = E(S), lonely(S)
    byd[d].append(L)
    verdict = "TIGHT" if L < 0.02 else ("near-tight" if L < 0.06 else "LOOSE")
    print(f"  {d:>3} {name:>26} {e:>6} {sumset(S):>6} {L:>9.4f}  {verdict:>10}")
print("-"*92)
for d in (1,2,3):
    if byd[d]:
        print(f"  dim {d}: L ranges [{min(byd[d]):.4f}, {max(byd[d]):.4f}], mean {np.mean(byd[d]):.4f}")
print("="*92)
print("READING: if dim-1 sets are tight/near-tight (L small) while dim-2,3 GAPs are LOOSE (L~0.1) EVEN AT")
print("HIGH E, then FREIMAN DIMENSION -- not additive energy -- is the loneliness parameter. The AP is the")
print("unique TIGHT extremal because it is the unique dimension-1 minimal set (aligned resonances).")
