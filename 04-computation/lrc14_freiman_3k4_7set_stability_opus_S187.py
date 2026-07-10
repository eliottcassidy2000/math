"""
lrc14_freiman_3k4_7set_stability_opus_S187.py  (opus-2026-07-09-S187)

THE FINISH-LINE PIECE (THM-675's missing STABILITY, = HYP-5682). The covering residual reduces to:
a 7-element majority-parity class with <= 12-13 distinct RESTRICTED pair sums {a_i+a_j : i<j} is NEAR-AP.
Freiman restricted-sumset floor: |A +^ A| >= 2k-3 = 11 for k=7. TARGET: characterize 7-sets with
|A+^A| in {11,12,13} (within 2 of the floor) -- confirm they are all contained in a SHORT AP (Freiman
3k-4 stability), giving the EXPLICIT finite perturbation family THM-675 needs.
Enumerate primitive 7-sets {0=a1<...<a7}, gcd(diffs)=1, spread<=B; classify by restricted-sumset size;
for the near-minimal ones, report the smallest containing AP length (Freiman: <= |A+^A|-k+2).
"""
from itertools import combinations
from math import gcd
from functools import reduce

def restricted_sumset(A):
    return len(set(A[i]+A[j] for i in range(len(A)) for j in range(i+1,len(A))))

def min_containing_AP_len(A):
    """smallest L such that A is contained in an AP of length L (some difference d)."""
    A=sorted(A); n=len(A); best=None
    span=A[-1]-A[0]
    # candidate differences d: divisors-ish of gaps; try all d dividing span with d>=1
    for d in range(1, span+1):
        if (A[-1]-A[0])%d!=0: 
            pass
        # A in AP a0 + d*Z iff all (a-A[0])%d==0
        if all((a-A[0])%d==0 for a in A):
            L=(A[-1]-A[0])//d + 1
            if best is None or L<best: best=L
    return best if best is not None else 10**9

B=16
buckets={}   # restricted-sumset size -> list of (set, minAPlen)
count=0
for combo in combinations(range(1,B+1), 6):   # a2..a7 in [1,B], a1=0
    A=(0,)+combo
    diffs=[A[i+1]-A[i] for i in range(6)]
    g=reduce(gcd,diffs)
    if g!=1: continue                          # primitive only (affine-normalized)
    count+=1
    r=restricted_sumset(A)
    if r<=15:
        buckets.setdefault(r,[]).append((A,min_containing_AP_len(A)))
print(f"enumerated {count} primitive 7-sets (a1=0, gcd(diffs)=1, spread<={B})")
print("="*80)
print(f"  {'|A+^A|':>7} {'#sets':>7} {'max minAP-len':>13} {'all in AP<=8?':>13}")
for r in sorted(buckets):
    sets=buckets[r]; maxlen=max(l for _,l in sets)
    allshort = all(l<=8 for _,l in sets)
    print(f"  {r:>7} {len(sets):>7} {maxlen:>13} {'YES' if allshort else 'NO':>13}")
print("-"*80)
print("NEAR-MINIMAL family (|A+^A| in {11,12,13}) -- the explicit perturbation shapes:")
for r in [11,12,13]:
    if r not in buckets: continue
    print(f"  |A+^A|={r}: {len(buckets[r])} primitive shapes, all in AP of length <= {max(l for _,l in buckets[r])}")
    for A,l in buckets[r][:6]:
        d=[A[i+1]-A[i] for i in range(6)]
        print(f"     {A}  minAP-len={l}  diffs={d}")
print("="*80)
print("READING: if every 7-set with |A+^A|<=13 sits in an AP of length <= 8 (Freiman 3k-4 stability),")
print("the near-minimal-burden majority classes are EXACTLY {7 points of a length<=8 AP} -- a finite,")
print("affine-normalized family => the finite check THM-675 needs to close the covering residual.")
