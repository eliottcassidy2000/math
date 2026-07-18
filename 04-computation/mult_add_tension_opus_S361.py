# opus-2026-07-17-S361 -- HYP-7440: THE MULTIPLICATIVE / ADDITIVE TENSION.
# MULTIPLICATIVE (THM-1035, kernel-pure): a counterexample must contain, for
#   each q in {8,...,14}, a speed divisible by q -- else the sieve gives 1/q.
# ADDITIVE: BONF5 fails on families with additive structure (the S331 failing
#   packet had 3-APs; the S332 certified packet was Sidon-like).  So the dense
#   core wants SMALL doubling.
# QUESTION: can a family do BOTH?  If not, the residual is empty.
from math import gcd
import random, itertools

MODULI = [8,9,10,11,12,13,14]
def blocks_all(V): return all(any(v % q == 0 for v in V) for q in MODULI)
def doubling(V):
    S = set(V); sums = {a+b for a in S for b in S}
    return len(sums) / len(S)          # |A+A|/|A|; Sidon ~ (k+1)/2, AP ~ 2

print("(1) DO COMPARABLE FAMILIES BLOCK ALL SEVEN MODULI?")
random.seed(361)
blk = tot = 0
for _ in range(4000):
    m = random.randint(23, 400)
    V = sorted(random.sample(range(m, 13*m), 13))
    tot += 1
    if blocks_all(V): blk += 1
print(f"    {tot} random comparable 13-families: block all seven in {blk}"
      f" ({100*blk/tot:.1f}%)")
print("    => the sieve already kills ~%.0f%% of comparable families outright." % (100*(tot-blk)/tot))

print()
print("(2) DO ARITHMETIC PROGRESSIONS BLOCK ALL SEVEN?  (small doubling ~2)")
apblk = aptot = 0
for _ in range(2000):
    a = random.randint(23, 500); d = random.randint(1, 200)
    V = [a + k*d for k in range(13)]
    if len(set(V)) < 13: continue
    aptot += 1
    if blocks_all(V): apblk += 1
print(f"    {aptot} random APs: block all seven in {apblk} ({100*apblk/max(aptot,1):.1f}%)")
print(f"    an AP has doubling {doubling([1+k*1 for k in range(13)]):.2f} (minimal)")

print()
print("(3) THE DECISIVE TEST: among families that BLOCK all seven, is the")
print("    doubling forced UP?  (tension would mean blockers cannot be additively structured)")
blk_d = []; non_d = []
for _ in range(6000):
    m = random.randint(23, 400)
    V = sorted(random.sample(range(m, 13*m), 13))
    d = doubling(V)
    (blk_d if blocks_all(V) else non_d).append(d)
blk_d.sort(); non_d.sort()
print(f"    blockers    : n={len(blk_d)}, median doubling {blk_d[len(blk_d)//2]:.2f}")
print(f"    non-blockers: n={len(non_d)}, median doubling {non_d[len(non_d)//2]:.2f}")
print()
print("(4) CAN A FAMILY HAVE BOTH minimal doubling AND block all seven?")
found = []
for _ in range(20000):
    a = random.randint(23, 3000); d = random.randint(1, 500)
    V = [a + k*d for k in range(13)]
    if len(set(V)) < 13: continue
    if blocks_all(V):
        found.append((a, d, doubling(V)))
        if len(found) >= 3: break
if found:
    print("    YES -- explicit APs (doubling 2.00, the additive minimum) blocking all seven:")
    for (a,d,db) in found:
        print(f"      a={a}, d={d}: {[a+k*d for k in range(13)][:4]}...  doubling {db:.2f}")
    print()
    print("    => NO TENSION.  The multiplicative sieve constraint and minimal")
    print("       additive doubling are SIMULTANEOUSLY satisfiable, by explicit")
    print("       arithmetic progressions.  The over-determination hope FAILS.")
else:
    print("    none found in 20000 tries -- tension may be real; investigate further")
