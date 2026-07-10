from fractions import Fraction as F
import random
random.seed(90)
def ni(x):
    r = x % 1
    return min(r, 1 - r)
# verify the descent lemma exactly: tau = (sigma+1)/2
viol = 0
for _ in range(200000):
    sigma = F(random.randrange(1, 5000), random.randrange(1, 5000))
    tau = (sigma + 1) / 2
    vp = random.randrange(1, 300)
    o = 2 * random.randrange(0, 150) + 1
    if ni(2 * vp * tau) != ni(vp * sigma): viol += 1
    if ni(o * tau) < ni(o * sigma) / 2: viol += 1
print(f"descent lemma exact checks: violations {viol} / 200k (evens exact, odds >= half)")
# the merge mechanics on monad's corner shape: doubling-rich even-rich families
def reduced(S):
    """one descent: R = {(v//2 if v even else v): threshold} taking max thresholds."""
    R = {}
    for v in S:
        if v % 2 == 0:
            w, th = v // 2, F(1, 14)
        else:
            w, th = v, F(1, 7)
        R[w] = max(R.get(w, F(0)), th)
    return R
# corner-style instances: long doubling chains + covering fillers
tests = [
    [3, 6, 12, 24, 5, 10, 20, 40, 7, 14, 28, 9, 11],     # two length-4 chains + one length-3 + odds
    [1, 2, 4, 8, 16, 3, 6, 12, 5, 10, 7, 9, 13],          # chains rooted 1,3,5 + odds
    [45, 90, 180, 360, 63, 126, 252, 33, 66, 132, 35, 70, 11],
]
for S in tests:
    R = reduced(S)
    merges = 13 - len(R)
    print(f"S={S}: after 1 descent: {len(R)} distinct (merges {merges}), "
          f"thresholds: {sorted(set(str(t) for t in R.values()))}")
