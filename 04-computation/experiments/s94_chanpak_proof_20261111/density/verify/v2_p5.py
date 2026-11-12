# Verifier (lane density, S94). P5: exact check of the counting inequality, and of what it does and does not obstruct.
import itertools
from math import comb, log2
ok = True
# (1) brute force over all odd classes of width w <= 12
for w in range(2, 13):
    byz = {}
    for mid in itertools.product((-1, 0, 1), repeat=w - 2):
        z = mid.count(0)
        byz[z] = byz.get(z, 0) + 4          # 4 choices of the two nonzero end digits
    tot = sum(byz.values())
    good = (tot == 4 * 3 ** (w - 2)) and all(byz[z] == 4 * comb(w - 2, z) * 2 ** (w - 2 - z) for z in byz)
    # number of pairs in the odd classes = 4^(w-1) (words differing in top and bottom bit)
    good &= sum(byz[z] * 2 ** z for z in byz) == 4 ** (w - 1)
    for m in range(0, w + 1):
        K = 2 ** m
        S = sum(byz[z] * min(2 ** z, K) for z in byz)
        rhs = 2 ** w * sum(comb(w - 2, z) for z in range(0, min(m, w - 2) + 1))
        if S < rhs: good = False
    ok &= good
    print("w=%2d odd classes %d = 4*3^(w-2): %s; class counts by z and the inequality for K = 2^m, m = 0..w: %s" % (w, tot, tot == 4 * 3 ** (w - 2), good))
# (2) exact large-w table: n = w + 2m (balanced cell t' = s = m).  log2 of
#   LB   = 2^w * sum_{z<=m} C(w-2,z)                    (the lane's lower bound)
#   SUM  = sum over odd classes of min(2^z, 2^m)         (exact)
#   CLm  = number of odd classes with at most m zeros     ("one solution per class with <= m zeros")
#   CL   = number of all odd classes = 4*3^(w-2)          ("one solution per class")
print()
print("n, m=t'=s, w | log2 LB - n | log2 SUM - n | log2 CLm - n | log2 CL - n")
def lg(x): return log2(x) if x > 0 else float('-inf')
for n in (200, 1000):
    for frac in (0.05, 0.10, 0.15, 0.17, 0.18, 0.19, 0.20, 0.22, 0.24, 0.25, 0.26):
        m = int(round(frac * n)); w = n - 2 * m
        LB = 2 ** w * sum(comb(w - 2, z) for z in range(0, min(m, w - 2) + 1))
        SUM = sum(4 * comb(w - 2, z) * 2 ** (w - 2 - z) * min(2 ** z, 2 ** m) for z in range(0, w - 1))
        CLm = sum(4 * comb(w - 2, z) * 2 ** (w - 2 - z) for z in range(0, min(m, w - 2) + 1))
        CL = 4 * 3 ** (w - 2)
        print("n=%4d m=%4d w=%4d | %+9.1f | %+9.1f | %+9.1f | %+9.1f" % (n, m, w, lg(LB) - n, lg(SUM) - n, lg(CLm) - n, lg(CL) - n))
print()
print("Reading: LB and SUM exceed 2^n for every m < n/4 (the lane's inequality).  CL = 4*3^(w-2) is below 2^n as soon as")
print("w < n log_3 2 = 0.6309 n, i.e. m > 0.1845 n: a bound of ONE solution per difference class is not obstructed for")
print("0.1845 n < t' < n/4.  The obstruction concerns min(2^z, 2^m) solutions per class (one per (P,c) near the line).")
print("ALL CHECKS PASS" if ok else "SOME CHECK FAILED")
