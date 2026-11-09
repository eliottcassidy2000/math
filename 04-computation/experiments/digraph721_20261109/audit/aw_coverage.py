"""AUDIT D1: coverage of the computational inputs of THM-5055 by independently verified witnesses.
Reads aw_verify output (one line per witness: tag set n k d=.. EVEN/ODD ... OK/BAD) and checks, using only the
lines marked OK (independently verified: d by inclusion-exclusion, even by a stabiliser-chain generating set,
flags from degrees):
  (i)   base of section 5: for every n <= 12 and every k in F(n), an even graph with d = k on <= n vertices
        (padding then gives exactly n vertices);
  (ii)  Lemma 7: for c in {1,3,5,7,11,13} and n in {8,9}, a GOOD graph on n vertices with d = 2^(n-1) + c;
  (iii) Lemma 8 base: for every odd j < 2^12, some n0 <= L(j) + 8 with GOOD graphs on n0 and n0 + 1 vertices with
        d = 2^(n0-1) - j and 2^n0 - j.
usage: python3 aw_coverage.py aw_verify_witnesses.out"""
import sys
from collections import defaultdict
def M(n): return 2**n - n if n % 2 else 2**n - n - 1
def inF(n, k):
    if k < 1 or k % 2 == 0 or k > M(n) or k in (3, 7, 13): return False
    if n % 2 == 1 and n >= 5 and k == 2**n - n - 2: return False
    return True
even_min_order = {}          # k -> least order of a verified even witness
good = set()                 # (n, k) verified good
for line in open(sys.argv[1]):
    p = line.split()
    if len(p) < 6 or p[0] in ("SUMMARY", "BF", "SKIP") or p[-1] != "OK": continue
    tag, st, n, k = p[0], p[1], int(p[2]), int(p[3])
    assert p[5] == "EVEN"
    even_min_order[k] = min(even_min_order.get(k, 99), n)
    if st == "G": good.add((n, k))
fail = 0
cnt_i = 0
for n in range(1, 13):
    for k in range(1, 2**n, 2):
        if inF(n, k):
            cnt_i += 1
            if even_min_order.get(k, 99) > n: fail += 1; print("MISSING base", n, k)
print("(i) F(n) for n <= 12: checked", cnt_i, "(n,k) pairs; missing", fail)
f2 = 0
for c in (1, 3, 5, 7, 11, 13):
    for n in (8, 9):
        if (n, 2**(n - 1) + c) not in good: f2 += 1; print("MISSING Lemma7", c, n)
print("(ii) Lemma 7 (12 good graphs): missing", f2)
f3 = 0; worst = 0
for j in range(1, 2**12, 2):
    L = j.bit_length(); ok = None
    for n0 in range(2, L + 9):
        if 2**(n0 - 1) > j and (n0, 2**(n0 - 1) - j) in good and (n0 + 1, 2**n0 - j) in good: ok = n0; break
    if ok is None: f3 += 1; print("MISSING seed", j)
    else: worst = max(worst, ok - L)
print("(iii) Lemma 8 base, odd j < 4096 (2048 values): missing", f3, "; max n0 - L(j) =", worst)
print("ALL COVERED" if fail == f2 == f3 == 0 else "COVERAGE FAILURE")
