"""HYP-9565 (A) heuristic (EXPLORATORY): random ribbon tournaments of span s (each pair with 1 <= j-i <= s reversed
independently with probability 1/2, all longer pairs forward) at n = 14: lambda, twin distance, module defect, H.
Their H is at most (s+3)^(2n) (Lemma B3), so under (S-inf) long ribbons are non-creators; the question is whether
the decomposability defect grows with s."""
import random, sys
from verify_escapers import lam, ham
from defect import defsub, twin
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
n = 14
for s in (2, 3, 4, 5, 6):
    rows = []
    for t in range(12):
        out = [0] * n
        for i in range(n):
            for j in range(i + 1, n):
                if j - i <= s and random.random() < 0.5:
                    out[j] |= 1 << i
                else:
                    out[i] |= 1 << j
        L = lam(out)
        if L == 0:
            continue
        d, M = defsub(out, n)
        rows.append((L, d, twin(out, n), min(L, d), ham(out, (1 << n) - 1)))
    if rows:
        print('s=%d strong samples=%d  lambda %s  module defect %s  twin %s  dec %s  H range [%d, %d]' % (
            s, len(rows), sorted(r[0] for r in rows), sorted(r[1] for r in rows), sorted(r[2] for r in rows),
            sorted(r[3] for r in rows), min(r[4] for r in rows), max(r[4] for r in rows)))
