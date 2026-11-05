#!/usr/bin/env python3
"""defects.py -- H of the tower {1}^n with letter defects: single and double defects; is the
response additive (numeration-like) or multiplicative?"""
import itertools
from ribbon import transfer_system, count_with_matrices, word_str
w = 2
sysw = transfer_system(w)
F = frozenset
one, two, both = F([1]), F([2]), F([1, 2])
def H(word): return count_with_matrices(w, list(word), sysw)
n = 13
tower = [one] * n
T = [H([one] * k) for k in range(0, n + 1)]
print("tribonacci T_k:", T)
for defect in (two, both):
    print(f"single defect {sorted(defect)} at position i (1-based), n={n}: D = H - T_n")
    row = []
    for i in range(3, n + 1):
        wd = tower[:]; wd[i - 1] = defect
        row.append((i, H(wd) - T[n]))
    print("  ", row)
# double defects: additivity test D(i,j) vs D(i)+D(j)
print("double defects {2} at i<j: H - T_n - D(i) - D(j):")
D = {}
for i in range(3, n + 1):
    wd = tower[:]; wd[i - 1] = two; D[i] = H(wd) - T[n]
for i in range(3, n + 1):
    line = []
    for j in range(i + 1, n + 1):
        wd = tower[:]; wd[i - 1] = two; wd[j - 1] = two
        line.append(H(wd) - T[n] - D[i] - D[j])
    print(f"  i={i}: {line}")
# tower length dependence of the single defect at fixed distance from the end
print("single {2}-defect at distance d from the end, as n grows: D(n, n-d)")
for d in range(0, 6):
    print(f"  d={d}:", [ (lambda wd: H(wd) - T[m])([one]*(m-d-1) + [two] + [one]*d) for m in range(d + 3, 14)])
