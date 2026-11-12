#!/usr/bin/env python3
"""Closed forms for B = transitive tournament on 0..m-1 with one back arc (j -> i), gap k = j - i - 1 >= 1.
With p = i (vertices below i), s = m - 1 - j (vertices above j):
  H(B) = 1 + 2^k,   W_ij = 2^k,
  U_B(y) = 2^(q+s+1) (3^k + 2^(k+1) - 1)                 y < i,  q = i - y - 1
  U_B(i) = 2^s (2^(k+1) - 1)
  U_B(y) = 2^s (2^(b+1) + 2^(a+1) 3^b - [a = 0])         i < y < j,  a = y - i - 1,  b = j - y - 1
  U_B(j) = 2^s (2^k + 1)
  U_B(y) = 2^(m-1-y) (2^k + 1)                           y > j
Checked against a subset DP for all m <= 11."""
from roots import make, root_data
ok = bad = 0
for m in range(3, 12):
    for i in range(m):
        for j in range(i + 2, m):
            k = j - i - 1; s = m - 1 - j
            H, U, W = root_data(make(m, {(j, i)}))
            exp = []
            for y in range(m):
                if y < i: q = i - y - 1; exp.append(2**(q + s + 1) * (3**k + 2**(k + 1) - 1))
                elif y == i: exp.append(2**s * (2**(k + 1) - 1))
                elif y < j: a = y - i - 1; b = j - y - 1; exp.append(2**s * (2**(b + 1) + 2**(a + 1) * 3**b - (1 if a == 0 else 0)))
                elif y == j: exp.append(2**s * (2**k + 1))
                else: exp.append(2**(m - 1 - y) * (2**k + 1))
            good = (H == 1 + 2**k) and (W[i][j] == 2**k) and (U == exp)
            ok += good; bad += (not good)
print("one back arc, 3 <= m <= 11: formulas hold in %d cases, fail in %d" % (ok, bad))
