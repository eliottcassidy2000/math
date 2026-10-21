# AUDIT a14 check: compare zero classes (a14_zero_classes_400.out) with the session's camion_verify.out (odd classes,
# p <= 199) and recompute the HYP-9358(b) Poisson statistics over all r in [4, p-1], 7 <= p <= 400.
import re
from math import exp, factorial
from collections import Counter
Z = {}
for l in open("a14_zero_classes_400.out"):
    a = list(map(int, l.split())); Z[a[0]] = a[1:]
print("universal {1,2,3} subset of Z_p for every p>=5:", all({1, 2, 3} <= set(Z[p]) for p in Z if p >= 5), " Z_3 =", Z[3])
sess = {}
for l in open("../camion_verify.out"):
    m = re.match(r"p=(\d+) eps=[+-]\d zero classes r mod \d+:([\d ]*)\s+\(", l)
    if m: sess[int(m.group(1))] = list(map(int, m.group(2).split()))
print("odd zero classes agree with camion_verify.out for all", len(sess), "primes <= 199:", all([r for r in Z[p] if r % 2 == 1] == sess[p] for p in sess))
print("even sporadic classes (r>=4) for p <= 199:", {p: [r for r in Z[p] if r % 2 == 0 and r >= 4] for p in Z if p <= 199 and any(r % 2 == 0 and r >= 4 for r in Z[p])})
ps = [p for p in Z if 7 <= p <= 400]
cnt = [len([r for r in Z[p] if r >= 4]) for p in ps]
lam = sum((p - 4) / p for p in ps) / len(ps)
C = Counter(cnt)
print("primes 7..400:", len(ps), " observed mean", sum(cnt) / len(ps), " model mean", lam)
print("k observed Poisson:", [(k, C[k], round(len(ps) * exp(-lam) * lam ** k / factorial(k), 1)) for k in range(6)])
