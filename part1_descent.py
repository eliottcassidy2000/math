"""Part 1: Descent/refinement. f(p^k) vs V=(1-2/n)^m."""
from lrc import f_pk, V_target
from fractions import Fraction

# Speed sets. Choose generic sets and some structured ones.
SETS = {
    "m=3 {1,2,3}": [1,2,3],
    "m=3 {1,3,4}": [1,3,4],
    "m=3 {2,3,7}": [2,3,7],
    "m=4 {1,2,3,4}": [1,2,3,4],
    "m=4 {1,3,4,7}": [1,3,4,7],
    "m=5 {1,2,3,4,5}": [1,2,3,4,5],
    "m=5 {1,3,4,5,9}": [1,3,4,5,9],
}
PRIMES = [2,3,5,7]
MAXQ = 4000

for name, speeds in SETS.items():
    m = len(speeds)
    V = V_target(m)
    print(f"\n=== {name}  V=(1-2/{m+1})^{m} = {V} = {float(V):.5f} ===")
    for p in PRIMES:
        row = []
        k = 1
        while p**k <= MAXQ:
            frac, cnt, tot = f_pk(speeds, p, k)
            row.append((k, p**k, frac, float(frac)))
            k += 1
        s = f"  p={p}: "
        s += "  ".join(f"k{k}(q={q}):{fl:.4f}" for k,q,fr,fl in row)
        print(s)
