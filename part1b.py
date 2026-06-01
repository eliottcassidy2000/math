"""Part 1b: focus on convergence to V for sets coprime to p, monotonicity, error."""
from lrc import f_pk, V_target
from fractions import Fraction

SETS = {
    "m=3 {1,2,3}": [1,2,3],
    "m=3 {2,3,7}": [2,3,7],
    "m=3 {1,3,4}": [1,3,4],
    "m=4 {1,3,4,7}": [1,3,4,7],
    "m=4 {2,3,5,7}": [2,3,5,7],
    "m=5 {1,3,4,5,9}": [1,3,4,5,9],
}
PRIMES = [2,3,5,7]
MAXQ = 4000

def coprime_set(speeds, p):
    return all(s % p != 0 for s in speeds)

for name, speeds in SETS.items():
    m = len(speeds)
    V = V_target(m)
    print(f"\n=== {name}  V={float(V):.5f} ===")
    for p in PRIMES:
        cop = coprime_set(speeds, p)
        vals = []
        k = 1
        while p**k <= MAXQ:
            frac, cnt, tot = f_pk(speeds, p, k)
            vals.append(float(frac))
            k += 1
        # error from V at deepest level, and monotone trend on |f-V|
        errs = [abs(v - float(V)) for v in vals]
        last = vals[-1]
        tag = "coprime" if cop else "p|some_speed"
        print(f"  p={p} [{tag}]: f-vals={['%.4f'%v for v in vals]}  last_err={errs[-1]:.4f}")
