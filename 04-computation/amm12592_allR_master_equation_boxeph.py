"""Master equation + dyadic parity theorem (exact checks).

With Delta_i = (p-q) + c_i (0<=i<=R-2), Delta_{R-1} = -1:
  (*)  <=>  sum_{i=0}^{R-2} x^i c_i = q^{R-1} - E_{R-1},   E_m := -1+x+...+x^m.
Parity: q^{R-1} - E_{R-1} == 0 coefficientwise mod 2  <=>  R is a power of 2.
"""
from math import comb
from amm12592_allR_family_toolbox_boxeph import *

def Em(m): return [-1] + [1]*m if m >= 0 else []

# (p-q) * sum_{i<=R-2} x^i == E_{R-2} + 2x^{R-1}
for R in range(2, 300):
    lhs = pmul([-1,2], [1]*(R-1))
    assert lhs == padd(Em(R-2), pshift([2], R-1)), R
print("(p-q)*geom == E_{R-2} + 2x^{R-1}: PASS for R=2..299")

# equivalence: sum x^i c_i = q^{R-1} - E_{R-1}
for R in (4, 8, 16, 32, 64):
    K = psub(qpow(R-1), Em(R-1))
    # reconstruct epoch sum from arbitrary even c's? just verify algebra:
    lhs = padd(padd(pmul([-1,2], [1]*(R-1)), K), pshift([-1], R-1))
    assert lhs == qpow(R-1), R
print("backbone + K + (-x^{R-1}) == q^{R-1}: PASS (master equation equivalence)")

# dyadic parity theorem
def parity_ok(R):
    K = psub(qpow(R-1), Em(R-1))
    return all(v % 2 == 0 for v in K)
good = [R for R in range(2, 600) if parity_ok(R)]
print(f"R with q^(R-1)-E_(R-1) even (pure-backbone parity-consistent): {good}")
print(f"== powers of 2: {good == [2**t for t in range(1, 10) if 2**t < 600]}")
