# Preimages of the aliquot map s(n) = sigma(n) - n for values v <= V (complete: a composite preimage of v has
# s(n) >= 1 + sqrt(n), so n <= (v-1)^2; primes give v = 1).  Untouchables, the twin-prime "double hits"
#   s(p^2) = p + 1 = s(2q) = q + 3  <=>  p = q + 2 (twin primes),
# and the unique self-hit 2q = q + 3 (q = 3, v = 6, the perfect number).  mac-mini-2026-10-25-pyth (2026-10-26).
import sys
from sympy import isprime, factorint
V = int(sys.argv[1]) if len(sys.argv) > 1 else 300
N = (V - 1) ** 2 + 1
import numpy as np
sig = np.zeros(N + 1, dtype=np.int64)
for d in range(1, N + 1): sig[d::d] += d
s = sig - np.arange(N + 1)
pre = {v: [] for v in range(V + 1)}
for n in range(2, N + 1):
    v = int(s[n])
    if v <= V and not isprime(n): pre[v].append(n)
pre[0] = [1]
def shape(n):
    f = factorint(n); return '*'.join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(f.items()))
print("v : composite preimages of v under s (v = 1 also has every prime)")
for v in range(0, 41):
    print(f"  {v:3d}: {pre[v]}" + ("   [untouchable]" if not pre[v] and v not in (1,) else "") + ("   [perfect]" if v in pre[v] else ""))
unt = [v for v in range(2, V + 1) if not pre[v] and v != 1]
print(f"untouchables <= {V}: {unt}")
print(f"odd untouchables <= {V}: {[v for v in unt if v % 2]};  prime untouchables: {[v for v in unt if isprime(v)]}")
# twin-prime double hits
dh = [v for v in range(2, V + 1) if any(int(round(m ** 0.5)) ** 2 == m and isprime(int(round(m ** 0.5))) for m in pre[v])
      and any(m % 2 == 0 and isprime(m // 2) for m in pre[v])]
pred = [p + 1 for p in range(5, V) if isprime(p) and isprime(p - 2)]
print(f"values with both a p^2 and a 2q preimage: {dh}")
print(f"  for v > 3 equals {{p + 1 : (p-2, p) twin primes}}: {[v for v in dh if v > 3] == [x for x in pred if x <= V]} (v = 3 also qualifies: 4 = 2^2 = 2*2)")
print(f"  self-hits (v = 2q in s^-1(v)): {[v for v in dh if v in pre[v]]}")
perfect = [n for n in range(2, N + 1) if s[n] == n]
print(f"perfect numbers <= {N}: {perfect}")
