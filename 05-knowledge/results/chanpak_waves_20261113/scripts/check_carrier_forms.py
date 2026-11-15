# Independent check of the Finite Carrier Conjecture on a range, using only the list of forms and the function p.
# Usage: python3 check_carrier_forms.py carrier_forms_orders3_5.txt XMAX [basis]
import sys
import numpy as np
forms=[]
only_basis = len(sys.argv) > 3 and sys.argv[3] == 'basis'
for line in open(sys.argv[1]):
    if line.startswith('#'): continue
    left, _, tag = line.split('|')
    if only_basis and 'basis' not in tag: continue
    forms.append(tuple(map(int, left.split())))
XMAX = int(sys.argv[2])
P = np.zeros(XMAX + 2, dtype=np.int64)            # p(u) for u <= XMAX+1
for u in range(2, XMAX + 2):
    P[u] = P[u >> 1] if (u & 1) else P[u >> 1] + (u >> 1) + ((u >> 1) & 1)
hit = np.zeros(XMAX + 1, dtype=bool)
for (a, q, s, e, t, f, n) in forms:
    X = 1
    while True:
        base = t * X + a * P[X] + n + f * (X & 1)
        if q * X + s + base > XMAX: break
        ymax = (XMAX - base) // (q * X + s) + 1
        Y = np.arange(1, ymax + 1, dtype=np.int64)
        F = (q * X + s) * Y + base + a * P[Y] + e * (Y & 1)
        F = F[F <= XMAX]
        hit[F] = True
        X += 1
odd = np.arange(333, XMAX + 1, 2)
miss = odd[~hit[odd]]
print(f"{len(forms)} forms, odd N in [333, {XMAX}]: misses {len(miss)}", (miss[:20].tolist() if len(miss) else ""))
