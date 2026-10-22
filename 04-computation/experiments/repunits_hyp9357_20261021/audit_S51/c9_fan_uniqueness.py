"""S51 audit c9: the late-added fan_uniqueness.out (commit 416b090462) and statement 6 of THM-4748, from the auditor's
own inventory (c4 cache): for H = 2^k - 1 (5 <= k <= 11), the strong classes, their least feedback arc numbers, the
ones with fa <= 3, and whether they are the Fibonacci fan F_k or its converse.  Prints ALL CHECKS PASSED or failures."""
import os, time
from collections import Counter
import numpy as np
import alib51 as L

FAIL = []
T0 = time.time()
SCR = r'C:\Users\Eliott\AppData\Local\Temp\claude\C--Users-Eliott-Documents-GitHub-ephrepos-math\e1fb236f-ec78-4211-aa52-c187130b32b0\scratchpad\audit_S51'


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


Z = np.load(os.path.join(SCR, 'inventory_2100.npz'))
orders = sorted(int(k[1:]) for k in Z.files if k.startswith('H'))
data = {n: (Z[f'H{n}'], Z[f'A{n}']) for n in orders}
sess = {5: (4, [6, 7], 1, True), 6: (2, [8], 2, False), 7: (4, [8, 9], 2, False), 8: (24, [9, 10, 11], 2, False),
        9: (2, [11], 2, False), 10: (8, [10, 11, 12], 2, False), 11: (2, [13], 2, False)}
mine = {}
for k in range(5, 12):
    h = 2 ** k - 1
    a, b = k // 2, k + 1 - k // 2
    fan = L.one_apex(k + 2, {a, b, a + b})
    cf = L.canon(k + 2, fan)
    cc = L.canon(k + 2, L.converse(fan))
    low = []; fas = Counter()
    for n in orders:
        H, A = data[n]
        for t in np.nonzero(H == h)[0]:
            f = int(L.fas(n, A[t]))
            fas[f] += 1
            if f <= 3:
                low.append((n, L.canon(n, A[t])))
    nf = sum(1 for n, c in low if n == k + 2 and (c[:2] == cf[:2] or c[:2] == cc[:2]))
    mine[k] = (len(low), sorted(set(n for n, c in low)), nf, cf[:2] == cc[:2])
    print(f'     k = {k}, H = {h}: strong classes {sum(fas.values())}, least-fa distribution {dict(sorted(fas.items()))}; with fa <= 3: {mine[k]}')
check(mine == sess, 'fan_uniqueness.out reproduced: number of strong classes with fa <= 3, their orders, fan/converse among them, fan self-converse')
check(all(min(int(L.fas(n, A[t])) for n in orders for H, A in [data[n]] for t in np.nonzero(H == 2 ** k - 1)[0]) == 3 for k in range(5, 12)),
      'least feedback arc number over ALL strong classes with H = 2^k - 1 is exactly 3 for 5 <= k <= 11')
print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f'({time.time() - T0:.0f}s)')
