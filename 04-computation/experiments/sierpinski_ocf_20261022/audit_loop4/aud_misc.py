# AUDIT misc Section-1 claims: (i) X_1, X_2 of loop_engine4 vs THM-4757 closed form / HYP-9367 engine at COMPLEX z
# (excess3_validate.out shows only z = 0.8); (ii) DP-optimal contraction cost of the worst excess-3 core at n = 7
# (claim: about 1e8 multiply-adds) and at n = 23.
import sys, time
import numpy as np
HERE = '/Users/e/Documents/GitHub/math/.claude/worktrees/gifted-chatterjee-298788/04-computation/experiments/sierpinski_ocf_20261022'
sys.path.insert(0, HERE); sys.path.insert(0, HERE + '/audit2')
import loop_engine4 as LE4
import loop_engine_fast as LEF
from second_order import X1_of
from common import regular_classes, code_to_A, skew, paley_A
for name, S in (('regular 7-tournament gentourng class 1', skew(code_to_A(regular_classes(7)[1]))), ('P_11', skew(paley_A(11)))):
    Sh = S / S.shape[0]
    for z in (0.7 * np.exp(0.4j), 0.5 + 0.6j):
        x1 = LE4.X_e(Sh, 1, z); x2 = LE4.X_e(Sh, 2, z)
        r1 = X1_of(Sh.astype(complex), z); r2 = LEF.XB(Sh, 2, complex(z))
        print(f"{name}, z = {z:.3f}: X_1 engine {x1:.13e} closed form {r1:.13e} (|diff| {abs(x1 - r1):.1e}); X_2 engine {x2:.13e} HYP-9367 engine {r2:.13e} (|diff| {abs(x2 - r2):.1e})")
for n in (7, 23):
    costs = []
    for i in range(len(LE4.load_cores(3))):
        cost, tree, opn, ops, out = LE4.network(3, i).plan(5, n, True)
        costs.append(cost)
    j = int(np.argmax(costs))
    print(f"n = {n}: DP-optimal contraction cost (multiply-adds, with W, 5-basis): max {max(costs):.3e} (core {j}, V = {LE4.load_cores(3)[j][0]}), median {np.median(costs):.2e}, total {sum(costs):.2e}")
