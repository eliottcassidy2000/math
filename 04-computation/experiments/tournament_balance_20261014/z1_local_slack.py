# Slack of the two corner conditions L4, U4 of the (Z1) certificate (all other conditions hold for every block):
#   L4: 1 + sigma/2 + 2(beta(v) + gamma(v) + alpha(v)) - (betabar(v) + gammabar(v) + alphabar(v) + s) >= 0
#   U4: 1 + sigma/2 + 2(betabar(v) + gammabar(v) + alphabar(v)) - (beta(v) + gamma(v) + alpha(v) + e) >= 0
# (bar = complementary part, X not containing v).  Minimum over all blocks of each size, all z-patterns, all v.
import sys
from itertools import product
from fractions import Fraction as Fr
from z1_cone_certificate import strong_blocks, params
M = int(sys.argv[1]) if len(sys.argv) > 1 else 6
for size in [1] + list(range(3, M + 1)):
    worst_L = worst_U = None
    for A, s in strong_blocks(size):
        if s != size: continue
        for zp in product((False, True), repeat=s):
            P, vp = params(A, s, zp)
            for v, V in enumerate(vp):
                bb, gb, ab = P['beta'] - V['beta'], P['gamma'] - V['gamma'], P['alpha'] - V['alpha']
                L4 = 1 + P['sigma'] / 2 + 2 * (V['beta'] + V['gamma'] + V['alpha']) - (bb + gb + ab + P['s'])
                U4 = 1 + P['sigma'] / 2 + 2 * (bb + gb + ab) - (V['beta'] + V['gamma'] + V['alpha'] + P['e'])
                if worst_L is None or L4 < worst_L[0]: worst_L = (L4, zp, v, P['sigma'])
                if worst_U is None or U4 < worst_U[0]: worst_U = (U4, zp, v, P['sigma'])
    print(f"size {size}: min L4 = {worst_L[0]} (sigma {worst_L[3]}, pattern {''.join('1' if x else '0' for x in worst_L[1])}, v {worst_L[2]}); "
          f"min U4 = {worst_U[0]} (sigma {worst_U[3]})", flush=True)
