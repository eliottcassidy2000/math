import sys
sys.argv = ['x', '2']
exec(open('trans_plus_one.py').read().split("M = int(sys.argv[1])")[0])
from itertools import product
from fractions import Fraction as F
for m in (4, 5, 6):
    for mid in product((0, 1), repeat=m - 2):
        w = (1,) + mid + (0,)
        H, inA, cnt = stats(w)
        pz = [F(inA[i], H) for i in range(m)]
        bal = [i + 1 for i in range(m) if F(1, 3) <= pz[i] <= F(2, 3)]
        print(''.join(map(str, w)), 'H =', H, ' P(u_i before z):', ' '.join(f"{float(p):.3f}" for p in pz), ' balanced i:', bal)
