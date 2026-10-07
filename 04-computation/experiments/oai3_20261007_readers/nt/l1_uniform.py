# uniform constant: for q in [Q_j, Q_{j+1}], L(1,chi) >= L_lower(Q_{j+1}) (monotone in q at fixed params),
# hence L(1,chi)*log q >= L_lower(Q_{j+1}) * log Q_j.  Odd chi (a=1).
import mpmath as mp
from l1_lower_bound import best
mp.mp.dps = 20
grid = [mp.mpf(10) ** (4 + j / 4) for j in range(0, 41)]   # 1e4 .. 1e14
vals = [best(q, 1)[0] for q in grid]
cmin = min(vals[j + 1] * mp.log(grid[j]) for j in range(len(grid) - 1))
print("uniform: L(1,chi) >= %s/log q for 1e4 <= q <= 1e14" % mp.nstr(cmin, 4))
cmin10 = min(vals[j + 1] * mp.log(grid[j]) for j in range(24, len(grid) - 1))
print("uniform: L(1,chi) >= %s/log q for 1e10 <= q <= 1e14" % mp.nstr(cmin10, 4))
for j in range(0, 41, 4): print("  q=1e%.2f  L>= %s  (%s/log q)" % (4 + j / 4, mp.nstr(vals[j], 5), mp.nstr(vals[j] * mp.log(grid[j]), 4)))
