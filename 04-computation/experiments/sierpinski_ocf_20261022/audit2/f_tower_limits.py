# Audit2 (D): two-loop limits for the ternary towers S_k = C_3^[k] to k = 11 (n = 177147), FFT predictor
# (e_predictor.two_loop_circ), with Richardson extrapolation in 1/n (ratio 3).  Output: f_tower_limits.out
import sys; sys.path.insert(0, '.')
from e_predictor import two_loop_circ
def tc(N):
    out = []
    for z in range(1, N):
        t = z
        while t % 3 == 0: t //= 3
        if t % 3 == 1: out.append(z)
    return out
vals = []
for k in range(6, 12):
    N = 3 ** k; t = two_loop_circ(N, tc(N)); vals.append((N, t['nTH'], t['y_pred'], t['nF8'], t['E0']))
    print(k, N, '%.7f %.7f %.7f %.7f' % (t['nTH'], t['y_pred'], t['nF8'], t['E0']))
for i in range(len(vals) - 1):
    print('Richardson nTheta', (3 * vals[i + 1][1] - vals[i][1]) / 2, ' y_pred', (3 * vals[i + 1][2] - vals[i][2]) / 2)
