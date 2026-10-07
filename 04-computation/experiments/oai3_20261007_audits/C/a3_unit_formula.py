# Check THM-4580 (3): Delta_(m-1) = 2^(1-m) Tr P_m(zeta_(2^m)), P_m = z^((10^m-1)/9) prod_{i<=m-4} (1-z^(9*10^i))/(1-z^(10^i))
import numpy as np, json
D = {int(k): v for k, v in json.load(open('/Users/e/Documents/GitHub/math-wt-chessboard-20261006/04-computation/experiments/oai3_20261007_readers/zeroless/zk_values.json'))['Delta'].items()}
ok = True
for m in range(1, 19):
    N = 1 << m
    t = np.arange(1, N, 2, dtype=np.int64)
    def z(e):  # zeta^(t*e)
        return np.exp(2j * np.pi * ((t * (e % N)) % N) / N)
    val = z((10**m - 1) // 9)
    for i in range(0, m - 3):
        val = val * (1 - z(9 * 10**i)) / (1 - z(10**i))
    tr = val.sum().real * 2.0**(1 - m)
    imag = val.sum().imag
    good = abs(tr - D[m - 1]) < 1e-6
    ok &= good
    print(f"m={m:2d}: 2^(1-m) Tr P_m = {tr: .6f} (imag {imag: .1e})   Delta_{m-1} = {D[m-1]}   {'OK' if good else 'MISMATCH'}   max|P_m| = {np.abs(val).max():.4g}  -> rho = {np.abs(val).max()**(1/m):.3f}")
print("ALL OK" if ok else "SOME MISMATCH")
