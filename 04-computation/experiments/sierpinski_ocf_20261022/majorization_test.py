# Does the squared skew spectrum (sigma_1^2 >= sigma_2^2 >= ...) of the rotational R_d majorise that of every regular
# d-tournament?  (Sufficient for Wormald's Conjecture 1 via Schur-convexity of rho.)  mac-mini-2026-10-22-sierp.
import subprocess, sys
import numpy as np
for d in [int(x) for x in sys.argv[1:]]:
    m = (d - 1) // 2
    gen = subprocess.run(['gentourng', f'-d{m}', f'-D{m}', str(d)], capture_output=True, text=True).stdout
    codes = [c for c in gen.split() if not c.startswith('>')]
    iu = np.triu_indices(d, 1)
    A = np.zeros((d, d))
    for i in range(d):
        for r in range(1, m + 1): A[i, (i + r) % d] = 1
    sR = np.sort(np.abs(np.linalg.eigvals(A - A.T).imag) ** 2)[::-1][0:2 * m:2]
    cR = np.cumsum(sR)
    viol = 0; worst = 0.0
    for start in range(0, len(codes), 20000):
        chunk = codes[start:start + 20000]
        B = np.zeros((len(chunk), d, d)); bits = np.array([[1 if ch == '1' else 0 for ch in c] for c in chunk])
        B[:, iu[0], iu[1]] = bits; B[:, iu[1], iu[0]] = 1 - bits
        s2 = np.sort(np.abs(np.linalg.eigvals(B - np.transpose(B, (0, 2, 1))).imag) ** 2, axis=1)[:, ::-1][:, 0:2 * m:2]
        c = np.cumsum(s2, axis=1)
        diff = c - cR[None, :]
        bad = np.any(diff > 1e-7, axis=1)
        viol += int(bad.sum()); worst = max(worst, float(diff.max()))
    print(f"d={d}: {len(codes)} regular classes; R_d majorises all: {viol == 0} (violations {viol}, max partial-sum excess {worst:.2e}); sigma^2(R_d) = {np.round(sR, 4).tolist()}", flush=True)
