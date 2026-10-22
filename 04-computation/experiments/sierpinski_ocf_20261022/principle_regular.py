# HYP-9364 beyond circulants: for every regular tournament on n vertices compare kappa = H 2^(n-1)/n! with the
# blow-up gain rho = prod_k (n^2 + s_k^2)/(n^2 - s_k^2).  Usage: principle_regular.py n HK_BINARY.  mac-mini-2026-10-22-sierp.
import subprocess, sys, math
import numpy as np
n = int(sys.argv[1]); HK = sys.argv[2]; m = (n - 1) // 2
gen = subprocess.run(['gentourng', f'-d{m}', f'-D{m}', str(n)], capture_output=True, text=True).stdout
codes = [c for c in gen.split() if not c.startswith('>')]
Hs = [int(x) for x in subprocess.run([HK], input=gen, capture_output=True, text=True).stdout.split()]
assert len(Hs) == len(codes)
iu = np.triu_indices(n, 1)
ratios = []; kap = []; rho = []
B = 20000
for start in range(0, len(codes), B):
    chunk = codes[start:start + B]
    A = np.zeros((len(chunk), n, n))
    bits = np.array([[1 if ch == '1' else 0 for ch in c] for c in chunk])
    A[:, iu[0], iu[1]] = bits
    A[:, iu[1], iu[0]] = 1 - bits
    S = A - np.transpose(A, (0, 2, 1))
    ev = np.linalg.eigvals(S)
    s2 = np.sort(np.abs(ev.imag), axis=1) ** 2
    s2 = s2[:, 1::2]
    r = np.prod((n * n + s2) / (n * n - s2), axis=1)
    rho.extend(r.tolist())
for H, r in zip(Hs, rho):
    k = H * 2 ** (n - 1) / math.factorial(n); kap.append(k); ratios.append(k / r)
ratios = np.array(ratios); kap = np.array(kap); rho = np.array(rho)
rk = lambda v: np.argsort(np.argsort(v))
sp = np.corrcoef(rk(kap), rk(rho))[0, 1]
print(f"n={n}: {len(codes)} regular classes; kappa/rho in [{ratios.min():.4f}, {ratios.max():.4f}], mean {ratios.mean():.4f}; "
      f"kappa in [{kap.min():.4f}, {kap.max():.4f}], rho in [{rho.min():.4f}, {rho.max():.4f}]; spearman(kappa, rho) = {sp:+.3f}; "
      f"argmax kappa = argmax rho: {int(np.argmax(kap)) == int(np.argmax(rho))}", flush=True)
