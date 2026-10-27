"""Maximum nuclear norm over all tournament classes of order n (gentourng), the maximizers, and the bound."""
import numpy as np, math, sys, subprocess
n = int(sys.argv[1])
bound = (n - 1) * math.sqrt(n) if n % 2 == 1 else n * math.sqrt(n - 1)
lines = subprocess.run(['/opt/homebrew/bin/gentourng', '-q', str(n)], capture_output=True, text=True).stdout.split('\n')
best = -1; arg = []; hist = {}
cnt = 0
for line in lines:
    if len(line) != n * (n - 1) // 2: continue
    cnt += 1
    S = np.zeros((n, n)); p = 0
    for i in range(n):
        for j in range(i + 1, n):
            S[i, j] = 1.0 if line[p] == '1' else -1.0; S[j, i] = -S[i, j]; p += 1
    nu = float(np.abs(np.linalg.eigvalsh(1j * S)).sum())
    sc = sorted(int((S[i] > 0).sum()) for i in range(n))
    if nu > best + 1e-9: best = nu; arg = [(line, sc)]
    elif abs(nu - best) < 1e-9: arg.append((line, sc))
print(f"n={n} classes={cnt} max nu={best:.6f} bound={bound:.6f} ratio={best/bound:.6f} maximizers={len(arg)} scores={[a[1] for a in arg[:4]]}")
