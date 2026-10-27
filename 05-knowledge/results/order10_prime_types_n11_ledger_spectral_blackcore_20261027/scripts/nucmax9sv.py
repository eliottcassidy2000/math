import numpy as np, math, subprocess, collections
n = 9
lines = subprocess.run(['/opt/homebrew/bin/gentourng', '-q', '9'], capture_output=True, text=True).stdout.split('\n')
best = -1; args = []
for line in lines:
    if len(line) != 36: continue
    S = np.zeros((n, n)); p = 0
    for i in range(n):
        for j in range(i + 1, n):
            S[i, j] = 1.0 if line[p] == '1' else -1.0; S[j, i] = -S[i, j]; p += 1
    ev = np.abs(np.linalg.eigvalsh(1j * S)); nu = float(ev.sum())
    if nu > best + 1e-9: best = nu; args = [(line, sorted(int((S[i] > 0).sum()) for i in range(n)), sorted(ev))]
    elif abs(nu - best) < 1e-9: args.append((line, sorted(int((S[i] > 0).sum()) for i in range(n)), sorted(ev)))
print("n=9 max nu", best, "count", len(args))
print("singular values of first maximizer:", [round(x, 6) for x in args[0][2]])
print("distinct singular-value multisets among maximizers:", len({tuple(round(x, 6) for x in a[2]) for a in args}))
print("score multisets:", collections.Counter(tuple(a[1]) for a in args).most_common(12))
print("8*sqrt(7)+... candidates:", 2*math.sqrt(7)+6*math.sqrt(10), 4*math.sqrt(7)+4*math.sqrt(10), 8*math.sqrt(8))
