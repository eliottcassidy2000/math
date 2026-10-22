# THM-086 check at n = 7 = 2*3+1: forward-edge polynomial F(T,x) mod 3 has nonzero coefficients only at x^0, x^3, x^6,
# each = H mod 3 (all 456 classes).  Steiner/Sierpinski collapse at orders 2*3^k+1.  mac-mini-2026-10-22-sierp.
import itertools, subprocess
lines = subprocess.run(['gentourng', '7'], capture_output=True, text=True).stdout.split()
n = 7; bad = 0; tot = 0
for code in lines:
    if code.startswith('>'): continue
    adj = [[0]*n for _ in range(n)]; p = 0
    for i in range(n):
        for j in range(i+1, n):
            if code[p] == '1': adj[i][j] = 1
            else: adj[j][i] = 1
            p += 1
    F = [0]*n
    for perm in itertools.permutations(range(n)):
        F[sum(adj[perm[t]][perm[t+1]] for t in range(n-1))] += 1
    H = F[n-1]
    for i in range(n):
        expect = H % 3 if i in (0, 3, 6) else 0
        if F[i] % 3 != expect: bad += 1
    tot += 1
print(f"n=7: {tot} classes, coefficient mismatches against H*(1 + x^3 + x^6) mod 3: {bad}")
