"""Cross-check of the C primality test (pair closure) against the definition (no module of size 2..n-1, all subsets)
on random tournaments and random g-invariant tournaments of order 9 (g of types 1^3 3^2, 3^3)."""
import random, subprocess, sys
random.seed(20261027)
n = 9
def brute_prime(out):
    for mask in range(1, 1 << n):
        k = bin(mask).count('1')
        if k < 2 or k > n - 1: continue
        ok = True
        for w in range(n):
            if mask >> w & 1: continue
            x = out[w] & mask
            if x and x != mask: ok = False; break
        if ok: return 0
    return 1
def arcs_of(out):
    return ''.join('1' if out[i] >> j & 1 else '0' for i in range(n) for j in range(i + 1, n))
tours = []
for t in range(3000):
    out = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.5: out[i] |= 1 << j
            else: out[j] |= 1 << i
    tours.append(out)
for t in range(3000):          # invariant under g = (3 4 5)(6 7 8) or (0 1 2)(3 4 5)(6 7 8)
    g = [0, 1, 2, 4, 5, 3, 7, 8, 6] if t % 2 else [1, 2, 0, 4, 5, 3, 7, 8, 6]
    out = [0] * n; done = set()
    for i in range(n):
        for j in range(n):
            if i == j or (i, j) in done or (j, i) in done: continue
            b = random.random() < 0.5
            a, c = i, j
            while (a, c) not in done:
                done.add((a, c))
                if b: out[a] |= 1 << c
                else: out[c] |= 1 << a
                a, c = g[a], g[c]
    tours.append(out)
inp = ''.join(f"x arcs={arcs_of(o)}\n" for o in tours)
res = subprocess.run(['./audit_reclass', str(n)], input=inp, capture_output=True, text=True).stdout.split('\n')
mism = 0; nprime = 0
for o, line in zip(tours, res):
    cp = int(line.split('prime=')[1]); bp = brute_prime(o)
    nprime += bp; mism += (cp != bp)
print(f"n={n}: {len(tours)} tournaments ({len(tours)//2} random, {len(tours)//2} invariant), brute-force prime {nprime}, mismatches {mism}")
