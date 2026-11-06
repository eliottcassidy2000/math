import numpy as np, cmath, math
def T(q, a):
    n = 2 * q; M = np.zeros((n, n), dtype=complex)
    for P in range(q):
        for l in (0, 1):
            i = 2 * P + l
            for b in (0, 1):
                add = (P + (1 if l == 1 else 0)) if b == 0 else 0
                j = 2 * ((2 * P + b) % q) + b
                M[j, i] += cmath.exp(2j * math.pi * a * add / q)
    return M
res = []
for q in range(3, 121):
    for a in range(1, q):
        if math.gcd(a, q) != 1: continue
        if 2 * a == q: continue
        r = max(abs(np.linalg.eigvals(T(q, a)))) / 2
        res.append((r, q, a))
res.sort(reverse=True)
print("worst decay rates |lambda|/2 over a/q, q <= 120 (theta not in (1/2)Z):")
for r, q, a in res[:25]: print(f"  {a}/{q}: {r:.4f}")
print("overall max:", res[0])
