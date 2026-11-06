import math
def bound(n, q):
    best = float('inf')
    for a in range(2, n):
        for L in range(1, n - a):
            l = n - a - L
            if l < 1: continue
            R = 2 * (2**L / q + 1) * math.sqrt((3 + 2**(3 - a) * q) * (3 + 2**(2 - l) * q))
            best = min(best, 2**(-L / 2) * math.sqrt(1 + R))
    return best
w1 = w2 = 0.0
for n in range(9, 91, 3):
    for e in [x / 4 for x in range(4, 4 * n)]:
        q = 2**e
        b = bound(n, q)
        if q <= 2**((n - 3) / 3): w1 = max(w1, b * math.sqrt(q))
        if 2**((n - 6) / 3) <= q <= 2**(n - 6): w2 = max(w2, b / (q / 2**n)**0.25)
print(f'max over tested (n, q) of bound * sqrt(q) for 2 <= q <= 2^((n-3)/3): {w1:.3f}  (claimed 6.5)')
print(f'max of bound / (q/2^n)^(1/4) for 2^((n-6)/3) <= q <= 2^(n-6): {w2:.3f}  (claimed 11.3)')
# minor-arc corollary: with Q = 2^((1-eps) n) and q > 2^(eps n): sup of bound * 2^(eps n / 4)
for eps in (0.1, 0.2, 1/3):
    worst = 0.0
    for n in range(30, 121, 6):
        for e in [x / 4 for x in range(int(4 * eps * n) + 1, int(4 * (1 - eps) * n))]:
            worst = max(worst, bound(n, 2**e) * 2**(eps * n / 4))
    print(f'eps={eps:.3f}: max of bound * 2^(eps n/4) over n in [30,120], 2^(eps n) < q <= 2^((1-eps) n): {worst:.3f}')
