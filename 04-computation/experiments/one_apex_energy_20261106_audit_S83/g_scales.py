# Independent audit of THM-4887, part 7 (upgrade check): (D1) at general m, Delta = 2^m, P = 2^(n-m-1) >= 8C.
# For m >= (n-1)/3 the P C^2 term dominates P^2, so (D1) reads E_n(2^m) <= O(n - 2m) 2^(n+m): linear smoothed energy,
# up to a factor O(n), at EVERY dyadic scale 2^m with n/3 <= m <= (n-4)/2.  Compare with exact E_n(2^m).
import math
import numpy as np

def level(n):
    cur = np.array([1], dtype=np.int64)
    for j in range(1, n):
        u = np.arange(1 << (j - 1), 1 << j, dtype=np.int64)
        nxt = np.empty(2 * len(cur), dtype=np.int64)
        nxt[0::2] = cur + u + (u & 1); nxt[1::2] = cur
        cur = nxt
    return cur

def D1(n, m):
    C = 2.0 ** m; P = 2.0 ** (n - m - 1); lB = math.log(P / C)
    assert P >= 8 * C
    return P * C + 10 * P * P + P * C * C * (35.74 + 12 * lB) + 3 * P * C * (1 + lB)

print(" n  m   E_n(2^m)/2^(n+m)   (D1)/2^(n+m)   (D1)/2^(n+m) formula 17.87+6 ln2 (n-2m-1)+2.5*2^(n-3m)+...")
for n in (18, 20, 22, 24):
    x = np.sort(level(n))
    for m in range(n // 3, (n - 4) // 2 + 1):
        D = 1 << m
        E = int((np.searchsorted(x, x + D, side="right") - np.searchsorted(x, x - D, side="left")).sum())
        b = D1(n, m)
        approx = 17.87 + 6 * math.log(2) * (n - 2 * m - 1) + 2.5 * 2.0 ** (n - 3 * m)
        print(f"{n:2d} {m:2d}   {E / 2 ** (n + m):10.4f}        {b / 2 ** (n + m):8.2f}        {approx:8.2f}")

# U1 constants: (D1)/2^(n+m) at the top admissible scale m = floor((n-4)/2), and the max over all
# floor(n/3) <= m <= floor((n-4)/2), for n = 12..400.
tops = []; allmax = 0.0
for n in range(12, 401):
    mt = (n - 4) // 2
    tops.append((n, D1(n, mt) / 2 ** (n + mt)))
    for m in range(n // 3, mt + 1):
        allmax = max(allmax, D1(n, m) / 2 ** (n + m) / max(1, (n - 2 * m)))
print("top scale m = floor((n-4)/2): (D1)/2^(n+m) for n = 12..21:", [round(v, 2) for _, v in tops[:10]])
print(f"  max over n = 12..400: {max(v for _, v in tops):.2f} (n = {max(tops, key=lambda t: t[1])[0]}); "
      f"max over n >= 20: {max(v for n, v in tops if n >= 20):.2f}; at n = 400 (even): "
      f"{tops[-1][1]:.2f}, at n = 399 (odd): {tops[-2][1]:.2f}")
print(f"  max over n = 12..400 and floor(n/3) <= m <= floor((n-4)/2) of (D1)/(2^(n+m) (n - 2m)) = {allmax:.2f}")
