# Refined (6-adic) classes of the negative cycle points: n = p (mod 3^a), backward one turn around p's cycle.
# m(n) = (2^L n - c)/3^a with T^L(m) = n, fixed point p; ratio 2^L/3^a < 1 for the negative cycles.
def T(x): return x // 2 if x % 2 == 0 else (3 * x + 1) // 2
def Tn(x, k):
    for _ in range(k): x = T(x)
    return x
import random
for p, L, a in [(-1, 1, 1), (-5, 3, 2), (-17, 11, 7)]:
    # c from T^L(p) = p with parity word of p:  (3^a p + c)/2^L = p
    c = p * 2**L - 3**a * p
    ok = True; worst = None
    for trial in range(2000):
        j = random.randint(1, 10**25)
        n = 3**a * j + p                    # n = p mod 3^a, positive
        m = 2**L * j + p                    # = (2^L n - c)/3^a
        assert (2**L * n - c) % 3**a == 0 and (2**L * n - c) // 3**a == m
        if not (0 < m < n and Tn(m, L) == n): ok = False; worst = n; break
    # small n too
    for j in range(1, 2000):
        n = 3**a * j + p; m = 2**L * j + p
        if not (0 < m < n and Tn(m, L) == n): ok = False; worst = n; break
    print(f"p={p}: class n = {p} mod 3^{a} = {p % 3**a} mod {3**a}: m(n) = (2^{L} n {'+' if -c>=0 else '-'} {abs(c)})/3^{a}, ratio 2^{L}/3^{a} = {2**L/3**a:.4f}, fixed point {p}; valid on all positive n tested: {ok}")
