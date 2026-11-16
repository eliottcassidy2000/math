#!/usr/bin/env python3
"""collatz_tree.py BOUND: the Collatz map in count coordinates h = 2x + 1: T(h) = 3h if h = 3 (mod 4), (h+1)/2 if h = 1 (mod 4).
Lists the odd h <= BOUND whose forward orbit meets {7, 21} (the backward tree of the forbidden pair), checks whether they
lie on 2-adic rays c*2^k + 1, and whether each (other than 7, 21) is a one-apex count 1 + p(u) (p the inversion statistic
of the apex packet: p(u) = p(w) if u odd else p(w) + w + (w & 1), w = u >> 1)."""
import sys
B = int(sys.argv[1])
def T(h): return 3 * h if h % 4 == 3 else (h + 1) // 2
M = 24; p = [0] * (1 << M)
for u in range(1, 1 << M):
    w = u >> 1; p[u] = p[w] if u & 1 else p[w] + w + (w & 1)
one = set(1 + v for v in p)
tree = []
for h in range(1, B + 1, 2):
    x = h; steps = 0
    while x > 21 and steps < 10000:
        x = T(x); steps += 1
        if x > 10**12: break
    if x in (7, 21): tree.append(h)
print(f"odd h <= {B} whose orbit meets {{7, 21}}: {len(tree)}; first 30: {tree[:30]}")
def oddpart(n):
    while n % 2 == 0: n //= 2
    return n
rays = {}
for h in tree:
    c = oddpart(h - 1) if h > 1 else 0
    rays.setdefault(c, []).append(h)
print("rays c*2^k + 1 (c = odd part of h - 1), first members:", {c: v[:5] for c, v in sorted(rays.items())[:12]}, "... number of distinct c:", len(rays))
miss = [h for h in tree if h not in (7, 21) and h not in one]
print(f"members not one-apex counts (excluding 7, 21): {len(miss)} {miss[:20]}")
print("one-apex membership of the rays 3*2^k+1, 5*2^k+1, 13*2^k+1 for k = 2..18:")
for c in (3, 5, 13):
    print(f"  c = {c}:", [(k, c * 2**k + 1 in one) for k in range(2, 19) if c * 2**k + 1 < (1 << 24)])
