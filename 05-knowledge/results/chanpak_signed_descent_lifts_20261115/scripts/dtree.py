"""The signed descent D on odd numbers (six_tutte_collatz_20261116/COUNTS.md) as a tree.
   D(h) = (h+1)/2 (h = 1 mod 4);  (3h+3)/4 (h = 3 mod 8);  (3h-1)/4 (h = 7 mod 8);  fixed points 1 and 3.
   Basins of the two non-counts 7 and 21, depths, and the slowest descents."""
import sys, math
X = int(sys.argv[1]) if len(sys.argv) > 1 else 10 ** 7
def D(h):
    if h % 4 == 1: return (h + 1) // 2
    if h % 8 == 3: return (3 * h + 3) // 4
    return (3 * h - 1) // 4
assert all(D(h) % 2 == 1 and D(h) < h for h in range(5, 20001, 2)) and D(1) == 1 and D(3) == 3
pre = lambda d: sorted({2 * d - 1} | ({(4 * d - 3) // 3} if d % 6 == 3 else set()) | ({(4 * d + 1) // 3} if d % 6 == 5 else set()))
print("preimages of 7:", [h for h in pre(7) if D(h) == 7], " of 21:", [h for h in pre(21) if D(h) == 21], " of 5:", pre(5), " of 3:", pre(3))
N = X // 2 + 1                       # index i <-> odd h = 2i+1
tag = bytearray(N)                   # 1: orbit passes through 7, 2: through 21, 0: neither
depth = [0] * N
tag[3] = 1; tag[10] = 2
for i in range(2, N):
    h = 2 * i + 1; d = D(h); j = (d - 1) // 2
    depth[i] = depth[j] + 1
    if i not in (3, 10): tag[i] = tag[j]
for B in [10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, X]:
    n = B // 2; c7 = sum(1 for i in range(2, n + 1) if tag[i] == 1); c21 = sum(1 for i in range(2, n + 1) if tag[i] == 2)
    tot = n - 1
    print("odd 5 <= h <= %8d: through 7: %.5f ; through 21: %.5f ; through 7 or 21: %.5f ; mean depth / ln h: %.4f" %
          (B, c7 / tot, c21 / tot, (c7 + c21) / tot, sum(depth[i] / math.log(2 * i + 1) for i in range(n // 2, n + 1)) / (n - n // 2 + 1)))
mu = 0.5 * math.log(2) + 0.5 * math.log(4 / 3)
print("heuristic: steps of size ln 2 (prob 1/2) and ln(4/3) (prob 1/2): depth ~ ln h / %.4f = %.4f ln h" % (mu, 1 / mu))
# first-level split below 5: everything >= 7 passes through 7 or 9
c9 = bytearray(N); c9[4] = 1
for i in range(5, N):
    c9[i] = c9[(D(2 * i + 1) - 1) // 2]
n = N - 1
print("through 9 (the other child of 5): %.5f ; through 7: %.5f" % (sum(c9[5:]) / (n - 4), sum(1 for i in range(5, n + 1) if tag[i] == 1) / (n - 4)))
# records of depth (slowest descents) and their classes
rec = []; best = -1
for i in range(1, N):
    if depth[i] > best: best = depth[i]; rec.append((2 * i + 1, depth[i]))
print("depth records (h, depth):", rec[:40])
print("records in the basin of 7 or 21:", [(h, dp, {1: 7, 2: 21}.get(tag[(h - 1) // 2], '-')) for h, dp in rec[:40]])
# the slowest branch: always the rising branches (h = 3 mod 4 at every step)
chain = [3]
for _ in range(14):
    d = chain[-1]; c = [h for h in pre(d) if h != 2 * d - 1 and h != d]
    if not c: break
    chain.append(c[0])
print("a chain of rising inverses from 3:", chain)
