#!/usr/bin/env python3
"""cover.py (M-5, M-6): the metagraph selector at order 8 -> 9.  Each class P of order 8 has the image I(P) of its
quadratic form H_P(x) = h + <B - e, x> - x^T M x over the cube (THM-5266 (A)): the counts of its 256 one-vertex extensions.
The new counts of order 9 (first born at 9) must be covered by these images (completeness of the extension grammar).
Greedy set cover by images; the best single class; whether the first creators of order 8 alone cover; where the
maximal-count classes rank."""
import numpy as np
from collections import defaultdict
P = np.load('profiles8.npz')
first = {}
for line in open('counts.txt'):
    n, h = map(int, line.split())
    if h not in first or n < first[h]: first[h] = n
new9 = {h for h, n in first.items() if n == 9}; print(f"new counts at order 9: {len(new9)}")
n = 8; A, h, e, B, M = P['A8'], P['h8'], P['e8'], P['B8'], P['M8']; k = len(h); N = 1 << n
X = np.array([[(m >> x) & 1 for x in range(n)] for m in range(N)], dtype=np.int64)
lin = X @ (B - e).T; quad = np.einsum('mi,kij,mj->mk', X, M, X); val = h[None, :] + lin - quad
images = [set(val[:, i].tolist()) & new9 for i in range(k)]
sizes = [len(s) for s in images]; best = int(np.argmax(sizes))
print(f"best single class: index {best}, H = {h[best]}, covers {sizes[best]} of {len(new9)} new counts ({sizes[best]/len(new9):.3f}); top five (H, covered): {[(int(h[i]), sizes[i]) for i in sorted(range(k), key=lambda i: -sizes[i])[:5]]}")
# first creators of order 8
fc8 = [i for i in range(k) if first.get(int(h[i])) == 8]; cov_fc = set().union(*[images[i] for i in fc8]) if fc8 else set()
print(f"(M-6) first creators of order 8: {len(fc8)}; their images cover {len(cov_fc)} of {len(new9)} new counts; missing: {sorted(new9 - cov_fc)[:20]}")
# greedy cover
rem = set(new9); chosen = []
while rem:
    i = max(range(k), key=lambda i: len(images[i] & rem)); g = len(images[i] & rem)
    if g == 0: break
    chosen.append((i, int(h[i]), g)); rem -= images[i]
print(f"(M-5) greedy cover of the {len(new9)} new counts by order-8 images: {len(chosen)} classes; first ten (index, H, newly covered): {chosen[:10]}; uncovered: {len(rem)}")
hmax = int(h.max()); top661 = [i for i in range(k) if h[i] == hmax]
print(f"classes of maximal count {hmax}: {len(top661)}; their image sizes among new counts: {sorted(sizes[i] for i in top661)}; rank of the best of them: {sorted(sizes, reverse=True).index(max(sizes[i] for i in top661)) + 1}")
