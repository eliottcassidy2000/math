import sys
from fractions import Fraction as Fr
from cutvertex_chain import brute, components
from check_minimiser import dec
s = sys.argv[1]; A = dec(s); n = len(A); H, cnt = brute(A)
for z in range(n):
    comps = components(A, [x for x in range(n) if x != z])
    if len(comps) > 1:
        info = []
        for c in comps:
            sub = [[A[u][v] for v in c] for u in c]
            h, _ = brute(sub) if len(c) > 1 else (1, None)
            info.append((c, h))
        print(f"cut vertex {z}: out-nbrs {[x for x in range(n) if A[z][x]]}; blocks (vertices, H): {info}")
pairs = sorted(((min(Fr(cnt[x][y], H), 1 - Fr(cnt[x][y], H)), x, y) for x in range(n) for y in range(x + 1, n)), reverse=True)[:6]
print("most balanced pairs:", [(x, y, str(p)) for p, x, y in pairs])
# which vertex deletions keep it strong, and delta of T - v
for v in range(n):
    rest = [x for x in range(n) if x != v]
    sub = [[A[a][b] for b in rest] for a in rest]
    if len(components(sub, range(n - 1))) == 1:
        h, c2 = brute(sub)
        d = max(min(Fr(c2[x][y], h), 1 - Fr(c2[x][y], h)) for x in range(n - 1) for y in range(x + 1, n - 1))
        print(f"  T - {v} strong: H = {h}, delta = {d} = {float(d):.4f}")
