"""aud_at_title_prefix.py -- independent audit code (uh_deletion_law_20261104_audit_S10).
THM-4874's title says "Every core block (AT_m, {m-1, m}) gives exactly three deletions"; (5) assumes a uniquely
Hamiltonian or one-vertex prefix. This checks that the prefix hypothesis is needed: with a rooted prefix (Q, z), z a sole
exit of Q and hc(Q) = 2, the one-arc join on (AT_m, {m-1, m}) has hc = 2(2^(m-3) + 1) and NO uniquely Hamiltonian deletion
(every block value is doubled)."""
from itertools import product
from aud_structure_check import hcycles

def at_rows(m):
    r = [0] * m
    for i in range(m):
        for j in range(i + 1, m):
            if (i, j) == (0, m - 1):
                r[m - 1] |= 1
            else:
                r[i] |= 1 << j
    return r

# find a rooted (Q, z) of order 5 with hc(Q) = 2 and z of out-degree 1
q = 5
pairs = [(i, j) for i in range(q) for j in range(i + 1, q)]
found = None
for bits in range(1 << len(pairs)):
    Q = [0] * q
    for b, (i, j) in enumerate(pairs):
        if (bits >> b) & 1:
            Q[i] |= 1 << j
        else:
            Q[j] |= 1 << i
    zs = [z for z in range(q) if bin(Q[z]).count("1") == 1]
    if zs and hcycles(Q, list(range(q))) == 2:
        found = (Q, zs[0])
        break
Q, z = found
l = Q[z].bit_length() - 1
print("prefix Q (out-nbr lists):", [[w for w in range(q) if (Q[v] >> w) & 1] for v in range(q)], "z =", z, "l =", l)
for m in (4, 5, 6):
    W = at_rows(m)
    R = [l] + [x for x in range(q) if x not in (z, l)]
    r = len(R)
    n = r + m + 1                      # chain p = 0
    u = n - 1
    rows = [0] * n
    idx = {x: i for i, x in enumerate(R)}
    for x in R:
        for y in R:
            if (Q[x] >> y) & 1:
                rows[idx[x]] |= 1 << idx[y]
    for i in range(r):
        for y in range(r, n):
            if not (i == 0 and y == u):
                rows[i] |= 1 << y
    rows[u] |= 1 << 0
    for i in range(m):
        for j in range(m):
            if (W[i] >> j) & 1:
                rows[r + i] |= 1 << (r + j)
        if i >= m - 2:
            rows[r + i] |= 1 << u
        else:
            rows[u] |= 1 << (r + i)
    h = hcycles(rows, list(range(n)))
    dels = [hcycles(rows, [x for x in range(n) if x != v]) for v in range(n)]
    print(f"m={m}: n={n}, min out-degree {min(bin(x).count('1') for x in rows)}, hc(T) = {h}, "
          f"#(hc(T-v) = 1) = {dels.count(1)}, block values {dels[r:r + m]}")
