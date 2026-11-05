"""aud_uh5_scores.py -- independent audit code (uh_deletion_law_20261104_audit_S10): the uniquely Hamiltonian
5-tournaments up to isomorphism, with score sequences and numbers of sole exits / sole entries, to check THM-4874's
description of W_5 ("the uniquely Hamiltonian 5-tournament with scores 1, 2, 2, 2, 3 and one sole exit")."""
from itertools import permutations

n = 5
pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]


def hc(rows):
    cnt = 0
    for p in permutations(range(1, n)):
        cyc = (0,) + p
        if all((rows[cyc[i]] >> cyc[(i + 1) % n]) & 1 for i in range(n)):
            cnt += 1
    return cnt


def key(rows):
    best = None
    for p in permutations(range(n)):
        code = tuple(sum(1 << p[w] for w in range(n) if (rows[v] >> w) & 1) for v in sorted(range(n), key=lambda x: p[x]))
        if best is None or code < best:
            best = code
    return best


seen = {}
for bits in range(1 << len(pairs)):
    rows = [0] * n
    for b, (i, j) in enumerate(pairs):
        if (bits >> b) & 1:
            rows[i] |= 1 << j
        else:
            rows[j] |= 1 << i
    if hc(rows) != 1:
        continue
    k = key(rows)
    if k not in seen:
        seen[k] = rows
print(f"uniquely Hamiltonian 5-tournament classes: {len(seen)}")
for k, rows in seen.items():
    sc = sorted(bin(r).count('1') for r in rows)
    se = sum(1 for r in rows if bin(r).count('1') == 1)
    sn = sum(1 for v in range(n) if sum(1 for w in range(n) if (rows[w] >> v) & 1) == 1)
    print(f"  scores {sc}, sole exits {se}, sole entries {sn}, out-nbr lists {[[w for w in range(n) if (r >> w) & 1] for r in rows]}")
