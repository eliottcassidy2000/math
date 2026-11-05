"""aud_at_block_core.py -- independent audit code (uh_deletion_law_20261104_audit_S10).
For m = 4..10: (AT_m, {m-1, m}) satisfies the block conditions of S75's Lemma B4 (all 2^m - 1 subsets Z), has type
2^(m-3) + 1, and in the smallest one-arc joins built on it (prefix {l}; chain p = 0, 1, 2) the minimal one-arc-leaving
set of T is exactly D (so the block is a terminal-core block, as THM-4874 (5) and its title say)."""
import sys
sys.setrecursionlimit(10000)
from aud_structure_check import block_conditions, end_count, leaving

def at_rows(m):          # vertices 0..m-1 stand for 1..m
    r = [0] * m
    for i in range(m):
        for j in range(i + 1, m):
            if (i, j) == (0, m - 1):
                r[m - 1] |= 1
            else:
                r[i] |= 1 << j
    return r

allok = True
for m in range(4, 11):
    W = at_rows(m)
    A = {m - 2, m - 1}
    bc = block_conditions(W, list(range(m)), A)
    typ = sum(end_count(W, list(range(m)), x) for x in A)
    line = f"m={m}: block conditions {bc}, type {typ} (2^(m-3)+1 = {2 ** (m - 3) + 1})"
    for p in (0, 1, 2):
        # T: l = 0, chain 1..p, W = p+1..p+m, u = p+m+1;  u -> w_2 (free arc chosen u -> w_2)
        n = 1 + p + m + 1
        l, u = 0, n - 1
        ch = list(range(1, 1 + p))
        Wv = list(range(1 + p, 1 + p + m))
        rows = [0] * n
        for y in range(1, n):
            if y != u:
                rows[l] |= 1 << y
        rows[u] |= 1 << l
        for a in range(p):
            for b in range(a + 1, p):
                rows[ch[a]] |= 1 << ch[b]
            for y in Wv:
                rows[ch[a]] |= 1 << y
        for i in range(m):
            for j in range(m):
                if (W[i] >> j) & 1:
                    rows[Wv[i]] |= 1 << Wv[j]
        for i in range(m):
            if i in A:
                rows[Wv[i]] |= 1 << u
            else:
                rows[u] |= 1 << Wv[i]
        for t in range(p):
            rows[u] |= 1 << ch[t]          # u -> w_1 (forced) and u -> w_2 (free arc, one choice)
        full = (1 << n) - 1
        one_arc = [D for D in range(1, full) if leaving(rows, n, D) == 1]
        minimal = [D for D in one_arc if not any((E & D) == E and E != D for E in one_arc)]
        Dm = full & ~1                      # D = everything but l
        ok = minimal == [Dm] and min(bin(r).count("1") for r in rows) >= 2
        allok &= ok and bc and typ == 2 ** (m - 3) + 1
        line += f"; p={p}: terminal core = D: {minimal == [Dm]}"
    print(line)
print("ALL OK:", allok)
