"""Helpers (digraph6 parsing, subset DP, Hamiltonian cycles, Phi) and the forward-chord analysis of the digit-1 ready tournaments through the entry attachment Y = Phi(T, a):
label Y's unique Hamiltonian cycle a = v1 -> v2 (= r) -> ... -> vn -> v1 and record the forward chords
(i, j), j >= i + 2, v_i -> v_j (all other chords are backward). Then tabulate H(T) against the chord set."""
import sys, itertools
from collections import defaultdict, Counter

def d6_to_rows(s):
    s = s.strip()
    assert s[0] == '&'
    s = s[1:]
    n = ord(s[0]) - 63
    bits = []
    for ch in s[1:]:
        v = ord(ch) - 63
        bits.extend((v >> k) & 1 for k in range(5, -1, -1))
    rows = [0] * n
    k = 0
    for i in range(n):
        for j in range(n):
            if bits[k]:
                rows[i] |= 1 << j
            k += 1
    return rows

def ham_paths_ending(rows, verts):
    verts = list(verts); idx = {v: k for k, v in enumerate(verts)}; m = len(verts)
    adj = [sum(1 << idx[u] for u in verts if rows[v] >> u & 1) for v in verts]
    full = (1 << m) - 1
    dp = [[0] * m for _ in range(1 << m)]
    for k in range(m):
        dp[1 << k][k] = 1
    for S in range(1, 1 << m):
        row = dp[S]
        for k in range(m):
            c = row[k]
            if not c: continue
            nxt = adj[k] & ~S
            while nxt:
                b = nxt & -nxt; nxt ^= b; j = b.bit_length() - 1
                dp[S | b][j] += c
    return {verts[k]: dp[full][k] for k in range(m)}

def ham_cycles(rows):
    """list of Hamiltonian cycles (as vertex sequences starting at 0)."""
    n = len(rows); res = []
    def rec(path, used):
        v = path[-1]
        if len(path) == n:
            if rows[v] >> 0 & 1: res.append(list(path))
            return
        m = rows[v] & ~used
        while m:
            b = m & -m; m ^= b; u = b.bit_length() - 1
            path.append(u); rec(path, used | b); path.pop()
    rec([0], 1)
    return res

def phi(rows, a):
    """reverse every arc not incident to a."""
    n = len(rows); out = [0] * n
    for i in range(n):
        for j in range(n):
            if i == j: continue
            if i == a or j == a:
                if rows[i] >> j & 1: out[i] |= 1 << j
            else:
                if rows[j] >> i & 1: out[i] |= 1 << j
    return out

def analyse(order, path):
    stats = defaultdict(list)
    chordsets = Counter()
    for line in open(path):
        parts = line.split()
        H = int(parts[-2]) if len(parts) > 2 else int(parts[0]); d6 = parts[-1]
        if not d6.startswith('&'): continue
        rows = d6_to_rows(d6); n = len(rows)
        for a in range(n):
            outs = [j for j in range(n) if rows[a] >> j & 1]
            if len(outs) != 1: continue
            r = outs[0]
            digit = ham_paths_ending(rows, [v for v in range(n) if v != a])[r]
            if digit != 1: continue
            Y = phi(rows, a)
            # relabel so that a = 0
            cyc = ham_cycles(Y) if a == 0 else None
            if a != 0:
                perm = [a] + [v for v in range(n) if v != a]
                inv = {v: k for k, v in enumerate(perm)}
                Y2 = [0] * n
                for i in range(n):
                    for j in range(n):
                        if Y[i] >> j & 1: Y2[inv[i]] |= 1 << inv[j]
                Y = Y2
                cyc = ham_cycles(Y)
            assert len(cyc) == 1, (H, d6, a, len(cyc))
            c = cyc[0]  # starts at a=0; c[1] should be r
            pos = {v: k for k, v in enumerate(c)}
            fwd = []
            for i in range(n):
                for j in range(i + 2, n):
                    if Y[c[i]] >> c[j] & 1:
                        fwd.append((i + 1, j + 1))
            key = tuple(sorted(fwd))
            stats[key].append(H)
            chordsets[key] += 1
    return stats

if __name__ == '__main__':
    order = int(sys.argv[1]); path = sys.argv[2]
    stats = analyse(order, path)
    print(f"order {order}: {len(stats)} distinct forward-chord sets over the digit-1 pairs")
    for key in sorted(stats, key=lambda k: (len(k), k)):
        hs = sorted(set(stats[key]))
        print(f"  chords {list(key)}: H values {hs} (pairs {len(stats[key])})")
