#!/usr/bin/env python3
"""width.py -- ribbon width of a tournament = the least w such that some vertex ordering has every
back arc of span <= w (bandwidth-type DP over subsets).  Classifies all tournaments of orders 3..8
by width, and locates the strong realisers of H = 47 (orders <= 8) and their widths."""
import sys, subprocess, itertools
sys.path.insert(0, '../chanpak_spine_20261104')
from ribbon import ham_count_bruteforce, is_strong

def gentourng(n):
    out = subprocess.run(['gentourng', '-q', str(n)], capture_output=True, text=True).stdout.split()
    res = []
    for line in out:
        bits = line.strip()
        adj = [0] * n
        k = 0
        for i in range(n):
            for j in range(i + 1, n):
                if bits[k] == '1':
                    adj[i] |= 1 << j
                else:
                    adj[j] |= 1 << i
                k += 1
        res.append(adj)
    return res

def width_le(adj, w):
    """is there an ordering with all back arcs of span <= w?  DP: state = (placed set, last w placed)."""
    n = len(adj)
    full = (1 << n) - 1
    # dominance: placing v after set S with window W (tuple of last <= w vertices): every u in S \ W must have u -> v
    from functools import lru_cache
    @lru_cache(maxsize=None)
    def go(S, W):
        if S == full:
            return True
        outside = S
        for u in W:
            outside &= ~(1 << u)
        for v in range(n):
            if S >> v & 1:
                continue
            # all u in outside must beat v
            ok = True
            o = outside
            while o:
                u = (o & -o).bit_length() - 1
                o &= o - 1
                if not (adj[u] >> v & 1):
                    ok = False; break
            if not ok:
                continue
            W2 = (W + (v,))[-w:] if w > 0 else ()
            if go(S | (1 << v), W2):
                return True
        return False
    return go(0, ())

def width(adj):
    n = len(adj)
    for w in range(0, n):
        if width_le(adj, w):
            return w
    return n

if __name__ == "__main__":
    for n in range(3, 9):
        ts = gentourng(n)
        from collections import Counter
        cnt = Counter(); cnt_strong = Counter()
        real47 = []
        for adj in ts:
            wd = width(adj)
            cnt[wd] += 1
            st = is_strong(adj)
            if st:
                cnt_strong[wd] += 1
            if n >= 6:
                h = ham_count_bruteforce(adj)
                if h == 47:
                    real47.append((wd, st, adj))
        print(f"order {n}: {len(ts)} classes; width distribution {dict(sorted(cnt.items()))}; strong: {dict(sorted(cnt_strong.items()))}")
        if real47:
            print(f"   H=47 realisers: {[(wd, st) for wd, st, adj in real47]}")
