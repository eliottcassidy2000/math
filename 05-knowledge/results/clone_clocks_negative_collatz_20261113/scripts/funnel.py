"""Nested one-start funnels: brute-force tournaments and the inherited six-state compiler.

Conventions (FUNNELS F3 / ENTRANCES of the repo):
  core U rooted at w; tail v_1..v_t with U -> tail, v_i -> v_j for i<j;
  new root r with r -> w, (U - w) -> r, and v_i -> r iff bit_i = 1 (r -> v_i iff bit_i = 0).
  A nonempty word ends in 1 (guard); the empty word is allowed on a nontrivial core.
Ready child: append z beaten only by the last root (root -> z, z -> everything else); delta_z = 1.
"""
from functools import lru_cache
import itertools, random

# ---------- brute force ----------
def hp_stats(adj, n):
    """Return (H, s[v], e[v]) for tournament given by out-neighbour bitmasks adj[0..n-1]."""
    N = 1 << n
    # f[S][v] = # paths covering S ending at v ; g[S][v] = # paths covering S starting at v
    f = [[0]*n for _ in range(N)]
    for v in range(n): f[1<<v][v] = 1
    for S in range(1, N):
        row = f[S]
        for v in range(n):
            c = row[v]
            if not c: continue
            out = adj[v] & ~S
            while out:
                u = (out & -out).bit_length() - 1
                out &= out - 1
                f[S | (1<<u)][u] += c
    full = N - 1
    e = [f[full][v] for v in range(n)]
    H = sum(e)
    # starts: reverse graph
    radj = [0]*n
    for v in range(n):
        for u in range(n):
            if adj[v] >> u & 1: radj[u] |= 1 << v
    g = [[0]*n for _ in range(N)]
    for v in range(n): g[1<<v][v] = 1
    for S in range(1, N):
        row = g[S]
        for v in range(n):
            c = row[v]
            if not c: continue
            out = radj[v] & ~S
            while out:
                u = (out & -out).bit_length() - 1
                out &= out - 1
                g[S | (1<<u)][u] += c
    s = [g[full][v] for v in range(n)]
    return H, s, e

def H_only(adj, n):
    N = 1 << n
    f = [[0]*n for _ in range(N)]
    for v in range(n): f[1<<v][v] = 1
    for S in range(1, N):
        row = f[S]
        for v in range(n):
            c = row[v]
            if not c: continue
            out = adj[v] & ~S
            while out:
                u = (out & -out).bit_length() - 1
                out &= out - 1
                f[S | (1<<u)][u] += c
    return sum(f[N-1])

def build_nested(words):
    """words: list of bit-strings applied to a rooted singleton. Returns (adj, n, root)."""
    adj = [0]; n = 1; root = 0
    for w in words:
        U = list(range(n))
        tail = list(range(n, n + len(w)))
        r = n + len(w)
        adj = adj + [0]*(len(w) + 1)
        # U -> tail
        for u in U:
            for v in tail: adj[u] |= 1 << v
        # tail transitive
        for i, v in enumerate(tail):
            for v2 in tail[i+1:]: adj[v] |= 1 << v2
        # root
        adj[r] |= 1 << root
        for u in U:
            if u != root: adj[u] |= 1 << r
        for i, v in enumerate(tail):
            if w[i] == '1': adj[v] |= 1 << r
            else: adj[r] |= 1 << v
        n = r + 1; root = r
    return adj, n, root

def append_ready(adj, n, root):
    """append z beaten only by root."""
    adj = adj + [0]
    z = n
    adj[root] |= 1 << z
    for u in range(n):
        if u != root: adj[z] |= 1 << u
    return adj, n + 1, z

def is_tournament(adj, n):
    for i in range(n):
        if adj[i] >> i & 1: return False
        for j in range(i+1, n):
            if ((adj[i] >> j) & 1) + ((adj[j] >> i) & 1) != 1: return False
    return True

# ---------- six-state compiler ----------
def word_polys(w, m):
    t = len(w)
    L = (m-1)**t + sum(m**i * (m-1)**(t-1-i) for i in range(t) if w[i] == '1')
    D = (m-1)**t + sum((m-1)**i * m**(t-1-i) for i in range(t) if w[i] == '0')
    # B via automaton: b in {E,G,X(bad)}, a in {0,1}
    st = {('E', 0): 1}
    for bit in w:
        ns = {}
        for (b, a), c in st.items():
            # ordinary bin
            ns[(b, a)] = ns.get((b, a), 0) + c*(m-1)
            # before bin
            nb = 'G' if bit == '1' else 'X'
            ns[(nb, a)] = ns.get((nb, a), 0) + c
            # after bin
            if a == 1 or bit == '0':
                ns[(b, 1)] = ns.get((b, 1), 0) + c
        st = ns
    B = sum(c for (b, a), c in st.items() if b != 'X')
    return L, D, B

def step6(state, w):
    h, e, s2, p2, e2, s3 = state
    t = len(w)
    L1, D1, B1 = word_polys(w, 1)
    L2, D2, B2 = word_polys(w, 2)
    L3, D3, B3 = word_polys(w, 3)
    z = 1 if t == 0 else 0
    return (L2*s2 + B1*h - D1*e,
            L1*h - z*e,
            2**t*s2 + D2*h,
            2*(L3*s3 + B2*p2 - D2*e2),
            L2*p2 - e2,
            3**t*s3 + D3*p2)

SINGLETON = (1, 1, 1, 2, 1, 1)

def ready_count(words):
    st = SINGLETON
    for w in words: st = step6(st, w)
    h, e, s2, p2, e2, s3 = st
    return h + e2 - 1, st

if __name__ == '__main__':
    rng = random.Random(1)
    bad = 0; tested = 0
    for trial in range(300):
        m = rng.randint(1, 4)
        words = []
        for i in range(m):
            t = rng.randint(0, 3)
            if t == 0 and i == 0: t = 1
            w = ''.join(rng.choice('01') for _ in range(t-1)) + ('1' if t else '')
            words.append(w)
        adj, n, root = build_nested(words)
        if n > 13: continue
        A2, n2, z = append_ready(adj, n, root)
        assert is_tournament(A2, n2)
        H, s, e = hp_stats(A2, n2)
        Hp, sp, ep = hp_stats(adj, n)
        pred, st = ready_count(words)
        tested += 1
        ok = (pred == H) and (sp[root] == 1) and (st[0] == Hp) and (st[1] == ep[root])
        if not ok:
            bad += 1; print('MISMATCH', words, H, pred, sp[root], st, Hp, ep[root])
    print('tested', tested, 'bad', bad)
