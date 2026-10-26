# Independent canonical labelling of vertex-coloured tournaments by individualisation-refinement (no pruning), written
# for the audit of HYP-9416.  A tournament is a tuple of out-neighbour bitmasks (bit w of out[v] set iff v -> w).
# canon(out, col) returns (code, aut, order): order lists the vertices in canonical position order, code encodes the
# relabelled upper triangle (bit index over pairs i<j, bit 1 iff order[i] -> order[j]), and aut counts the leaves that
# realise the minimum code, which equals the order of the colour-preserving automorphism group (the action of Aut on
# the leaves of the unpruned search tree is free and transitive on leaves with equal code).
# Initial colours must be integers 0..k-1; the canonical order lists colour-0 vertices first, then colour 1, etc.
from collections import Counter

def refine(out, n, col):
    k = max(col) + 1
    while True:
        masks = [0] * k
        for v in range(n): masks[col[v]] |= 1 << v
        sig = [(col[v],) + tuple((out[v] & m).bit_count() for m in masks) for v in range(n)]
        u = sorted(set(sig))
        if len(u) == k: return col
        d = {s: i for i, s in enumerate(u)}
        col = [d[s] for s in sig]
        k = len(u)

def leaf_code(out, n, col):
    order = [0] * n
    for v in range(n): order[col[v]] = v
    code = 0; b = 0
    for i in range(n - 1):
        row = out[order[i]]
        for j in range(i + 1, n):
            if (row >> order[j]) & 1: code |= 1 << b
            b += 1
    return code, order

def _search(out, n, col):
    col = refine(out, n, col)
    k = max(col) + 1
    if k == n:
        code, order = leaf_code(out, n, col)
        return code, 1, order
    cnt = Counter(col)
    c = min(x for x in range(k) if cnt[x] > 1)
    best = None; bcount = 0; border = None
    for v in range(n):
        if col[v] != c: continue
        raw = [2 * col[u] + (1 if (col[u] == c and u != v) else 0) for u in range(n)]
        rank = {x: i for i, x in enumerate(sorted(set(raw)))}
        code, m, order = _search(out, n, [rank[x] for x in raw])
        if best is None or code < best: best, bcount, border = code, m, order
        elif code == best: bcount += m
    return best, bcount, border

def canon(out, col=None):
    n = len(out)
    if col is None: col = [0] * n
    else:
        rank = {x: i for i, x in enumerate(sorted(set(col)))}
        col = [rank[x] for x in col]
    return _search(tuple(out), n, list(col))

def relabel(out, order):
    """new vertex i = old vertex order[i]"""
    n = len(out); pos = [0] * n
    for i, v in enumerate(order): pos[v] = i
    new = [0] * n
    for i, v in enumerate(order):
        row = out[v]; r = 0
        while row:
            w = row & -row; r |= 1 << pos[w.bit_length() - 1]; row ^= w
        new[i] = r
    return tuple(new)

def free_key(out):
    return (len(out), canon(out)[0])

def rooted_key(out, s, t):
    n = len(out); col = [2] * n; col[s] = 0; col[t] = 1
    return (n, 'st', canon(out, col)[0])

def onesided_key(out, r):
    n = len(out); col = [1] * n; col[r] = 0
    return (n, 'r', canon(out, col)[0])

def from_code(n, code):
    out = [0] * n; b = 0
    for i in range(n - 1):
        for j in range(i + 1, n):
            if (code >> b) & 1: out[i] |= 1 << j
            else: out[j] |= 1 << i
            b += 1
    return tuple(out)

def all_tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for m in range(1 << len(pairs)):
        out = [0] * n
        for b, (i, j) in enumerate(pairs):
            if (m >> b) & 1: out[i] |= 1 << j
            else: out[j] |= 1 << i
        yield tuple(out)

def is_strong(out):
    n = len(out)
    sc = sorted(o.bit_count() for o in out)
    s = 0
    for k in range(1, n):
        s += sc[k - 1]
        if s == k * (k - 1) // 2: return False     # Landau: a proper initial segment is a dominated set
    return True

def ham_paths(out):
    n = len(out); full = (1 << n) - 1; res = []
    def rec(path, used, last):
        if used == full:
            res.append(tuple(path)); return
        m = out[last] & ~used
        while m:
            w = m & -m; m ^= w; i = w.bit_length() - 1
            path.append(i); rec(path, used | w, i); path.pop()
    for v in range(n): rec([v], 1 << v, v)
    return res

def partner(out, P):
    """C^op xor E(P): path arcs kept, every other arc reversed"""
    n = len(out); full = (1 << n) - 1
    op = [full & ~(1 << u) & ~out[u] for u in range(n)]
    for a, b in zip(P, P[1:]):
        op[a] |= 1 << b; op[b] &= ~(1 << a)
    return tuple(op)

if __name__ == '__main__':
    # Self-test 1: exhaustive class counts from ALL labelled tournaments, and the orbit-stabiliser identity
    #   sum over classes n!/|Aut| = 2^(n choose 2).
    import itertools, random, time
    from math import factorial
    for n in range(1, 7):
        t0 = time.time()
        cls = {}
        for out in all_tournaments(n):
            code, aut, order = canon(out)
            cls.setdefault(code, aut)
        tot = sum(factorial(n) // a for a in cls.values())
        print(f"n = {n}: {len(cls)} classes from all {2 ** (n * (n - 1) // 2)} labelled tournaments; "
              f"sum n!/|Aut| = {tot}  [{time.time() - t0:.1f}s]")
        assert tot == 2 ** (n * (n - 1) // 2)
    # Self-test 2: brute-force comparison of rooted canonical forms at n = 5 (all labelled tournaments, all (s, t)):
    # the partition into rooted classes must equal the partition by brute-force minimum over all relabellings fixing
    # (s, t) in position (0, 1).
    n = 5
    def brute_rooted(out, s, t):
        rest = [v for v in range(n) if v not in (s, t)]
        best = None
        for p in itertools.permutations(rest):
            order = [s, t] + list(p)
            c = leaf_code(out, n, [order.index(v) for v in range(n)])[0]
            if best is None or c < best: best = c
        return best
    rnd = random.Random(1)
    allt = list(all_tournaments(n)); sample = rnd.sample(allt, 300)
    pairs_mine = {}; pairs_brute = {}
    for out in sample:
        for s in range(n):
            for t in range(n):
                if s == t: continue
                pairs_mine[(out, s, t)] = rooted_key(out, s, t)
                pairs_brute[(out, s, t)] = brute_rooted(out, s, t)
    m2b = {}; b2m = {}
    for key in pairs_mine:
        m2b.setdefault(pairs_mine[key], set()).add(pairs_brute[key])
        b2m.setdefault(pairs_brute[key], set()).add(pairs_mine[key])
    ok = all(len(x) == 1 for x in m2b.values()) and all(len(x) == 1 for x in b2m.values())
    print(f"rooted canonical partition equals the brute-force partition on {len(pairs_mine)} rooted labelled "
          f"tournaments (n = 5; {len(m2b)} classes): {ok}")
    assert ok
    # Self-test 2b: n = 7 representatives from nauty's gentourng are pairwise non-isomorphic for this canonical form,
    # and sum 7!/|Aut| = 2^21 (completeness and correctness of the automorphism counts).
    import subprocess
    reps = []
    for s in subprocess.run(['gentourng', '-q', '7'], capture_output=True, text=True).stdout.split():
        out = [0] * 7; k = 0
        for i in range(7):
            for j in range(i + 1, 7):
                if s[k] == '1': out[i] |= 1 << j
                else: out[j] |= 1 << i
                k += 1
        reps.append(tuple(out))
    cs = [canon(o) for o in reps]
    assert len({c[0] for c in cs}) == len(reps) == 456
    assert sum(factorial(7) // c[1] for c in cs) == 2 ** 21
    print("n = 7: 456 gentourng classes pairwise distinct; sum 7!/|Aut| = 2^21: ok")
    # Self-test 3: invariance under random relabelling at n = 9, 10, 12 with colours
    for n in (9, 10, 12):
        for _ in range(300):
            out = from_code(n, rnd.getrandbits(n * (n - 1) // 2))
            col = [rnd.randrange(3) for _ in range(n)]
            perm = list(range(n)); rnd.shuffle(perm)       # new vertex i = old perm[i]
            out2 = relabel(out, perm); col2 = [col[perm[i]] for i in range(n)]
            assert canon(out, col)[0] == canon(out2, col2)[0]
    print("relabelling invariance at n = 9, 10, 12 (random coloured tournaments): ok")
