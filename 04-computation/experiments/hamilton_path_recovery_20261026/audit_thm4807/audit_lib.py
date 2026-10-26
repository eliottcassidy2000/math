"""Independent audit library for THM-4807 (auditor's own code).

Imports nothing from the lane's scripts.  Tournaments are lists of out-neighbour
bitmasks: out[i] has bit j iff i -> j.  Vertex words are tuples.

Two independent ways to get Hamilton-path data are provided:
  * explicit DFS enumeration of all Hamilton words (all_hps), and
  * subset dynamic programming (hp_tables): F[S] = set of possible last vertices,
    B[S] = set of possible first vertices of a Hamilton path of T[S].
The census scripts cross-check the two.
"""

import random


# ---------------------------------------------------------------- basics

def bits_iter(x):
    while x:
        low = x & -x
        yield low.bit_length() - 1
        x ^= low


def popcount(x):
    return bin(x).count("1")


def pair_list(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def from_upper_string(n, s):
    """gentourng format: upper triangle row by row, '1' at (i,j), i<j, means i -> j."""
    assert len(s) == n * (n - 1) // 2, (n, s)
    out = [0] * n
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if s[k] == "1":
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
            k += 1
    return out


def from_mask(n, mask):
    """bit k (LSB first) of mask <-> k-th lexicographic pair (i<j); 1 means i -> j."""
    out = [0] * n
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (mask >> k) & 1:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
            k += 1
    return out


def to_mask(n, out):
    m = 0
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (out[i] >> j) & 1:
                m |= 1 << k
            k += 1
    return m


def check_tournament(n, out):
    for i in range(n):
        if (out[i] >> i) & 1:
            return False
        for j in range(i + 1, n):
            a = (out[i] >> j) & 1
            b = (out[j] >> i) & 1
            if a + b != 1:
                return False
    return True


def in_masks(n, out):
    inn = [0] * n
    for i in range(n):
        for j in bits_iter(out[i]):
            inn[j] |= 1 << i
    return inn


def reverse_arcs(n, out, arcs):
    """Return a copy of T with each arc (a,b) (meaning a->b in T) reversed."""
    o = out[:]
    for a, b in arcs:
        assert (o[a] >> b) & 1, "arc not present"
        o[a] &= ~(1 << b)
        o[b] |= 1 << a
    return o


def relabel(n, out, perm):
    """perm[i] = new name of old vertex i."""
    o = [0] * n
    for i in range(n):
        for j in bits_iter(out[i]):
            o[perm[i]] |= 1 << perm[j]
    return o


def is_legal(out, word):
    for a, b in zip(word, word[1:]):
        if not (out[a] >> b) & 1:
            return False
    return True


def is_ham_word(n, out, word):
    return len(word) == n and len(set(word)) == n and is_legal(out, word)


# ---------------------------------------------------------------- explicit Hamilton words

def all_hps(n, out, allowed=None):
    """All Hamilton words of T[allowed] by depth-first search (explicit)."""
    if allowed is None:
        allowed = (1 << n) - 1
    res = []
    path = []

    def rec(last, used):
        if used == allowed:
            res.append(tuple(path))
            return
        x = out[last] & allowed & ~used
        while x:
            low = x & -x
            w = low.bit_length() - 1
            x ^= low
            path.append(w)
            rec(w, used | low)
            path.pop()

    for s in bits_iter(allowed):
        path.append(s)
        rec(s, 1 << s)
        path.pop()
    return res


# ---------------------------------------------------------------- subset DP

def hp_tables(n, out, inn=None):
    """F[S]: vertices v in S that END some Hamilton path of T[S];
       B[S]: vertices v in S that START some Hamilton path of T[S]."""
    if inn is None:
        inn = in_masks(n, out)
    N = 1 << n
    F = [0] * N
    B = [0] * N
    for v in range(n):
        F[1 << v] = 1 << v
        B[1 << v] = 1 << v
    for S in range(1, N):
        if S & (S - 1) == 0:
            continue
        f = 0
        b = 0
        x = S
        while x:
            low = x & -x
            v = low.bit_length() - 1
            x ^= low
            S2 = S ^ low
            if F[S2] & inn[v]:
                f |= low
            if B[S2] & out[v]:
                b |= low
        F[S] = f
        B[S] = b
    return F, B


def traced_from_tables(n, out, F, B):
    """tr[u] = mask of v such that u -> v is consecutive in some Hamilton word."""
    full = (1 << n) - 1
    tr = [0] * n
    for S in range(1, full):
        f = F[S]
        if not f:
            continue
        b = B[full ^ S]
        if not b:
            continue
        for u in bits_iter(f):
            tr[u] |= out[u] & b
    return tr


def traced_from_words(n, words):
    tr = [0] * n
    for w in words:
        for a, b in zip(w, w[1:]):
            tr[a] |= 1 << b
    return tr


def path_from_tables(out, inn, F, B, S, start=None, end=None):
    """Some Hamilton path of T[S] with the given start OR the given end (not both)."""
    assert start is None or end is None
    if start is None and end is None:
        start = next(bits_iter(B[S]))
    if start is not None:
        assert (B[S] >> start) & 1
        path = [start]
        rem = S & ~(1 << start)
        cur = start
        while rem:
            cand = out[cur] & B[rem]
            assert cand
            w = (cand & -cand).bit_length() - 1
            path.append(w)
            rem &= ~(1 << w)
            cur = w
        return path
    assert (F[S] >> end) & 1
    path = [end]
    rem = S & ~(1 << end)
    cur = end
    while rem:
        cand = inn[cur] & F[rem]
        assert cand
        w = (cand & -cand).bit_length() - 1
        path.append(w)
        rem &= ~(1 << w)
        cur = w
    path.reverse()
    return path


def hp_count_matrix_row(n, out, s, allowed=None):
    """cnt[t] = number of Hamilton paths of T[allowed] from s to t (DP with fixed start)."""
    if allowed is None:
        allowed = (1 << n) - 1
    N = 1 << n
    dp = {1 << s: {s: 1}}
    order = sorted([S for S in range(N) if (S & allowed) == S and (S >> s) & 1], key=popcount)
    table = {}
    table[1 << s] = {s: 1}
    for S in order:
        if S not in table:
            continue
        row = table[S]
        for last, c in row.items():
            x = out[last] & allowed & ~S
            for w in bits_iter(x):
                S2 = S | (1 << w)
                r2 = table.setdefault(S2, {})
                r2[w] = r2.get(w, 0) + c
    return table.get(allowed, {})


# ---------------------------------------------------------------- strong components

def reach_set(out, v, allowed):
    seen = 1 << v
    frontier = seen
    while frontier:
        nxt = 0
        x = frontier
        while x:
            low = x & -x
            w = low.bit_length() - 1
            x ^= low
            nxt |= out[w]
        nxt &= allowed & ~seen
        seen |= nxt
        frontier = nxt
    return seen


def scc_ordered(n, out, allowed=None):
    """Strong components of the digraph out restricted to `allowed`, ordered so that
    the reach set strictly decreases (valid when the condensation is a chain, e.g. any
    tournament).  Returns (list of component masks, flag 'condensation is a chain')."""
    if allowed is None:
        allowed = (1 << n) - 1
    verts = list(bits_iter(allowed))
    R = {v: reach_set(out, v, allowed) for v in verts}
    comps = []
    assigned = 0
    for v in verts:
        if (assigned >> v) & 1:
            continue
        c = 0
        for w in verts:
            if (R[v] >> w) & 1 and (R[w] >> v) & 1:
                c |= 1 << w
        comps.append(c)
        assigned |= c
    rep = {c: (c & -c).bit_length() - 1 for c in comps}
    comps.sort(key=lambda c: -popcount(R[rep[c]]))
    chain = True
    for t in range(len(comps) - 1):
        # C_t must reach C_{t+1} and the reach sets must be nested and strictly decreasing
        if not (R[rep[comps[t]]] >> rep[comps[t + 1]]) & 1:
            chain = False
        if R[rep[comps[t]]] != comps[t] | R[rep[comps[t + 1]]]:
            chain = False
    return comps, chain


def is_strong(n, out, allowed=None):
    if allowed is None:
        allowed = (1 << n) - 1
    if allowed == 0:
        return True
    v = (allowed & -allowed).bit_length() - 1
    if reach_set(out, v, allowed) != allowed:
        return False
    inn = in_masks(n, out)
    return reach_set(inn, v, allowed) == allowed


def comp_index(n, comps):
    idx = [-1] * n
    for t, c in enumerate(comps):
        for v in bits_iter(c):
            idx[v] = t
    return idx


# ---------------------------------------------------------------- decoders

def decode_endpoint_from_words(n, words):
    """THM-4807 (e), computed only from the list of Hamilton words.
    Returns (out, info) where info records failures of the 'exactly one direction'."""
    info = {"bad_pairs": 0}
    traced = traced_from_words(n, words)
    cuts = []
    for m in range(n + 1):
        sets = {frozenset(w[:m]) for w in words}
        if len(sets) == 1:
            cuts.append(m)
    w0 = words[0]
    blocks = [frozenset(w0[cuts[t]:cuts[t + 1]]) for t in range(len(cuts) - 1)]
    bidx = {}
    for t, b in enumerate(blocks):
        for x in b:
            bidx[x] = t
    out = traced[:]
    for x in range(n):
        for y in range(n):
            if x != y and bidx[x] < bidx[y]:
                out[x] |= 1 << y
    for B in blocks:
        if len(B) < 3:
            continue
        ends = set()
        for w in words:
            r = [x for x in w if x in B]
            ends.add((r[0], r[-1]))
        bl = sorted(B)
        for a in range(len(bl)):
            for b in range(a + 1, len(bl)):
                x, y = bl[a], bl[b]
                if (traced[x] >> y) & 1 or (traced[y] >> x) & 1:
                    continue
                fx = (x, y) in ends
                fy = (y, x) in ends
                if fx == fy:
                    info["bad_pairs"] += 1
                    continue
                if fx:
                    out[x] |= 1 << y
                else:
                    out[y] |= 1 << x
    info["cuts"] = cuts
    info["blocks"] = blocks
    return out, info


def decode_G1(n, D):
    """guard_reconstruction G1 two-step-count decoder, from the trace support D only."""
    info = {"bad_pairs": 0, "stats": []}
    Din = in_masks(n, D)
    comps, chain = scc_ordered(n, D)
    info["chain"] = chain
    idx = comp_index(n, comps)
    out = D[:]
    for x in range(n):
        for y in range(n):
            if x != y and idx[x] < idx[y]:
                out[x] |= 1 << y
    for C in comps:
        cl = list(bits_iter(C))
        for a in range(len(cl)):
            for b in range(a + 1, len(cl)):
                x, y = cl[a], cl[b]
                if (D[x] >> y) & 1 or (D[y] >> x) & 1:
                    continue
                cxy = popcount(D[x] & Din[y] & C)
                cyx = popcount(D[y] & Din[x] & C)
                info["stats"].append((x, y, cxy, cyx))
                if (cxy >= 2) == (cyx >= 2):
                    info["bad_pairs"] += 1
                    continue
                if cxy >= 2:
                    out[x] |= 1 << y
                else:
                    out[y] |= 1 << x
    return out, info


def decode_G4(n, D):
    """Session70 G4/G5 deleted-vertex one-way reachability decoder, from D only."""
    info = {"bad_pairs": 0, "eligible": []}
    comps, chain = scc_ordered(n, D)
    idx = comp_index(n, comps)
    out = D[:]
    for x in range(n):
        for y in range(n):
            if x != y and idx[x] < idx[y]:
                out[x] |= 1 << y
    for C in comps:
        cl = list(bits_iter(C))
        for a in range(len(cl)):
            for b in range(a + 1, len(cl)):
                x, y = cl[a], cl[b]
                if (D[x] >> y) & 1 or (D[y] >> x) & 1:
                    continue
                dirs = set()
                elig = 0
                for v in cl:
                    if v == x or v == y:
                        continue
                    A = C & ~(1 << v)
                    xy = (reach_set(D, x, A) >> y) & 1
                    yx = (reach_set(D, y, A) >> x) & 1
                    if xy != yx:
                        elig += 1
                        dirs.add((x, y) if xy else (y, x))
                info["eligible"].append(elig)
                if len(dirs) != 1:
                    info["bad_pairs"] += 1
                    continue
                (p, q), = dirs
                out[p] |= 1 << q
    return out, info


# ---------------------------------------------------------------- random generators

def random_tournament(n, rng):
    out = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
    return out


def random_strong_on(vertices, out, rng, tries=1000):
    """Place a random strong subtournament (or singleton) on the given vertex list,
    writing into out (all pairs among `vertices` overwritten)."""
    vs = list(vertices)
    m = len(vs)
    assert m == 1 or m >= 3
    if m == 1:
        return
    for _ in range(tries):
        # random Hamiltonian cycle plus random chords -> always strong
        rng.shuffle(vs)
        for a in range(m):
            for b in range(a + 1, m):
                x, y = vs[a], vs[b]
                out[x] &= ~(1 << y)
                out[y] &= ~(1 << x)
                if b == a + 1:
                    out[x] |= 1 << y
                elif a == 0 and b == m - 1:
                    out[y] |= 1 << x
                else:
                    if rng.random() < 0.5:
                        out[x] |= 1 << y
                    else:
                        out[y] |= 1 << x
        return


def ordered_join(sizes, rng):
    """Random tournament C_1 => C_2 => ... with random strong components of the given sizes
    (each size 1 or >= 3).  Returns (n, out, comps as lists)."""
    n = sum(sizes)
    out = [0] * n
    comps = []
    start = 0
    for s in sizes:
        comps.append(list(range(start, start + s)))
        start += s
    for a in range(len(comps)):
        for b in range(a + 1, len(comps)):
            for x in comps[a]:
                for y in comps[b]:
                    out[x] |= 1 << y
    for c in comps:
        random_strong_on(c, out, rng)
    return n, out, comps


def random_size(rng, lo=1, hi=5):
    while True:
        s = rng.randint(lo, hi)
        if s != 2:
            return s
