"""aud_common.py (AUDIT, independent code) - small tournament toolkit written for the audit of the S23 packet.
Nothing here imports the author's scripts."""
import itertools, subprocess, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))


def decode(code, m):
    """adjacency (list of sets of out-neighbours); bit b (LSB first) for the pair (i<j) in lexicographic order, 1 iff i -> j."""
    out = [set() for _ in range(m)]
    b = 0
    for i in range(m):
        for j in range(i + 1, m):
            if code >> b & 1:
                out[i].add(j)
            else:
                out[j].add(i)
            b += 1
    return out


def masks(out):
    return [sum(1 << u for u in s) for s in out]


def ham_py(out):
    """number of Hamiltonian paths, pure Python subset DP (for small n)."""
    n = len(out)
    om = masks(out)
    f = [dict() for _ in range(1 << n)]
    for v in range(n):
        f[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, full + 1):
        fm = f[mask]
        for v, c in fm.items():
            nb = om[v] & ~mask
            while nb:
                low = nb & -nb
                u = low.bit_length() - 1
                nb ^= low
                d = f[mask | low]
                d[u] = d.get(u, 0) + c
    return sum(f[full].values())


def ham_brute(out):
    """number of Hamiltonian paths by enumerating permutations (independent of the DP; n <= 8)."""
    n = len(out)
    c = 0
    for p in itertools.permutations(range(n)):
        ok = True
        for i in range(n - 1):
            if p[i + 1] not in out[p[i]]:
                ok = False
                break
        if ok:
            c += 1
    return c


def P_tuple(out, kmax):
    """P_k(T) for k = 1..kmax: number of ordered k-tuples of disjoint NONEMPTY paths covering T.
    DP over (mask, last vertex, number of paths so far)."""
    n = len(out)
    om = masks(out)
    full = (1 << n) - 1
    f = {}
    for v in range(n):
        f[(1 << v, v, 1)] = 1
    res = [0] * (kmax + 1)
    # process masks in increasing order of popcount
    layers = [dict() for _ in range(n + 1)]
    for key, c in f.items():
        layers[1][key] = c
    for size in range(1, n + 1):
        for (mask, v, k), c in layers[size].items():
            if mask == full:
                if k <= kmax:
                    res[k] += c
                continue
            rest = full & ~mask
            r = rest
            while r:
                low = r & -r
                u = low.bit_length() - 1
                r ^= low
                # continue the current path along an arc v -> u
                if om[v] >> u & 1:
                    key = (mask | low, u, k)
                    layers[size + 1][key] = layers[size + 1].get(key, 0) + c
                # start a new path at u
                if k + 1 <= kmax:
                    key = (mask | low, u, k + 1)
                    layers[size + 1][key] = layers[size + 1].get(key, 0) + c
    return res  # res[k], k >= 1


def P_tuple_brute(out, k):
    """brute force: permutations of the vertices with k-1 of the n-1 gaps chosen as breaks."""
    n = len(out)
    c = 0
    for p in itertools.permutations(range(n)):
        arcs = [p[i + 1] in out[p[i]] for i in range(n - 1)]
        # gaps that are not arcs must be breaks; gaps that are arcs may be breaks or not
        forced = sum(1 for a in arcs if not a)
        free = sum(1 for a in arcs if a)
        need = k - 1 - forced
        if 0 <= need <= free:
            # choose 'need' of the free gaps as breaks
            from math import comb
            c += comb(free, need)
    return c


def substitute(K, ports, modules):
    """K: list of out-sets; ports: list of vertices of K; modules: list of tournaments (out-sets) put at the ports.
    Returns the out-sets of K[...]."""
    m = len(K)
    blocks = []
    start = 0
    where = {}
    for v in range(m):
        if v in ports:
            mod = modules[ports.index(v)]
            sz = len(mod)
        else:
            mod = [set()]
            sz = 1
        blocks.append((start, sz, mod))
        start += sz
    N = start
    out = [set() for _ in range(N)]
    for v in range(m):
        s, sz, mod = blocks[v]
        for a in range(sz):
            for b in mod[a]:
                out[s + a].add(s + b)
        for u in K[v]:
            s2, sz2, _ = blocks[u]
            for a in range(sz):
                for b in range(sz2):
                    out[s + a].add(s2 + b)
    return out


def transitive(a):
    return [set(range(i + 1, a)) for i in range(a)]


def random_tournament(n, rng):
    out = [set() for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5:
                out[i].add(j)
            else:
                out[j].add(i)
    return out


def ham_c(list_of_out):
    """batch Hamiltonian path counts through the audit's own C program hp.exe."""
    exe = os.path.join(HERE, "hp.exe")
    lines = []
    for out in list_of_out:
        om = masks(out)
        lines.append("%d %s" % (len(out), " ".join(str(x) for x in om)))
    r = subprocess.run([exe], input="\n".join(lines) + "\n", capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stderr)
    return [int(x) for x in r.stdout.split()]


def read_classes(m):
    path = os.path.join(HERE, "..", "data", "cls_%d.txt" % m)
    res = []
    with open(path) as fh:
        for l in fh:
            t = l.split()
            if len(t) >= 2:
                res.append((int(t[0]), int(t[1])))
    return res


def visit_matrix_dfs(K, x, y):
    """N[(i, j)] = number of walks in K visiting every vertex other than x, y exactly once, x exactly i times and
    y exactly j times.  Plain depth-first enumeration of walks (no DP)."""
    m = len(K)
    others = [v for v in range(m) if v != x and v != y]
    fullmask = sum(1 << v for v in others)
    N = {}

    def rec(v, mask, i, j):
        if mask == fullmask:
            N[(i, j)] = N.get((i, j), 0) + 1
        for u in K[v]:
            if u == x:
                rec(u, mask, i + 1, j)
            elif u == y:
                rec(u, mask, i, j + 1)
            elif not (mask >> u & 1):
                rec(u, mask | (1 << u), i, j)
    import sys
    sys.setrecursionlimit(10000)
    for v in range(m):
        if v == x:
            rec(v, 0, 1, 0)
        elif v == y:
            rec(v, 0, 0, 1)
        else:
            rec(v, 1 << v, 0, 0)
    return N
