"""aud_route_sets.py - what are the 38 shortest routes, as arc SETS?

A shortest route of length 4 from the labelled A to a labelled copy B' of B reverses exactly the 4 arcs of the
symmetric difference A xor B', in some order.  So routes = (copy B' at Hamming distance 4) x (orderings of its 4 arcs
whose 3 intermediates are strong).  This script finds the copies by relabelling B over all 720 permutations fixing the
root (pure Python, independent of the C programs), and counts strong orderings and below-start orderings per copy,
in both directions.
"""
import itertools

N = 7
A = [8, 73, 3, 4, 15, 31, 61]
B = [2, 84, 1, 7, 13, 31, 61]


def adj_of(rows):
    return [[(rows[a] >> b) & 1 for b in range(N)] for a in range(N)]


def rows_of(adj):
    return [sum(1 << b for b in range(N) if adj[a][b]) for a in range(N)]


def relabel(rows, p):  # p: old -> new
    adj = adj_of(rows)
    new = [[0] * N for _ in range(N)]
    for a in range(N):
        for b in range(N):
            new[p[a]][p[b]] = adj[a][b]
    return rows_of(new)


def strong(rows):
    full = (1 << N) - 1
    f = 1
    for _ in range(N):
        for x in range(N):
            if (f >> x) & 1:
                f |= rows[x]
    b = 1
    for _ in range(N):
        for y in range(N):
            if rows[y] & b:
                b |= 1 << y
    return f == full and b == full


def profile(rows):
    full = (1 << N) - 1
    end = [[0] * N for _ in range(1 << N)]
    for x in range(N):
        end[1 << x][x] = 1
    for S in range(1, full + 1):
        for x in range(N):
            c = end[S][x]
            if c:
                for y in range(N):
                    if not (S >> y) & 1 and (rows[x] >> y) & 1:
                        end[S | (1 << y)][y] += c
    HS = [1] + [sum(end[S]) for S in range(1, full + 1)]
    h = HS[full]
    F = sum(end[S][0] * HS[full ^ S] for S in range(1, full + 1, 2))
    # s by brute force DFS from 0
    s = 0
    stack = [(0, 1)]
    while stack:
        cur, vis = stack.pop()
        if vis == full:
            s += 1
            continue
        for y in range(N):
            if not (vis >> y) & 1 and (rows[cur] >> y) & 1:
                stack.append((y, vis | (1 << y)))
    return h, s, F


def flip(rows, a, b):
    r = rows[:]
    if (r[a] >> b) & 1:
        r[a] &= ~(1 << b)
        r[b] |= 1 << a
    else:
        r[b] &= ~(1 << a)
        r[a] |= 1 << b
    return r


def diff_arcs(r1, r2):
    return [(a, b) for a in range(N) for b in range(a + 1, N) if ((r1[a] >> b) & 1) != ((r2[a] >> b) & 1)]


for src_name, src, dst in (("A -> B", A, B), ("B -> A", B, A)):
    h0, s0, F0 = profile(src)
    copies = {}
    for perm in itertools.permutations(range(1, N)):
        p = [0] + list(perm)
        c = tuple(relabel(dst, p))
        copies.setdefault(c, 0)
    dmin = min(len(diff_arcs(src, list(c))) for c in copies)
    near = [list(c) for c in copies if len(diff_arcs(src, list(c))) == dmin]
    print("%s: %d labelled copies of the target (root fixed); minimum Hamming distance %d attained by %d copies" % (src_name, len(copies), dmin, len(near)))
    tot = tot_s1 = tot_below = 0
    for c in near:
        arcs = diff_arcs(src, c)
        nstrong = ns1 = nbelow = 0
        for order in itertools.permutations(arcs):
            cur = src
            ok = True
            s1 = True
            below = True
            for (a, b) in order[:-1]:
                cur = flip(cur, a, b)
                if not strong(cur):
                    ok = False
                    break
                h, s, F = profile(cur)
                if s != 1:
                    s1 = False
                if not (F <= F0 and h + 2 * F < h0 + 2 * F0):
                    below = False
            if ok:
                nstrong += 1
                ns1 += s1
                nbelow += below
        tot += nstrong
        tot_s1 += ns1
        tot_below += nbelow
        print("   copy with arc set %s: %d of 24 orderings stay strong; keep s = 1: %d; keep wrappers below the start: %d" % (arcs, nstrong, ns1, nbelow))
    print("   total shortest routes %d; keep s = 1: %d; below the start: %d" % (tot, tot_s1, tot_below))
