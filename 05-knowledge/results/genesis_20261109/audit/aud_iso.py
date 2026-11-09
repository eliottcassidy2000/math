"""AUDIT G1: exact tournament isomorphism by colour refinement + backtracking (own code, no canon.c).

A tournament is a list out[v] of out-neighbour bitmasks.  Refined colours are label-independent hashes
(iterated: own colour, sorted out-neighbour colours, sorted in-neighbour colours), so they can be compared
across tournaments; colours only prune the search, the final test is exact arc-by-arc consistency.
"""
import hashlib


def decode(code, n):
    out = [0] * n
    p = 0
    for i in range(n):
        for j in range(i + 1, n):
            if code[p] == '1':
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
            p += 1
    assert p == len(code)
    return out


def encode(out):
    n = len(out)
    return ''.join('1' if out[i] >> j & 1 else '0' for i in range(n) for j in range(i + 1, n))


def converse(out):
    n = len(out)
    full = (1 << n) - 1
    return [full & ~out[v] & ~(1 << v) for v in range(n)]


def _h(*parts):
    return hashlib.sha1(repr(parts).encode()).hexdigest()[:16]


def refine(out):
    n = len(out)
    # initial colour: score and number of 3-cycles through v
    col = []
    for v in range(n):
        c3 = 0
        for a in range(n):
            if out[v] >> a & 1:
                for b in range(n):
                    if out[a] >> b & 1 and out[b] >> v & 1:
                        c3 += 1
        col.append(_h(bin(out[v]).count('1'), c3))
    for _ in range(n):
        new = []
        for v in range(n):
            o = sorted(col[w] for w in range(n) if out[v] >> w & 1)
            i = sorted(col[w] for w in range(n) if w != v and not out[v] >> w & 1)
            new.append(_h(col[v], tuple(o), tuple(i)))
        if len(set(new)) == len(set(col)):
            col = new
            break
        col = new
    return col


def invariant(out):
    return tuple(sorted(refine(out)))


def _search(o1, c1, o2, c2, count_all=False):
    n = len(o1)
    from collections import Counter
    if Counter(c1) != Counter(c2):
        return 0
    cls = Counter(c1)
    order = sorted(range(n), key=lambda v: (cls[c1[v]], c1[v]))
    cand = {v: [w for w in range(n) if c2[w] == c1[v]] for v in range(n)}
    f = [-1] * n
    used = [False] * n
    found = [0]

    def rec(k):
        if k == n:
            found[0] += 1
            return not count_all
        v = order[k]
        for w in cand[v]:
            if used[w]:
                continue
            ok = True
            for t in range(k):
                u = order[t]
                if (o1[v] >> u & 1) != (o2[w] >> f[u] & 1):
                    ok = False
                    break
            if not ok:
                continue
            f[v] = w
            used[w] = True
            if rec(k + 1):
                return True
            used[w] = False
            f[v] = -1
        return False

    rec(0)
    return found[0]


def isomorphic(o1, o2, c1=None, c2=None):
    if len(o1) != len(o2):
        return False
    c1 = c1 or refine(o1)
    c2 = c2 or refine(o2)
    return _search(o1, c1, o2, c2) > 0


def aut_count(o):
    c = refine(o)
    return _search(o, c, o, c, count_all=True)


def classes(outs):
    """Partition a list of tournaments into isomorphism classes; returns list of lists of indices."""
    cols = [refine(o) for o in outs]
    inv = [tuple(sorted(c)) for c in cols]
    buckets = {}
    for i, k in enumerate(inv):
        buckets.setdefault(k, []).append(i)
    res = []
    for k, idx in buckets.items():
        reps = []
        for i in idx:
            for r in reps:
                if isomorphic(outs[r[0]], outs[i], cols[r[0]], cols[i]):
                    r.append(i)
                    break
            else:
                reps.append([i])
        res.extend(reps)
    return res


if __name__ == '__main__':
    # self-test: random relabellings are isomorphic; single arc flips of rigid ones mostly are not
    import random
    random.seed(5)
    bad = 0
    for t in range(200):
        n = random.randint(3, 12)
        out = [0] * n
        for i in range(n):
            for j in range(i + 1, n):
                if random.random() < 0.5:
                    out[i] |= 1 << j
                else:
                    out[j] |= 1 << i
        p = list(range(n))
        random.shuffle(p)
        o2 = [0] * n
        for v in range(n):
            for w in range(n):
                if out[v] >> w & 1:
                    o2[p[v]] |= 1 << p[w]
        if not isomorphic(out, o2):
            bad += 1
        # brute force aut count for n <= 7
        if n <= 7:
            import itertools
            a = sum(1 for q in itertools.permutations(range(n))
                    if all(((out[v] >> w) & 1) == ((out[q[v]] >> q[w]) & 1) for v in range(n) for w in range(n)))
            if a != aut_count(out):
                bad += 1
    print('iso self-test failures:', bad)
    # brute-force cross-check of the iso test on random pairs with equal score multisets (n <= 7)
    import itertools
    tested = 0
    for t in range(4000):
        n = random.randint(4, 7)
        def rt():
            o = [0] * n
            for i in range(n):
                for j in range(i + 1, n):
                    if random.random() < 0.5:
                        o[i] |= 1 << j
                    else:
                        o[j] |= 1 << i
            return o
        a, b = rt(), rt()
        if sorted(bin(x).count('1') for x in a) != sorted(bin(x).count('1') for x in b):
            continue
        tested += 1
        bf = any(all(((a[v] >> w) & 1) == ((b[q[v]] >> q[w]) & 1) for v in range(n) for w in range(n))
                 for q in itertools.permutations(range(n)))
        if bf != isomorphic(a, b):
            bad += 1
    print('brute-force iso pairs tested:', tested, ' total failures:', bad)
