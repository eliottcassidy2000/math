"""AUDIT G1: structural census of the 193 order-12 escaper classes (own code; no session code).

(a) cuts (A, B), |A|, |B| >= 2, with exactly 3 arcs from B back to A, both T[A], T[B] first creators
    (FC of order k: H not in S(k - 1)); per class: sizes, H values, whether the 3 back arcs form a matching.
(b) some T - v strong and non-prime whose maximal modules (found as the classes of the relation
    "closure({u, w}) != V", a different algorithm from the session's subset enumeration) and quotient are all FCs
    (a single vertex counts as the FC of 1, as in the session).
(c) the "order-11 FC with a doubled vertex and one reversal" description, tested directly.
"""
import os
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
S = {n: set(int(x) for x in open(os.path.join(HERE, 'S_%d.txt' % n))) for n in range(1, 12)}
S[0] = set()


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
    return out


def ham(out, vs):
    """Hamiltonian paths of the subtournament on vertex list vs (own DP, extension at the front)."""
    k = len(vs)
    if k <= 1:
        return 1
    loc = [0] * k
    for a in range(k):
        for b in range(k):
            if a != b and out[vs[a]] >> vs[b] & 1:
                loc[a] |= 1 << b
    # f[mask][a] = number of paths covering mask that START at a (prepend vertices)
    f = [None] * (1 << k)
    for a in range(k):
        f[1 << a] = {a: 1}
    for m in range(1, 1 << k):
        cur = f[m]
        if not cur:
            continue
        for a, c in cur.items():
            # prepend b with b -> a, b not in m
            for b in range(k):
                if not m >> b & 1 and loc[b] >> a & 1:
                    d = f[m | 1 << b]
                    if d is None:
                        d = f[m | 1 << b] = {}
                    d[b] = d.get(b, 0) + c
        if m != (1 << k) - 1:
            f[m] = None
    return sum(f[(1 << k) - 1].values())


def is_fc(out, vs):
    k = len(vs)
    if k == 1:
        return True
    return ham(out, vs) not in S[k - 1]


def bits(m):
    return [v for v in range(12) if m >> v & 1]


def cut_census(out, n):
    full = (1 << n) - 1
    res = []
    for B in range(1, full):
        k = bin(B).count('1')
        if k < 2 or k > n - 2:
            continue
        back = [(b, a) for b in bits(B) for a in bits(full & ~B) if out[b] >> a & 1]   # arcs B -> A
        if len(back) != 3:
            continue
        A = full & ~B
        av, bv = bits(A), bits(B)
        fa, fb = is_fc(out, av), is_fc(out, bv)
        tails = set(x for x, y in back)
        heads = set(y for x, y in back)
        matching = len(tails) == 3 and len(heads) == 3
        res.append((len(av), len(bv), fa, fb, ham(out, av), ham(out, bv), matching))
    return res


def closure(out, n, verts, u, w):
    M = {u, w}
    ch = True
    while ch:
        ch = False
        for x in verts:
            if x in M:
                continue
            o = sum(1 for m in M if out[x] >> m & 1)
            if 0 < o < len(M):
                M.add(x)
                ch = True
    return M


def strong(out, verts):
    vs = set(verts)
    for start in (verts[0],):
        for rev in (False, True):
            seen = {start}
            st = [start]
            while st:
                v = st.pop()
                for w in vs:
                    if w not in seen and ((out[w] >> v & 1) if rev else (out[v] >> w & 1)):
                        seen.add(w)
                        st.append(w)
            if seen != vs:
                return False
    return True


def subst_fc(out, verts):
    """verts: vertex list of an induced subtournament U. True iff U is strong, non-prime, and its maximal modules
    and quotient are first creators."""
    if not strong(out, verts):
        return False
    V = set(verts)
    # relation u ~ w iff closure({u,w}) != V ; for strong U its classes are the maximal modules
    parent = {v: v for v in verts}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    nonprime = False
    for i, u in enumerate(verts):
        for w in verts[i + 1:]:
            if closure(out, len(out), verts, u, w) != V:
                nonprime = True
                parent[find(u)] = find(w)
    if not nonprime:
        return False
    classes = {}
    for v in verts:
        classes.setdefault(find(v), []).append(v)
    parts = list(classes.values())
    # sanity: each class is a module
    for p in parts:
        ps = set(p)
        for x in verts:
            if x in ps:
                continue
            o = sum(1 for m in p if out[x] >> m & 1)
            assert o in (0, len(p)), 'class is not a module'
    reps = [p[0] for p in parts]
    if not is_fc(out, reps):            # the quotient is the subtournament on one representative per class
        return False
    return all(is_fc(out, p) for p in parts)


def doubled(out, n):
    """pairs with exactly one distinguisher y; reverse y-z (z in {u,w}); quotient = result minus w; FC test."""
    hits = 0
    for u in range(n):
        for w in range(u + 1, n):
            d = [y for y in range(n) if y not in (u, w) and ((out[y] >> u) ^ (out[y] >> w)) & 1]
            if len(d) != 1:
                continue
            y = d[0]
            for z in (u, w):
                o = out[:]
                o[y] ^= 1 << z
                o[z] ^= 1 << y
                q = [v for v in range(n) if v != w]
                if ham(o, q) not in S[n - 2]:
                    hits += 1
    return hits


def main():
    rows = [l.split() for l in open(os.path.join(HERE, 'classes193_props.txt'))]
    tab = Counter()
    pairs = Counter()
    sizes = Counter()
    fcfc_per_class = Counter()
    special = {}
    dbl = 0
    for r in rows:
        code, H = r[0], int(r[1][2:])
        out = decode(code, 12)
        cuts = cut_census(out, 12)
        good = [c for c in cuts if c[2] and c[3]]
        a = bool(good)
        fcfc_per_class[len(set((c[0], c[1], c[4], c[5]) for c in good))] += 1
        for c in set((c[0], c[1], c[4], c[5], c[6]) for c in good):
            pairs[(c[2], c[3])] += 1
            sizes[(c[0], c[1])] += 1
        b = any(subst_fc(out, [x for x in range(12) if x != v]) for v in range(12))
        tab[(a, b)] += 1
        if H in (48537, 48849, 33021, 57643):
            special[H] = good
        dbl += doubled(out, 12) > 0
    print('(re-connected sum of two FCs via a 3-back-arc cut, FC-substitution plus a vertex) -> classes:')
    for k in sorted(tab):
        print('  ', k, tab[k])
    print('classes with an FC-FC 3-cut:', sum(v for k, v in tab.items() if k[0]))
    print('classes with an FC-substitution deletion:', sum(v for k, v in tab.items() if k[1]))
    print('side sizes of FC-FC 3-cuts:', dict(sizes))
    print('FC-FC cut H pairs (A earlier, B later) -> classes:', pairs.most_common())
    print('distinct FC-FC cuts per class:', dict(fcfc_per_class))
    for H, g in sorted(special.items()):
        print('  H=%d FC-FC cuts (|A|,|B|,fcA,fcB,H(A),H(B),back arcs a matching): %s' % (H, sorted(set(g))))
    print('classes where the doubled-vertex-plus-one-reversal quotient is an order-11 FC:', dbl)


if __name__ == '__main__':
    main()
