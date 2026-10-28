"""Referee core utilities (independent audit of kind-pasteur S8, THM-4835/4836, HYP-9435).

Written from scratch by the referee. Imports nothing from the author's code.

Conventions
- A permutation of [n] is a tuple p with p[i] the image of i.  compose(p, q) = p o q.
- An oriented carrier is (n, arcs) with arcs a list of ordered pairs, no 2-cycles, len(arcs) = n-1.
- Groups are given as explicit element lists (small groups only) and wrapped in PermGroup,
  which indexes elements and keeps multiplication / conjugation tables.
"""
import itertools
import math
from collections import Counter
from fractions import Fraction


# ---------------------------------------------------------------- permutations
def compose(p, q):
    return tuple(p[i] for i in q)


def inverse(p):
    inv = [0] * len(p)
    for i, j in enumerate(p):
        inv[j] = i
    return tuple(inv)


def cycles(p):
    n = len(p)
    seen = [False] * n
    res = []
    for i in range(n):
        if not seen[i]:
            c = []
            j = i
            while not seen[j]:
                seen[j] = True
                c.append(j)
                j = p[j]
            res.append(c)
    return res


def perm_order(p):
    o = 1
    for c in cycles(p):
        o = o * len(c) // math.gcd(o, len(c))
    return o


# ---------------------------------------------------- automorphisms (backtracking)
def automorphisms(n, arcs):
    """All permutations of [n] mapping the arc set onto itself (orientation preserving)."""
    adj = [[0] * n for _ in range(n)]
    for a, b in arcs:
        adj[a][b] = 1
    outd = [sum(adj[v]) for v in range(n)]
    ind = [sum(adj[u][v] for u in range(n)) for v in range(n)]

    def inv(v):
        return (outd[v], ind[v],
                tuple(sorted((outd[w], ind[w]) for w in range(n) if adj[v][w])),
                tuple(sorted((outd[w], ind[w]) for w in range(n) if adj[w][v])))

    invs = [inv(v) for v in range(n)]
    # assignment order: BFS through the underlying graph so later vertices are constrained
    und = [[w for w in range(n) if adj[v][w] or adj[w][v]] for v in range(n)]
    order = []
    placed = [False] * n
    for s in sorted(range(n), key=lambda v: (-len(und[v]), v)):
        if placed[s]:
            continue
        queue = [s]
        placed[s] = True
        while queue:
            v = queue.pop(0)
            order.append(v)
            for w in und[v]:
                if not placed[w]:
                    placed[w] = True
                    queue.append(w)
    img = [-1] * n
    used = [False] * n
    res = []

    def bt(k):
        if k == n:
            res.append(tuple(img))
            return
        v = order[k]
        for w in range(n):
            if used[w] or invs[w] != invs[v]:
                continue
            ok = True
            for t in range(k):
                u = order[t]
                if adj[u][v] != adj[img[u]][w] or adj[v][u] != adj[w][img[u]]:
                    ok = False
                    break
            if ok:
                img[v] = w
                used[w] = True
                bt(k + 1)
                used[w] = False
                img[v] = -1

    bt(0)
    return res


# ------------------------------------------------------------------ group class
class PermGroup:
    def __init__(self, elems):
        self.elems = sorted(set(elems))
        self.deg = len(self.elems[0])
        self.N = len(self.elems)
        self.idx = {g: i for i, g in enumerate(self.elems)}
        self.e = self.idx[tuple(range(self.deg))]
        E = self.elems
        self.mul = [[self.idx[compose(a, b)] for b in E] for a in E]
        self.inv = [self.idx[inverse(a)] for a in E]
        self.conj = [[self.mul[self.mul[c][s]][self.inv[c]] for s in range(self.N)]
                     for c in range(self.N)]
        self.order_of = [perm_order(g) for g in E]
        self._pc = {}

    def gen(self, gens):
        seen = {self.e}
        frontier = [self.e]
        while frontier:
            new = []
            for x in frontier:
                for g in gens:
                    y = self.mul[g][x]
                    if y not in seen:
                        seen.add(y)
                        new.append(y)
            frontier = new
        return frozenset(seen)

    def all_subgroups(self):
        """dict: subgroup (frozenset of element indices) -> list of generator indices."""
        cyc = {}
        for g in range(self.N):
            S = self.gen([g])
            if S not in cyc:
                cyc[S] = g
        subs = {frozenset([self.e]): []}
        for S, g in cyc.items():
            if S not in subs:
                subs[S] = [g]
        frontier = list(subs)
        while frontier:
            new = []
            for H in frontier:
                for C, g in cyc.items():
                    if g in H:
                        continue
                    J = self.gen(subs[H] + [g])
                    if J not in subs:
                        subs[J] = subs[H] + [g]
                        new.append(J)
            frontier = new
        return subs

    def subgroup_classes(self, subs):
        """list of (representative, set of conjugates)."""
        seen = set()
        classes = []
        for H in sorted(subs, key=lambda s: (len(s), sorted(s))):
            if H in seen:
                continue
            orb = {frozenset(self.conj[c][h] for h in H) for c in range(self.N)}
            seen |= orb
            classes.append((H, orb))
        return classes

    def pair_class(self, S, K):
        """Canonical A-conjugacy key of a pair (S, K <= S) of subgroups."""
        key0 = (S, K)
        if key0 in self._pc:
            return self._pc[key0]
        best = None
        for c in range(self.N):
            cand = (tuple(sorted(self.conj[c][s] for s in S)),
                    tuple(sorted(self.conj[c][k] for k in K)))
            if best is None or cand < best:
                best = cand
        self._pc[key0] = best
        return best

    def coset_action(self, H):
        """Transitive A-set A/H: returns list act[g] = tuple of images of the cosets."""
        cos = []
        cidx = {}
        for a in range(self.N):
            c = frozenset(self.mul[a][h] for h in H)
            if c not in cidx:
                cidx[c] = len(cos)
                cos.append(a)
        act = []
        for g in range(self.N):
            act.append(tuple(cidx[frozenset(self.mul[self.mul[g][a]][h] for h in H)]
                             for a in cos))
        return act


def disjoint_union(actions):
    """Concatenate A-set actions (lists over group elements of image tuples)."""
    if not actions:
        return None, 0
    N = len(actions[0])
    out = []
    for g in range(N):
        img = []
        off = 0
        for act in actions:
            img.extend(off + y for y in act[g])
            off += len(act[0])
        out.append(tuple(img))
    return out, len(out[0])


# ------------------------------------------------ signed pair sets and marks
def signed_pair_types(G, act, npts, excluded=frozenset()):
    """Multiset of (A x Z2)-orbit types of ordered pairs (Z2 = swap) on the unordered pairs
    of the A-set (act, npts) minus `excluded` (set of (y,z) with y<z).  The orbit of the
    ordered pair (y,z) has stabiliser = graph of eps: S -> Z2 with S = Stab{y,z}; its type is the
    A-conjugacy class of (S, ker eps)."""
    seen = set()
    types = Counter()
    for y in range(npts):
        for z in range(y + 1, npts):
            if (y, z) in excluded or (y, z) in seen:
                continue
            S = []
            K = []
            for g in range(G.N):
                a, b = act[g][y], act[g][z]
                key = (a, b) if a < b else (b, a)
                seen.add(key)
                if key == (y, z):
                    S.append(g)
                    if a == y:
                        K.append(g)
            types[G.pair_class(frozenset(S), frozenset(K))] += 1
    return types


def ordered_pair_orbits(act_gens, npts, excluded):
    """Union-find orbits of a group (given by generator actions) on ordered non-excluded pairs.
    Returns (number of unordered-pair orbits, reversal_flag)."""
    parent = {}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    op = [(y, z) for y in range(npts) for z in range(npts)
          if y != z and (min(y, z), max(y, z)) not in excluded]
    for p in op:
        parent[p] = p
    for act in act_gens:
        for (y, z) in op:
            union((y, z), (act[y], act[z]))
    reversed_ = any(find((y, z)) == find((z, y)) for (y, z) in op)
    norb = len({find(p) for p in op})
    if reversed_:
        return None, True
    return norb // 2, False


def fix_count_exponent(act_g, npts, excluded):
    """For a single element: None if it fixes no 'tournament' on the pair set, else the number
    of <g>-orbits on unordered pairs (so it fixes 2^that many)."""
    return ordered_pair_orbits([act_g], npts, excluded)[0]


def tournament_count(G, act, npts, excluded):
    """|{orientations of the non-excluded pairs}/A| by Burnside over elements."""
    tot = 0
    cache = {}
    for g in range(G.N):
        ex = fix_count_exponent(act[g], npts, excluded)
        if ex is not None:
            tot += 2 ** ex
    q = Fraction(tot, G.N)
    assert q.denominator == 1
    return int(q)


def tournament_marks(G, act, npts, excluded, reps, subs):
    """Marks of the A-set of orientations at subgroup representatives: tuple of exponents
    (None = mark 0)."""
    out = []
    for H in reps:
        gens = subs[H] if subs[H] else []
        acts = [act[g] for g in gens]
        if not acts:
            acts = [act[G.e]]
        out.append(ordered_pair_orbits(acts, npts, excluded)[0])
    return tuple(out)


# ------------------------------------------------------------- carrier helpers
def carrier_data(n, arcs):
    E = [tuple(sorted(a)) for a in arcs]
    Eset = set(E)
    assert len(Eset) == len(arcs), "2-cycle or repeated arc"
    return E, Eset


def chi_induced(Vset, Eset):
    e = sum(1 for (a, b) in Eset if a in Vset and b in Vset)
    return len(Vset) - e, e


def fixed_set(G, K):
    return frozenset(v for v in range(G.deg) if all(G.elems[k][v] == v for k in K))
