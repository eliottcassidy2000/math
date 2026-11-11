"""cd_extremal.py - extremal questions for sign tournaments on PG(2,2) and PG(3,2) (S23).

1. Fano plane: distribution of the number of even zero-sum 4-sets (line complements) over all 2^21 tournaments.
2. PG(3,2): the 135 maximal sets of 8 planes that can be odd at once: types, solution dimensions, compatibility with
   'all lines cyclic'; a representative of each type: even 4-sets, cyclic lines, scores, H, double regularity.
3. PG(3,2): least number of even 4-sets over all tournaments: lower bound from planes, upper bound by annealing.
4. Split forms of dimension 16: even 4-sets for every choice of the parameters."""
import itertools, random, sys
from cd_more import F2System, planes, quads, c3_of, table
from cd_signs import ham_paths


def popc(x):
    return bin(x).count("1")


def fano_distribution():
    S = F2System(3)
    Q = quads(3)
    eqs = []
    for q in Q:
        eqs.append(S.quad_rows(q))
    dist = {}
    for x in range(1 << 21):
        bad = 0
        for (r1, r2) in eqs:
            if (popc(x & r1[0]) & 1) != r1[1] or (popc(x & r2[0]) & 1) != r2[1]:
                bad += 1
        dist[bad] = dist.get(bad, 0) + 1
    return dist


def solve_full(rows, nv):
    """Gaussian elimination; returns (consistent, rank, one solution as int or None)."""
    piv = {}
    ok = True
    for mask, const in rows:
        while mask:
            hb = mask.bit_length() - 1
            if hb in piv:
                pm, pc = piv[hb]
                mask ^= pm
                const ^= pc
            else:
                piv[hb] = (mask, const)
                break
        else:
            if const:
                ok = False
    if not ok:
        return False, len(piv), None
    x = 0
    for hb in sorted(piv):                 # back substitution, free variables = 0
        mask, const = piv[hb]
        val = const ^ (popc(x & (mask ^ (1 << hb))) & 1)
        if val:
            x |= 1 << hb
    return True, len(piv), x


def describe(S, x, P, Q4, label):
    N = 16
    def arc(u, v):
        if u < v:
            return not (x >> S.var[(u, v)] & 1)
        return bool(x >> S.var[(v, u)] & 1)
    V = range(1, N)
    lines = sorted({tuple(sorted((a, b, a ^ b))) for a in V for b in V if a != b})
    cyc = sum(1 for L in lines if c3_of(arc, L) == 1)
    ev = sum(1 for q in Q4 if c3_of(arc, q) % 2 == 0)
    scores = sorted(sum(1 for j in V if j != i and arc(i, j)) for i in V)
    common = set()
    for i in V:
        for j in V:
            if i < j:
                common.add(sum(1 for k in V if k != i and k != j and arc(i, k) and arc(j, k)))
    per_plane = sorted(sum(1 for q in Q4 if set(q) <= set(pl) and c3_of(arc, q) % 2 == 0) for pl in P)
    adj = [[1 if (i != j and arc(i + 1, j + 1)) else 0 for j in range(N - 1)] for i in range(N - 1)]
    print("   %s: cyclic lines %d/35, even 4-sets %d, per plane %s, scores %s, common-out counts %s, H = %d"
          % (label, cyc, ev, per_plane, sorted(set(scores)), sorted(common), ham_paths(adj)))


if __name__ == "__main__":
    random.seed(23)
    if "fano" in sys.argv:
        print("Fano plane, number of even 4-sets -> number of tournaments:", sorted(fano_distribution().items()))
    S = F2System(4)
    nv = len(S.var)
    P = planes(4)
    Q4 = quads(4)
    rows_by_plane = []
    for pl in P:
        s = set(pl)
        rr = []
        for q in Q4:
            if s >= set(q):
                rr += S.quad_rows(q)
        rows_by_plane.append(rr)
    line_rows = []
    V = range(1, 16)
    lines = sorted({tuple(sorted((a, b, a ^ b))) for a in V for b in V if a != b})
    for (a, b, c) in lines:
        for (p, q) in (((a, b), (b, c)), ((b, c), (c, a))):
            m1, c1 = S.term(*p)
            m2, c2 = S.term(*q)
            line_rows.append((m1 ^ m2, c1 ^ c2))
    types = {}
    for mask in range(1, 1 << 15):
        if popc(mask) != 8:
            continue
        rr = []
        for i in range(15):
            if mask >> i & 1:
                rr += rows_by_plane[i]
        ok, rk, x = solve_full(rr, nv)
        if not ok:
            continue
        sel = [set(P[i]) for i in range(15) if mask >> i & 1]
        rest = [set(P[i]) for i in range(15) if not (mask >> i & 1)]
        # points on no selected plane / on all unselected planes
        missing = [p for p in V if all(p not in s for s in sel)]
        common_rest = [p for p in V if all(p in s for s in rest)]
        t = "A: the 8 planes off a point" if missing else ("B: a plane and the 7 planes through a point off it" if True else "?")
        if not missing:
            # verify type B: one selected plane P0 and a point p not on it with the other seven through p
            okB = False
            for s0 in sel:
                others = [s for s in sel if s is not s0]
                pts = [p for p in V if all(p in s for s in others)]
                if len(pts) == 1 and pts[0] not in s0:
                    okB = True
            if not okB:
                t = "other"
        ok2, rk2, x2 = solve_full(rr + line_rows, nv)
        key = (t, nv - rk, ok2, (nv - rk2) if ok2 else None)
        types.setdefault(key, []).append((mask, x, x2))
    for key, lst in sorted(types.items()):
        print("maximal odd plane sets of type [%s]: %d sets; solution dimension %d; with all lines cyclic: consistent %s, dimension %s"
              % (key[0], len(lst), key[1], key[2], key[3]))
        mask, x, x2 = lst[0]
        describe(S, x, P, Q4, "a solution (free variables 0)")
        if x2 is not None:
            describe(S, x2, P, Q4, "a solution with all lines cyclic")
    # 3. least number of even 4-sets: annealing
    quad_eq = [S.quad_rows(q) for q in Q4]
    by_var = [[] for _ in range(nv)]
    for qi, (r1, r2) in enumerate(quad_eq):
        for v in range(nv):
            if (r1[0] | r2[0]) >> v & 1:
                by_var[v].append(qi)

    def badq(x, qi):
        r1, r2 = quad_eq[qi]
        return (popc(x & r1[0]) & 1) != r1[1] or (popc(x & r2[0]) & 1) != r2[1]
    import math
    best_all = 10 ** 9
    hist = {}
    for rep in range(int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 40):
        x = random.getrandbits(nv)
        cur = sum(badq(x, qi) for qi in range(len(Q4)))
        best = cur
        T = 1.5
        for step in range(60000):
            v = random.randrange(nv)
            before = sum(badq(x, qi) for qi in by_var[v])
            y = x ^ (1 << v)
            after = sum(badq(y, qi) for qi in by_var[v])
            d = after - before
            if d <= 0 or random.random() < math.exp(-d / T):
                x = y
                cur += d
                if cur < best:
                    best = cur
            T = max(0.05, T * 0.99992)
        hist[best] = hist.get(best, 0) + 1
        best_all = min(best_all, best)
    print("annealing, least number of even 4-sets found per run:", sorted(hist.items()), " overall", best_all)
    # 4. split forms of dimension 16
    for gs in itertools.product((-1, 1), repeat=4):
        g = table(list(gs))
        arc = lambda i, j: g[i][j] == 1
        ev = sum(1 for q in Q4 if c3_of(arc, q) % 2 == 0)
        neg = sum(1 for i in range(1, 16) if g[i][i] == -1)
        print("   parameters %s: units of square -1: %d, even 4-sets %d (pairings %d)" % (gs, neg, ev, 2 * ev))
