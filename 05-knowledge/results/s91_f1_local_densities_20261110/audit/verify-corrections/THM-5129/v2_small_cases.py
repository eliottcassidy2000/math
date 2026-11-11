#!/usr/bin/env python3
"""Verifier of the corrections to THM-5129 (2026-11-11): independent small-case checks of the statements and
of the proof steps added or reworded after the audit. Own code; reads only slope_lists_le73.txt (to compare).

(1) all labelled tournaments of order <= 5 and a random sample of orders 6, 7: head law, tail law (t = 1, 2, 3;
    with the orientation stated in the Setting), E = U_odd - U_even, U >= D, sink doubling, source padding,
    the deletion identity H(T) - H(T - v) = 2 sum over odd cycles Z through v of H(T - V(Z)), the insertion
    bound, and the product law (A1) with literal strong components.
(2) strong floors f(3..7) by exhaustion (labelled), f(n) >= f(n-1) + 2, the +2 step on every strong
    tournament of order 4..6: some vertex deletion stays strong and loses at least 2 paths.
(3) (A2), both inclusions, for every attained offset D <= 23: slopes predicted from the strong labelled
    tournaments of order <= 6 against (a) the census of all rooted labelled tournaments of order <= 6 and
    (b) explicit ordered joins realising every listed pair (odd part, least exponent); comparison with the
    session's slope_lists_le73.txt for D <= 23.
(4) (A3) on those lists; the bound j >= s - 3; the witnesses of the last section; the order-10 base of
    (55, 352); the four head triples of (C); Pi_2 of the strong tournament of order 4.
"""
import io
import itertools
import random
import sys
from math import log2

ROOT = "C:/tmp/wt-s91/"
LISTS = ROOT + ("04-computation/experiments/s91_f1_local_densities_20261110/rooted_spectrum/"
                "slope_lists_le73.txt")


def tables(out, n):
    """out[i] = bitmask of the vertices beaten by i. Returns h, ends, starts indexed by subset."""
    N = 1 << n
    ends = [[0] * n for _ in range(N)]
    starts = [[0] * n for _ in range(N)]
    h = [0] * N
    h[0] = 1
    inn = [0] * n
    for i in range(n):
        for j in range(n):
            if out[i] >> j & 1:
                inn[j] |= 1 << i
    for S in range(1, N):
        tot = 0
        for v in range(n):
            if not S >> v & 1:
                continue
            R = S ^ (1 << v)
            if R == 0:
                ends[S][v] = 1
                starts[S][v] = 1
            else:
                e = 0
                m = inn[v] & R
                while m:
                    u = (m & -m).bit_length() - 1
                    e += ends[R][u]
                    m &= m - 1
                ends[S][v] = e
                s = 0
                m = out[v] & R
                while m:
                    u = (m & -m).bit_length() - 1
                    s += starts[R][u]
                    m &= m - 1
                starts[S][v] = s
            tot += ends[S][v]
        h[S] = tot
    return h, ends, starts


def popc(x):
    return bin(x).count("1")


def stats(out, n, y, tab=None):
    """D, E, U_even, U_odd, Pi_2 of (B, y)."""
    h, ends, starts = tab if tab else tables(out, n)
    full = (1 << n) - 1
    D = h[full]
    E = ends[full][y]
    Ue = Uo = 0
    for S in range(1, 1 << n):
        if S >> y & 1 and starts[S][y]:
            c = starts[S][y] * h[full ^ S]
            if popc(S) % 2 == 0:
                Ue += c
            else:
                Uo += c
    return D, E, Ue, Uo


def pi2(out, n, tab=None):
    h = (tab if tab else tables(out, n))[0]
    full = (1 << n) - 1
    return sum(h[S] * h[full ^ S] for S in range(1 << n))


def H(out, n):
    return tables(out, n)[0][(1 << n) - 1]


def check_tournament(out, n):
    for i in range(n):
        assert not out[i] >> i & 1
        for j in range(i + 1, n):
            assert (out[i] >> j & 1) != (out[j] >> i & 1), (i, j)


def close_head(out, n, y):
    """B + z: z -> y, every other vertex -> z."""
    o = list(out) + [1 << y]
    for v in range(n):
        if v != y:
            o[v] |= 1 << n
    return o, n + 1


def close_tail(out, n, y, t):
    """B => TT_t => z, z -> y, every other vertex of B -> z."""
    o = list(out)
    m = n + t + 1
    z = n + t
    for v in range(n):
        for k in range(t):
            o[v] |= 1 << (n + k)
        if v != y:
            o[v] |= 1 << z
    for k in range(t):
        row = 1 << z
        for l in range(k + 1, t):
            row |= 1 << (n + l)
        o.append(row)
    o.append(1 << y)
    return o, m


def join(parts):
    """ordered join X_1 => X_2 => ...; parts = list of (out, n). Returns (out, n, offsets)."""
    tot = sum(p[1] for p in parts)
    o = []
    offs = []
    base = 0
    for (po, pn) in parts:
        offs.append(base)
        later = ((1 << tot) - 1) ^ ((1 << (base + pn)) - 1)
        for i in range(pn):
            o.append((po[i] << base) | later)
        base += pn
    return o, tot, offs


def strong_components(out, n):
    """list of vertex lists in the order C_1 => C_2 => ..."""
    reach = [(1 << i) for i in range(n)]
    changed = True
    while changed:
        changed = False
        for i in range(n):
            r = reach[i]
            m = r
            while m:
                u = (m & -m).bit_length() - 1
                r |= out[u] | reach[u]
                m &= m - 1
            if r != reach[i]:
                reach[i] = r
                changed = True
    comps = {}
    for i in range(n):
        key = 0
        for j in range(n):
            if reach[i] >> j & 1 and reach[j] >> i & 1:
                key |= 1 << j
        comps[key] = popc(reach[i])
    order = sorted(comps, key=lambda k: -comps[k])
    return [[j for j in range(n) if k >> j & 1] for k in order]


def sub(out, verts):
    idx = {v: i for i, v in enumerate(verts)}
    o = []
    for v in verts:
        row = 0
        for w in verts:
            if out[v] >> w & 1:
                row |= 1 << idx[w]
        o.append(row)
    return o, len(verts)


def all_tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in range(1 << len(pairs)):
        o = [0] * n
        for k, (i, j) in enumerate(pairs):
            if bits >> k & 1:
                o[i] |= 1 << j
            else:
                o[j] |= 1 << i
        yield o


def random_tournament(n, rng):
    o = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5:
                o[i] |= 1 << j
            else:
                o[j] |= 1 << i
    return o


def odd_cycle_sum(out, n, v, tab):
    """sum over odd directed cycles Z through v of H(T - V(Z))."""
    h = tab[0]
    full = (1 << n) - 1
    # paths starting at v: f[S][end]
    f = {(1 << v): {v: 1}}
    total = 0
    layer = {(1 << v): {v: 1}}
    for size in range(1, n):
        nxt = {}
        for S, d in layer.items():
            for e, c in d.items():
                m = out[e] & ~S & full
                while m:
                    u = (m & -m).bit_length() - 1
                    nxt.setdefault(S | (1 << u), {}).setdefault(u, 0)
                    nxt[S | (1 << u)][u] += c
                    m &= m - 1
        layer = nxt
        if (size + 1) % 2 == 1 and size + 1 >= 3:
            for S, d in layer.items():
                cyc = sum(c for e, c in d.items() if out[e] >> v & 1)
                if cyc:
                    total += cyc * h[full ^ S]
    return total


def main():
    ok = True
    rng = random.Random(20261111)

    # ---------------- (1) laws
    cnt = fails = 0
    nonstrong = 0
    for n in range(1, 8):
        if n <= 5:
            gen = all_tournaments(n)
        else:
            gen = (random_tournament(n, rng) for _ in range(400 if n == 6 else 120))
        for out in gen:
            tab = tables(out, n)
            h = tab[0]
            full = (1 << n) - 1
            comps = strong_components(out, n)
            compdata = []
            if len(comps) > 1:
                nonstrong += 1
                for C in comps:
                    so, sn = sub(out, C)
                    st = tables(so, sn)
                    compdata.append((C, so, sn, st))
                prodH = 1
                for (C, so, sn, st) in compdata:
                    prodH *= st[0][(1 << sn) - 1]
                if prodH != h[full]:
                    fails += 1
            for v in range(n):
                # deletion identity and insertion bound
                if n >= 2:
                    rest = [w for w in range(n) if w != v]
                    hv = h[full ^ (1 << v)]
                    if h[full] - hv != 2 * odd_cycle_sum(out, n, v, tab) or h[full] < hv:
                        fails += 1
            for y in range(n):
                cnt += 1
                D, E, Ue, Uo = stats(out, n, y, tab)
                U = Ue + Uo
                if E != Uo - Ue or U < D:
                    fails += 1
                ho, hn = close_head(out, n, y)
                Hh = H(ho, hn)
                if Hh != D - E + U or Hh != D + 2 * Ue:
                    fails += 1
                for t in (1, 2, 3):
                    if n + t + 1 > 9:
                        continue
                    to, tn = close_tail(out, n, y, t)
                    if H(to, tn) != D + 2 ** t * U:
                        fails += 1
                # appended sink doubles U, keeps D; prepended source keeps both
                so, sn, _ = join([(out, n), ([0], 1)])
                D2, E2, a, b = stats(so, sn, y)
                if (D2, a + b) != (D, 2 * U):
                    fails += 1
                so, sn, offs = join([([0], 1), (out, n)])
                D3, E3, a, b = stats(so, sn, y + 1)
                if (D3, a + b) != (D, U):
                    fails += 1
                # product law with literal components
                if len(comps) > 1:
                    k = next(i for i, (C, so, sn, st) in enumerate(compdata) if y in C)
                    pred = 1
                    for i, (C, so, sn, st) in enumerate(compdata):
                        if i < k:
                            pred *= st[0][(1 << sn) - 1]
                        elif i == k:
                            d, e, a, b = stats(so, sn, C.index(y), st)
                            pred *= a + b
                        else:
                            pred *= pi2(so, sn, st)
                    if pred != U:
                        fails += 1
    print("(1) rooted tournaments checked: %d (all labelled of order <= 5, samples of orders 6, 7; %d non-strong "
          "tournaments): head law, tail law t = 1, 2, 3, E = U_odd - U_even, U >= D, sink doubling, source "
          "padding, deletion identity, insertion bound, product law (A1): failures %d" % (cnt, nonstrong, fails))
    ok &= fails == 0

    # ---------------- (2) strong floors and the +2 step
    strong = {}       # n -> list of out
    f = {1: 1}
    plus2_fail = 0
    for n in range(3, 7):
        lst = []
        for out in all_tournaments(n):
            if len(strong_components(out, n)) == 1:
                lst.append(out)
        strong[n] = lst
        hs = [H(o, n) for o in lst]
        f[n] = min(hs)
        if n >= 4:
            for o, hh in zip(lst, hs):
                good = False
                for v in range(n):
                    rest = [w for w in range(n) if w != v]
                    so, sn = sub(o, rest)
                    if len(strong_components(so, sn)) == 1 and hh - H(so, sn) >= 2 and H(so, sn) >= f[n - 1]:
                        good = True
                        break
                if not good:
                    plus2_fail += 1
    print("(2) strong labelled tournaments of order 3..6: %s; least H: %s; strong tournaments of order 4..6 "
          "without a vertex v with T - v strong and H(T) >= H(T - v) + 2: %d"
          % ([len(strong[n]) for n in range(3, 7)], [f[n] for n in range(3, 7)], plus2_fail))
    ok &= plus2_fail == 0 and [f[n] for n in range(3, 7)] == [3, 5, 9, 15]
    busch = lambda n: {1: 5 ** ((n - 1) // 3), 0: 3 * 5 ** ((n - 3) // 3), 2: 9 * 5 ** ((n - 5) // 3)}[n % 3]
    fl = [busch(n) for n in range(3, 17)]
    print("    Busch formula n = 3..16:", fl, "; steps >= 2:", all(b - a >= 2 for a, b in zip(fl, fl[1:])),
          "; formula at n = 2 would be 9 * 5^(-1) (not an integer)")
    ok &= fl[:10] == [3, 5, 9, 15, 25, 45, 75, 125, 225, 375]
    print("    Moon 1972 lower bound 6^((n-1)/4) at n = 3..12:", [round(6 ** ((n - 1) / 4), 1) for n in range(3, 13)])

    # ---------------- (3) (A2) for D <= 23
    # strong data: H -> {U values}, {Pi_2 values}, representatives
    rootrep = {}   # (H, U) -> (out, n, y)
    sufrep = {}    # (H, Pi2) -> (out, n)
    hrep = {1: ([0], 1)}  # H -> some tournament with H paths (strong or point)
    for n in range(3, 7):
        for o in strong[n]:
            tab = tables(o, n)
            hh = tab[0][(1 << n) - 1]
            if hh > 23:
                continue
            hrep.setdefault(hh, (o, n))
            p = pi2(o, n, tab)
            sufrep.setdefault((hh, p), (o, n))
            for y in range(n):
                D, E, a, b = stats(o, n, y, tab)
                rootrep.setdefault((hh, a + b), (o, n, y))
    strongH = sorted(set(k[0] for k in rootrep))
    print("(3) path counts <= 23 of strong tournaments (order <= 6, complete since f(7) = 25):", strongH)
    # attained numbers <= 23 = products of strong H values, with a representative
    att = dict(hrep)
    changed = True
    while changed:
        changed = False
        for a in list(att):
            for s in strongH:
                if a * s <= 23 and a * s not in att:
                    o, n, _ = join([att[a], hrep[s]])
                    att[a * s] = (o, n)
                    changed = True
    print("    attained odd numbers <= 23:", sorted(att), "(missing:",
          [d for d in range(1, 24, 2) if d not in att], ")")
    # Suf(d): ordered factorisations into strong H values -> products of Pi_2, with representatives
    suf = {1: {1: []}}
    for d in range(3, 24, 2):
        res = {}
        for s in strongH:
            if d % s == 0 and d // s in suf:
                for (hh, p), rep in sufrep.items():
                    if hh != s:
                        continue
                    for w, reps in suf[d // s].items():
                        res.setdefault(p * w, [rep] + reps)
        if res:
            suf[d] = res
    lists = {}
    realised = bad = 0
    for D in sorted(att):
        prim = {}   # slope value -> construction
        for d1 in sorted(att):
            if D % d1:
                continue
            rem = D // d1
            for hh in [1] + strongH:
                if rem % hh or rem // hh not in suf:
                    continue
                roots = [(1, ([0], 1, 0))] if hh == 1 else [(u, rep) for (x, u), rep in rootrep.items() if x == hh]
                for u, rrep in roots:
                    for w, reps in suf[rem // hh].items():
                        prim.setdefault(d1 * u * w, (d1, rrep, reps))
        odd = {}
        for val, cons in prim.items():
            j = 0
            v = val
            while v % 2 == 0:
                v //= 2
                j += 1
            if v not in odd or j < odd[v][0]:
                odd[v] = (j, cons)
        lists[D] = {d: j for d, (j, cons) in odd.items()}
        # (b) realise every (odd part, least exponent) by an explicit ordered join and recount
        for d, (j, (d1, rrep, reps)) in odd.items():
            parts = []
            if d1 != 1:
                parts.append(att[d1])
            parts.append((rrep[0], rrep[1]))
            ypos = len(parts) - 1
            parts.extend(reps)
            o, n, offs = join(parts)
            if n <= 13:
                y = offs[ypos] + rrep[2]
                check_tournament(o, n)
                DD, EE, a, b = stats(o, n, y)
                realised += 1
                if (DD, a + b) != (D, d << j):
                    bad += 1
    print("    explicit ordered joins realising (D, d 2^j) for every listed odd part with its least exponent: "
          "%d built and recounted, mismatches %d" % (realised, bad))
    ok &= bad == 0
    # (a) census of all rooted labelled tournaments of order <= 6 with D <= 23
    census = set()
    for n in range(1, 7):
        for out in all_tournaments(n):
            tab = tables(out, n)
            if tab[0][(1 << n) - 1] > 23:
                continue
            for y in range(n):
                D, E, a, b = stats(out, n, y, tab)
                census.add((D, a + b))
    outside = []
    for (D, U) in census:
        v, j = U, 0
        while v % 2 == 0:
            v //= 2
            j += 1
        if v not in lists[D] or j < lists[D][v]:
            outside.append((D, U))
    print("    census of all rooted labelled tournaments of order <= 6: %d pairs (D, U) with D <= 23; outside "
          "Slopes_0(D) * 2^j: %d" % (len(census), len(outside)))
    ok &= not outside
    # comparison with the session lists
    sess = {}
    with io.open(LISTS, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            head, tail = line.split(": ", 1)
            sess[int(head)] = {int(t.split(":")[0]): int(t.split(":")[1]) for t in tail.split()}
    same = all(lists[D] == sess[D] for D in lists)
    print("    lists for the %d attained offsets <= 23 identical to slope_lists_le73.txt: %s" % (len(lists), same))
    ok &= same

    # ---------------- (4) (A3), bounds, witnesses
    a3 = True
    for D, lst in sess.items():
        for X in (1, 2, 5, 17, 100, 1000, 12345):
            c = 0
            for d, j in lst.items():
                v = d << j
                while v <= X:
                    c += 1
                    v <<= 1
            if c > len(lst) * (1 + log2(X)) + 1e-9:
                a3 = False
        for a in (1, 3, 7, 20, 64, 333):
            c = 0
            for d, j in lst.items():
                v = d << j
                while v < 2 * a:
                    if v >= a:
                        c += 1
                    v <<= 1
            if c > len(lst):
                a3 = False
    print("(4) (A3) on the 35 lists (X = 1..12345, a = 1..333):", a3)
    ok &= a3
    for s in range(2, 13):
        j = 0
        while 9 * 2 ** j < 2 ** s + 1:
            j += 1
        assert j == max(0, s - 3)
        j = 0
        while 2 ** j < 2 ** s + 1:
            j += 1
        assert j == s + 1
    print("    least j with 9 * 2^j >= 2^s + 1 is max(0, s - 3), least j with 2^j >= 2^s + 1 is s + 1 (s = 2..12): True")
    wit = [(7, (444, 273, 258, 6, 268, 286, 319, 382, 8), 8, (129, 144)),
           (8, (82, 36, 129, 135, 142, 157, 190, 3), 7, (257, 288)),
           (9, (402, 360, 3, 261, 270, 285, 317, 382, 4), 8, (513, 576))]
    for s, B, y, pair in wit:
        n = len(B)
        check_tournament(list(B), n)
        D, E, a, b = stats(list(B), n, y)
        U = a + b
        strong_ = len(strong_components(list(B), n)) == 1
        vals = []
        for t in (1, 2):
            to, tn = close_tail(list(B), n, y, t)
            vals.append(H(to, tn))
        good = (D, U) == pair and D == 2 ** s + 1 and U == 9 * 2 ** (s - 3) and vals == [D + 2 * U, D + 4 * U]
        print("    s = %d: order %d, strong %s, (D, U) = (%d, %d), 9 * 2^(s-3) = %d, H closed (t = 1, 2) = %s : %s"
              % (s, n, strong_, D, U, 9 * 2 ** (s - 3), vals, "ok" if good else "FAIL"))
        ok &= good
    # the order-10 base of (55, 352)
    c5 = next(o for o in strong[5] if H(o, 5) == 11)
    c4 = strong[4][0]
    o, n, offs = join([(c5, 5), ([0], 1), (c4, 4)])
    D, E, a, b = stats(o, n, 5)
    print("    C(H = 11, order 5) => {y} => C'(H = 5, order 4): order %d, (D, U) = (%d, %d); Pi_2(C') = %d"
          % (n, D, a + b, pi2(c4, 4)))
    ok &= (D, a + b) == (55, 352) and pi2(c4, 4) == 32 and n == 10
    # the four head triples (audit's out-masks, s9_misc.out (3))
    heads = [((2, 8, 3, 5), 3, (5, 1, 7), 11), ((2, 8, 3, 5), 2, (5, 1, 9), 13),
             ((2, 32, 3, 7, 15, 29), 5, (17, 1, 31), 47), ((12, 9, 10, 48, 7, 23), 0, (33, 7, 41), 67)]
    for B, y, trip, N in heads:
        n = len(B)
        check_tournament(list(B), n)
        D, E, a, b = stats(list(B), n, y)
        ho, hn = close_head(list(B), n, y)
        good = (D, E, a + b) == trip and H(ho, hn) == N == D - E + a + b and E < D
        print("    head %d = %d - %d + %d at a root of out-degree %d of a base of order %d: %s"
              % (N, D, E, a + b, popc(B[y]), n, "ok" if good else "FAIL"))
        ok &= good
    print("ALL CHECKS PASSED" if ok else "SOME CHECK FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
