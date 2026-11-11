#!/usr/bin/env python3
"""Auditor numbers-vs-outputs: independent recomputation of the numbers of THM-5129 (A)-(D) and of the claims that
have no committed output line (number of strong classes per offset, correlation 0.68, 'a base of order 10' at
offset 55, uniqueness of the level 7, 8, 9 witnesses).

Own kernels (written for the audit; they share only the data file reps_le8.npz and its pair coding with the
session's scripts: for i < j the bit j(j-1)/2 + i of the code is 1 iff i -> j).

Usage: python3 verify_supply.py PATH/reps_le8.npz
"""
import sys
from collections import defaultdict
from math import factorial

import numpy as np
from numba import njit

DMAX = 73


@njit(cache=False)
def tables(n, out, hS, st, en):
    """hS[S] = H(T[S]); st[S, v] = paths of T[S] starting at v; en[S, v] = paths of T[S] ending at v."""
    full = (1 << n) - 1
    hS[0] = 1
    for S in range(1, full + 1):
        tot = 0
        for v in range(n):
            st[S, v] = 0
            en[S, v] = 0
            if (S >> v) & 1:
                R = S ^ (1 << v)
                if R == 0:
                    st[S, v] = 1
                    en[S, v] = 1
                else:
                    a = 0
                    b = 0
                    for w in range(n):
                        if (R >> w) & 1:
                            if (out[v] >> w) & 1:      # v -> w: path v, then a path of R starting at w
                                a += st[R, w]
                            else:                      # w -> v: a path of R ending at w, then v
                                b += en[R, w]
                    st[S, v] = a
                    en[S, v] = b
                tot += st[S, v]
        hS[S] = tot


@njit(cache=False)
def decode(code, n, out):
    for i in range(n):
        out[i] = 0
    for j in range(1, n):
        base = j * (j - 1) // 2
        for i in range(j):
            if (code >> (base + i)) & 1:
                out[i] |= (1 << j)
            else:
                out[j] |= (1 << i)


@njit(cache=False)
def is_strong_scores(n, out):
    """Landau-type criterion: T is strong iff the k smallest scores sum to more than C(k, 2) for every k < n."""
    sc = np.zeros(n, np.int64)
    for v in range(n):
        c = 0
        m = out[v]
        while m:
            c += 1
            m &= m - 1
        sc[v] = c
    sc = np.sort(sc)
    tot = 0
    for k in range(1, n):
        tot += sc[k - 1]
        if tot <= k * (k - 1) // 2:
            return False
    return True


@njit(cache=False)
def class_data(n, reps, res):
    """res[r] = (H, Pi_2, strong, U(0..n-1), E(0..n-1)) for every class representative."""
    full = (1 << n) - 1
    hS = np.zeros(1 << n, np.int64)
    st = np.zeros((1 << n, n), np.int64)
    en = np.zeros((1 << n, n), np.int64)
    out = np.zeros(n, np.int64)
    for r in range(reps.shape[0]):
        decode(reps[r], n, out)
        tables(n, out, hS, st, en)
        res[r, 0] = hS[full]
        p = 0
        for S in range(full + 1):
            p += hS[S] * hS[full ^ S]
        res[r, 1] = p
        res[r, 2] = 1 if is_strong_scores(n, out) else 0
        for y in range(n):
            u = 0
            for S in range(1, full + 1):
                if (S >> y) & 1:
                    u += st[S, y] * hS[full ^ S]
            res[r, 3 + y] = u
            res[r, 3 + n + y] = en[full, y]


@njit(cache=False)
def product_law_check(n, reps):
    """(A1) on every root of every class: returns the number of failures and the number of rooted tournaments."""
    full = (1 << n) - 1
    hS = np.zeros(1 << n, np.int64)
    st = np.zeros((1 << n, n), np.int64)
    en = np.zeros((1 << n, n), np.int64)
    out = np.zeros(n, np.int64)
    comp = np.zeros(n, np.int64)          # component mask of each vertex
    reach = np.zeros(n, np.int64)
    bad = 0
    cnt = 0
    for r in range(reps.shape[0]):
        decode(reps[r], n, out)
        tables(n, out, hS, st, en)
        # reachability closure
        for v in range(n):
            reach[v] = out[v] | (1 << v)
        changed = True
        while changed:
            changed = False
            for v in range(n):
                m = reach[v]
                new = m
                for w in range(n):
                    if (m >> w) & 1:
                        new |= reach[w]
                if new != m:
                    reach[v] = new
                    changed = True
        for v in range(n):
            c = 0
            for w in range(n):
                if ((reach[v] >> w) & 1) and ((reach[w] >> v) & 1):
                    c |= (1 << w)
            comp[v] = c
        D = hS[full]
        for y in range(n):
            cnt += 1
            u_direct = 0
            for S in range(1, full + 1):
                if (S >> y) & 1:
                    u_direct += st[S, y] * hS[full ^ S]
            # product law
            prodD = 1
            prodU = 1
            seen = 0
            for v in range(n):
                c = comp[v]
                if seen & c:
                    continue
                seen |= c
                hc = hS[c]
                prodD *= hc
                if c == comp[y]:
                    # U of the root component
                    u = 0
                    sub = c
                    while True:
                        if (sub >> y) & 1:
                            u += st[sub, y] * hS[c ^ sub]
                        if sub == 0:
                            break
                        sub = (sub - 1) & c
                    prodU *= u
                elif (reach[v] >> y) & 1:
                    # component before the root (it reaches y)
                    prodU *= hc
                else:
                    p2 = 0
                    sub = c
                    while True:
                        p2 += hS[sub] * hS[c ^ sub]
                        if sub == 0:
                            break
                        sub = (sub - 1) & c
                    prodU *= p2
            if prodD != D or prodU != u_direct:
                bad += 1
    return bad, cnt


@njit(cache=False)
def aut_order(n, out):
    """|Aut(T)| by backtracking over vertex images."""
    perm = np.zeros(n, np.int64)
    used = np.zeros(n, np.int64)
    count = 0
    # iterative backtracking
    pos = 0
    perm[0] = -1
    while pos >= 0:
        perm[pos] += 1
        while perm[pos] < n:
            c = perm[pos]
            ok = used[c] == 0
            if ok:
                for i in range(pos):
                    a = (out[i] >> pos) & 1
                    b = (out[perm[i]] >> c) & 1
                    if a != b:
                        ok = False
                        break
            if ok:
                break
            perm[pos] += 1
        if perm[pos] >= n:
            pos -= 1
            if pos >= 0:
                used[perm[pos]] = 0
            continue
        if pos == n - 1:
            count += 1
            continue
        used[perm[pos]] = 1
        pos += 1
        perm[pos] = -1
    return count


@njit(cache=False)
def mass(n, reps):
    out = np.zeros(n, np.int64)
    tot = 0
    f = 1
    for i in range(2, n + 1):
        f *= i
    for r in range(reps.shape[0]):
        decode(reps[r], n, out)
        tot += f // aut_order(n, out)
    return tot


@njit(cache=False)
def scan_targets(n, reps, targets, hits, maxhits):
    """rooted tournaments of order n + 1 = (class representative of order n, cut), root = the new vertex z, which
    beats exactly the vertices of `cut`.  Records every (r, cut) whose (H(B + z), U_{B+z}(z)) is a target."""
    full = (1 << n) - 1
    hS = np.zeros(1 << n, np.int64)
    st = np.zeros((1 << n, n), np.int64)
    en = np.zeros((1 << n, n), np.int64)
    out = np.zeros(n, np.int64)
    Q = np.zeros((n, n), np.int64)
    ym = np.zeros(n, np.int64)
    nh = 0
    for r in range(reps.shape[0]):
        decode(reps[r], n, out)
        tables(n, out, hS, st, en)
        for a in range(n):
            ym[a] = 0
            for b in range(n):
                Q[a, b] = 0
        for S in range(1, full + 1):
            C = full ^ S
            for a in range(n):
                if (S >> a) & 1:
                    ym[a] += st[S, a] * hS[C]            # covers (P, R): R starts at a and covers S
                    if en[S, a]:
                        for b in range(n):
                            if (C >> b) & 1:
                                Q[a, b] += en[S, a] * st[C, b]   # path of S ending at a, path of C starting at b
        HB = hS[full]
        for cut in range(1 << n):
            H = 0
            U = HB
            for b in range(n):
                if (cut >> b) & 1:          # z -> b
                    H += st[full, b]        # z first
                    U += ym[b]
                else:                       # b -> z
                    H += en[full, b]        # z last
                    for c in range(n):
                        if (cut >> c) & 1:
                            H += Q[b, c]    # ... b, z, c ...
            for k in range(targets.shape[0]):
                if targets[k, 0] == H and targets[k, 1] == U:
                    if nh < maxhits:
                        hits[nh, 0] = k
                        hits[nh, 1] = r
                        hits[nh, 2] = cut
                    nh += 1
    return nh


def literal_U(out, n, y):
    """U from the definition, pure Python (enumeration of all directed paths from y)."""
    def ham(mask):
        vs = [v for v in range(n) if (mask >> v) & 1]
        if not vs:
            return 1
        cnt = 0
        stack = [(v, 1 << v) for v in vs]
        while stack:
            v, used = stack.pop()
            if used == mask:
                cnt += 1
                continue
            m = out[v] & mask & ~used
            while m:
                w = (m & -m).bit_length() - 1
                m &= m - 1
                stack.append((w, used | (1 << w)))
        return cnt
    full = (1 << n) - 1
    total = 0
    stack = [(y, 1 << y)]
    while stack:
        v, used = stack.pop()
        total += ham(full ^ used)
        m = out[v] & ~used & full
        while m:
            w = (m & -m).bit_length() - 1
            m &= m - 1
            stack.append((w, used | (1 << w)))
    return total


def py_decode(code, n):
    out = [0] * n
    for j in range(1, n):
        base = j * (j - 1) // 2
        for i in range(j):
            if (code >> (base + i)) & 1:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
    return out


def odd_part(u):
    j = (u & -u).bit_length() - 1
    return u >> j, j


def main():
    z = np.load(sys.argv[1])
    reps = {n: z["n%d" % n].astype(np.int64) for n in range(1, 9)}
    print("class counts n = 1..8:", [len(reps[n]) for n in range(1, 9)])
    # mass formula: the representatives are a complete system iff sum n!/|Aut| = 2^C(n,2) (given the count)
    for n in range(3, 9):
        m = int(mass(n, reps[n]))
        print(f"  mass formula n = {n}: sum n!/|Aut| = {m}  (2^C(n,2) = {1 << (n * (n - 1) // 2)})  {'ok' if m == 1 << (n * (n - 1) // 2) else 'MISMATCH'}")
    data = {}
    for n in range(1, 9):
        res = np.zeros((len(reps[n]), 3 + 2 * n), np.int64)
        class_data(n, reps[n], res)
        data[n] = res
    # literal check of U on every root of every class of order <= 6
    lit_bad = 0
    lit_cnt = 0
    for n in range(1, 7):
        for r, code in enumerate(reps[n].tolist()):
            out = py_decode(code, n)
            for y in range(n):
                lit_cnt += 1
                if literal_U(out, n, y) != int(data[n][r, 3 + y]):
                    lit_bad += 1
    print(f"U (DP) = U (literal enumeration of the qualifying paths) on {lit_cnt} rooted tournaments of order <= 6: failures {lit_bad}")
    # product law (A1) on every root of every class of order <= 8
    tot_bad = 0
    tot_cnt = 0
    for n in range(1, 9):
        bad, cnt = product_law_check(n, reps[n])
        tot_bad += bad
        tot_cnt += cnt
    print(f"product law (A1) on every root of every class of order <= 8: {tot_cnt} rooted tournaments, failures {tot_bad}")
    # U >= D, U = E mod 2
    viol = 0
    for n in range(1, 9):
        res = data[n]
        for y in range(n):
            viol += int(np.sum(res[:, 3 + y] < res[:, 0])) + int(np.sum((res[:, 3 + y] - res[:, 3 + n + y]) % 2 != 0))
    print("violations of U >= D or U = E (mod 2), orders <= 8:", viol)
    # strong data
    print("(1) order: classes, strong classes, least H strong, #(H, Pi_2), #(H, U)")
    root = defaultdict(dict)     # h -> {u: least order}
    pis = defaultdict(dict)      # h -> {pi2: least order}
    nstrong = defaultdict(int)   # h -> number of strong classes
    strong_orders = defaultdict(set)
    for n in range(3, 9):
        res = data[n]
        st = res[res[:, 2] == 1]
        pp = set()
        pu = set()
        for row in st.tolist():
            h = row[0]
            nstrong[h] += 1
            strong_orders[h].add(n)
            pp.add((h, row[1]))
            if row[1] not in pis[h] or pis[h][row[1]] > n:
                pis[h][row[1]] = n
            for y in range(n):
                u = row[3 + y]
                pu.add((h, u))
                if u not in root[h] or root[h][u] > n:
                    root[h][u] = n
        print(f"    n = {n}: {len(res)} classes, {len(st)} strong, least H = {int(st[:, 0].min())}, {len(pp)} pairs (H, Pi_2), {len(pu)} pairs (H, U)")
    strong_vals = sorted(h for h in root if h <= DMAX)
    # attained numbers (products of strong values) with least order; suffix sets with least order
    INF = 10 ** 9
    pre = {1: 0}
    suf = {1: {1: 0}}
    for d in range(3, DMAX + 1, 2):
        best = INF
        s = {}
        for h in strong_vals:
            if d % h == 0:
                e = d // h
                if e in pre:
                    hmin = min(strong_orders[h])
                    best = min(best, pre[e] + hmin)
                if e in suf:
                    for w, o in suf[e].items():
                        for p, n in pis[h].items():
                            v = p * w
                            if v not in s or s[v] > o + n:
                                s[v] = o + n
        if best < INF:
            pre[d] = best
        if s:
            suf[d] = s
    missing_attained = sorted(set(range(1, DMAX + 1, 2)) - set(pre))
    print("    attained odd numbers <= 73: all except", missing_attained)
    slopes0 = {}
    for D in sorted(pre):
        s = {}
        for d1, o1 in pre.items():
            if D % d1:
                continue
            rest = D // d1
            for h in [1] + strong_vals:
                if rest % h:
                    continue
                d3 = rest // h
                if d3 not in suf:
                    continue
                us = {1: 1} if h == 1 else root[h]
                for u, on in us.items():
                    for w, o3 in suf[d3].items():
                        v = d1 * u * w
                        o = o1 + on + o3
                        if v not in s or s[v] > o:
                            s[v] = o
        slopes0[D] = s
    prim = {}
    for D, s in slopes0.items():
        m = {}
        for u in s:
            dl, j = odd_part(u)
            if dl not in m or j < m[dl]:
                m[dl] = j
        prim[D] = m

    def min_order(D, dl, j):
        """least order of a base with (H, U) = (D, dl 2^j)"""
        best = INF
        for u, o in slopes0[D].items():
            d2, i = odd_part(u)
            if d2 == dl and i <= j:
                best = min(best, o + (j - i))
        return best

    def is_slope(D, u):
        if D not in prim:
            return False
        dl, j = odd_part(u)
        return dl in prim[D] and j >= prim[D][dl]

    print("(3) odd parts at the Fermat-type offsets (d:j), with the least order of a base in brackets when it is >= 9")
    for s in range(1, 7):
        D = (1 << s) + 1
        m = prim[D]
        cells = []
        for dl, j in sorted(m.items()):
            o = min_order(D, dl, j)
            cells.append(f"{dl}:{j}" + (f"[{o}]" if o >= 9 else ""))
        print(f"    D = {D}: {len(m)} odd parts, {sum(1 for j in m.values() if j <= s)} with j <= s")
        print("        " + " ".join(cells))
    # primitive slopes (odd part, least exponent) whose least base has order >= 10, for all D <= 73
    late = []
    hist = defaultdict(int)
    for D in sorted(prim):
        for dl, j in sorted(prim[D].items()):
            o = min_order(D, dl, j)
            hist[o] += 1
            if o >= 10:
                late.append((D, dl, j, o))
    print("(2') primitive slopes (D <= 73, odd part with its least exponent) by least order of a base:", dict(sorted(hist.items())))
    print("     those needing a base of order >= 10:", late)
    # which primitive slopes need order exactly 9 (so the order <= 9 census is needed to see them)
    # direct census of order <= 8 against the product law
    direct = defaultdict(set)
    for n in range(1, 9):
        res = data[n]
        for row in res.tolist():
            if row[0] <= DMAX:
                for y in range(n):
                    direct[row[0]].add(row[3 + y])
    bad = [(d, u) for d in direct for u in direct[d] if not is_slope(d, u)]
    print(f"     direct census of order <= 8: {sum(len(v) for v in direct.values())} pairs (D, U) with D <= 73; not in the product-law set: {len(bad)}")
    notseen = []
    for D in sorted(prim):
        for dl, j in sorted(prim[D].items()):
            o = min_order(D, dl, j)
            if o <= 8 and (dl << j) not in direct[D]:
                notseen.append((D, dl, j, o))
    print("     primitive slopes with predicted least order <= 8 that are absent from the direct census of order <= 8:", notseen)
    # (4)
    NMAX = 225
    r = {}
    for N in range(3, NMAX, 2):
        r[N] = sum(1 for D in range(1, N // 3 + 1, 2) if is_slope(D, (N - D) // 2))
    zero = [N for N in r if r[N] == 0]
    print("(4) odd N in [3, 225) with r(N) = 0:", zero)
    print("    r(N), N = 69..131:", " ".join(str(r[N]) for N in range(69, 133, 2)))
    vals = [r[N] for N in range(129, NMAX, 2)]
    print(f"    least r on [129, 225): {min(vals)} (at N = {[N for N in range(129, NMAX, 2) if r[N] == min(vals)]}); mean {np.mean(vals):.4f}")
    # (5)
    print("(5) windows [h, 3h/2]")
    share = {}
    cells = []
    for h in sorted(prim):
        hi = (3 * h) // 2
        c = sum(1 for u in range(h, hi + 1) if is_slope(h, u))
        share[h] = (c, hi - h + 1)
        cells.append(f"{h}:{c}/{hi - h + 1}")
    print("    " + " ".join(cells))
    print("    whole window present:", [h for h in share if share[h][0] == share[h][1]])
    hs = [h for h in sorted(prim) if 37 <= h <= 73]
    fr = [share[h][0] / share[h][1] for h in hs]
    print(f"    37 <= h <= 73: mean {np.mean(fr):.4f}, min {min(fr):.4f} at h = {hs[int(np.argmin(fr))]}, max {max(fr):.4f} at h = {hs[int(np.argmax(fr))]}")
    # strong classes per h and correlation
    print("(6) number of strong classes with H = h (all orders; complete for h <= 73 since f(9) = 75), covered share of [h, 3h/2], |Delta(h)|")
    for h in sorted(prim):
        if h >= 3:
            print(f"    h = {h:2d}: strong classes {nstrong.get(h, 0):3d} (orders {sorted(strong_orders.get(h, []))}), share {share[h][0]}/{share[h][1]} = {share[h][0] / share[h][1]:.3f}, odd parts {len(prim[h])}")
    for lo in (23, 15, 3, 37):
        hs = [h for h in sorted(prim) if lo <= h <= 73]
        a = np.array([share[h][0] / share[h][1] for h in hs])
        b = np.array([nstrong.get(h, 0) for h in hs], dtype=float)
        c1 = float(np.corrcoef(a, b)[0, 1])
        ra = np.argsort(np.argsort(a)).astype(float)
        rb = np.argsort(np.argsort(b)).astype(float)
        c2 = float(np.corrcoef(ra, rb)[0, 1])
        c3 = float(np.corrcoef(a, np.log(np.maximum(b, 1)))[0, 1])
        print(f"    {lo} <= h <= 73 ({len(hs)} attained h): Pearson correlation(share, #strong classes) = {c1:.4f}; rank correlation (no tie correction) {c2:.4f}; with log(#strong) {c3:.4f}")
    two = [h for h in sorted(prim) if h >= 3 and nstrong.get(h, 0) == 2]
    one = [h for h in sorted(prim) if h >= 3 and nstrong.get(h, 0) == 1]
    none = [h for h in sorted(prim) if h >= 3 and nstrong.get(h, 0) == 0]
    print("    attained h <= 73 with exactly two strong classes:", two)
    print("    with exactly one:", one, " with none (products only):", none)
    order_by_share = sorted((share[h][0] / share[h][1], h) for h in sorted(prim) if h >= 23)
    print("    23 <= h <= 73 sorted by covered share (poorest first):", [(h, round(s, 3), nstrong.get(h, 0)) for s, h in order_by_share[:10]])
    # uniqueness of the witnesses at levels 7, 8, 9
    print("(7) rooted tournaments (class representative of order m - 1, cut; root = new vertex) with a given (D, U)")
    T = np.array([[129, 144], [257, 288], [513, 576], [129, 288], [65, 72], [65, 288], [5, 9], [9, 9], [17, 18], [33, 36]], dtype=np.int64)
    for n in range(3, 9):
        hits = np.zeros((4096, 3), np.int64)
        nh = int(scan_targets(n, reps[n], T, hits, 4096))
        per = defaultdict(list)
        for k, rr, cut in hits[:min(nh, 4096)].tolist():
            per[k].append((rr, cut))
        for k in sorted(per):
            lst = per[k]
            # isomorphism classes of rooted tournaments: same base class and cuts in one Aut(base) orbit
            bases = defaultdict(list)
            for rr, cut in lst:
                bases[rr].append(cut)
            ncls = 0
            for rr, cuts in bases.items():
                out = py_decode(int(reps[n][rr]), n)
                auts = automorphisms(out, n)
                seen = set()
                for cut in cuts:
                    if cut in seen:
                        continue
                    ncls += 1
                    for p in auts:
                        img = 0
                        for v in range(n):
                            if (cut >> v) & 1:
                                img |= 1 << p[v]
                        seen.add(img)
            print(f"    order {n + 1}: (D, U) = ({T[k, 0]}, {T[k, 1]}): {len(lst)} (class, cut) pairs on {len(bases)} base classes; rooted tournaments up to isomorphism: {ncls}")


def automorphisms(out, n):
    import itertools
    res = []
    score = [bin(out[v]).count("1") for v in range(n)]
    for p in itertools.permutations(range(n)):
        if any(score[v] != score[p[v]] for v in range(n)):
            continue
        ok = True
        for i in range(n):
            for j in range(i + 1, n):
                if ((out[i] >> j) & 1) != ((out[p[i]] >> p[j]) & 1):
                    ok = False
                    break
            if not ok:
                break
        if ok:
            res.append(p)
    return res


if __name__ == "__main__":
    main()
