#!/usr/bin/env python3
"""procgen_tourn_20260924, part B: tournaments on Collatz orbit windows, both sheets.

Window tournament W_b(x; m): vertices x_0..x_{m-1} (x_{j+1} = T_b(x_j), all distinct); path arcs
x_j -> x_{j+1}; every other pair oriented by numerical order (smaller -> larger).  Generic tournament
G_w: the same with the order replaced by the word prediction x_k > x_j iff 3^a > 2^(k-j) (a = odd steps).

Sections (every check raises on failure):
  B0  sheet-blindness of map+order tournaments on a FIXED vertex set V (T_+ and T_- agree on V c Z_{>=3}).
  B1  generic windows: G_{w^R} = r(G_w)^op exactly (word reversal = converse), H/c3/scores per word,
      the Redei tower digit H mod 4 = 1 + 2*alpha_1 (THM-466) on generic windows m <= 9.
  B2  actual windows versus generic: crossing pairs, their time direction (plus: forward, minus: backward;
      the direction lemma), crossing census by window length and start range, both sheets.
  B3  distributions of score sequences, c3, forward-arc count phi and H on complete residue systems:
      identical on the two sheets when no window crosses; the differences sit exactly on crossings.
  B4  Syracuse (odd-iterate) windows through the known gate crossings (minus 165->163 at lag 12, plus
      lag-17 crossings of the order-laws note): H, c3, phi of actual versus generic tournaments.
  B5  typing table of the statistics under negation nu, time reversal tau, converse kappa, word reversal.
Run: python3 04-computation/experiments/procgen_tourn_20260924_windows.py   (about 1-3 minutes, < 400 MB)
"""
import itertools
import math
import sys
import time

import numpy as np


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)


# ---------------------------------------------------------------- tournament utilities (adjacency bitmasks)
def ham_count(out):
    """Number of directed Hamiltonian paths; out[v] = bitmask of out-neighbours. Vectorised Held-Karp."""
    n = len(out)
    if n == 1:
        return 1
    N = 1 << n
    dp = np.zeros((N, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    allm = np.arange(N, dtype=np.int64)
    pc = np.zeros(N, dtype=np.int64)
    for v in range(n):
        pc += (allm >> v) & 1
    layers = [allm[pc == k] for k in range(n + 1)]
    for k in range(1, n):
        L = layers[k]
        for v in range(n):
            sel = L[((L >> v) & 1) == 1]
            vals = dp[sel, v]
            nz = vals != 0
            if not nz.any():
                continue
            sel, vals = sel[nz], vals[nz]
            ov = out[v]
            for u in range(n):
                if not (ov >> u) & 1:
                    continue
                m2 = ((sel >> u) & 1) == 0
                if m2.any():
                    dp[sel[m2] | (1 << u), u] += vals[m2]
    return int(dp[N - 1].sum())


def odd_cycle_count(out):
    """Number of directed cycles of odd length >= 3 (alpha_1 of THM-002/THM-466)."""
    n = len(out)
    total = 0
    for s in range(n):
        # paths starting at s using only vertices > s
        dp = {(1 << s, s): 1}
        frontier = dict(dp)
        for length in range(1, n - s):
            nxt = {}
            for (mask, v), c in frontier.items():
                ov = out[v]
                for u in range(s + 1, n):
                    if (ov >> u) & 1 and not (mask >> u) & 1:
                        key = (mask | (1 << u), u)
                        nxt[key] = nxt.get(key, 0) + c
            for (mask, v), c in nxt.items():
                if (out[v] >> s) & 1 and (length + 1) % 2 == 1 and length + 1 >= 3:
                    total += c
            frontier = nxt
    return total


def scores_of(out):
    return [bin(o).count("1") for o in out]


def c3_of(out):
    n = len(out)
    return n * (n - 1) * (n - 2) // 6 - sum(s * (s - 1) // 2 for s in scores_of(out))


def converse(out):
    n = len(out)
    res = [0] * n
    for v in range(n):
        for u in range(n):
            if (out[v] >> u) & 1:
                res[u] |= 1 << v
    return res


def relabel(out, perm):
    n = len(out)
    res = [0] * n
    for v in range(n):
        for u in range(n):
            if (out[v] >> u) & 1:
                res[perm[v]] |= 1 << perm[u]
    return res


def tournament_from_up(m, up):
    """up[(j,k)] for j<k non-adjacent: True iff arc j->k.  Path arcs j->j+1."""
    out = [0] * m
    for j in range(m - 1):
        out[j] |= 1 << (j + 1)
    for j in range(m):
        for k in range(j + 2, m):
            if up[(j, k)]:
                out[j] |= 1 << k
            else:
                out[k] |= 1 << j
    return out


def is_tournament(out):
    n = len(out)
    for v in range(n):
        if (out[v] >> v) & 1:
            return False
        for u in range(v + 1, n):
            if ((out[v] >> u) & 1) + ((out[u] >> v) & 1) != 1:
                return False
    return True


# ---------------------------------------------------------------- maps and words
def T(x, b):
    return x // 2 if x % 2 == 0 else (3 * x + b) // 2


GROW = {}  # GROW[(a, p)] = 3^a > 2^p (exact)


def grow(a, p):
    key = (a, p)
    if key not in GROW:
        GROW[key] = 3 ** a > 2 ** p
    return GROW[key]


def generic_up(word, m):
    """word = parities e_0..e_{m-2}; returns up[(j,k)] = word prediction x_k > x_j."""
    pre = [0]
    for e in word:
        pre.append(pre[-1] + e)
    return {(j, k): grow(pre[k] - pre[j], k - j) for j in range(m) for k in range(j + 2, m)}


# ---------------------------------------------------------------- B0
def section_B0():
    print("=" * 100)
    print("B0  a fixed vertex set V: map arcs {x, T_b(x)} inside V + ascending order on the other pairs")
    rng = np.random.default_rng(20260924)
    trials = 0
    for _ in range(3000):
        size = int(rng.integers(4, 12))
        V = sorted(set(int(v) for v in rng.integers(3, 60, size=size)))
        # add some images to create map arcs
        extra = [T(v, 1) for v in V[:3]] + [T(v, -1) for v in V[3:6]] + [2 * v for v in V[:2]]
        V = sorted(set(V + [e for e in extra if e >= 3]))
        res = []
        for b in (1, -1):
            arcs = set()
            maps = set()
            for x in V:
                y = T(x, b)
                if y in V and y != x:
                    arcs.add((x, y))
                    maps.add(frozenset((x, y)))
            for x, y in itertools.combinations(V, 2):
                if frozenset((x, y)) not in maps:
                    arcs.add((x, y) if x < y else (y, x))
            res.append(frozenset(arcs))
        check(res[0] == res[1], "map+order tournament on a fixed V is sheet-independent: %r" % V)
        trials += 1
    print("  %d random vertex sets V c [3, ~90]: the tournaments of T_+ and T_- coincide on every V." % trials)
    print("  Reason (PROVED): on positive x every odd step goes up and every even step halves, so each map")
    print("  arc agrees with the order except the halving arcs {x, x/2}, which are the same on both sheets.")
    print("  => a map+order tournament can see the sheet only through the choice of the vertex set.")


# ---------------------------------------------------------------- B1
def section_B1():
    print("=" * 100)
    print("B1  generic window tournaments G_w (the order is the word's multiplicative skeleton)")
    Hcache = {}
    summary = []
    for m in range(3, 13):
        words = list(itertools.product((0, 1), repeat=m - 1))
        stats_H = {}
        rev_ok = 0
        for w in words:
            G = tournament_from_up(m, generic_up(w, m))
            check(is_tournament(G), "tournament")
            wr = tuple(reversed(w))
            GR = tournament_from_up(m, generic_up(wr, m))
            r = [m - 1 - j for j in range(m)]
            check(converse(relabel(G, r)) == GR, "EXACT: G_{w^R} = r(G_w)^op (m=%d, w=%r)" % (m, w))
            rev_ok += 1
            key = tuple(G)
            if key not in Hcache:
                Hcache[key] = ham_count(G) if m <= 12 else None
            h = Hcache[key]
            check(h % 2 == 1, "Redei")
            stats_H[h] = stats_H.get(h, 0) + 1
            if m <= 9:
                a1 = odd_cycle_count(G)
                check((h - 1 - 2 * a1) % 4 == 0, "THM-466 digit: H = 1 + 2 alpha_1 mod 4")
        Hs = sorted(stats_H)
        mean = sum(h * c for h, c in stats_H.items()) / len(words)
        summary.append((m, len(words), rev_ok, Hs[0], Hs[-1], mean, len(stats_H)))
    print("  m  words  G_{w^R}=r(G_w)^op  min H  max H   mean H over words  #distinct H")
    for m, nw, ok, hmin, hmax, mean, nd in summary:
        print("  %2d %6d %8d / %-6d %6d %6d %14.3f %10d" % (m, nw, ok, nw, hmin, hmax, mean, nd))
    print("  PROVED: y_j = a_j log3 - j log2 (skeleton); for the reversed word y^R_j = y_{m-1} - y_{m-1-j}, so")
    print("  the reversed word has the time-reversed, order-reversed skeleton: G_{w^R} = r(G_w)^op, i.e. WORD")
    print("  REVERSAL = CONVERSE for generic windows.  Hence H(G_{w^R}) = H(G_w) and c3 agree (converse-even).")
    print("  Redei (H odd) and the THM-466 digit H = 1 + 2 alpha_1 (mod 4) hold on every generic window, m <= 9.")
    # growth of H: log H / m for the max and mean
    print("  growth: max_w log2 H(G_w) / m = %s" % ", ".join("%d:%.3f" % (m, math.log2(hmax) / m) for m, _, _, _, hmax, _, _ in summary))


# ---------------------------------------------------------------- vectorised windows
def windows(b, xs0, m):
    """Return int64 array (N, m) of T_b-orbit windows and a mask of windows with m distinct values."""
    N = len(xs0)
    W = np.zeros((N, m), dtype=np.int64)
    W[:, 0] = xs0
    for j in range(1, m):
        x = W[:, j - 1]
        W[:, j] = np.where(x % 2 == 0, x // 2, (3 * x + b) // 2)
    # distinctness: sort each row and compare neighbours
    S = np.sort(W, axis=1)
    ok = np.all(S[:, 1:] != S[:, :-1], axis=1)
    check(np.all(np.abs(W) < 2 ** 62), "no overflow")
    return W, ok


def window_stats(b, xs0, m, want_H=False):
    W, ok = windows(b, xs0, m)
    W = W[ok]
    N = W.shape[0]
    par = (W[:, :-1] & 1).astype(np.int16)
    pre = np.zeros((N, m), dtype=np.int16)
    pre[:, 1:] = np.cumsum(par, axis=1, dtype=np.int16)
    score = np.zeros((N, m), dtype=np.int16)
    score[:, :-1] += 1  # path arcs
    phi = np.zeros(N, dtype=np.int16)
    ncross = np.zeros(N, dtype=np.int16)
    bad_dir = 0
    max_cross_start = 0
    for j in range(m):
        for k in range(j + 2, m):
            up = W[:, k] > W[:, j]
            a = pre[:, k] - pre[:, j]
            gen = np.array([grow(aa, k - j) for aa in range(k - j + 1)], dtype=bool)[a]
            cross = up != gen
            if cross.any():
                ncross += cross
                max_cross_start = max(max_cross_start, int(W[cross, j].max()))
                # direction lemma: plus crossings are forward (up), minus crossings backward (down)
                if b == 1:
                    bad_dir += int(np.sum(cross & ~up))
                else:
                    bad_dir += int(np.sum(cross & up))
            score[:, j] += up
            score[:, k] += ~up
            phi += up
    sc = score.astype(np.int64)
    c3 = m * (m - 1) * (m - 2) // 6 - np.sum(sc * (sc - 1) // 2, axis=1)
    return dict(W=W, score=score, phi=phi.astype(np.int64), ncross=ncross, c3=c3, bad_dir=bad_dir,
                nterm=int(np.sum(~ok)), N=N, max_cross_start=max_cross_start)


def out_from_row(Wrow):
    m = len(Wrow)
    out = [0] * m
    for j in range(m - 1):
        out[j] |= 1 << (j + 1)
    for j in range(m):
        for k in range(j + 2, m):
            if Wrow[k] > Wrow[j]:
                out[j] |= 1 << k
            else:
                out[k] |= 1 << j
    return out


# ---------------------------------------------------------------- B2
def section_B2():
    print("=" * 100)
    print("B2  actual windows versus the word prediction: gate crossings and their direction")
    print("    (start x over a range; windows that revisit a value (a cycle) are excluded as terminal)")
    print("  range            m   sheet  windows  terminal  with-crossing  crossing-pairs  max x_j at a crossing  dir.viol.")
    rows = []
    for (lo, hi) in ((3, 3 + 2 ** 16), (10 ** 6, 10 ** 6 + 2 ** 16)):
        for m in (8, 12, 16, 20, 24, 32, 40, 48, 64):
            xs0 = np.arange(lo, hi, dtype=np.int64)
            for b in (1, -1):
                st = window_stats(b, xs0, m)
                wc = int(np.sum(st["ncross"] > 0))
                cp = int(np.sum(st["ncross"]))
                check(st["bad_dir"] == 0, "direction lemma (plus crossings forward, minus backward)")
                rows.append((lo, m, b, st["N"], st["nterm"], wc, cp, st["max_cross_start"]))
                print("  [%7d, +2^16) %3d   %+d   %7d  %8d  %13d  %14d  %21d  %d"
                      % (lo, m, b, st["N"], st["nterm"], wc, cp, st["max_cross_start"], st["bad_dir"]))
    big = [r for r in rows if r[0] == 10 ** 6 and r[1] <= 24]
    check(all(r[5] == 0 for r in big), "no crossing for x ~ 10^6 and windows of <= 24 values")
    tot = {b: sum(r[6] for r in rows if r[2] == b) for b in (1, -1)}
    print("  totals over all rows: %d plus and %d minus (window, crossing pair) incidences, 0 direction violations"
          % (tot[1], tot[-1]))
    mx = max(r[7] for r in rows)
    print("  every crossing pair, in every row, starts at a value x_j <= %d (a window from 10^6 crosses only" % mx)
    print("  after descending into the small-value region where the gates are).")
    print("  Direction lemma (order-laws note 3.1, T-map form): on the plus sheet every crossing pair is a decay")
    print("  pair that goes UP (a forward-in-time arc); on the minus sheet a growth pair that goes DOWN: 0 violations.")
    print("  For starts near 10^6 no window of <= 24 values crosses: those tournaments are word functions.")


# ---------------------------------------------------------------- B3
HCACHE = {}


def section_B3():
    print("=" * 100)
    print("B3  statistic distributions on complete residue systems (each word of length m-1 equally often)")
    print("  x runs over [x0, x0 + K*2^(m-1)); the word of a window is a function of x mod 2^(m-1) (Terras).")
    for (x0, K) in ((10 ** 6, 4), (3, 4)):
        for m in (6, 8, 10, 12, 14):
            xs0 = np.arange(x0, x0 + K * 2 ** (m - 1), dtype=np.int64)
            if m > 12:
                pass  # H not computed (DP cost); scores, c3 and phi still compared
            per = {}
            for b in (1, -1):
                st = window_stats(b, xs0, m)
                scoreseq = [tuple(sorted(r)) for r in st["score"].tolist()]
                d_sc, d_c3, d_phi = {}, {}, {}
                for s, c, p in zip(scoreseq, st["c3"].tolist(), st["phi"].tolist()):
                    d_sc[s] = d_sc.get(s, 0) + 1
                    d_c3[c] = d_c3.get(c, 0) + 1
                    d_phi[p] = d_phi.get(p, 0) + 1
                # H via a cache over distinct tournaments, shared by the two sheets (m <= 12)
                d_H = {}
                seen = set()
                for row in st["W"]:
                    out = tuple(out_from_row(row.tolist()))
                    seen.add(out)
                    if m <= 12:
                        if out not in HCACHE:
                            HCACHE[out] = ham_count(list(out))
                        h = HCACHE[out]
                        d_H[h] = d_H.get(h, 0) + 1
                per[b] = dict(sc=d_sc, c3=d_c3, phi=d_phi, H=d_H, N=st["N"], term=st["nterm"],
                              cross=int(np.sum(st["ncross"] > 0)), phisum=int(st["phi"].sum()),
                              ntour=len(seen))
            def tv(d1, d2, n1, n2):
                keys = set(d1) | set(d2)
                return 0.5 * sum(abs(d1.get(k, 0) / n1 - d2.get(k, 0) / n2) for k in keys)
            p, q = per[1], per[-1]
            tvs = [tv(p[s], q[s], p["N"], q["N"]) for s in ("sc", "c3", "phi", "H")]
            print("  x0=%-7d m=%2d  windows %5d/%5d (terminal %d/%d)  crossing windows %4d/%4d  distinct tournaments %4d/%4d"
                  % (x0, m, p["N"], q["N"], p["term"], q["term"], p["cross"], q["cross"], p["ntour"], q["ntour"]))
            print("      TV distance plus vs minus: scores %.4f  c3 %.4f  phi %.4f  H %s   mean phi %+.4f / %+.4f"
                  % (*tvs[:3], ("%.4f" % tvs[3]) if m <= 12 else "n/a", p["phisum"] / p["N"], q["phisum"] / q["N"]))
            if x0 == 10 ** 6:
                check(p["cross"] == 0 and q["cross"] == 0, "no crossings at 10^6")
                check(all(t == 0 for t in tvs), "identical distributions without crossings")
    print("  => on complete residue systems far from the gates the four statistics have IDENTICAL distributions")
    print("     on the two sheets (TV = 0 exactly).  Near the root (x0 = 3) the differences come from terminal")
    print("     windows (the cycles differ) and from crossings, whose direction is the sign law.")


# ---------------------------------------------------------------- B4
def syracuse_window(x, b, w):
    xs = [x]
    ks = []
    for _ in range(w - 1):
        y = 3 * xs[-1] + b
        k = (y & -y).bit_length() - 1
        ks.append(k)
        xs.append(y >> k)
    return xs, ks


def syracuse_generic_up(ks, w):
    K = [0]
    for k in ks:
        K.append(K[-1] + k)
    return {(i, j): grow(j - i, K[j] - K[i]) for i in range(w) for j in range(i + 2, w)}


def section_B4():
    print("=" * 100)
    print("B4  what a crossing does to H, c3 and phi (actual tournament A versus its word's generic G)")
    signs = {}
    print("  (i) T-map windows from [3, 3 + 2^16) that contain a crossing:")
    print("      sheet  m  windows  dH>0  dH=0  dH<0   dc3>0  dc3=0  dc3<0   dphi = +-#crossings")
    for m in (12, 14, 16):
        for b in (1, -1):
            xs0 = np.arange(3, 3 + 2 ** 16, dtype=np.int64)
            W, ok = windows(b, xs0, m)
            W = W[ok]
            cnt = dict(hp=0, h0=0, hn=0, cp=0, c0=0, cn=0, n=0)
            for row in W.tolist():
                word = [v & 1 for v in row[:-1]]
                gen = generic_up(word, m)
                act = {pq: row[pq[1]] > row[pq[0]] for pq in gen}
                cr = [pq for pq in gen if gen[pq] != act[pq]]
                if not cr:
                    continue
                A, G = tournament_from_up(m, act), tournament_from_up(m, gen)
                dH = ham_count(A) - ham_count(G)
                dc3 = c3_of(A) - c3_of(G)
                dphi = sum(act.values()) - sum(gen.values())
                check(dH % 2 == 0, "Redei: an arc flip changes H by an even number")
                check(dphi == (len(cr) if b == 1 else -len(cr)), "phi moves by +#crossings (plus), -#crossings (minus)")
                cnt["n"] += 1
                cnt["hp" if dH > 0 else "h0" if dH == 0 else "hn"] += 1
                cnt["cp" if dc3 > 0 else "c0" if dc3 == 0 else "cn"] += 1
                signs.setdefault(b, set()).add((dH > 0) - (dH < 0))
            print("       %+d  %2d  %7d  %4d  %4d  %4d   %5d  %5d  %5d   ok" % (
                b, m, cnt["n"], cnt["hp"], cnt["h0"], cnt["hn"], cnt["cp"], cnt["c0"], cnt["cn"]))
    print("  (ii) Syracuse (odd-iterate) windows through the gate crossings of the order-laws note section 4")
    print("       (minus: 165, 309, 549 at lag 12, complete list; plus: the 12 lag-17 windows, complete list):")
    print("      sheet  start  values  crossing pairs  H(actual)  H(generic)        dH   c3 act/gen  phi act/gen")
    cases = [(-1, x, 13) for x in (165, 309, 549)] + \
            [(1, x, 18) for x in (165, 171, 231, 257, 259, 387, 389, 391, 437, 581, 587, 589)]
    for b, x, w in cases:
        xs, ks = syracuse_window(x, b, w)
        check(len(set(xs)) == w, "distinct")
        act = {(i, j): xs[j] > xs[i] for i in range(w) for j in range(i + 2, w)}
        gen = syracuse_generic_up(ks, w)
        cross = sorted(pq for pq in act if act[pq] != gen[pq])
        check(len(cross) >= 1, "window contains a crossing")
        for pq in cross:
            check(act[pq] == (b == 1), "direction lemma on Syracuse windows")
        A = tournament_from_up(w, act)
        G = tournament_from_up(w, gen)
        hA, hG = ham_count(A), ham_count(G)
        check((hA - hG) % 2 == 0, "Redei parity of dH")
        signs.setdefault(b, set()).add((hA > hG) - (hA < hG))
        print("       %+d  %5d   %4d   %-14s %10d  %10d  %+9d    %3d/%-3d     %3d/%-3d"
              % (b, x, w, str(cross), hA, hG, hA - hG, c3_of(A), c3_of(G), sum(act.values()), sum(gen.values())))
    check(signs[1] >= {1, -1} and signs[-1] >= {1, -1}, "dH takes both signs on each sheet")
    print("  => a crossing flips exactly its own arcs: phi (forward order arcs) moves by +1 per plus crossing and")
    print("     -1 per minus crossing (the sign law), while dH is even (Redei) and takes BOTH signs on BOTH")
    print("     sheets (most T-window plus crossings raise H, the Syracuse lag-17 ones lower it or leave it):")
    print("     H and c3 do not carry the direction of the sign law; phi does, by construction.")


# ---------------------------------------------------------------- B5
def section_B5():
    print("=" * 100)
    print("B5  typing of window statistics (generic windows, all words, m = 7 and m = 9)")
    for m in (7, 9):
        ok = dict(H_conv=True, c3_conv=True, sc_conv=True, H_nu_tau=True, c3_nu_tau=True, H_rev=True)
        words = list(itertools.product((0, 1), repeat=m - 1))
        for w in words:
            up = generic_up(w, m)
            G = tournament_from_up(m, up)
            Gc = converse(G)
            # nu: order reversed, path kept
            Gnu = tournament_from_up(m, {p: not v for p, v in up.items()})
            # tau: path reversed, order kept
            Gtau = [0] * m
            for j in range(m - 1):
                Gtau[j + 1] |= 1 << j
            for (j, k), v in up.items():
                if v:
                    Gtau[j] |= 1 << k
                else:
                    Gtau[k] |= 1 << j
            hG, hc, hnu, htau = ham_count(G), ham_count(Gc), ham_count(Gnu), ham_count(Gtau)
            ok["H_conv"] &= (hG == hc)
            ok["c3_conv"] &= (c3_of(G) == c3_of(Gc))
            ok["sc_conv"] &= (sorted(scores_of(Gc)) == sorted(m - 1 - s for s in scores_of(G)))
            ok["H_nu_tau"] &= (hnu == htau)
            ok["c3_nu_tau"] &= (c3_of(Gnu) == c3_of(Gtau))
            ok["H_rev"] &= (hG == ham_count(tournament_from_up(m, generic_up(tuple(reversed(w)), m))))
            npairs = len(up)
            phiG = sum(up.values())
            check(sum(1 for (j, k) in up if (Gnu[j] >> k) & 1) == npairs - phiG, "phi o nu = npairs - phi")
            check(sum(1 for (j, k) in up if (Gtau[j] >> k) & 1) == phiG, "phi o tau = phi")
            check(sum(1 for (j, k) in up if (Gc[j] >> k) & 1) == npairs - phiG, "phi o kappa = npairs - phi")
        check(ok["H_conv"] and ok["c3_conv"] and ok["sc_conv"] and ok["H_nu_tau"] and ok["c3_nu_tau"] and ok["H_rev"],
              "typing identities")
        print("  m=%d: H(G^op)=H(G), c3(G^op)=c3(G), scores(G^op)=m-1-scores(G), H o nu = H o tau, c3 o nu = c3 o tau,"
              " phi o nu = phi o kappa = #pairs - phi, phi o tau = phi, H(G_{w^R}) = H(G_w): all words" % m)
    print("  side separation by H, word by word: words w with H(nu G_w) = H(G_w) (nu G_w = the same word on the")
    print("  other side of 0; H(nu G) = H(tau G) since the two are converse):")
    exc = []
    mod4 = set()
    for m in range(3, 12):
        eq = []
        for w in itertools.product((0, 1), repeat=m - 1):
            up = generic_up(w, m)
            h1 = ham_count(tournament_from_up(m, up))
            h2 = ham_count(tournament_from_up(m, {pq: not v for pq, v in up.items()}))
            mod4.add((h1 % 4, h2 % 4))
            if h1 == h2:
                eq.append("".join(map(str, w)))
        exc += eq
        print("    m=%2d: %4d words, equal H on both sides: %s" % (m, 2 ** (m - 1), eq if eq else "none"))
    check(exc == ["101", "01110"], "only 101 and 01110 have side-equal H for m <= 11")
    check(mod4 == {(1, 1), (1, 3), (3, 1), (3, 3)}, "all four (H mod 4, H(nu G) mod 4) pairs occur")
    print("  => for m <= 11 the Hamiltonian-path count of a generic window differs between the two sides of 0 for")
    print("     every word except 101 and 01110 (FINITE-EXACT; not a mod-4 law: all four pairs of H mod 4 occur).")
    print("  Typing (PROVED, from converse = nu o tau on window tournaments):")
    print("    * every window statistic uses the order, so it is side-aware in principle (nu changes it);")
    print("    * H, c3, the OCF data alpha_k are converse-EVEN: f o nu = f o tau (negation = time reversal for them);")
    print("    * scores are converse-ODD; phi (forward order arcs on time labels) is nu-odd, tau-even, kappa-odd;")
    print("    * on a fixed side all of them are word functions away from the gates (B2, B3), hence sheet-blind;")
    print("      they separate the sheets only through crossings, whose direction is the sign law (B4).")


def main():
    t0 = time.time()
    section_B0()
    section_B1()
    section_B2()
    section_B3()
    section_B4()
    section_B5()
    print("=" * 100)
    print("PART B: ALL CHECKS PASSED")
    print("(time %.0f s)" % (time.time() - t0), file=sys.stderr)


if __name__ == "__main__":
    main()
