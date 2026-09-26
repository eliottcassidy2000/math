#!/usr/bin/env python3
"""procgen_landing_20260926_run.py -- single runner of the `landing` lane (session
collatz-procgen-20260922, 2026-09-26).  Every printed claim is a check(...) that raises on
failure; the transcript is written to 05-knowledge/results/procgen_landing_20260926.out.

Sections
  A  constants, Theorem 1 (a*(mu) and saturation of the recursion)
  B  exhaustive worst case (C scan, all y <= 2^(k+1)-1): Proposition U is attained exactly
  C  Sturmian hover-then-halve integers (constructed hostile segments, k <= 400)
  D  Proposition H (a residue class of heavy dippers with a uniform margin)
  E  Lemma A (depth averaging) on actual segments
  F  delay and path records <= 10^12 (OEIS A006877, A006884)
  G  the collatz-landing probe segments; climb-then-drop; a W-bit hover
  H  2-adic hostile word and its integer realisation
Usage: python3 -u procgen_landing_20260926_run.py
"""
import hashlib
import math
import os
import random
import resource
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
import procgen_landing_20260926_lib as P  # noqa: E402

SCR = os.path.join(ROOT, "scratch", "procgen_landing")
OUT = os.path.join(ROOT, "05-knowledge", "results", "procgen_landing_20260926.out")
FALLBACK = [os.path.join(ROOT, "scratch", "procgen_family27", "oeis_cache")]
LINES = []
T0 = time.time()


def emit(s=""):
    print(s, flush=True)
    LINES.append(s)


def check(cond, msg):
    if not cond:
        raise AssertionError("CHECK FAILED: " + msg)
    emit("ok: " + msg)


def sha(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


# ----------------------------------------------------------------------------------------
def section_A():
    emit("== A. constants and Theorem 1 ==")
    emit("alpha = %.7f, h* = %.7f, lambda* = %.7f, lambda*/h* = %.6f" % (P.ALPHA, P.H_STAR, P.LAMBDA_STAR, P.LAMBDA_STAR / P.H_STAR))
    check(abs(P.H_STAR - 0.9499555) < 5e-7 and abs(P.LAMBDA_STAR - 0.488077) < 5e-7, "h* = 0.9499555 and lambda* = 0.488077 (THM-4499)")
    for mu, want in ((1.0, -0.986211), (0.5, -1.243105), (0.0, -1.5)):
        emit("a*(%.2f) = mu lambda*/h* - 3/2 = %.6f" % (mu, P.a_star(mu)))
        check(abs(P.a_star(mu) - want) < 1e-6, "a*(%.1f) = %.6f" % (mu, want))
    check(abs(P.a_star(0.5) - (P.LAMBDA_STAR / (2 * P.H_STAR) - 1.5)) < 1e-12,
          "multiplicity L^(1/2) gives lambda*/(2h*) - 3/2 (THM-4499 remark confirmed)")
    worst = {}
    for mu in (0.0, 0.25, 0.5, 0.75, 1.0):
        m = min(P.saturation_margin(mu, L, D) for L in (16, 64, 256, 1024, 2 ** 14, 2 ** 20)
                for D in list(range(1, min(L // 2, 400))) + [L // 2 - 1])
        worst[mu] = m
        emit("saturation mu=%.2f: min over L in {2^4..2^20}, 1 <= D < L/2 of log2(RHS/psi) = %.4f" % (mu, m))
    check(all(v >= 0 for v in worst.values()),
          "Theorem 1(ii): psi(X) = X^h* L^(a*(mu)) satisfies the recursion (R_mu) at every depth D, K=C=1, E=C^(-lambda*/h*)")
    f = [(2 ** ((W - 1) * P.H_STAR) / W, W) for W in range(1, 40)]
    emit("Lemma A: effective multiplicity factor 2^((W-1)h*)/W of the best of W depths: " + ", ".join("W=%d: %.4f" % (W, v) for v, W in f[:4]))
    check(min(f)[1] == 2 and abs(min(f)[0] - 0.965907) < 1e-6,
          "depth averaging lowers the effective multiplicity by at most the factor 0.965907 (W = 2): mu unchanged")


# ----------------------------------------------------------------------------------------
def parse_scan(path):
    rows = []
    for line in open(path):
        d = {}
        for kv in line.split():
            key, val = kv.split("=", 1)
            d[key] = val
        rows.append({"b": int(d["b"]), "k": int(d["k"]), "D2": int(d["D2"]), "max": int(d["max"]),
                     "nmax": int(d["nmax"]), "argy": int(d["argy"]), "viol": int(d["viol"]),
                     "skip": int(d["skip"]), "H": [int(x) for x in d["H"].split(",")] if d["H"] else []})
    return rows


def section_B():
    emit("== B. exhaustive worst case (C scan) ==")
    src = os.path.join(HERE, "procgen_landing_20260926_scan.c")
    exe = os.path.join(SCR, "scan")
    subprocess.run(["cc", "-O2", "-o", exe, src, "-lm"], check=True)
    jobs = [("1", "10", "25", "2,3,4,5,6,7,8,10,12"), ("-1", "10", "25", "2,3,4,5,6,7,8,10,12"),
            ("5", "10", "22", "2,3,4,6,8,12"), ("-5", "10", "22", "2,3,4,6,8,12")]
    rows = []
    for pair in (jobs[:2], jobs[2:]):                     # at most two processes at a time
        procs = []
        for j in pair:
            outp = os.path.join(SCR, "scan_b%s.txt" % j[0].replace("-", "m"))
            procs.append((outp, subprocess.Popen([exe] + list(j), stdout=open(outp, "w"))))
        for outp, pr in procs:
            check(pr.wait() == 0, "scan %s finished" % os.path.basename(outp))
            rows += parse_scan(outp)
    check(len(rows) == 144 + 144 + 78 + 78, "444 (b, k, D) cells scanned (b = 1, -1: k = 10..25; b = 5, -5: k = 10..22)")
    check(all(r["viol"] == 0 for r in rows),
          "no violation of Lemma S (shell), Lemma O (odd separation), halving into the landing index, or Proposition U in any cell")
    exact, total = {}, {}
    for r in rows:
        D = r["D2"] / 2
        tight = math.ceil((r["k"] - D) / P.ALPHA - 1e-12)
        check_ub = P.ub_mult(r["k"], D, r["b"])
        if r["max"] > check_ub:
            raise AssertionError("bound exceeded %r" % r)
        total[r["b"]] = total.get(r["b"], 0) + 1
        exact[r["b"]] = exact.get(r["b"], 0) + (r["max"] == tight)
    for b in (1, 5, -1, -5):
        emit("b=%2d: cells with max multiplicity == ceil((k-D)/alpha): %d of %d" % (b, exact[b], total[b]))
    check(exact[1] == total[1] and exact[5] == total[5] and exact[-1] == total[-1],
          "b = 1, 5, -1: the worst case over all integers equals ceil((k-D)/alpha) in every cell (Proposition U is sharp)")
    emit("table b=1 (X = 2^(k+1)-1, all y <= X): max multiplicity / ceil((k-D)/alpha) / k")
    emit("   k |" + "".join("  D=%-4s" % (d / 2) for d in (2, 4, 6, 8, 10, 12)))
    for k in (15, 20, 25):
        cells = []
        for d in (2, 4, 6, 8, 10, 12):
            r = [x for x in rows if x["b"] == 1 and x["k"] == k and x["D2"] == d][0]
            cells.append("%3d/%-3d" % (r["max"], math.ceil((k - d / 2) / P.ALPHA - 1e-12)))
        emit("  %2d |" % k + " ".join(cells) + "   (k = %d)" % k)
    r = [x for x in rows if x["b"] == 1 and x["k"] == 25 and x["D2"] == 4][0]
    X = 2 ** 26 - 1
    emit("heavy counts H_M = #{y <= X : multiplicity counted from y >= M}, b=1, k=25, D=2:")
    emit("   M: " + " ".join("%6d" % M for M in range(1, len(r["H"]) + 1)))
    emit("   log2(H_M/X): " + " ".join("%6.2f" % math.log2(h / X) for h in r["H"]))
    dec = [math.log2(r["H"][i] / r["H"][i + 1]) for i in range(len(r["H"]) - 1)]
    emit("   bits lost per unit M: " + " ".join("%5.2f" % x for x in dec))
    check(all(x < 1.1 for x in dec[:8]),
          "heavy counts lose < 1.1 bits per unit M for M <= 9, well inside Proposition H's lower-bound rate 2 alpha = 3.17")
    return rows


# ----------------------------------------------------------------------------------------
def section_C():
    emit("== C. constructed hostile segments: Sturmian hover, then floor(D)+1 halvings ==")
    emit("   b     k      D   m(j)  ceil((k-D)/alpha)  m/k    word length   y (first dipper)")
    gaps = []
    for b in (1, -1):
        for k in (30, 40, 64, 100, 200, 400):
            for D in (P.depth(1.05 * math.log2(k)), P.depth(0.1 * k), P.depth(0.2 * k)):
                X = (1 << (k + 1)) - 1
                m, y, a, n = P.best_hostile(k, D, b, X, tries=24 if k >= 200 else 48)
                ys = P.segment(y, n + k + 1, b)
                land, _, _ = P.landing_map(ys, X, D, k)
                ds = land[n]
                bad = P.check_structure(ys, n, ds, D, b)
                ub = P.ub_mult(k, D, b)
                tight = math.ceil((k - float(D)) / P.ALPHA - 1e-12)
                gaps.append(tight - m)
                check(not bad and len(ds) == m and max(ys[i] for i in ds) <= X and len(set(ys)) == len(ys),
                      "b=%d k=%d D=%s: actual segment of distinct integers, landing index %d has %d dippers <= X, Lemmas S/O hold" % (b, k, D, n, m))
                emit("  %2d  %4d  %6.3f  %4d  %6d  %17.3f   %4d   %s" % (b, k, float(D), m, tight, m / k, n, str(y)[:24] + ("..." if len(str(y)) > 24 else "")))
                check(m <= ub, "within Proposition U")
    check(max(gaps) <= 2, "constructed integers reach ceil((k-D)/alpha) - 2 or better in all 36 cases (max gap %d)" % max(gaps))
    emit("gap histogram (ceil((k-D)/alpha) - m): %s" % sorted((g, gaps.count(g)) for g in set(gaps)))


def section_D():
    emit("== D. Proposition H: residue classes of heavy dippers with margin 1/8 ==")
    rng = random.Random(20260926)
    emit("   O     D     l = word length   safe visits   samples   min m   log2 class density")
    for O, D, k in ((10, P.depth(3), 40), (20, P.depth(4), 60), (40, P.depth(5.5), 120), (80, P.depth(7), 200),
                    (120, P.depth(9.25), 260)):
        X = (1 << (k + 1)) - 1
        word, times, safe = P.margin_word(O, D)
        rr = [r for O2 in (1,) for r in P.check_margin_class(O, D, 1, X, 40, rng)]
        l = len(word)
        check(len(rr) >= 30 and all(r[3] for r in rr),
              "O=%d D=%s: every sampled y in the class (48 l < y <= X/4) has all %d safe visits as dippers at index %d" % (O, D, len(safe), l))
        check(2 * len(safe) >= O + 2, "at least ceil(O/2)+1 safe visits (one of any two consecutive r has ||r alpha + D|| >= 1/8)")
        check(l <= 2 * P.ALPHA * ((len(safe) - 1)) + float(D) + 3 + 2 * P.ALPHA,
              "class length l = %d <= 2 alpha (M-1) + D + 3 + 2 alpha with M = #safe" % l)
        emit("  %3d  %6.3f   %5d            %4d          %4d     %4d    %8.2f" % (O, float(D), l, len(safe), len(rr), min(r[1] for r in rr), -l))
    emit("consequence (Proposition H): #H_M(X) >= X 2^-(2 alpha (M-1) + D + 5) - O(1); an orbit-blind split needs")
    for L in (100, 1000, 10 ** 4, 10 ** 5):
        D = 1.05 * math.log2(L)
        M0 = ((1 - P.H_STAR) * L - D - 5) / (2 * P.ALPHA) + 1
        emit("   L = %6d: heavy threshold M >= ((1-h*)L - D - 5)/(2 alpha) + 1 = %9.1f  (= %.4f L)" % (L, M0, M0 / L))
    check(((1 - P.H_STAR) * 10 ** 5 - 1.05 * math.log2(10 ** 5) - 5) / (2 * P.ALPHA) / 10 ** 5 > 0.015,
          "the forced threshold is linear in L (slope (1-h*)/(2 alpha) = %.4f)" % ((1 - P.H_STAR) / (2 * P.ALPHA)))


def depth_count(ys, i, j):
    """Integers D >= 1 with landing_D(i) = j (Lemma A); returns the list."""
    if i == j - 1:
        hi = math.log2(ys[i] / ys[j])
        return [D for D in range(1, 4) if D < hi - 1e-12 and P.below(ys[j], ys[i], D)]
    mu = min(ys[i + 1:j])
    lo, hi = math.log2(ys[i] / mu), math.log2(ys[i] / ys[j])
    cands = range(max(1, math.floor(lo) - 1), math.ceil(hi) + 2)
    return [D for D in cands if P.below(ys[j], ys[i], D) and not any(P.below(m, ys[i], D) for m in ys[i + 1:j])]


def section_E(orbits):
    emit("== E. Lemma A (depth averaging): each pair (i, j) serves at most one integer depth ==")
    pairs = worst = land_pts = 0
    for ys in orbits:
        k = 40
        for j in range(k, len(ys)):
            if ys[j - 1] != 2 * ys[j]:
                continue
            land_pts += 1
            tot = 0
            for i in range(j - k, j):
                c = len(depth_count(ys, i, j))
                pairs += 1
                worst = max(worst, c)
                tot += c
            if tot > k - 1:
                raise AssertionError("Lemma A sum exceeded at j=%d" % j)
    check(worst <= 1, "Lemma A on %d windows (k = 40) of the record orbits: %d pairs, each serves <= 1 integer depth, so sum_D m_D(j) <= k - 1" % (land_pts, pairs))


# ----------------------------------------------------------------------------------------
def count_below(ys, logs, X, D):
    """N(X 2^-D) = #{i : y_i <= X 2^-D} (float with exact fallback)."""
    thr = math.log2(X) - float(D)
    n = 0
    for v, lv in zip(ys, logs):
        if lv < thr - 1e-7:
            n += 1
        elif lv <= thr + 1e-7 and not P.below(X, v, D):      # v 2^D <= X  <=>  not X < v 2^D
            n += 1
    return n


def section_F():
    emit("== F. delay records (A006877) and path records (A006884) <= 10^12, T_1-orbits to 1 ==")
    cache = os.path.join(SCR, "oeis")
    p1, dl = P.load_bfile("b006877.txt", cache, FALLBACK)
    p2, pr = P.load_bfile("b006884.txt", cache, FALLBACK)
    dl = [x for x in dl if x <= 10 ** 12]
    pr = [x for x in pr if x <= 10 ** 12]
    emit("b006877.txt sha256 %s (%d delay records <= 10^12, largest %d)" % (sha(p1), len(dl), dl[-1]))
    emit("b006884.txt sha256 %s (%d path records <= 10^12, largest %d)" % (sha(p2), len(pr), pr[-1]))
    check(len(dl) == 90 and dl[-1] == 989345275647 and len(pr) == 61 and pr[-1] == 871673828443, "record lists as cached")
    starts = sorted(set(dl) | set(pr))
    orbits = [P.orbit_to_one(n) for n in starts]
    check(all(len(set(o)) == len(o) for o in orbits), "%d distinct record starts; every orbit segment up to 1 has distinct terms" % len(starts))
    emit("   L  depth       D   landing  dippers  mean-mult  max-mult (record)       max/ceil((k-D)/a)  #D/N(X2^-D)")
    worst_ratio = 0.0
    viol = 0
    for L in (20, 30, 40, 50, 60):
        X = 1 << L
        for name, D in (("theta_X", P.depth(1.05 * math.log2(L))), ("0.1L", P.depth(0.1 * L)), ("0.2L", P.depth(0.2 * L))):
            nl = nd_ = nb = 0
            mx, arg = 0, None
            for n, ys in zip(starts, orbits):
                lg = [math.log2(v) for v in ys]
                land, _, _ = P.landing_map(ys, X, D, L, lg)
                for j, ds in land.items():
                    viol += len(P.check_structure(ys, j, ds, D, 1))
                    if len(ds) > mx:
                        mx, arg = len(ds), n
                nl += len(land)
                nd_ += sum(len(v) for v in land.values())
                nb += count_below(ys, lg, X, D)
            ub = P.ub_mult(L, D, 1)
            worst_ratio = max(worst_ratio, mx / ub)
            emit("  %2d  %-7s %6.3f  %7d  %7d   %7.3f   %4d (%13d)   %8.3f          %7.3f" % (L, name, float(D), nl, nd_, nd_ / max(nl, 1), mx, arg or 0, mx / ub, nd_ / max(nb, 1)))
    check(viol == 0, "Lemmas S, O and the halving into the landing index hold at every landing point of every record orbit")
    check(worst_ratio <= 1.0, "every record-orbit multiplicity is within Proposition U (largest ratio to ceil((k-D)/alpha): %.3f)" % worst_ratio)
    return orbits


def section_G():
    emit("== G. the collatz-landing probe segments; climb-then-drop; a W-bit hover ==")
    emit("   L  start        D      max-mult  dippers(index range)  initial odd run  max dippers from the run per landing point")
    for L in (20, 30, 40, 60, 80):
        D = P.depth(1.05 * math.log2(L))
        for name, y0 in (("2^L-1", 2 ** L - 1), ("2^(L-1)+1", 2 ** (L - 1) + 1)):
            ys = P.orbit_to_one(y0)[:40 * L + 1]
            w = P.parity_word(y0, len(ys) - 1)
            run = 0
            while run < len(w) and w[run] == 1:
                run += 1
            land, _, _ = P.landing_map(ys, 1 << L, D, L)
            j, ds = max(land.items(), key=lambda kv: (len(kv[1]), -kv[0]))
            from_run = max(sum(1 for i in dd if i < run) for dd in land.values())
            emit("  %2d  %-10s %6.3f   %4d     %4d..%-4d               %4d              %d" % (L, name, float(D), len(ds), ds[0], ds[-1], run, from_run))
            check(from_run <= 2, "L=%d %s: the initial odd run contributes <= 2 dippers to any landing point" % (L, name))
    # a genuine climb-then-drop: 40 odd steps, then 35 halvings, scale 2^100
    L, D = 100, P.depth(1.05 * math.log2(100))
    word = [1] * 40 + [0] * 35
    y = P.realize(word, 1, 1 << 75, (1 << 76) - 1)
    ys = P.segment(y, len(word) + L + 1)
    check(P.parity_word(y, len(word)) == word, "climb-then-drop integer y = %d follows 1^40 0^35" % y)
    land, _, _ = P.landing_map(ys, 1 << L, D, L)
    per = [sum(1 for i in dd if i < 40) for dd in land.values()]
    emit("climb-then-drop (40 odd steps then 35 halvings, X = 2^100, D = %.3f): climb dippers %d spread over %d landing points, max %d per landing point" % (float(D), sum(per), sum(1 for p in per if p), max(per)))
    check(max(per) <= 2 and sum(per) >= 35, "climb points share landing points at most in pairs (Lemma S: one dyadic shell holds <= 2 points of an odd run)")
    # a W-bit hover (W = 4) of 60 steps, then D + W + 2 halvings
    rng = random.Random(4)
    W, m = 4.0, 60
    pos, hw = 0.0, []
    for _ in range(m):
        c = rng.randrange(2)
        nxt = pos + (P.ALPHA - 1 if c else -1)
        if not (-W <= nxt <= 0):
            c = 1 - c
            nxt = pos + (P.ALPHA - 1 if c else -1)
        hw.append(c)
        pos = nxt
    D = P.depth(1.05 * math.log2(100))
    word = hw + [0] * (int(W) + math.floor(float(D)) + 3)
    y = P.realize(word, 1, 1 << 90, (1 << 96) - 1)
    ys = P.segment(y, len(word) + L + 1)
    land, _, _ = P.landing_map(ys, 1 << L, D, L)
    hov = {j: [i for i in dd if i < m] for j, dd in land.items()}
    hov = {j: v for j, v in hov.items() if v}
    tot = sum(len(v) for v in hov.values())
    emit("W-bit hover (W = 4, 60 steps, then %d halvings, X = 2^100): %d hover dippers on %d landing points, multiplicities %s" % (len(word) - m, tot, len(hov), sorted(len(v) for v in hov.values())))
    check(len(hov) <= int(W) + 1 and tot / len(hov) >= tot / (W + 1),
          "the hover's dippers split among <= ceil(W)+1 landing points (one per unit shell); average >= (hover dippers)/(W+1)")


def cycle_word(kd, Dd, R):
    base, a, O = P.aligned_hostile_word(kd, Dd)
    word, mains = [], []
    for _ in range(R):
        word += base
        mains.append(len(word))                              # the designed landing index
        pos = P.walk_from_word(word)[-1]
        while pos[0] * P.ALPHA - pos[1] < -0.2:            # climb back to the cycle's start level
            word.append(1)
            pos = (pos[0] + 1, pos[1] + 1)
    return word, mains


def section_H():
    emit("== H. a 2-adic hostile word (repeated hover-crash-climb cycles) and its integer realisation ==")
    R = 20
    ratios = {}
    for kd in (40, 80):
        Dd = P.depth(4)
        word, mains = cycle_word(kd, Dd, R + 1)             # one extra cycle: no truncated window
        mains = mains[:R]
        walk = P.walk_from_word(word)
        top = max(o * P.ALPHA - t for o, t in walk)
        land, nd, e = P.landing_map_walk(walk, 0.0, top + 0.01, Dd, kd)
        ms = [len(v) for v in land.values()]
        nbelow = sum(1 for o, t in walk if o * P.ALPHA - t <= top + 0.01 - float(Dd))
        ub = P.ub_mult(kd, Dd, 1)
        main_m = [len(land.get(j, [])) for j in mains]
        ratios[kd] = sum(ms) / nbelow
        emit("pure walk, %d cycles, word length %d, design k = %d, D = %s: landing points %d, mean multiplicity %.2f (= %.3f k), designed landing points: min %d max %d, ceil((k-D)/alpha) = %d, #D/N(X2^-D) = %.2f" % (R, len(word), kd, Dd, len(ms), sum(ms) / len(ms), sum(ms) / len(ms) / kd, min(main_m), max(main_m), ub, ratios[kd]))
        check(max(ms) <= ub and min(main_m) >= ub - 1,
              "k=%d: every cycle's designed landing point has multiplicity >= ceil((k-D)/alpha) - 1 (2-adic model, all %d cycles)" % (kd, R))
        check(sum(ms) / len(ms) >= 0.25 * kd, "k=%d: mean multiplicity over all landing points >= 0.25 k" % kd)
        if kd == 40:
            N = len(word)
            r = P.terras_inverse(word, 1)
            y = r if r >= (1 << (N - 1)) else r + (1 << N)
            check(P.parity_word(y, N) == word, "the integer y = r + 2^N u (%d bits) follows the whole word (Terras)" % y.bit_length())
            ys = P.segment(y, 3 * N)
            Lact = max(v.bit_length() for v in ys[:N + 1])
            Dact = P.depth(1.05 * math.log2(Lact))
            landa, _, _ = P.landing_map(ys, 1 << Lact, Dact, Lact)
            in_word = [len(v) for j, v in landa.items() if j <= N]
            mx = max(in_word) if in_word else 0
            emit("integer realisation: scale 2^%d, k = %d, D = %s: landing points inside the word %d, max multiplicity %d, max/k = %.3f" % (Lact, Lact, Dact, len(in_word), mx, mx / Lact))
            check(mx <= 0.1 * Lact,
                  "at the integer's own scale the cycles are invisible (crash depth 5 < D): max multiplicity/k <= 0.1")
    check(ratios[80] >= 1.6 * ratios[40], "the averaged ratio #D/N(X2^-D) grows linearly with the design window (%.2f -> %.2f)" % (ratios[40], ratios[80]))


# ----------------------------------------------------------------------------------------
def main():
    os.makedirs(SCR, exist_ok=True)
    emit("procgen_landing_20260926_run.py -- landing multiplicity of the THM-4476 / THM-4499 recursion")
    emit("python %s" % sys.version.split()[0])
    section_A()
    section_B()
    section_C()
    section_D()
    orbits = section_F()
    section_E(orbits)
    section_G()
    section_H()
    emit("== reproduction ==")
    for name in ("procgen_landing_20260926_lib.py", "procgen_landing_20260926_run.py", "procgen_landing_20260926_scan.c"):
        emit("sha256 %s  %s" % (sha(os.path.join(HERE, name)), name))
    ru_self = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    ru_child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    scale = 1 if sys.platform == "darwin" else 1024              # bytes on macOS, KiB on Linux
    emit("peak RSS: runner %.1f MB, largest child %.1f MB" % (ru_self * scale / 2 ** 20, ru_child * scale / 2 ** 20))
    check(ru_self * scale < 700 * 2 ** 20 and ru_child * scale < 700 * 2 ** 20, "every process under 700 MB RSS")
    emit("wall time %.1f s" % (time.time() - T0))
    emit("ALL CHECKS PASSED")
    with open(OUT, "w") as f:
        f.write("\n".join(LINES) + "\n")


if __name__ == "__main__":
    main()
