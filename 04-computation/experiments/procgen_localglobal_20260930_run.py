#!/usr/bin/env python3
"""procgen_localglobal_20260930_run.py -- runner of the local/global lane (LRC(14) vs Collatz and the primes).

Re-checks every finite claim of 05-knowledge/results/procgen_localglobal_20260930_lrc_collatz_primes.md.
stdout only; every check raises on failure; the last line is ALL CHECKS PASSED.
  A1  LRC prime layer: AP rows, exhaustive tight sets (k<=7), deep wells, LRC(14) named rows,
      consecutive rows, random rows k = 8..14, S-smooth runner sets (kappa(S) = 1/P(S))
  A2  Collatz prime layer: gate census p <= 40 (local counts at every prime power of D, primitive
      words, gluing, classification), identities (zero-carry, repetition/cyclotomic, the a = p-2 mod-5
      family), convergent gates p <= 1054, the coverage regime |D| <= C+1 (p <= 5000), 5x+1 control
  B   duality: Poisson/relation identity, the single-relation lemma, mod-ell local badness = short
      relations (k = 3..6), Collatz objects miss the primes 3 / 5, the information budget
Run: nice python3 -u procgen_localglobal_20260930_run.py > ../../05-knowledge/results/procgen_localglobal_20260930.out
"""
import itertools
import math
import os
import random
import sys
import time
from collections import Counter
from fractions import Fraction

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from procgen_localglobal_20260930_lrc import (PRIMES, M_exact, covers, divisor_complete,  # noqa: E402
                                             least_lonely_den, least_lonely_prime, lonely_count,
                                             lonely_primes, prime_cover_depth, smooth_upto)
from procgen_localglobal_20260930_gates import (arch_class, brute, classify, convergent_gates,  # noqa: E402
                                               count_residues_dp, count_zero_dp, count_zero_mod,
                                               count_zero_prim,
                                               density_dp, gate, gate_row, lam_loc, lyndon_count,
                                               mult_order)
from procgen_localglobal_20260930_duality import (bad_classes, bad_classes_fast,  # noqa: E402
                                                 lonely_count_brute, lonely_count_fourier,
                                                 min_relation_l1, on_tight_line, short_relation_rate,
                                                 single_relation_blocks, tight_line_prediction)

T0 = time.time()
NCHECK = 0


def check(cond, msg):
    global NCHECK
    if not cond:
        raise AssertionError(msg)
    NCHECK += 1


def hdr(s):
    print()
    print("=" * 100)
    print(s, f"   [t={time.time() - T0:.0f}s]")
    print("=" * 100, flush=True)


def phi6(n):
    return n * n - n + 1


def least_uncovered(V):
    q = 2
    while covers(V, q):
        q += 1
    return q


def least_uncovered_prime(V):
    return next(q for q in PRIMES if not covers(V, q))


def joint_failures(V, lld):
    return [q for q in range(2, lld) if not covers(V, q)]


def spearman(x, y):
    x = np.argsort(np.argsort(x, kind="stable"), kind="stable")
    y = np.argsort(np.argsort(y, kind="stable"), kind="stable")
    return float(np.corrcoef(x, y)[0, 1])


# =====================================================================================================
def part_A1():
    hdr("A1.1  AP rows {1..k}: M = 1/(k+1), lonely only at j/(k+1); a prime-denominator lonely time exists iff k+1 is prime")
    for k in range(7, 15):
        V = list(range(1, k + 1))
        M, arg = M_exact(V)
        lld = least_lonely_den(V, 500)
        lp = lonely_primes(V, 500)
        want = [Fraction(j, k + 1) for j in range(1, (k + 1) // 2 + 1) if math.gcd(j, k + 1) == 1]
        print(f"  k={k:2d}  M={M}  maximisers={[str(t) for t in arg]}  lld={lld}  lonely primes<=500: {lp}")
        check(M == Fraction(1, k + 1) and arg == want and lld == k + 1, "AP")
        check(lp == ([k + 1] if (k + 1) in PRIMES else []), "AP primes")

    hdr("A1.2  exhaustive primitive sets (max speed <= B): tight sets, divisor-completeness, hardest covering sets")
    expect_tight = {3: [(1, 2, 3)], 4: [(1, 2, 3, 4), (1, 3, 4, 7)], 5: [(1, 2, 3, 4, 5), (1, 3, 4, 5, 9)],
                    6: [(1, 2, 3, 4, 5, 6)], 7: [(1, 2, 3, 4, 5, 6, 7), (1, 2, 3, 4, 5, 7, 12), (1, 4, 5, 6, 7, 11, 13)]}
    for k, B in [(3, 16), (4, 18), (5, 20), (6, 21), (7, 21)]:
        tight, best_cov, n, sig_dc, sig_nd = [], None, 0, [], []
        for S in itertools.combinations(range(1, B + 1), k):
            if math.gcd(*S) != 1:
                continue
            n += 1
            M = _M_fast(S)
            sig = (k + 1) * M - 1
            dc = divisor_complete(S)
            (sig_dc if dc else sig_nd).append(float(sig))
            if sig == 0:
                tight.append(S)
            if dc and (best_cov is None or M < best_cov[0]):
                best_cov = (M, S)
        print(f"  k={k} B={B}: {n} primitive sets; tight = {tight}; DC sets {len(sig_dc)} (median sigma "
              f"{np.median(sig_dc):.3f}) vs non-DC {len(sig_nd)} (median {np.median(sig_nd):.3f}); hardest primitive "
              f"covering M = {best_cov[0]} at {best_cov[1]}")
        check(tight == expect_tight[k], f"tight k={k}")
        check(all(not divisor_complete(S) for S in tight), "tight sets are not DC")
        check(best_cov[0] == Fraction(2, 2 * k + 1), "hardest covering = 2/(2k+1)")
        for S in tight:
            check(least_lonely_den(S, 100) == k + 1 == least_uncovered(S) and joint_failures(S, k + 1) == [],
                  "tight sets: lld = least uncovered = k+1, no joint failures")

    hdr("A1.3  deep wells {1..k-1, k(k+1)}, the LRC(14) named rows, joint (non-divisibility) failures")
    for k in range(4, 15):
        V = list(range(1, k)) + [k * (k + 1)]
        M, arg = M_exact(V)
        lld = least_lonely_den(V, 500)
        llp = least_lonely_prime(V, 5000)
        jf = joint_failures(V, lld)
        print(f"  k={k:2d}  M={M} (= (k+1)/Phi6(k+1): {M == Fraction(k + 1, phi6(k + 1))})  maximiser={arg}  lld={lld}  "
              f"least uncovered={least_uncovered(V)}  llp={llp}  joint failures={len(jf)}")
        check(M == Fraction(k + 1, phi6(k + 1)) and arg == [Fraction(k + 1, phi6(k + 1))], "deep well M")
        check(lld == 2 * k + 1 and lonely_count(V, 2 * k + 1) > 0, "deep well lonely at 2/(2k+1)")
        t = Fraction(2, 2 * k + 1)
        check(all((k + 1) * min((t * v) % 1, 1 - (t * v) % 1) >= 1 for v in V), "t = 2/(2k+1) lonely (exact)")
    V = list(range(1, 13)) + [182]
    lo = [Fraction(j, 10 ** 6) for j in range(70000, 78000)]
    good = [x for x in lo if all(14 * min((x * v) % 1, 1 - (x * v) % 1) >= 1 for v in V)]
    print(f"  deep well {{1..12,182}}: lonely grid points in [0.07,0.078] (step 1e-6): {len(good)}, span "
          f"[{float(min(good)):.6f}, {float(max(good)):.6f}]; t = 2/27 lonely: "
          f"{all(14 * min((Fraction(2, 27) * v) % 1, 1 - (Fraction(2, 27) * v) % 1) >= 1 for v in V)}")
    check(float(min(good)) < 0.0719 and float(max(good)) > 0.0765, "deep well lonely interval")
    rows = {"1..12,182": V, "1..12,5460": list(range(1, 13)) + [5460],
            "26*(1..12),339": [26 * i for i in range(1, 13)] + [339]}
    expM = {"1..12,182": Fraction(14, 183), "1..12,5460": Fraction(420, 5461), "26*(1..12),339": Fraction(1, 13)}
    for name, W in rows.items():
        M, arg = M_exact(W)
        lld = least_lonely_den(W, 3000)
        llp = least_lonely_prime(W, 20000)
        jf = joint_failures(W, lld)
        print(f"  {name:16s} M={M}  maximisers ({len(arg)}): {[str(t) for t in arg[:3]]}  DC={divisor_complete(W)}  lld={lld}  "
              f"least uncovered={least_uncovered(W)}  llp={llp}  joint failures {jf}")
        check(M == expM[name] and lld == 27 and llp == 41 and divisor_complete(W), f"named row {name}")
    hard = [(1, 3, 4, 5), (1, 3, 4, 5, 18), (1, 2, 5, 6, 7, 8), (1, 4, 5, 6, 7, 11, 16)]
    for S in hard:
        lld = least_lonely_den(S, 200)
        jf = joint_failures(S, lld)
        print(f"  hardest covering {S}: lld={lld} = 2k+1, least uncovered={least_uncovered(S)}, joint failures {jf}")
        check(lld == 2 * len(S) + 1 and len(jf) >= 3, "hardest covering rows: joint failures")

    hdr("A1.4  consecutive rows {w, ..., w+12} (k = 13), w <= 400")
    res = []
    for w in range(1, 401):
        V = list(range(w, w + 13))
        M, _ = M_exact(V)
        res.append((14 * M - 1, w, divisor_complete(V)))
    res.sort()
    dcmin = min(r for r in res if r[2])
    print(f"  hardest: {[(str(s), w, dc) for s, w, dc in res[:5]]};  DC rows: {sum(r[2] for r in res)}/400, "
          f"min sigma over DC rows = {dcmin[0]} at w = {dcmin[1]}")
    check(res[0][1] == 1 and res[0][0] == 0 and not res[0][2] and dcmin[0] == Fraction(3, 4) and dcmin[1] == 2,
          "consecutive rows")

    hdr("A1.5  random primitive rows, k = 8..14 (200 rows per cell, seed 20260930)")
    rng = random.Random(20260930)
    print("  k  Vmax  DC%   P(lld=u)  P(llp=u_p)  joint/row  rho(depth,sigma)  rho(depth,lld)  med sigma DC/nonDC")
    cells = []
    for k in range(8, 15):
        for Vmax in (2 * k + 2, 60, 1000):
            recs = []
            while len(recs) < 200:
                V = sorted(rng.sample(range(1, Vmax + 1), k))
                if math.gcd(*V) != 1:
                    continue
                M, _ = M_exact(V)
                lld = least_lonely_den(V, 3000)
                llp = least_lonely_prime(V, 20000)
                recs.append((float((k + 1) * M - 1), divisor_complete(V), prime_cover_depth(V), lld,
                             lld == least_uncovered(V), llp == least_uncovered_prime(V), len(joint_failures(V, lld))))
            a = np.array(recs, dtype=float)
            dc = a[:, 1] == 1
            r1, r2 = spearman(a[:, 2], a[:, 0]), spearman(a[:, 2], a[:, 3])
            m_dc = np.median(a[dc, 0]) if dc.any() else float("nan")
            m_nd = np.median(a[~dc, 0]) if (~dc).any() else float("nan")
            print(f"  {k:2d} {Vmax:5d} {100 * dc.mean():4.0f}  {a[:, 4].mean():8.3f}  {a[:, 5].mean():9.3f}  "
                  f"{a[:, 6].mean():8.3f}  {r1:+15.2f}  {r2:+14.2f}   {m_dc:.3f}/{m_nd:.3f}")
            cells.append((k, Vmax, a[:, 4].mean(), a[:, 5].mean(), a[:, 6].mean(), r1, r2))
    for k, Vmax, p_u, p_up, jr, r1, r2 in cells:
        check(p_u >= (0.9 if Vmax >= 60 else 0.8) and p_up >= 0.65 and jr <= 0.4, "random rows: lld law")
        check(abs(r1) <= 0.4 and r2 > 0.2, "random rows: depth vs sigma weak, depth vs lld positive")

    hdr("A1.6  S-smooth runner sets: kappa(S) = 1/P(S), P(S) = least prime outside S (S = all primes < P)")
    for S, P in [((2,), 3), ((2, 3), 5), ((2, 3, 5), 7), ((2, 3, 5, 7), 11), ((2, 3, 5, 7, 11), 13)]:
        M, arg = M_exact(list(range(1, P)))
        sm = smooth_upto(10 ** 6, S)
        ok = all(v % P != 0 for v in sm) and set(range(1, P)) <= set(sm)
        print(f"  S={S}: M({{1..{P - 1}}}) = {M}, maximisers {[str(t) for t in arg]}; {len(sm)} S-smooth numbers <= 1e6, "
              f"none divisible by {P}: kappa(S) = 1/{P}")
        check(M == Fraction(1, P) and ok, "kappa")
    M, arg = M_exact([1, 2, 3, 4])
    check(arg == [Fraction(1, 5), Fraction(2, 5)], "kappa(2,3) extremisers")


def _M_fast(V):
    V = np.array(V, dtype=np.int64)
    k = len(V)
    D = np.unique((V[:, None] + V[None, :])[np.triu_indices(k)])
    c = np.concatenate([np.arange(1, d // 2 + 1) for d in D])
    d = np.concatenate([np.full(d // 2, d) for d in D])
    r = (c[:, None] * V[None, :]) % d[:, None]
    num = np.minimum(r, d[:, None] - r).min(axis=1)
    i = int(np.argmax(num / d))
    return Fraction(int(num[i]), int(d[i]))


# =====================================================================================================
def census(q, pmax):
    rows = []
    for p in range(2, pmax + 1):
        for a in range(1, p):
            r = gate_row(p, a, q)
            r["cls"] = classify(r)
            rows.append(r)
    return rows


def part_A2():
    hdr("A2.1  gate census, q = 3, all (p,a) with p <= 40 (exact meet-in-the-middle counts, primitive words)")
    rows = census(3, 40)
    cnt = Counter(r["cls"] for r in rows)
    print("  classes:", dict(sorted(cnt.items())))
    glob = [(r["p"], r["a"], r["D"], r["N_Dp"]) for r in rows if r["cls"] == "GLOBAL"]
    print("  primitive cycles found (p, a, D, #words):", glob)
    check([g[:2] for g in glob] == [(2, 1), (3, 2), (11, 7)], "census: exactly the known primitive cycles")
    check(cnt["ARCH"] == sum(1 for r in rows if r["D"] > 0 and 2 ** r["p"] > 4 ** r["a"]), "ARCH = dyadic a < p/2")
    check(all(r["arch"] != "none" for r in rows if r["D"] < 0), "3-adic side has no perigee obstruction (q = 3)")
    rep = [r for r in rows if r["cls"] == "REPEAT"]
    check(all(r["N_D"] in (2, 3, 11) for r in rep) and len(rep) == 33, "repeats of {1,2}, {-5,-7,-10}, {-17,...}")
    side = Counter((r["cls"], "dyadic" if r["D"] > 0 else "3-adic") for r in rows)
    print("  by side:", dict(sorted(side.items())))
    H = [r for r in rows if r["cls"].startswith("HASSE")]
    mx = max(H, key=lam_loc)
    print(f"  HASSE gates: {len(H)}; largest local-product prediction {lam_loc(mx):.3f} primitive necklaces at "
          f"(p,a)=({mx['p']},{mx['a']}), D={mx['D']}={mx['fac']}, local={mx['localp']}")
    check(lam_loc(mx) < 1, "no HASSE gate predicts a cycle")
    # local densities, zero classes by prime type
    stat = {}
    for r in rows:
        p, a, L = r["p"], r["a"], r["L"]
        if r["arch"] == "none" or abs(r["D"]) == 1 or L == 0 or r["cls"] == "GLOBAL":
            continue
        g = math.gcd(p, a)
        for qq, d in r["primes"].items():
            if d["e"] != 1:
                continue
            kind = "primitive"
            if any(g % dd == 0 and gate(p // dd, a // dd) % qq == 0 for dd in range(2, g + 1)):
                kind = "sub-gate"
            elif any(g % dd == 0 and (r["D"] // gate(p // dd, a // dd)) % qq == 0 for dd in range(2, g + 1)):
                kind = "cofactor"
            s = stat.setdefault(kind, dict(n=0, exp=0.0, obs=0, ratio=[], z=[]))
            lam = L / qq
            s["n"] += 1
            s["exp"] += math.exp(-lam)
            s["obs"] += d["Nqp"] == 0
            if lam >= 5:
                s["ratio"].append(d["Nqp"] / (p * L / qq))
    print("  prime type   pairs  zero-class empty (obs / necklace-Poisson)   N/(Cp/q) for L/q>=5: n, mean, sd")
    for kind, s in sorted(stat.items()):
        rr = np.array(s["ratio"])
        print(f"  {kind:10s} {s['n']:6d}   {s['obs']:5d} / {s['exp']:7.1f}                      "
              f"{len(rr):4d}  {rr.mean():.4f}  {rr.std():.4f}")
        check(abs(s["obs"] - s["exp"]) <= 3 * math.sqrt(s["exp"]) + 3 and abs(rr.mean() - 1) < 0.05,
              f"local zeros/densities consistent with 1/q ({kind})")
    # gluing on the bulk
    rat = [Nm / pred for r in rows for (_, m, Nm, pred) in r["glue"]
           if pred >= 20 and r["arch"] != "none" and 0.5 <= r["a"] / r["p"] <= 0.8]
    rat_all = [Nm / pred for r in rows for (_, m, Nm, pred) in r["glue"] if pred >= 20]
    print(f"  CRT gluing N_m / (Cp prod nu) (prediction >= 20): bulk 0.5<=a/p<=0.8 eligible n={len(rat)} mean "
          f"{np.mean(rat):.3f} sd {np.std(rat):.3f};  all gates n={len(rat_all)} mean {np.mean(rat_all):.3f} sd {np.std(rat_all):.3f}")
    check(abs(np.mean(rat) - 1) < 0.15, "CRT gluing independent in the bulk")
    # all-words vs primitive: the repetition effect
    zs_all = []
    for r in rows:
        for qq, d in r["primes"].items():
            mu = r["C"] / qq
            if d["e"] == 1 and mu >= 200:
                zs_all.append((d["Nq"] - mu) / math.sqrt(r["p"] * mu))
    print(f"  all words (repetitions included): max z-score of N_q against C/q = {max(zs_all):.1f}  (cyclotomic "
          f"cofactor effect; e.g. (40,22): 2^40-3^22 = (2^20-3^11)(2^20+3^11) and every square word u u is 0 mod 2^20+3^11)")
    check(max(zs_all) > 20, "repetition effect visible on all words")
    surprising = [(r["p"], r["a"], r["fac"], r["localp"]) for r in rows if r["cls"].endswith("!")]
    print("  '!' gates (a local zero with necklace-Poisson expectation >= 3):")
    for s in surprising:
        print("     ", s)
    for side, sg in (("dyadic", 1), ("3-adic", -1)):
        el = [r for r in rows if r["D"] * sg > 0 and abs(r["D"]) > 1 and r["arch"] != "none" and r["cls"] != "GLOBAL"]
        su = sum(r["L"] / abs(r["D"]) for r in el)
        sl = sum(lam_loc(r) for r in el)
        print(f"  {side}: {len(el)} eligible gates without a primitive cycle: uniform prediction sum L/|D| = {su:.3f} "
              f"primitive necklaces; local-product prediction sum = {sl:.3f}; observed 0")
        check(sl < su and sl < 1, "local data lower the prediction below 1")
    for (p, a, qq) in [(36, 32, 263), (40, 35, 1931)]:
        h = np.array(count_residues_dp(p, a, qq), dtype=np.int64)
        empty = [int(x) for x in np.nonzero(h == 0)[0]]
        print(f"  sub-gate anomaly ({p},{a}) mod {qq}: C = {int(h.sum())}, mean per class {h.sum() / qq:.1f}, "
              f"empty classes {empty} (gate D({p // math.gcd(p, a)},{a // math.gcd(p, a)}) = {gate(p // math.gcd(p, a), a // math.gcd(p, a))})")
        check(empty == [0], "the only empty residue class is 0")
    cov = [(r["p"], r["a"], r["D"], r["cls"]) for r in rows if 1 < abs(r["D"]) <= r["C"] + 1]
    print("  coverage regime |D| <= C+1 at p <= 40:", cov)
    check([c[:2] for c in cov] == [(4, 2), (5, 3), (8, 5), (11, 7), (16, 10), (19, 12), (27, 17), (38, 24)], "coverage p<=40")
    ordtab = []
    for r in rows:
        if r["cls"] in ("HASSE", "PRIME-GATE") and 0.58 <= r["a"] / r["p"] <= 0.7 and r["p"] >= 24:
            ordtab.append((r["p"], r["a"], r["D"], {qq: (d["o2"], d["o3"], d["Nqp"]) for qq, d in r["primes"].items()}))
    print("  near-critical HASSE / PRIME-GATE gates, p >= 24: (p, a, D, {q: (ord_q 2, ord_q 3, primitive N_q)})")
    for t in ordtab[:14]:
        print("     ", t)

    hdr("A2.2  identities: zero-carry, repetition (cyclotomic cofactor), the a = p-2 family mod 5")
    n = 0
    for p in range(2, 13):
        for a in range(1, p):
            D = gate(p, a)
            for S in itertools.combinations(range(p), a):
                Z = [s for s in range(p) if s not in S]
                c = sum(3 ** (a - 1 - i) * 2 ** s for i, s in enumerate(S))
                check(c == -D + sum(3 ** (a + j - z) * 2 ** z for j, z in enumerate(Z)), "zero-carry")
                check(c % 3 == pow(2, S[-1], 3), "c_w = 2^{s_{a-1}} mod 3")
                n += 1
    print(f"  c_w = -(2^p - 3^a) + sum_j 3^(a+j-1-z_j) 2^(z_j)  (zeros z_1<...<z_z) and c_w = 2^(s_(a-1)) (mod 3): {n} words")
    n = 0
    for p in range(2, 13):
        for a in range(1, p):
            g = math.gcd(p, a)
            for d in range(2, g + 1):
                if g % d:
                    continue
                for S in itertools.combinations(range(p // d), a // d):
                    u = [1 if i in S else 0 for i in range(p // d)]
                    w = u * d
                    cu = sum(3 ** (a // d - 1 - i) * 2 ** s for i, s in enumerate(S))
                    Sw = [i for i, b in enumerate(w) if b]
                    cw = sum(3 ** (a - 1 - i) * 2 ** s for i, s in enumerate(Sw))
                    check(cw * gate(p // d, a // d) == cu * gate(p, a), "repetition identity")
                    n += 1
    print(f"  c_(u^d) = c_u * (2^p - 3^a)/(2^(p/d) - 3^(a/d)) for every d-fold repetition: {n} cases")
    fam = []
    for p in range(5, 62, 2):
        a = p - 2
        if gate(p, a) % 5 == 0:
            fam.append((p, count_zero_dp(p, a, 5)))
    print(f"  a = p-2, 5 | D (exactly the odd p): N_5 = {set(x[1] for x in fam)} for p in {fam[0][0]}..{fam[-1][0]} "
          f"({len(fam)} gates); 2/3 = -1 mod 5, so c_w = 3^a((-1)^z1 + 3(-1)^z2) != 0 mod 5")
    check(all(x[1] == 0 for x in fam) and len(fam) == 29, "a=p-2 family")

    hdr("A2.3  convergent / semiconvergent gates of log_2 3 up to p = 1054")
    import sympy
    G = [(5, 3), (11, 7), (27, 17), (46, 29), (65, 41), (84, 53), (149, 94), (233, 147), (317, 200), (401, 253),
         (485, 306), (569, 359), (1054, 665)]
    check(set(convergent_gates(1100)) <= set(G) | {(2, 1), (3, 2), (8, 5), (19, 12)}, "convergent list")
    for p, a in G:
        D = gate(p, a)
        M = abs(D)
        amp = max(2 ** p, 3 ** a) / M
        l2 = (math.lgamma(p + 1) - math.lgamma(a + 1) - math.lgamma(p - a + 1)) / math.log(2) - math.log2(M)
        f = sympy.factorint(M, limit=10 ** 6)
        comp = [x for x in f if not sympy.isprime(x)]
        small = {x: e for x, e in f.items() if x < 10 ** 6}
        dens = {}
        for x in sorted(small):
            if x <= 20000:
                dens[x] = x * density_dp(p, a, x) - 1
        big = [len(str(x)) for x in f if x >= 10 ** 6]
        print(f"  ({p},{a}) {'dyadic' if D > 0 else '3-adic'} |D|~10^{len(str(M)) - 1} amp={amp:.3g} "
              f"log2 C - log2|D| = {l2:+.2f}  small primes {small}  other factors (digits) {big} "
              f"{'(composite cofactor left)' if comp else '(fully factored)'}  "
              f"ords {[(x, mult_order(2, x), mult_order(3, x)) for x in sorted(small)]}  q*nu_q - 1: "
              f"{ {x: float(f'{v:.2e}') for x, v in dens.items()} }", flush=True)
        if (p, a) not in ((5, 3), (11, 7), (27, 17)):
            check(all(abs(v) < 1e-5 for v in dens.values()), "no local bias at small gate primes (large p)")

    hdr("A2.4  the coverage regime |2^p - 3^a| <= C(p,a)+1 (where the LRC trivial lemma speaks), p <= 5000")
    out = []
    P3 = [1]
    for _ in range(3400):
        P3.append(P3[-1] * 3)
    x = math.log(2) / math.log(3)
    for p in range(2, 5001):
        for a in sorted({math.floor(p * x), math.ceil(p * x)}):
            if 1 <= a <= p:
                M = abs((1 << p) - P3[a])
                if 1 < M <= math.comb(p, a) + 1:
                    out.append((p, a))
    print("  gates:", out)
    check(out == [(4, 2), (5, 3), (8, 5), (11, 7), (16, 10), (19, 12), (27, 17), (38, 24), (46, 29), (84, 53)],
          "coverage regime list")

    hdr("A2.5  control: 5x+1 (q = 5), all (p,a) with p <= 22")
    rows5 = census(5, 22)
    cnt5 = Counter(r["cls"] for r in rows5)
    glob5 = [(r["p"], r["a"], r["D"], r["N_Dp"]) for r in rows5 if r["cls"] == "GLOBAL"]
    print("  classes:", dict(sorted(cnt5.items())))
    print("  primitive cycles (p, a, D, #words):", glob5)
    check((5, 2, 7, 5) in glob5 and (7, 3, 3, 14) in glob5, "5x+1 sporadic cycles found")


# =====================================================================================================
def part_B():
    hdr("B1  Poisson / relation-lattice identity  N_lone(V,q) = q * sum_{xi in Lambda_q(V)} prod ghat(xi_i)")
    for V, q in [((1, 2, 3), 7), ((1, 3, 4), 11), ((2, 5, 7), 13), ((1, 2, 3, 4), 7), ((1, 3, 4, 7), 5), ((1, 2, 5), 9)]:
        a, b = lonely_count_brute(V, q), lonely_count_fourier(V, q)
        print(f"  V={V} q={q}: brute {a}, Fourier over Lambda_q {b.real:+.10f}{b.imag:+.1e}i")
        check(abs(a - b) < 1e-8, "Poisson identity")

    hdr("B2  single-relation lemma: one relation xi of Lambda_q(V) excludes every lonely m/q by itself iff ||xi||_1 = 1")
    for k in range(3, 9):
        n = 0
        for xi in itertools.product(range(-3, 4), repeat=min(k, 4)):
            xi = tuple(xi) + (0,) * (k - len(xi))
            if not any(xi):
                continue
            l1 = sum(abs(t) for t in xi)
            check(single_relation_blocks(xi, k) == (l1 == 1), "single relation lemma")
            n += 1
        print(f"  k={k}: checked {n} relation vectors (entries in [-3,3] on 4 coordinates)")

    hdr("B3  mod-ell local badness: every bad residue vector (no lonely m/ell) carries a relation of l1-norm <= 4")
    for k, Lmax in [(3, 61), (4, 41), (5, 29), (6, 23)]:
        for ell in [p for p in PRIMES if k + 2 <= p <= Lmax]:
            bad, beta = bad_classes(k, ell)
            norms = [min_relation_l1(u, ell, 2 * k + 7) for u in bad]
            mx = max(norms) if norms else 0
            hist = dict(sorted(Counter(norms).items()))
            print(f"  k={k} ell={ell:2d}: bad classes {len(bad):4d}  beta={beta:.5f}  ell*beta={ell * beta:6.3f}  "
                  f"min-relation norms {hist}")
            check(None not in norms and mx <= 4, "short relations")

    print("  base rate: fraction of ALL vectors carrying a relation of norm <= 3 / <= 4 (the test is informative only "
          "where this is < 1):")
    for k, ell in [(3, 37), (3, 47), (3, 61), (4, 41), (5, 29), (6, 23)]:
        r3, r4 = short_relation_rate(k, ell, 3), short_relation_rate(k, ell, 4)
        print(f"     k={k} ell={ell}: {r3:.4f} / {r4:.4f}")
        if k == 3:
            check(r3 < 0.6, "k=3 base rate informative")

    hdr("B5  local stabilization (tight-line principle): for large primes ell the bad vectors are exactly the "
        "multiples of the reductions of the tight sets")
    for k, ells, exc in [(3, [31] + [p for p in PRIMES if 47 <= p <= 199], [29, 37, 41, 43]),
                         (4, [59, 79, 83, 97, 107, 109, 113, 127] + [p for p in PRIMES if 137 <= p <= 199],
                          [61, 67, 71, 73, 89, 101, 103, 131]),
                         (5, [101, 151], [53, 59, 61, 67, 71, 89, 127])]:
        for ell in exc + ells:
            bad, beta = bad_classes_fast(k, ell)
            off = [u for u in bad if not on_tight_line(u, ell, k)]
            ratio = beta / tight_line_prediction(k, ell)
            tag = "exception" if ell in exc else "stabilised"
            print(f"  k={k} ell={ell:3d} [{tag:10s}]: bad classes {len(bad):4d}, density / tight-line prediction = "
                  f"{ratio:.4f}, off the tight lines {len(off):4d} {off[:2]}")
            if ell in exc:
                check(len(off) > 0, "listed exception")
            else:
                check(len(off) == 0 and abs(ratio - 1) < 1e-12, "stabilised: bad set = tight lines")

    hdr("B4  Collatz objects sit in LRC's trivial (non-covering) class; the information budget")
    for p in range(2, 15):
        for a in range(1, p):
            check(all(c % 3 for c in brute(p, a)), "no carry divisible by 3")
    cyc = [[1, 2], [5, 7, 10], [17, 25, 37, 55, 82, 41, 61, 91, 136, 68, 34]]
    check(all(x % 3 for c in cyc for x in c), "known cycles avoid 3")
    check(all(v % 5 for v in smooth_upto(10 ** 6, (2, 3))), "{2,3}-runners avoid 5")
    print("  every carry c_w is prime to 3 (c_w = 2^(s_(a-1)) mod 3) => the carry set is lonely at t = 1/3 (level 1/3);")
    print("  cycle elements are prime to 3 (known cycles checked; general: a multiple of 3 has no odd-step preimage);")
    print("  the runner set {2^s 3^r} is prime to 5 => lonely at t = 1/5 (kappa(2,3) = 1/5 exactly, A1.6).")
    th = math.log(2) / math.log(3)
    H = -(th * math.log2(th) + (1 - th) * math.log2(1 - th))
    print(f"  information budget at the critical ratio a/p = log_3 2 = {th:.6f}: log2 C(p,a)/p -> H = {H:.5f} bits, "
          f"log2|D|/p -> 1 (Baker/Rhin): deficit {1 - H:.5f} bits per step")
    check(abs(H - 0.94996) < 1e-4, "entropy")
    for p, a in [(84, 53), (485, 306), (1054, 665)]:
        lc = (math.lgamma(p + 1) - math.lgamma(a + 1) - math.lgamma(p - a + 1)) / math.log(2)
        print(f"  ({p},{a}): log2 C = {lc:.1f}, log2|D| = {math.log2(abs(gate(p, a))):.1f}")


if __name__ == "__main__":
    part_A1()
    part_A2()
    part_B()
    print()
    print(f"{NCHECK} checks, total time {time.time() - T0:.0f}s")
    print("ALL CHECKS PASSED")
