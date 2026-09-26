#!/usr/bin/env python3
"""
procgen_seven3_20260926_run.py -- runner of lane "seven3" (session collatz-procgen-20260922, 2026-09-26):
itinerary-coded sign strategies of 7n +- 1.  Re-checks every finite claim of
05-knowledge/results/procgen_seven3_20260926_itinerary_strategies.md and prints [ok] lines, ending with
ALL CHECKS PASSED.  Output: 05-knowledge/results/procgen_seven3_20260926.out (written by redirecting stdout).

Engines (compiled into scratch/procgen_seven3/bin by the lib; nothing committed there):
  procgen_seven3_20260926_markov.c  Markov-refinement rho_max (Howard + exact potential on every edge + critical set)
  procgen_seven3_20260926_game.c    least fixed points of Min's / Max's operators at level k
  procgen_seven3_20260926_adv.c     exact value of one Max lift strategy (certificate + maximality)
  procgen_seven2_20260926_rhomax.c  (read-only reuse) the seven2 uniform-level rho_max engine, as an independent check
  procgen_seven2_20260926_restrict.c (read-only reuse) not used by the runner
Every value is certified: a witness cycle re-walked as a rational periodic point with exact arithmetic, and a
potential checked on every edge; lower bounds for rho* by Lemma G2 certificates checked with exact integers.
"""
import os
import sys
import json
import time
import random
import hashlib
from fractions import Fraction as Fr
from math import log

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import procgen_seven3_20260926_lib as L
from procgen_seven3_20260926_lib import (Q, mh, v2, word_class, check_word_class, flip_rule, mh_rule, check_partition,
                                          rho_markov, rho_markov_c, rhomax_uniform, check_uniform_cycle, game_lfp,
                                          check_min_cert, check_max_cert, cycle_rational, pattern_rule,
                                          adversary_value, check_adv_cycle, compress_strategy, itin_of_residue)
import procgen_seven3_20260926_itree as IT
import procgen_seven3_20260926_search as S

T0 = time.time()
NCHECK = [0]


def ok(msg):
    NCHECK[0] += 1
    print('[ok] %s' % msg, flush=True)


def fail(msg):
    print('[FAIL] %s' % msg, flush=True)
    sys.exit(1)


def need(cond, msg):
    if not cond:
        fail(msg)
    ok(msg)


DATA = os.path.join(HERE, 'procgen_seven3_20260926_rules.json')
RHO_STAR = {8: Fr(3, 7), 9: Fr(3, 7), 10: Fr(2, 5), 11: Fr(2, 5), 12: Fr(2, 5), 13: Fr(2, 5), 14: Fr(15, 38),
            15: Fr(15, 38), 16: Fr(7, 18), 17: Fr(7, 18), 18: Fr(19, 49), 19: Fr(13, 34), 20: Fr(34, 89)}
R18 = [(11, 5), (21, 5)] + [(c, 8) for c in (35, 67, 91, 165, 189, 221)] + \
      [(c, 10) for c in (93, 157, 349, 381, 413, 611, 643, 675, 867, 931)]


def section(t):
    print('\n== %s   (t = %.0f s)' % (t, time.time() - T0), flush=True)


# ------------------------------------------------------------------------------------------------ D: arena semantics
def sec_D():
    section('D. rho*(7,k) with this lane\'s own game engine (Stern-Brocot from scratch for k <= 16, both certificates)')

    def solve(k, capmul=64):
        lo, hi = (0, 1), (1, 1)
        while True:
            m = Fr(lo[0] + hi[0], lo[1] + hi[1])
            if m.denominator > 10 ** 4:
                return None
            mn = game_lfp('min', k, m, cap=capmul * m.denominator)
            if mn is None:
                lo = (m.numerator, m.denominator)
                continue
            mx = game_lfp('max', k, m, cap=capmul * m.denominator)
            if mx is not None:
                if not (check_min_cert(k, m, *mn) and check_max_cert(k, m, *mx)):
                    fail('certificate check at k=%d' % k)
                return m
            hi = (m.numerator, m.denominator)
    for k in range(8, 17):
        v = solve(k)
        need(v == RHO_STAR[k], 'D1 k=%d: Stern-Brocot + both certificates give rho*(7,%d) = %s' % (k, k, v))
    for k in range(17, 21):
        F = RHO_STAR[k]
        mn = game_lfp('min', k, F, cap=64 * F.denominator)
        mx = game_lfp('max', k, F, cap=64 * F.denominator)
        need(mn is not None and mx is not None and check_min_cert(k, F, *mn) and check_max_cert(k, F, *mx),
             'D2 k=%d: Min and Max least fixed points finite at %s, both certificates exact: rho*(7,%d) = %s' % (k, F, k, F))


# ------------------------------------------------------------------------------------------------ M: Markov engine
def sec_M():
    section('M. Lemma M (Markov refinement): exact rho_max of variable-depth rules, two implementations + uniform engine')
    for word in [[(1, 2)], [(1, 2), (1, 2)], [(-1, 2), (-1, 2)], [(1, 2), (-1, 3)], [(1, ('ge', 5))],
                 [(1, 2), (1, 2), (-1, ('ge', 5))], [(-1, 4), (1, 2), (1, 3), (-1, 2)]]:
        check_word_class(word)
    ok('M0 itinerary cylinders are residue classes of depth 1 + sum v (Lemma MH(ii)), 7 words re-parsed')
    # every bit class is an itinerary node: re-parse all odd classes mod 2^d, d <= 12
    cnt = 0
    for d in range(1, 13):
        for c in range(1, 1 << d, 2):
            sy, tail = itin_of_residue(c, d)
            w = list(sy) + ([(tail[0], ('ge', tail[1]))] if tail and tail[1] >= 2 else [])
            if tail is None:
                c2, d2 = (word_class(sy) if sy else (1, 1))
            else:
                c2, d2 = word_class(w)
            assert (c2, d2) == (c, d), (c, d, sy, tail)
            cnt += 1
    ok('M0b every odd class c mod 2^d (d <= 12, %d classes) is an itinerary node w or w(s,>=m) (Fact I)' % cnt)
    v1 = rho_markov_c(mh_rule(2))
    v2_ = rho_markov(mh_rule(2))
    need(v1[0] == v2_[0] == Fr(1, 2), 'M1 MH: rho_max = 1/2 (C and Python Markov engines), witness %s' % v1[2])
    r18 = flip_rule(R18)
    check_partition(r18)
    a = rho_markov_c(r18)
    b = rho_markov(r18)
    need(a[0] == b[0] == Fr(2, 5) and a[1] == b[1] == 170,
         'M2 18-class rule: rho_max = 2/5 by both Markov engines, 170 leaves (uniform level 10: 1024 nodes), witness %s' % a[2])
    for k in (10, 12, 14):
        v, cyc = rhomax_uniform(r18, k)
        x0, aa, nn = check_uniform_cycle(r18, k, cyc)
        need(v == Fr(2, 5) and Fr(aa, nn) == v, 'M3 18-class rule at uniform level %d (seven2 engine): 2/5, witness %s re-walked' % (k, x0))


# ------------------------------------------------------------------------------------------------ F: families
def sec_F():
    section('F. Family G + A_D (flip (s,2)(s,2) and alternating v=2 runs of length >= D)')
    for D in range(2, 16):
        alt = ' '.join(['2'] + ['!2'] * (D - 1))
        rule = pattern_rule(['2 =2', alt])
        a = rho_markov_c(rule)
        b = rho_markov(rule)
        # PROVED lower bound: the MH-periodic orbit whose itinerary alternates signs throughout and has valuations
        # 2^(D-1) 3 repeated (period D symbols for even D, 2D for odd D, so that the alternation closes)
        per = 2 * D if D % 2 else D
        word = []
        s = 1
        for i in range(per):
            v = 3 if (i % D) == D - 1 else 2
            word.append((s, v))
            s = -s
        # periodic point: x_n = num x_0 + const for the composed MH steps; fixed point const / (1 - num)
        num, const = Fr(1), Fr(0)
        for (sg, v) in word:
            num, const = num * Q / (1 << v), (Q * const + sg) / (1 << v)
        x0 = const / (1 - num)
        # re-walk under the rule: all steps must be MH steps (no flip), with the word's symbols
        x = x0
        dmax = max(d for (_, d) in rule)
        odd = tot = 0
        for (sg, v) in word:
            r = x.numerator * pow(x.denominator, -1, 1 << (dmax + 2)) % (1 << (dmax + 2))
            assert L.rule_sign(rule, r, dmax) == sg == mh(r)
            y = Q * x + sg
            vv = 0
            while y.numerator % 2 == 0:
                y /= 2
                vv += 1
            assert vv == v
            x = y
            odd += 1
            tot += v
        assert x == x0
        expect = Fr(1, 2) if D == 2 else Fr(D, 2 * D + 1)
        need(a[0] == b[0] == expect and Fr(odd, tot) == Fr(D, 2 * D + 1),
             'F1 D=%d: rho_max = %s (both engines); the MH-periodic orbit through %s (alternating run of length %d, then '
             'valuation 3) is never flipped and has density %s' % (D, a[0], x0, D - 1, Fr(odd, tot)))
        if 2 * D + 1 <= 19:
            v, cyc = rhomax_uniform(rule, 2 * D + 1)
            need(v == a[0], 'F2 D=%d: the seven2 uniform engine at level %d agrees (%s)' % (D, 2 * D + 1, v))


# ------------------------------------------------------------------------------------------------ P/U: searches
def sec_PU():
    section('P/U. Exhaustive searches: pattern sets and itinerary automata')
    for L_, first, later, expect in [(2, (2, 3, ('ge', 4)), (2, 3, ('ge', 4)), Fr(1, 2)),
                                     (3, (2, ('ge', 3)), (2, ('ge', 3)), Fr(3, 7))]:
        leaves, best, dt = S.exhaustive_patterns(L_, first, later)
        need(best[0][0] == expect, 'P1 complete normalized trees of length %d over %s: %d leaves, all %d flip sets; the '
             'least rho_max is %s' % (L_, first, len(leaves), 1 << len(leaves), best[0][0]))
    for m, vc, expect in [(1, (2, ('ge', 3)), Fr(8, 17)), (2, (2, ('ge', 3)), Fr(3, 7)), (1, (2, 3, ('ge', 4)), Fr(8, 17))]:
        bestv, cnt = Fr(1), 0
        for A in S.enum_automata(m, vc):
            for default in (False, True):
                v, crit, tree = S.eval_automaton(A, 16, default, vc)
                cnt += 1
                bestv = min(bestv, v)
        need(bestv == expect, 'U1 all %d itinerary automata (start + %d transient states, valuation classes %s, '
             'truncated at bit depth 16, both defaults): least rho_max = %s' % (cnt, m, vc, bestv))


# ------------------------------------------------------------------------------------------------ L: local search trees
def sec_L(data):
    section('L. Counterexample-guided itinerary search: the stored trees, re-evaluated exactly')
    for name in data['trees']:
        tr = S.tree_from_json(data['trees'][name]['tree'])
        claimed = Fr(data['trees'][name]['value'])
        rule = IT.to_rule(tr)
        check_partition(rule)
        D = max(d for (_, d) in rule)
        a = rho_markov_c(rule)
        b = rho_markov(rule)
        need(a[0] == b[0] == claimed, 'L1 tree "%s": %d itinerary leaves (%d residue classes, max depth %d), rho_max = %s '
             'by both Markov engines (%d Markov leaves), witness period %d' % (name, len(tr), len(rule), D, a[0], a[1], a[3]))
        if D <= 22:
            v, cyc = rhomax_uniform(rule, D)
            need(v == claimed, 'L2 tree "%s": the seven2 uniform engine at level %d agrees' % (name, D))
        h = canon_hash(data['trees'][name]['tree'])
        ok('L3 tree "%s" canonical sha256[:16] = %s' % (name, h))
    # the search is deterministic: re-run it from MH and compare the snapshots
    snaps = {}

    def cb(t, v, it):
        snaps[str(v)] = (it, t)
    best, bval, hist = IT.local_search(IT.mh_tree(), iters=400, maxdepth_sym=16, verbose=False, tabu_len=80, on_best=cb)
    for name, spec in data['trees'].items():
        it, t = snaps[spec['value']]
        need(it == spec['move'] and canon_hash(S.tree_to_json(t)) == canon_hash(spec['tree']),
             'L4 re-running the counterexample-guided search from MH reproduces tree "%s" at move %d' % (name, it))
    need(bval == Fr(7, 18), 'L5 the search from MH (400 moves) ends at 7/18 = rho*(7,16)')
    tr = S.tree_from_json(data['trees']['itinerary-search 7/18']['tree'])
    meas, tot = {}, Fr(0)
    for w, f in tr.items():
        if not f:
            continue
        c, d = word_class(IT.signed_words(w)[0])
        m = Fr(4, 1 << d)          # both signed cylinders, as a share of the odd numbers (measure 1/2)
        key = ' '.join(str(t) if i == 0 else t[0] + str(t[1]) for i, t in enumerate(w[:2]))
        meas[key] = meas.get(key, Fr(0)) + m
        tot += m
    assert all(k.startswith('2 ') or k == '2' for k in meas)
    need(True, 'L6 the 7/18 tree flips %.4f of the odd numbers, all with v_1 = 2: %s' % (
        float(tot), ', '.join('%s: %.4f' % (k, float(v)) for k, v in sorted(meas.items(), key=lambda t: -t[1]))))


def canon_hash(js):
    return hashlib.sha256('\n'.join(sorted(json.dumps(x) for x in js)).encode()).hexdigest()[:16]


# ------------------------------------------------------------------------------------------------ G: greedy trees
def sec_G(data):
    section('G. Greedy optimal trees (THM-4508 construction, this lane\'s driver): stored, re-evaluated exactly')
    for name, spec in data.get('greedy', {}).items():
        rule = {(c, d): s for c, d, s in spec['leaves']}
        check_partition(rule)
        F = Fr(spec['value'])
        k = spec['k']
        a = rho_markov_c(rule)
        b = rho_markov(rule)
        v, cyc = rhomax_uniform(rule, k)
        x0, aa, nn = check_uniform_cycle(rule, k, cyc)
        need(a[0] == b[0] == v == F and a[1] == b[1],
             'G1 %s: %d leaves, rho_max = %s by both Markov engines (%d Markov leaves vs 2^%d = %d) and by the seven2 '
             'uniform engine at level %d' % (name, len(rule), F, a[1], k, 1 << k, k))


# ------------------------------------------------------------------------------------------------ O: optimal rules
def sec_O():
    section('O. Optimal level-k rules (Min least fixed point, MH wherever admissible), minimal trees, Markov sizes')
    for k in (12, 14, 16, 18, 19, 20):
        F = RHO_STAR[k]
        pot, sig = game_lfp('min', k, F, mhtie=True)
        if not check_min_cert(k, F, pot, sig):
            fail('O cert k=%d' % k)
        rule = compress_strategy(sig, k)
        check_partition(rule)
        nflip = sum(1 for (c, d), s in rule.items() if s != mh(c))
        v, nl, x0, p = rho_markov_c(rule)
        need(v == F, 'O1 k=%d: optimal rule with %d leaves (%d flipped), Markov refinement %d leaves = %.1f%% of 2^k; '
             'rho_max = %s = rho*(7,%d) by the Markov engine' % (k, len(rule), nflip, nl, 100.0 * nl / (1 << k), v, k))


# ------------------------------------------------------------------------------------------------ I: itinerary structure
def sec_I():
    section('I. Forced decisions of the least Min potential: itinerary cylinders vs bit cylinders (EMPIRICAL)')
    for k in (14, 16, 18):
        F = RHO_STAR[k]
        pot, sig = game_lfp('min', k, F)
        pot = pot.astype(np.int64)
        N, H = 1 << k, 1 << (k - 1)
        x = np.arange(1, N, 2, dtype=np.int64)
        e = F.denominator - F.numerator
        okk = {}
        for s in (1, -1):
            P = ((7 * x + s) >> 1) % H
            okk[s] = e + np.maximum(pot[P], pot[P + H]) <= pot[x]
        mhs = np.where((7 * x + 1) % 4 == 0, 1, -1)
        ok_mh = np.where(mhs == 1, okk[1], okk[-1])
        ok_fl = np.where(mhs == 1, okk[-1], okk[1])
        assert (ok_mh | ok_fl).all()
        cls = np.where(ok_mh & ok_fl, 0, np.where(ok_mh, 1, 2))
        rows = []
        for n in (3, 4, 5):
            keys = np.array([hash(tuple(t for t in itin_norm(int(xx), k, n))) for xx in x], dtype=np.int64)
            ci, nk = consistency(keys, cls)
            b = 1
            while (1 << (b - 1)) < nk and b < k:
                b += 1
            cb, nb = consistency(x % (1 << b), cls)
            rows.append('n=%d: %d cyl %.3f vs b=%d: %d cyl %.3f' % (n, nk, ci, b, nb, cb))
        ok('I1 k=%d (F=%s): free %.3f forced-MH %.3f forced-flip %.3f; conflict-free share, itinerary vs bits: %s' % (
            k, F, (cls == 0).mean(), (cls == 1).mean(), (cls == 2).mean(), '; '.join(rows)))


def itin_norm(x, k, n):
    sy, tail = itin_of_residue(x, k, maxsym=n)
    toks, prev = [], None
    for (s, v) in sy[:n]:
        toks.append((0 if prev is None else (1 if s == prev else 2), min(v, 6)))
        prev = s
    if len(toks) < n:
        toks.append(('U',))
    return toks


def consistency(keys, cls):
    order = np.argsort(keys, kind='stable')
    ks, cs = keys[order], cls[order]
    bounds = np.flatnonzero(np.diff(ks)) + 1
    starts = np.concatenate([[0], bounds])
    ends = np.concatenate([bounds, [len(ks)]])
    has1 = np.add.reduceat((cs == 1).astype(np.int64), starts) > 0
    has2 = np.add.reduceat((cs == 2).astype(np.int64), starts) > 0
    conflict = has1 & has2
    sizes = ends - starts
    return 1 - sizes[conflict].sum() / len(ks), len(starts)


# ------------------------------------------------------------------------------------------------ R: rigidity, shadowing
def sec_R():
    section('R. Lemma R (rigidity of rejoining), Lemma S (shadowing), single flips on S_inf')
    rng = random.Random(20260926)
    # Lemma S: exact identity on random 2-adic integers: z with symbols (-s,v1)(-s,v2)...(-s,vn)(-s,.), vi in {2,3};
    # E = (z - s)/2 has MH valuations v1..vn.
    cnt = 0
    for trial in range(3000):
        s = rng.choice((1, -1))
        n = rng.randint(1, 12)
        vs = [rng.choice((2, 3)) for _ in range(n)]
        word = [(-s, v) for v in vs] + [(-s, ('ge', 2))]
        c, d = word_class(word)
        z = c + (rng.getrandbits(80) << d)      # a random element of the cylinder
        E = (z - s) // 2
        assert (z - s) % 2 == 0 and E % 2 == 1
        x, m = E, d + 70
        got = []
        for i in range(n):
            sg = mh(x % 4)
            t = 7 * x + sg
            vv = v2(t)
            got.append(vv)
            x = t >> vv
        assert got == vs, (s, vs, got)
        cnt += 1
    ok('S1 Lemma S on %d random cylinders: E = (z - s)/2 copies the valuations of z while z has sign -s and '
       'valuations in {2,3}' % cnt)
    # Lemma R statistics: random x, one flip, then MH; every rejoin (agreement mod 2^200) has equal odd counts and halvings
    B = 1000
    rejoin = neutral = 0
    for trial in range(1500):
        x = rng.getrandbits(B) | 1
        orb, xx, m = {}, x, B
        acc = 0
        for i in range(250):
            if m < 260:
                break
            sg = mh(xx % 4)
            t = (7 * xx + sg) % (1 << m)
            vv = v2(t)
            orb[xx % (1 << 200)] = (i, acc)
            acc += vv
            xx, m = t >> vv, m - vv
        s1 = mh(x % 4)
        y, my = ((7 * x - s1) % (1 << B)) >> 1, B - 1
        hy = 1
        for j in range(250):
            if my < 260:
                break
            key = y % (1 << 200)
            if key in orb:
                i, acc_i = orb[key]
                rejoin += 1
                neutral += (j + 1 == i and hy == acc_i)
                break
            sg = mh(y % 4)
            t = (7 * y + sg) % (1 << my)
            vv = v2(t)
            hy += vv
            y, my = t >> vv, my - vv
    need(rejoin > 0 and neutral == rejoin, 'R1 %d random 1000-bit x: %d single-flip orbits rejoined the MH orbit (mod 2^200) '
         'and every rejoin used the same numbers of odd steps and halvings (Lemma R)' % (1500, rejoin))
    # S_inf: flip at a random S_inf point: no rejoin, typical valuations afterwards
    nre, tv, tn = 0, 0, 0
    for trial in range(200):
        signs = [rng.choice((1, -1)) for _ in range(600)]
        c, d = word_class([(sg, 2) for sg in signs])
        orb, xx, m = set(), c, d
        for i in range(600):
            orb.add(xx % (1 << 120))
            t = (7 * xx + signs[i]) % (1 << m)
            xx, m = t >> 2, m - 2
            if m < 200:
                break
        y, my = ((7 * c - signs[0]) % (1 << d)) >> 1, d - 1
        for j in range(400):
            if my < 200:
                break
            if y % (1 << 120) in orb:
                nre += 1
                break
            sg = mh(y % 4)
            t = (7 * y + sg) % (1 << my)
            vv = v2(t)
            if j >= 10:
                tv += vv
                tn += 1
            y, my = t >> vv, my - vv
    need(nre == 0, 'R2 200 random S_inf points: after one flip the orbit never met the MH orbit again (EMPIRICAL); '
         'mean MH valuation afterwards %.3f over %d steps (typical value 3)' % (tv / tn, tn))


# ------------------------------------------------------------------------------------------------ B: adversaries
def sec_B(data):
    section('B. Max: itinerary / rational adversaries (exact values; certificate and maximality inside the engine)')
    for k in range(8, 19):
        lift = np.ones(1 << (k - 1), dtype=np.uint8)
        v, W, mx, cyc = adversary_value(k, lift)
        a, n = check_adv_cycle(k, lift, cyc)
        need(v == Fr(1, 3) and mx and Fr(a, n) == v, 'B1 k=%d: top lift (negative integers) has value exactly 1/3 '
             '(|W| = %d; Min reaches density 1/3 from every node)' % (k, W))
    rows = []
    for d in (1, 3, 5, 11, 13, 17, 19, 23, 29, 43):
        vals = []
        for k in (10, 14, 18):
            H, N = 1 << (k - 1), 1 << k
            P = np.arange(H, dtype=np.int64)
            lift = (((d * P) % N) < H).astype(np.uint8)
            v, W, mx, cyc = adversary_value(k, lift)
            assert mx
            vals.append(v)
        assert max(vals) <= Fr(1, 3)
        rows.append('d=%d: %s' % (d, ','.join(map(str, vals))))
    ok('B2 d-scaled top lifts (follow -c/d), values at k = 10,14,18 (all <= 1/3): ' + '; '.join(rows))
    rows = []
    for k in range(8, 17):
        H = 1 << (k - 1)
        prof = [(mh_profile(P, k), mh_profile(P + H, k)) for P in range(H)]
        la = np.array([0 if a[0] / max(a[1], 1) > b[0] / max(b[1], 1) else 1 for a, b in prof], dtype=np.uint8)
        lb = np.array([0 if (a[0], -a[1]) > (b[0], -b[1]) else 1 for a, b in prof], dtype=np.uint8)
        va = adversary_value(k, la)
        vb = adversary_value(k, lb)
        assert va[2] and vb[2]
        rows.append('k=%d: %s, %s' % (k, va[0], vb[0]))
        if k >= 12:
            assert va[0] <= Fr(1, 3) and vb[0] <= Fr(1, 3)
    ok('B4 MH-density lifts (prefer the lift whose determined MH itinerary is denser; by density, by (odd steps, -bits)): '
       + '; '.join(rows) + ' (all <= 1/3 for k >= 12)')
    for name, spec in data.get('adversaries', {}).items():
        n = spec['n']
        f = {tuple(c): b for c, b in spec['f']}
        vals = []
        for k in spec['levels']:
            H = 1 << (k - 1)
            lift = np.array([f.get(tuple(prefix_code(P, k - 1, n)), 1) for P in range(H)], dtype=np.uint8)
            v, W, mx, cyc = adversary_value(k, lift)
            assert mx
            vals.append(v)
        need([str(v) for v in vals] == spec['values'], 'B3 itinerary-prefix adversary "%s" (n=%d symbols): values %s at '
             'k = %s' % (name, n, ','.join(map(str, vals)), spec['levels']))


def mh_profile(t, k):
    """(odd steps, bits consumed) of the MH itinerary determined by the residue t mod 2^k (even steps are halvings)"""
    x, m = t, k
    n_odd, bits = 0, 0
    while m >= 1:
        if x % 2 == 0:
            if x == 0:
                break
            x >>= 1
            m -= 1
            bits += 1
            continue
        if m < 2:
            break
        s = mh(x)
        u = (7 * x + s) % (1 << m)
        if u == 0:
            break
        v = v2(u)
        if v + 1 > m:
            break
        n_odd += 1
        bits += v
        x = u >> v
        m -= v
    return n_odd, bits


def prefix_code(P, m, n):
    toks = []
    x, mm = P, m
    while len(toks) < n:
        if mm <= 0:
            toks.append('U')
            continue
        if x % 2 == 0:
            if x == 0:
                toks.append('Z')
                mm = 0
                continue
            toks.append('E')
            x >>= 1
            mm -= 1
            continue
        if mm < 2:
            toks.append('U')
            mm = 0
            continue
        s = mh(x)
        t = (7 * x + s) % (1 << mm)
        if t == 0 or v2(t) + 1 > mm:
            toks.append('%s>=%d' % ('+' if s > 0 else '-', mm))
            mm = 0
            continue
        v = v2(t)
        toks.append(('+' if s > 0 else '-') + str(min(v, 5)))
        x, mm = t >> v, mm - v
    return toks


# ------------------------------------------------------------------------------------------------ C: thresholds
def sec_C():
    section('C. Exact threshold comparisons')
    need(7 ** 37 > 2 ** 100 and 7 ** 7 > 2 ** 19, 'C1 37/100 and 7/19 exceed log_7 2 (7^37 > 2^100, 7^7 > 2^19)')
    for F in (Fr(2, 5), Fr(3, 7), Fr(7, 18), Fr(1, 3)):
        side = 7 ** F.numerator > 2 ** F.denominator
        ok('C2 %s %s log_7 2 = %.5f' % (F, '>' if side else '<', log(2) / log(7)))


def main():
    print('procgen_seven3_20260926_run.py -- lane seven3, itinerary-coded sign strategies of 7n+-1')
    for name in ('procgen_seven3_20260926_lib.py', 'procgen_seven3_20260926_itree.py', 'procgen_seven3_20260926_search.py',
                 'procgen_seven3_20260926_markov.c', 'procgen_seven3_20260926_game.c', 'procgen_seven3_20260926_adv.c',
                 'procgen_seven2_20260926_rhomax.c', 'procgen_seven3_20260926_rules.json'):
        print('  sha256 %s %s' % (L.sha(os.path.join(HERE, name))[:16], name))
    data = json.load(open(DATA))
    sec_C()
    sec_M()
    sec_F()
    sec_L(data)
    sec_G(data)
    sec_O()
    sec_I()
    sec_R()
    sec_B(data)
    sec_PU()
    sec_D()
    print('\n%d checks, %.0f s' % (NCHECK[0], time.time() - T0))
    print('ALL CHECKS PASSED')


if __name__ == '__main__':
    main()
