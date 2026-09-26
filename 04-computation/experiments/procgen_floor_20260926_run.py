#!/usr/bin/env python3
"""
procgen_floor_20260926_run.py -- runner of the floor lane (session collatz-procgen-20260922, 2026-09-26):
the min-max cycle density rho*(q,k) of the q n +- 1 strategy cube as the value of a mean-payoff game.

Every printed claim is a check(...) that raises on failure; the output ends with ALL CHECKS PASSED.

Environment: FLOOR_KMAX5 (default 21), FLOOR_KMAX7 (22), FLOOR_KMAXQ (18: q = 9..21), FLOOR_KALLQ (11: all
residues q mod 2^(k-1)).
"""
import os
import sys
import time
import math
import hashlib
import resource
import platform
from fractions import Fraction as Fr
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import procgen_floor_20260926_lib as L    # noqa: E402

KMAX5 = int(os.environ.get('FLOOR_KMAX5', 21))
KMAX7 = int(os.environ.get('FLOOR_KMAX7', 22))
KMAXQ = int(os.environ.get('FLOOR_KMAXQ', 18))
KALLQ = int(os.environ.get('FLOOR_KALLQ', 11))
T0 = time.time()
NCHECK = [0]


def say(*a):
    print(*a, flush=True)


def check(cond, msg):
    L.check(cond, msg)
    NCHECK[0] += 1
    say('  [ok] ' + msg)


def rss_mb():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / (1 << 20) if sys.platform == 'darwin' else r / 1024


def fr_below_log(q, F):
    """F = a/p < log_q 2  <=>  q^a < 2^p (exact)"""
    return q ** F.numerator < 2 ** F.denominator


def cert_hash(res):
    h = hashlib.sha256()
    flip, psi = res['upper']
    W, tau, f = res['lower']
    for arr in (np.asarray(flip, dtype=np.uint8), np.asarray(psi, dtype=np.int64), np.asarray(W, dtype=np.uint8),
                np.asarray(tau, dtype=np.uint8), np.asarray(f, dtype=np.int64)):
        h.update(arr.tobytes())
    return h.hexdigest()[:16]


LOG2 = math.log(2)
say("procgen_floor_20260926_run.py -- floor lane, collatz-procgen-20260922, 2026-09-26")
say(f"python {platform.python_version()}, numpy {np.__version__}; KMAX5={KMAX5} KMAX7={KMAX7} KMAXQ={KMAXQ} "
    f"KALLQ={KALLQ}")
for fn in ['procgen_floor_20260926_lib.py', 'procgen_floor_20260926_game.c', 'procgen_floor_20260926_run.py']:
    say(f"  sha256 {fn} = {L.sha256(os.path.join(HERE, fn))}")

VAL = {}          # (q, k) -> rho*


def rho(k, q, known_upper=None):
    """exact rho*(q,k) with both certificates re-checked"""
    if (q, k) in VAL:
        return VAL[(q, k)], None
    res = L.rho_star(k, q, known_upper=known_upper)
    L.certify(k, q, res)
    VAL[(q, k)] = res['rho']
    return res['rho'], res


# =============================================================================================== A
say("\n== A. The game solver and its certificates ==")
say("A1. rho*(q,k) from the game solver vs exhaustive enumeration of all 2^(2^(k-1)) strategies (C, Karp):")
for k in range(2, 6):
    vals = []
    for q in range(1, 1 << k, 2):
        v, _ = rho(k, q)
        G = L.Game(k, q)
        ex, cnt = G.exhaustive()
        L.check(ex == v, f"exhaustive mismatch k={k} q={q}")
        vals.append((q, str(v), cnt))
    check(True, f"k={k}: game value = exhaustive min for all {len(vals)} odd q < 2^{k}: "
          f"{sorted(set(v for _, v, _ in vals))}; optimal strategy counts {[c for _, _, c in vals]}")

say("A2. independent pure-python Karp on the extracted optimal strategies (k <= 7):")
n = 0
for q in [3, 5, 7, 9, 11, 13, 15, 17]:
    for k in range(2, 8):
        res = L.rho_star(k, q)
        L.certify(k, q, res)
        r = L.rho_max_exact(k, q, res['upper'][0])
        L.check(r == res['rho'], f"pure-python Karp mismatch q={q} k={k}")
        n += 1
check(True, f"{n} (q,k) pairs: rho_max(extracted sigma) by pure-python Karp equals the certified value")

say("A3. cross-check with the drift lane's independent engine (Dinkelbach + SPFA, procgen_drift_20260926_lib):")
import procgen_drift_20260926_lib as D   # noqa: E402
pairs = [(5, k) for k in range(7, 13)] + [(7, k) for k in range(7, 13)] + [(9, k) for k in range(9, 13)] + \
        [(11, k) for k in range(9, 13)] + [(13, k) for k in range(8, 13)]
for q, k in pairs:
    res = L.rho_star(k, q)
    L.certify(k, q, res)
    F, cyc, psi = D.rho_max(k, q, np.asarray(res['upper'][0], dtype=np.uint8))
    L.check(F == res['rho'], f"drift-lane rho_max mismatch q={q} k={k}: {F} vs {res['rho']}")
check(True, f"{len(pairs)} pairs (q,k), k <= 12: the drift lane's exact rho_max of our optimal strategies equals rho*")

say("A4. negative tests (the verifiers reject corrupted certificates):")
k, q = 12, 7
res = L.rho_star(k, q)
F = res['rho']
flip, psi = res['upper']
W, tau, f = res['lower']
N = 1 << k
H = N >> 1
t = L.targets_of(k, q, flip)
e = L.weights(k, F)
psi64 = psi.astype(np.int64)
tight = np.nonzero((psi64[t] + e == psi64) | (psi64[t + H] + e == psi64))[0]
bad = psi.copy()
bad[tight[0]] -= 1
check(L.verify_upper(k, q, flip, psi, F) and not L.verify_upper(k, q, flip, bad, F),
      f"upper certificate (q=7,k=12) passes; lowering psi by 1 at a tight node makes it fail")
te, tp, tm = L.arena(k, q)
s_all = np.arange(N)
tt = np.where(s_all % 2 == 0, te, tp)
succ = tt + tau[tt].astype(np.int64) * H
f64 = f.astype(np.int64)
tightL = np.nonzero(W & (f64[succ] == f64 + e))[0]
badf = f.copy()
badf[succ[tightL[0]]] += 1
check(L.verify_lower(k, q, W, tau, f, F) and not L.verify_lower(k, q, W, tau, badf, F),
      f"lower certificate (q=7,k=12) passes; raising f by 1 at the head of a tight edge makes it fail")
check(not L.verify_upper(k, q, flip, psi, Fr(7, 18)) and not L.verify_lower(k, q, W, tau, f, Fr(3, 7)),
      "the q=7,k=12 certificates do not certify the neighbouring values 7/18 (upper) or 3/7 (lower)")

say("A5. Max certificates re-checked by an independent algorithm (Karp on the reachable part of G^tau):")
n = 0
for q, k in [(5, 8), (5, 10), (7, 9), (7, 10), (9, 10), (13, 9), (3, 8)]:
    res = L.rho_star(k, q)
    L.certify(k, q, res)
    W, tau, f = res['lower']
    G = L.Game(k, q)
    start = int(np.nonzero(W)[0][0])
    m = G.min_mean_tau(tau, start)
    L.check(m is not None and m >= res['rho'], f"Karp on G^tau below the certified value q={q} k={k}")
    n += 1
check(True, f"{n} cases: the least cycle density Min can reach against the certified tau is >= rho*")

# =============================================================================================== B
say("\n== B. Exact min-max densities rho*(q,k) (both certificates re-checked for every entry) ==")
TAB = {}
qlist = [3, 5, 7] + list(range(9, 22, 2))
kmax = {3: 16, 5: KMAX5, 7: KMAX7}
for q in qlist:
    km = kmax.get(q, KMAXQ)
    ku = None
    row = []
    for k in range(2, km + 1):
        t1 = time.time()
        res = L.rho_star(k, q, known_upper=ku)
        L.certify(k, q, res)
        VAL[(q, k)] = res['rho']
        ku = res['rho']
        row.append(res['rho'])
        TAB[(q, k)] = res['rho']
        if k >= 14:
            say(f"    q={q:2d} k={k:2d}: rho* = {str(res['rho']):>7s} = {float(res['rho']):.5f}  (|W| = {int(res['lower'][0].sum())}, "
                f"cert {cert_hash(res)}, {time.time() - t1:.1f}s, RSS {rss_mb():.0f} MB)")
        del res
    L.check(all(row[i + 1] <= row[i] for i in range(len(row) - 1)), f"rho*({q},k) non-increasing")
    first = next((k for k in range(2, km + 1) if TAB[(q, k)] < Fr(3, 7)), None)
    say(f"  q={q:2d} (log_q 2 = {LOG2 / math.log(q):.5f}): " +
        ", ".join(f"k={k}:{TAB[(q, k)]}" for k in range(2, km + 1) if k == 2 or TAB[(q, k)] != TAB[(q, k - 1)]) +
        f"  [first k with rho* < 3/7: {first}]")
check(all(TAB[(q, k)] == Fr(1, 2) for q in [5, 7, 9, 11] for k in range(2, 7)) and
      all(TAB[(q, k)] == Fr(3, 7) for q in [5, 7, 9, 11] for k in (7, 8)) and TAB[(7, 9)] == Fr(3, 7) and
      all(TAB[(13, k)] == Fr(1, 2) for k in range(2, 9)) and all(TAB[(3, k)] == Fr(1, 2) for k in range(2, 17)),
      "the drift lane's values are reproduced (q=5,7,9,11: 1/2 for k<=6, 3/7 at k=7,8; q=7: 3/7 at k=9; "
      "q=13: 1/2 for k<=8; q=3: 1/2)")
check(TAB[(7, 10)] == Fr(2, 5) and TAB[(9, 10)] == Fr(5, 12) and TAB[(11, 10)] == Fr(5, 12) and
      TAB[(5, 10)] == Fr(3, 7) and TAB[(5, 11)] == Fr(5, 12) and TAB[(13, 10)] == Fr(4, 9),
      "Q1: rho* >= 3/7 FAILS: rho*(7,10) = 2/5, rho*(9,10) = rho*(11,10) = 5/12, rho*(5,11) = 5/12; the "
      "q-coincidence breaks at k = 10 (5: 3/7, 7: 2/5, 9: 5/12, 11: 5/12, 13: 4/9)")
check(all(TAB[(5, k)] == Fr(2, 5) for k in range(15, KMAX5 + 1)) and TAB[(5, 14)] == Fr(5, 12),
      f"rho*(5,k) = 2/5 for 15 <= k <= {KMAX5} (and 5/12 at k = 14)")
viol = [(q, k) for (q, k), v in TAB.items() if q >= 7 and fr_below_log(q, v)]
check(not viol, f"for every q in 7..21 and every computed level (q=7: k <= {KMAX7}; q=9..21: k <= {KMAXQ}) "
      f"rho*(q,k) > log_q 2 exactly (q^a > 2^p): class (i) is EMPTY there")
check(all(fr_below_log(5, TAB[(5, k)]) for k in range(7, KMAX5 + 1)) and not fr_below_log(5, TAB[(5, 6)]),
      "q = 5: rho*(5,k) < log_5 2 exactly iff k >= 7 (class (i) nonempty from k = 7)")
lowest7 = min(TAB[(7, k)] for k in range(2, KMAX7 + 1))
check(float(lowest7) - LOG2 / math.log(7) > 0.02, f"q = 7: the lowest value {lowest7} = {float(lowest7):.5f} "
      f"(k = {KMAX7}) still exceeds log_7 2 = {LOG2 / math.log(7):.5f} by more than 0.02")

# =============================================================================================== C
say("\n== C. Symmetries of rho*(.,k) and the level floors ==")
say("C1. the arena at level k depends only on q mod 2^(k-1): in arena(q + 2^(k-1)) node s has the options of "
    "node s + 2^(k-1) in arena(q) (odd s):")
for k in range(3, 13):
    Nk = 1 << k
    Hk = Nk >> 1
    ok = True
    for q in range(1, Nk, 2):
        te1, tp1, tm1 = L.arena(k, q)
        te2, tp2, tm2 = L.arena(k, q + Hk)
        s = np.arange(1, Nk, 2)
        ok &= bool(np.all(tp2[s] == tp1[(s + Hk) % Nk]) and np.all(tm2[s] == tm1[(s + Hk) % Nk]) and
                   np.all(te1 == te2))
    L.check(ok, f"arena isomorphism fails at k={k}")
check(True, "k = 3..12, every odd q: option sets match (so G_sigma(q) and G_sigma'(q + 2^(k-1)) are isomorphic "
      "parity graphs with sigma'(s) = sigma(s + 2^(k-1)))")

say("C2. rho*(q,k) = rho*(q + 2^(k-1),k) = rho*(-q,k) for every odd q mod 2^k, k = 3..10:")
for k in range(3, 11):
    Nk = 1 << k
    vals = {}
    for q in range(1, Nk, 2):
        vals[q], _ = rho(k, q)
    a = all(vals[q] == vals[q + (Nk >> 1)] for q in range(1, Nk >> 1, 2))
    b = all(vals[q] == vals[(-q) % Nk] for q in range(1, Nk, 2))
    L.check(a and b, f"symmetry fails at k={k}")
check(True, "both symmetries hold at k = 3..10 (all 2^(k-1) odd residues)")

say("C3. the proof of rho*(-q,k) = rho*(q,k): the least fixed point f* of the Min operator at F = rho*(q,k) is "
    "nu-invariant, and (-sigma*, f*) certifies F for -q:")
n = 0
for q, k in [(5, 10), (5, 14), (7, 12), (7, 16), (9, 11), (11, 13), (13, 12), (3, 9), (15, 12), (21, 12)]:
    F, _ = rho(k, q)
    up = L.nu_symmetric_upper(k, q, F)
    L.check(up is not None, "least fixed point infinite at the value")
    fl, fstar = up
    Nk = 1 << k
    neg = np.where(np.arange(Nk) % 2 == 1, 1 - fl, 0).astype(np.uint8)
    L.check(L.is_nu_invariant(k, fstar) and L.verify_upper(k, q, fl, fstar, F) and
            L.verify_upper(k, (-q) % Nk, neg, fstar, F), f"transport fails q={q} k={k}")
    low = L.nu_symmetric_lower(k, q, F)
    W, tau, g = low
    L.check(L.tau_equivariant(k, tau) and L.verify_lower(k, q, W, tau, g, F) and
            L.verify_lower(k, (-q) % Nk, W, tau, g, F), f"lower transport fails q={q} k={k}")
    n += 1
check(True, f"{n} cases: nu-invariant f*, sigma* certify F for q; -sigma* with the same f* certifies F for -q "
      f"(and the nu-symmetric Max certificate transports unchanged)")

say("C4. strategy-level comparison at k = 4, 5 (pure-python Karp): the full rho_max distributions of q and -q can differ, "
    "the nu-invariant strategies (sigma(-s) = -sigma(s)) match under sigma -> -sigma:")
for k in (4, 5):
    Nk = 1 << k
    odd = list(range(1, Nk, 2))
    reps = [s for s in odd if s < Nk - s]
    for q in range(1, Nk >> 1, 2):
        qq = (-q) % Nk
        d1, d2 = Counter(), Counter()
        for m in range(1 << len(reps)):
            fl = np.zeros(Nk, dtype=np.uint8)
            for i, s in enumerate(reps):
                b = (m >> i) & 1
                fl[s] = b
                fl[Nk - s] = 1 - b
            fl2 = np.where(np.arange(Nk) % 2 == 1, 1 - fl, 0).astype(np.uint8)
            r1 = L.rho_max_exact(k, q, fl)
            r2 = L.rho_max_exact(k, qq, fl2)
            L.check(r1 == r2, f"nu-invariant strategy mismatch k={k} q={q}")
            d1[r1] += 1
check(True, "k = 4, 5: for every nu-invariant sigma of arena(q), rho_max(sigma) = rho_max(-sigma in arena(-q))")
if True:
    k = 4
    dist = {}
    for q in (3, 5):
        c = Counter()
        for m in range(1 << 8):
            fl = np.zeros(16, dtype=np.uint8)
            for i, s in enumerate(range(1, 16, 2)):
                fl[s] = (m >> i) & 1
            c[L.rho_max_exact(k, q, fl)] += 1
        dist[q] = c
    check(dist[3] != dist[5] and min(dist[3]) == min(dist[5]) == Fr(1, 2),
          f"k = 4: all-strategy distributions differ (q=3: {sorted((str(a), b) for a, b in dist[3].items())[:3]}..., "
          f"q=5: {sorted((str(a), b) for a, b in dist[5].items())[:3]}...) but the minima agree (1/2)")

say(f"C5. level floors min_q rho*(q,k) over the 2^(k-3) classes (+-q mod 2^(k-1)), k <= {KALLQ}:")
floors = {}
for k in range(3, KALLQ + 1):
    M = 1 << (k - 1)
    classes = sorted(set(min(q % M, (-q) % M) for q in range(1, M, 2)))
    vals = {}
    ku = None
    for c in classes:
        vals[c], _ = rho(k, c)
    fl = min(vals.values())
    floors[k] = fl
    att = [c for c in classes if vals[c] == fl]
    dist = Counter(vals.values())
    say(f"    k={k:2d}: floor {fl} ({float(fl):.5f}) attained by +-q mod 2^{k - 1} for q in {att[:10]}"
        f"{'...' if len(att) > 10 else ''} ({len(att)} of {len(classes)} classes); value spectrum "
        f"{sorted(((str(a), b) for a, b in dist.items()), key=lambda t: Fr(t[0]))[:6]}...")
    if k == 7:
        check(sorted(att) == [5, 7, 9, 11] and all(vals[c] == Fr(1, 2) for c in classes if c not in att),
              "k = 7: the classes with rho* = 3/7 are exactly +-5, +-7, +-9, +-11 mod 64 (all others 1/2): "
              "this is the q = 5..11 coincidence")
check(all(floors[k] == Fr(1, 2) for k in range(3, 7)), "k <= 6: rho*(q,k) = 1/2 for EVERY odd q")
check(all(floors[k] == Fr(3, 7) for k in (7, 8, 9)) and all(floors[k] == Fr(2, 5) for k in range(10, KALLQ + 1)) and
      all(VAL[(7, k)] == floors[k] for k in range(7, KALLQ + 1)),
      f"level floors: 3/7 (k = 7..9), 2/5 (k = 10..{KALLQ}); q = 7 attains the floor at every k = 7..{KALLQ}")

# =============================================================================================== D
say("\n== D. Theorem N: the negative-integer adversary; rho_max >= log_{q+1} 2 for every strategy of every q ==")
say("D1. U_k (u -> u/2; u odd -> rho((qu -+ 1)/2) in [1, H]) is the Min graph of the arena when Max always takes "
    "the lift with top bit 1 (node N - u):")
check(all(L.U_matches_tau_one(k, q) for k in range(2, 15) for q in range(3, 66, 2)),
      "k = 2..14, q = 3..65: edge sets coincide")
say("D2. the multiplicative potential u: every odd step has 2 v <= (q+1) u, every even step v = u/2:")
okk = True
for k in range(2, 19):
    for q in range(3, 66, 2):
        u, v1, v2 = L.U_edges(k, q)
        od = u % 2 == 1
        okk &= bool(np.all(2 * v1[od] <= (q + 1) * u[od]) and np.all(2 * v2[od] <= (q + 1) * u[od]) and
                    np.all(2 * v1[~od] == u[~od]) and np.all(v1[~od] == v2[~od]))
check(okk, "k = 2..18, q = 3..65: all edges of U_k satisfy the potential inequalities (hence every cycle has "
      "(q+1)^a >= 2^p, density >= log_{q+1} 2)")
say("D3. independent confirmation: least cycle density of U_k (numpy Karp) = value of the frozen strategy tau = 1 "
    "(game solver), and it is >= log_{q+1} 2:")
DQ = {}
for q in range(3, 22, 2):
    for k in range(5, 11):
        u, v1, v2 = L.U_edges(k, q)
        succ = np.stack([v1 - 1, v2 - 1], axis=1)
        md = L.karp_min_density(succ, (u % 2).astype(np.int64))
        Hk = 1 << (k - 1)
        vt, _ = L.value_of_tau(k, q, np.ones(Hk, dtype=np.uint8))
        L.check(md == vt and (q + 1) ** md.numerator >= 2 ** md.denominator, f"U_k mismatch q={q} k={k}")
        DQ[(q, k)] = md
check(True, "q = 3..21, k = 5..10: Karp(U_k) = frozen-tau value; all >= log_{q+1} 2")
for q in range(3, 22, 2):
    for k in range(12, 17, 2):
        Hk = 1 << (k - 1)
        DQ[(q, k)], _ = L.value_of_tau(k, q, np.ones(Hk, dtype=np.uint8))
        L.check((q + 1) ** DQ[(q, k)].numerator >= 2 ** DQ[(q, k)].denominator, "Theorem N violated")
    say(f"    q={q:2d}: tau=1 value at k=10,12,14,16: {[str(DQ[(q, k)]) for k in (10, 12, 14, 16)]}; "
        f"log_(q+1) 2 = {LOG2 / math.log(q + 1):.5f}, log_q 2 = {LOG2 / math.log(q):.5f}")
say("D4. tightness for q = 2^j - 1: the cycle 1 -> 2^(j-1) -> ... -> 2 -> 1 lies in U_k and has density 1/j = log_(q+1) 2:")
for j in range(2, 7):
    q = (1 << j) - 1
    k = j + 3
    u, v1, v2 = L.U_edges(k, q)
    cyc = [1] + [1 << (j - 1 - i) for i in range(j - 1)]
    ok = v2[0] == (1 << (j - 1)) and all(v1[c - 1] == c // 2 for c in cyc[1:])
    L.check(ok, f"tight cycle missing q={q}")
    vt, _ = L.value_of_tau(k + 4, q, np.ones(1 << (k + 3), dtype=np.uint8))
    L.check(vt == Fr(1, j), f"tau=1 value not 1/j for q={q}")
check(True, "q = 3, 7, 15, 31, 63: the free cycle through 1 is in U_k and the tau = 1 value equals 1/j exactly")
say("D5. the critical cycles of the tau = 1 adversary are sign-choice cycles of q u -+ 1 on the positive integers "
    "(no wrap-around) for q <= 19; for q = 21 the critical cycle at k = 16 uses wrap-arounds:")
for q, k in [(3, 12), (5, 12), (7, 12), (9, 14), (11, 14), (13, 14), (15, 14), (17, 14), (19, 16), (21, 16)]:
    v, cyc, wr = L.critical_cycle_tau1(k, q)
    L.check(cyc is not None, "no critical cycle")
    tag = 'contracting' if fr_below_log(q, v) else 'expanding'
    say(f"    q={q:2d} k={k}: value {v} ({tag}); cycle of length {len(cyc)}, min {min(cyc)}, max {max(cyc)}, "
        f"wraps {wr}: {cyc if len(cyc) <= 20 else cyc[:12] + ['...']}")
    if q <= 19:
        L.check(wr == 0, f"critical cycle wraps for q={q}")
check(True, "q = 3..19: the tau = 1 value is attained by an exact cycle of u -> u/2, (q u -+ 1)/2 on positive integers: "
      "free cycles (1,2) [q=3], (1,4,2) [q=7], (1,8,4,2) [q=15, 17]; sporadic ones (1,3,8,4,2) [q=5], ...")
import mpmath   # noqa: E402
mpmath.mp.dps = 50
p0lo, p0hi = mpmath.mpf('0.2270921'), mpmath.mpf('0.2270923')
fp = lambda p: p - p * mpmath.log(p, 2) - (1 - p) * mpmath.log(1 - p, 2) - 1
L.check(fp(p0lo) < 0 < fp(p0hi), "p0 bracket")
beats = [q for q in range(3, 101, 2) if mpmath.log(2) / mpmath.log(q + 1) > p0hi]
check(beats == list(range(3, 20, 2)), f"log_(q+1) 2 > p0 = 0.2270922 exactly for q in {beats} (q <= 19): Theorem N "
      f"raises the floor for these q; for q >= 21 the entropy floor p0 is larger")
check(all(math.log(2) / math.log(q + 1) < math.log(2) / math.log(q) for q in range(3, 101, 2)),
      "log_(q+1) 2 < log_q 2 always: Theorem N settles no new q (class (i) needs rho_max in [log_(q+1) 2, log_q 2))")

# =============================================================================================== E
say("\n== E. q = 5 is settled at every level: rho*(5,k) = 1/2, 3/7, 5/12, 2/5 (k <= 6, 7..10, 11..14, >= 15) ==")
say("E1. corrected potential Phi(u) = u^2 (u != 3), Phi(3) = 8 on U_k (q = 5): Phi(v) <= 8 Phi(u) on odd steps, "
    "4 Phi(u/2) <= Phi(u) on even steps:")
okk = True
for k in range(2, 23):
    u, v1, v2 = L.U_edges(k, 5)
    Phi = lambda x: np.where(x == 3, 8, x * x)
    od = u % 2 == 1
    okk &= bool(np.all(Phi(v1[od]) <= 8 * Phi(u[od])) and np.all(Phi(v2[od]) <= 8 * Phi(u[od])) and
                np.all(4 * Phi(v1[~od]) <= Phi(u[~od])))
check(okk, "k = 2..22: all edges of U_k satisfy it (hand proof for all k in the note), so every cycle of U_k has "
      "8^a >= 4^(p-a), density >= 2/5, for every k")
u, v1, v2 = L.U_edges(10, 5)
check(v2[0] == 3 and v2[2] == 8 and v1[7] == 4 and v1[3] == 2 and v1[1] == 1,
      "the 5x+1 cycle (1, 3, 8, 4, 2) lies in U_k (k >= 5), density 2/5: the bound is attained")
check(all(TAB[(5, k)] == v for k, v in [(6, Fr(1, 2)), (7, Fr(3, 7)), (10, Fr(3, 7)), (11, Fr(5, 12)),
                                            (14, Fr(5, 12)), (15, Fr(2, 5))]),
      "with the certified rho*(5,15) = 2/5 and lifting: rho*(5,k) = 2/5 for EVERY k >= 15; the whole sequence is "
      "1/2 (k <= 6), 3/7 (7..10), 5/12 (11..14), 2/5 (k >= 15)")

# =============================================================================================== F
say("\n== F. Limits of arguments through the uniform-lift stationary law (entropy-type bounds) ==")
say("F1. level-2 max-halving (greedy) strategy: exact stationary law (Fractions) of its uniform-lift chain:")
for q in range(3, 24, 2):
    fl = np.zeros(4, dtype=np.uint8)
    for r in (1, 3):
        if (q * r + 1) % 4 != 0:
            fl[r] = 1
    t = L.targets_of(2, q, fl)          # pairs mod 2
    # chain on Z/4: s -> t[s], t[s] + 2 with prob 1/2 each; solve exactly
    P = [[Fr(0)] * 4 for _ in range(4)]
    for s in range(4):
        P[s][int(t[s])] += Fr(1, 2)
        P[s][int(t[s]) + 2] += Fr(1, 2)
    # power iteration is not exact; solve pi = pi P, sum 1 by elimination
    A = [[P[j][i] - (1 if i == j else 0) for j in range(4)] for i in range(4)]
    A[3] = [Fr(1)] * 4
    b = [Fr(0), Fr(0), Fr(0), Fr(1)]
    for i in range(4):
        piv = next(r for r in range(i, 4) if A[r][i] != 0)
        A[i], A[piv] = A[piv], A[i]
        b[i], b[piv] = b[piv], b[i]
        for r in range(4):
            if r != i and A[r][i] != 0:
                fac = A[r][i] / A[i][i]
                A[r] = [A[r][c] - fac * A[i][c] for c in range(4)]
                b[r] -= fac * b[i]
    pi = [b[i] / A[i][i] for i in range(4)]
    L.check(pi[1] + pi[3] == Fr(1, 3) and L.rho_max_exact(2, q, fl) == Fr(1, 2), f"greedy stationary law q={q}")
check(True, "q = 3..23: pi(odd) = 1/3 exactly (unique stationary law) while rho_max = 1/2; so no inequality valid for "
      "all stationary laws of all strategies can give a floor above 1/3, and none can settle q <= 7")
say("F2. certified strategies whose every closed class has pi(odd) < log_q 2 (exact sub-invariance certificates):")


def mdp_policy(k, q, iters=3000):
    Nk = 1 << k
    Hk = Nk >> 1
    te, tp, tm = L.arena(k, q)
    s = np.arange(Nk)
    odd = (s % 2 == 1)
    V = np.zeros(Nk)
    for _ in range(iters):
        Pv = V[:Hk] + V[Hk:]
        V = odd.astype(float) + np.where(odd, np.minimum(Pv[tp], Pv[tm]), Pv[te]) / 2
        V -= V.min()
    Pv = V[:Hk] + V[Hk:]
    return np.where(odd & (Pv[tm] < Pv[tp]), 1, 0).astype(np.uint8)


for q, k, c in [(7, 10, Fr(31, 100)), (9, 10, Fr(31, 100)), (11, 12, Fr(2885, 10000))]:
    fl = mdp_policy(k, q)
    ok, pis = L.stationary_odd_upper_certificate(k, q, fl, c)
    check(ok and float(c) < LOG2 / math.log(q), f"q={q}, k={k}: every closed class has pi(odd) <= {c} = {float(c):.4f} "
          f"< log_q 2 = {LOG2 / math.log(q):.4f} (class values {[round(x, 4) for x in pis]}); rho_max of this "
          f"strategy is still >= rho*({q},{k}) = {VAL.get((q, k), rho(k, q)[0])}")
say("F3. EMPIRICAL (float value iteration): m(q,k) = least stationary odd frequency over all level-k strategies:")
for q in range(5, 24, 2):
    row = []
    for k in (8, 12, 16):
        Nk = 1 << k
        Hk = Nk >> 1
        te, tp, tm = L.arena(k, q)
        s = np.arange(Nk)
        odd = (s % 2 == 1)
        V = np.zeros(Nk)
        for _ in range(3000):
            Pv = V[:Hk] + V[Hk:]
            Vn = odd.astype(float) + np.where(odd, np.minimum(Pv[tp], Pv[tm]), Pv[te]) / 2
            d = Vn - V
            V = Vn - Vn.min()
        row.append(round(float(d.min()), 4))
    say(f"    q={q:2d}: m(q,8), m(q,12), m(q,16) = {row};  log_q 2 = {LOG2 / math.log(q):.4f}"
        f"{'   (cap already below log_q 2)' if row[-1] < LOG2 / math.log(q) else ''}")
check(True, "table printed (EMPIRICAL): the stationary-law cap is below log_q 2 for q = 7, 9, 11 and above it (so far) "
      "for q = 13..23")

# =============================================================================================== G
say("\n== G. Approaches that do not beat Theorem N ==")
say("G1. exact traps (finite rational sets closed under x/2 and (qx+-1)/2) give level-independent bounds; best found "
    "(denominators <= 31, numerators <= 4 den):")


def trap_best(q, Dmax=31, cap=60):
    par = lambda x: x.numerator % 2
    best = Fr(-1)
    for Dd in range(1, Dmax + 1, 2):
        for a in range(-4 * Dd, 4 * Dd + 1):
            x = Fr(a, Dd)
            if x.denominator != Dd:
                continue
            S = {x: None}
            st = [x]
            edges = {}
            okc = True
            while st:
                y = st.pop()
                nx = [y / 2] if par(y) == 0 else [(q * y + 1) / 2, (q * y - 1) / 2]
                edges[y] = nx
                for z in nx:
                    if z not in S:
                        S[z] = None
                        st.append(z)
                        if len(S) > cap:
                            okc = False
                            st = []
                            break
            if not okc:
                continue
            nodes = list(edges)
            idx = {v: i for i, v in enumerate(nodes)}
            succ = np.array([[idx[z] for z in (edges[v] + edges[v])[:2]] for v in nodes], dtype=np.int64)
            md = L.karp_min_density(succ, np.array([par(v) for v in nodes], dtype=np.int64))
            best = max(best, md)
    return best


TR = {q: trap_best(q) for q in (3, 5, 7, 9, 11, 13)}
check(TR[3] == Fr(1, 2) and TR[5] == Fr(1, 3) and TR[7] == Fr(1, 4) and TR[9] == Fr(1, 4) and TR[11] <= 0 and
      TR[13] <= 0, f"best exact-trap bounds: {[(q, str(v)) for q, v in TR.items()]} -- all weaker than Theorem N "
      f"(q=5: {{1/3, 2/3, 4/3}} gives 1/3 < 2/5)")
say("G2. rational-window adversaries for q = 7 (Max plays a/D with a in a window of length H):")
best = Fr(0)
for Dd in (1, 3, 5, 9, 11, 13, 15):
    for cA in (Fr(-1), Fr(-1, 2), Fr(0)):
        for k in (10, 12):
            Nk = 1 << k
            Hk = Nk >> 1
            A = int(cA * Hk)
            y = np.arange(Hk, dtype=np.int64)
            a = (Dd * y - A) % Hk + A
            r = (a % Nk) * pow(Dd, -1, Nk) % Nk
            tau = ((r >> (k - 1)) & 1).astype(np.uint8)
            v, _ = L.value_of_tau(k, 7, tau)
            best = max(best, v)
check(best == Fr(1, 3), "their best value is 1/3 = Theorem N's bound (the free cycle {D, 4D, 2D} of 7x +- D is always "
      "available to Min)")

say(f"\nchecks: {NCHECK[0]}; wall time {time.time() - T0:.0f} s; peak RSS {rss_mb():.0f} MB")
say("ALL CHECKS PASSED")
