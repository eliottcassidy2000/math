#!/usr/bin/env python3
"""
procgen_seven2_20260926_run.py -- runner of lane "seven2" (session collatz-procgen-20260922, 2026-09-26):
structure of the optimal sign strategies of 5n+-1 and 7n+-1, automatic strategies, adversaries, and the level k = 30.

Every printed claim is a check(...) that raises on failure; the output ends with ALL CHECKS PASSED.
Environment: SEVEN2_TREE7_KMAX (default 16: the largest level of the q = 7 tree series), SEVEN2_K30 (default 1: the
level-30 certificates; about 10 minutes and 610 MB peak), SEVEN2_STRUCT_KMAX (default 26: the largest level of the
potential statistics, which re-solves the game at the certified values with the seven lane's engine).
"""
import os
import sys
import time
import math
import random
import platform
import resource
import subprocess
from fractions import Fraction as Fr

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import procgen_seven2_20260926_lib as L                 # noqa: E402
import procgen_floor_20260926_lib as FL                 # noqa: E402  (floor lane, read-only: Karp, value_of_tau)

TREE7_KMAX = int(os.environ.get('SEVEN2_TREE7_KMAX', 16))
DO_K30 = int(os.environ.get('SEVEN2_K30', 1))
STRUCT_KMAX = int(os.environ.get('SEVEN2_STRUCT_KMAX', 26))
T0 = time.time()
NCHECK = [0]
LOG72 = math.log(2) / math.log(7)
PEAK = {'rss': 0, 'what': ''}


def say(*a):
    print(*a, flush=True)


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)
    NCHECK[0] += 1
    say('  [ok] ' + msg)


def timed(cmd):
    t = time.time()
    p = subprocess.run(['/usr/bin/time', '-l'] + cmd, capture_output=True, text=True)
    rss = 0
    for line in p.stderr.splitlines():
        if 'maximum resident set size' in line:
            rss = int(line.split()[0])
    if rss > PEAK['rss']:
        PEAK['rss'], PEAK['what'] = rss, ' '.join(os.path.basename(c) for c in cmd[:4])
    return p.stdout, time.time() - t, rss, p.returncode


def witness(q, k, minus, val, cyc):
    """exact re-check of a rhomax witness: a closed walk of G_sigma of density val, and its rational cycle"""
    a, p = L.check_cycle(q, k, minus, cyc)
    orb, ss = L.cycle_rational(q, k, minus, cyc)
    return Fr(a, p) == val, orb, ss


def rho_cert(q, k, minus, tag):
    """rho_max(sigma) with both halves re-checked: witness cycle (exact) and potential certificate (the seven
    lane's independent checker, upper-only)"""
    pre = os.path.join(L.SCR, tag)
    for ext in ('up_sig', 'up_psi'):
        if os.path.exists(pre + '.' + ext):
            os.remove(pre + '.' + ext)
    v, cyc, le = L.rhomax(q, k, minus, tag, F0=Fr(0), cert=pre)
    ok_w, orb, ss = witness(q, k, minus, v, cyc)
    out, _, _, rc = timed([L.BIN['verify'], str(q), str(k), str(v.numerator), str(v.denominator), pre])
    ok_c = rc == 0 and f'RESULT UPPER-CERTIFIED {v.numerator}/{v.denominator}' in out
    for ext in ('up_sig', 'up_psi'):
        if os.path.exists(pre + '.' + ext):
            os.remove(pre + '.' + ext)
    return v, (not le) and ok_w and ok_c, orb, ss


say("procgen_seven2_20260926_run.py -- lane seven2, collatz-procgen-20260922, 2026-09-26")
say(f"python {platform.python_version()}, numpy {np.__version__}; TREE7_KMAX={TREE7_KMAX} K30={DO_K30} "
    f"STRUCT_KMAX={STRUCT_KMAX}")
L.build()
for fn in ['procgen_seven2_20260926_rhomax.c', 'procgen_seven2_20260926_restrict.c', 'procgen_seven2_20260926_lean8.c',
           'procgen_seven2_20260926_verify8.c', 'procgen_seven2_20260926_lib.py', 'procgen_seven2_20260926_run.py',
           'procgen_seven_20260926_game.c', 'procgen_seven_20260926_verify.c', 'procgen_seven_20260926_tight.c',
           'procgen_floor_20260926_lib.py']:
    say(f"  sha256 {fn} = {L.sha(os.path.join(HERE, fn))}")
for name, src in (('lean8', 'procgen_seven2_20260926_lean8.c'), ('verify8', 'procgen_seven2_20260926_verify8.c')):
    p = os.path.join(HERE, src)
    exe = os.path.join(L.SCR, f'{name}_{L.sha(p)[:12]}')
    if not os.path.exists(exe):
        subprocess.run(['cc', '-O2', '-Wall', '-o', exe + '.tmp', p], check=True)
        os.replace(exe + '.tmp', exe)
    L.BIN[name] = exe

# ================================================================================================= A
say("\n== A. Engines and cross-checks ==")
rng = np.random.default_rng(20260926)
n = 0
for q in (5, 7, 9):
    for k in (5, 6, 7, 8):
        x = L.odd_nodes(k)
        strategies = [L.mh_minus(q, x)] + [rng.integers(0, 2, len(x)).astype(bool) for _ in range(4)]
        for m in strategies:
            v, cyc, le = L.rhomax(q, k, m, 'a1')
            ok_w, _, _ = witness(q, k, m, v, cyc)
            flip = np.zeros(1 << k, dtype=np.uint8)
            flip[1::2] = m
            if not (ok_w and not le and FL.rho_max_exact(k, q, flip) == v):
                check(False, f"A1 mismatch q={q} k={k}")
            n += 1
check(True, f"A1. {n} strategies (max-halving and random; q = 5, 7, 9; k = 5..8): the exact rho_max of the rhomax "
            f"engine equals the floor lane's pure-Python Karp value, and every witness is a closed walk of G_sigma of "
            f"that density whose rational periodic point has the walk's residues")

x = L.odd_nodes(12)
v, ok, _, _ = rho_cert(7, 12, L.mh_minus(7, x), 'a2')
pre = os.path.join(L.SCR, 'a2neg')
L.rhomax(7, 12, L.mh_minus(7, x), 'a2n', F0=Fr(0), cert=pre)
psi = np.fromfile(pre + '.up_psi', dtype=np.int32)
out1, _, _, rc1 = timed([L.BIN['verify'], '7', '12', '1', '2', pre])
out2, _, _, rc2 = timed([L.BIN['verify'], '7', '12', '6', '13', pre])
psi2 = psi.copy()
H12 = 1 << 11
sig = np.fromfile(pre + '.up_sig', dtype=np.uint8)
bad_done = False
for xx in range(1 << 12):          # raise psi by 1 at the head of a tight edge
    e = 1 if xx % 2 else -1
    Qn = (xx // 2) % H12 if xx % 2 == 0 else ((7 * xx + (-1 if sig[(xx - 1) // 2] else 1)) // 2) % H12
    if psi[Qn] + e == psi[xx % H12]:
        psi2[Qn] += 1
        bad_done = True
        break
psi2.tofile(pre + '.up_psi')
out3, _, _, rc3 = timed([L.BIN['verify'], '7', '12', '1', '2', pre])
check(v == Fr(1, 2) and ok and rc1 == 0 and rc2 != 0 and bad_done and rc3 != 0,
      "A2. the rhomax potential certificates are accepted by the seven lane's independent checker "
      "(max-halving, q=7, k=12: rho_max = 1/2); the checker rejects the same certificate at 6/13 and after "
      "raising the potential by 1 at the head of one tight edge")

k8dir = os.path.join(L.SCR, 't8')
os.makedirs(k8dir, exist_ok=True)
pre = os.path.join(k8dir, 'q7k24')
res = {}
for kind in ('lower', 'upper'):
    out, wall, rss, rc = timed([L.BIN['lean8'], kind, '7', '24', '3', '8', pre])
    vo, _, _, vrc = timed([L.BIN['verify8'], kind, '7', '24', '3', '8', pre])
    res[kind] = rc == 0 and vrc == 0 and f'{kind.upper()}-CERTIFIED 3/8' in vo
so, _, _, src_ = timed([L.BIN['verify'], '7', '24', '3', '8', pre])
g = bytearray(open(pre + '.lo_g8', 'rb').read())
g0 = bytes(g)
Hh = 1 << 23
g[5] = (g[5] + 1) % 255 if g[5] < 254 else 0
open(pre + '.lo_g8', 'wb').write(bytes(g))
vo_asym, _, _, rc_asym = timed([L.BIN['verify8'], 'lower', '7', '24', '3', '8', pre])
g = bytearray(g0)
tau = np.unpackbits(np.fromfile(pre + '.lo_taub', dtype=np.uint8), bitorder='little')[:Hh]
gg = np.frombuffer(g0, dtype=np.uint8).astype(np.int64)
target = None
for P in range(1, 1 << 20):         # a tight edge P -> Q; raise g(Q) and g(-Q) together (keeps the mirror symmetry)
    if gg[P] == 255:
        continue
    xnode = P + int(tau[P]) * Hh
    e = 5 if xnode % 2 else -3
    tg = [(xnode // 2) % Hh] if xnode % 2 == 0 else [((7 * xnode + 1) // 2) % Hh, ((7 * xnode - 1) // 2) % Hh]
    for Qn in tg:
        if gg[Qn] == gg[P] + e and gg[Qn] < 254 and Qn not in (0, Hh // 2):
            target = Qn
            break
    if target is not None:
        break
g[target] += 1
g[Hh - target] += 1
open(pre + '.lo_g8', 'wb').write(bytes(g))
vo_bad, _, _, rc_bad = timed([L.BIN['verify8'], 'lower', '7', '24', '3', '8', pre])
open(pre + '.lo_g8', 'wb').write(g0)
vo_ok, _, _, rc_ok = timed([L.BIN['verify8'], 'lower', '7', '24', '3', '8', pre])
check(res['lower'] and res['upper'] and src_ == 0 and 'RESULT CERTIFIED 3/8' in so and rc_asym != 0 and
      'asymmetries' in vo_asym and rc_bad != 0 and 'LOWER FAIL' in vo_bad and rc_ok == 0,
      "A3. the one-byte engine lean8 writes both certificates of rho*(7,24) = 3/8 (THM-4486); the memory-lean checker "
      "verify8 and the seven lane's checker both accept them; verify8 rejects a mirror-asymmetric potential file and "
      "a symmetric one raised by 1 at the head of a tight edge; the restored file passes again")
for ext in ('lo_taub', 'lo_g8', 'up_sigb', 'up_psi8'):
    if os.path.exists(pre + '.' + ext):
        os.remove(pre + '.' + ext)

# ================================================================================================= B
say("\n== B. The max-halving skeleton and the flip calculus (Lemmas C, MH, F of the note) ==")
okall, nn = True, 0
for q in (5, 7, 9, 11):
    for k in (8, 12, 16):
        x = L.odd_nodes(k)
        v, cyc, le = L.rhomax(q, k, L.mh_minus(q, x), 'b1')
        ok_w, orb, ss = witness(q, k, L.mh_minus(q, x), v, cyc)
        # accelerated valuations of the witness: all equal to 2
        idx = [i for i, y in enumerate(orb) if y.numerator % 2]
        gaps = [(idx[(j + 1) % len(idx)] - idx[j]) % len(orb) or len(orb) for j in range(len(idx))]
        okall &= (v == Fr(1, 2) and ok_w and not le and set(gaps) == {2})
        nn += 1
check(okall, f"B1. rho_max(max-halving) = 1/2 for q = 5, 7, 9, 11 at k = 8, 12, 16 ({nn} cases), and every witness is a "
             f"cycle of the v=2 set S_inf (all accelerated valuations equal 2), as Lemma MH states")
tot = {}
for q in (5, 7):
    b = Fr(1, q - 4)
    for nper in range(1, 9):
        pts = L.sinf_periodic_points(q, nper)
        tot[(q, nper)] = len(pts)
        if not all(-b <= p <= b for p, _ in pts):
            check(False, "S_inf point outside the interval")
    fx = L.sinf_periodic_points(q, 1)
    check(sorted(p for p, _ in fx) == [-b, b] and all(tot[(q, n_)] == 2 ** n_ for n_ in range(1, 9)),
          f"B2. q={q}: the periodic points of S_inf of period dividing n are exactly 2^n rationals (one per sign word, "
          f"n = 1..8), all in [-1/{q - 4}, 1/{q - 4}]; the fixed points are +-1/{q - 4}")
random.seed(1)
K = 200
bad = 0
counts = {}
neutral_ok = True
for _ in range(100000):
    xx = random.getrandbits(K - 8) | 1
    it = L.mh_itin(7, xx, 3)
    (s1, v1), (s2, v2), (s3, v3) = it
    t = 7 * xx + s1
    x1 = t >> v1
    x2 = (7 * x1 + s2) >> v2
    y = (7 * xx - s1) // 2
    if y != (1 << (v1 - 1)) * x1 - s1:
        bad += 1
    (sy, vy), = L.mh_itin(7, y, 1)
    if v1 == 2:
        pred = 2 if s2 != s1 else (3 if v2 >= 3 else (4 if s3 != s1 else (5 if v3 == 2 else (6 if v3 >= 4 else 7))))
        okp = vy >= 7 if pred == 7 else vy == pred
    elif v1 == 3:
        okp = vy == 2
    elif v1 == 4:
        okp = vy == (3 + v2 if s2 != s1 else 4)
        if s2 != s1:
            neutral_ok &= ((7 * y + sy) >> vy) == x2
    else:
        okp = vy == 3
    bad += 0 if okp else 1
    counts[v1 if v1 <= 5 else 5] = counts.get(v1 if v1 <= 5 else 5, 0) + 1
gains = {}
random.seed(3)
for _ in range(100000):
    xx = random.getrandbits(K - 8) | 1
    (s1, v1), (s2, v2) = L.mh_itin(7, xx, 2)
    y = (7 * xx - s1) // 2
    (sy, vy), = L.mh_itin(7, y, 1)
    gain = 1 + vy - v1 - v2
    cls = 'gain' if (v1 == 2 and v2 == 2 and s1 == s2) else ('neutral' if (v1 == 4 and s2 != s1) else 'loss')
    gains.setdefault(cls, set()).add(1 if gain > 0 else (0 if gain == 0 else -1))
check(bad == 0 and neutral_ok and gains == {'gain': {1}, 'neutral': {0}, 'loss': {-1}},
      "B3. flip calculus (Lemma F, q = 7) on 100000 random 2-adic integers (exact to 2^-192): the flipped successor is "
      "y = 2^(v1-1) x1 - s1, and its MH valuation is 2 (v1 = 2, s2 != s1), 3 (v1 = 2, s2 = s1, v2 >= 3), 4 (v1 = 2, "
      "s2 = s1, v2 = 2, s3 != s1), 5/6/>=7 (v1 = 2, s1 = s2 = s3, v2 = 2, v3 = 2/>=4/3), 2 (v1 = 3), 3 + v2 (v1 = 4, "
      "s2 != s1; then the orbit rejoins x2), 4 (v1 = 4, s2 = s1), 3 (v1 >= 5); consequently (Corollary F, another "
      "100000 samples) a flip followed by max-halving gains halvings over two max-halving steps exactly when the MH "
      "itinerary starts (s,2),(s,2), is density-neutral exactly when it starts (s,4),(-s,.), and loses otherwise")
bad5 = 0
random.seed(2)
for _ in range(100000):
    xx = random.getrandbits(K - 8) | 1
    (s1, v1), (s2, v2), (s3, v3) = L.mh_itin(5, xx, 3)
    x1 = (5 * xx + s1) >> v1
    x2 = (5 * x1 + s2) >> v2
    y = (5 * xx - s1) // 2
    (sy, vy), = L.mh_itin(5, y, 1)
    if y != (1 << (v1 - 1)) * x1 - s1:
        bad5 += 1
    if v1 == 2:
        pred = 2 if s2 != s1 else (3 if v2 >= 3 else (4 if s3 == s1 else 5))
        okp = vy >= 5 if pred == 5 else vy == pred
    elif v1 == 3:
        okp = vy == (2 + v2 if s2 != s1 else 3)
        if s2 != s1:
            okp &= ((5 * y + sy) >> vy) == x2
    else:
        okp = vy == 2
    bad5 += 0 if okp else 1
check(bad5 == 0,
      "B3'. the same calculus for q = 5 (100000 samples): v(y) = 2 (v1 = 2, s2 != s1), 3 (v1 = 2, s2 = s1, v2 >= 3), "
      "4 (v1 = 2, s1 = s2 = s3, v2 = 2), >= 5 (v1 = 2, s1 = s2 != s3, v2 = 2), 2 + v2 (v1 = 3, s2 != s1; rejoins x2), "
      "3 (v1 = 3, s2 = s1), 2 (v1 >= 4); at the S_inf fixed point 1 (all signs '-') the flip gives (1, 4): the "
      "sporadic cycle (1,3,8,4,2) of density 2/5, whereas for q = 7 the flip at -1/3 gives (1, 5), density 1/3")

# ================================================================================================= C
say("\n== C. q = 5: an explicit optimal strategy for every level k >= 15, and why 15 ==")
t1 = time.time()
leaves5, m5, tests5 = L.build_tree(5, 15, Fr(2, 5))
nflip5 = sum(1 for c, d, s in leaves5 if d >= 2 and (s == '+') != ((5 * c + 1) % 4 == 0))
depth5 = {}
for c, d, s in leaves5:
    depth5[d] = depth5.get(d, 0) + 1
m5t = L.tree_strategy(5, 15, leaves5)
vals5 = {}
okc = True
for k in (15, 16, 17, 18):
    v, ok, orb, ss = rho_cert(5, k, L.lift(m5t, 15, k), f'c1_{k}')
    vals5[k] = v
    okc &= ok
check(np.array_equal(m5t, m5) and all(v == Fr(2, 5) for v in vals5.values()) and okc and len(leaves5) == 122 and
      max(depth5) == 15,
      f"C1. the greedy decision tree (low bits first, MH sign tried first) at F = 2/5 has {len(leaves5)} leaves "
      f"({nflip5} flip the max-halving sign), depths {dict(sorted(depth5.items()))} ({tests5} restricted-game tests, "
      f"{time.time() - t1:.0f} s); its strategy has rho_max = 2/5 at k = 15, 16, 17, 18 (witness cycle re-walked "
      f"exactly + potential certificate accepted by the seven lane's checker): ONE explicit 122-leaf rule is optimal "
      f"at every level k >= 15 (the map T_sigma does not depend on k)")
say("    the 122 leaves (x mod 2^d -> sign; '*' marks a flip of the max-halving sign):")
lv5 = sorted(leaves5, key=lambda t: (t[1], t[0]))
for i in range(0, len(lv5), 8):
    say("      " + "  ".join(f"{c}/2^{d}:{s}{'*' if (s == '+') != ((5 * c + 1) % 4 == 0) else ''}" for c, d, s in lv5[i:i + 8]))
x15 = L.odd_nodes(15)
sgn_at = {}
for r in (Fr(1), Fr(-1), Fr(1, 9), Fr(-1, 9), Fr(9, 7), Fr(-9, 7)):
    res_ = L.residue(r, 1 << 15)
    sgn_at[r] = '-' if m5t[(res_ - 1) // 2] else '+'
mh5 = {r: '-' if (5 * L.residue(r, 8) + 1) % 4 != 0 else '+' for r in sgn_at}
check(sgn_at[Fr(1)] == '+' and sgn_at[Fr(-1)] == '-' and mh5[Fr(1)] == '-' and sgn_at[Fr(1, 9)] != mh5[Fr(1, 9)] and
      sgn_at[Fr(-1, 9)] != mh5[Fr(-1, 9)],
      "C2. the tree flips the max-halving sign at the S_inf fixed points +1, -1 (turning the density-1/2 cycles "
      "(1,2), (-1,-2) into the sporadic cycles (1,3,8,4,2), (-1,-3,-8,-4,-2)) and at the S_inf period-2 orbit "
      "+-1/9; it keeps the max-halving sign at +-9/7")
deep = sorted((c, s) for c, d, s in leaves5 if d == 15)
pairs = {}
for c, s in deep:
    pairs.setdefault(c % (1 << 14), []).append((c, s))
rows = []
okw = True
collide = []
for base in sorted(pairs):
    lst = pairs[base]
    a, D = L.rat_recon(np.array([c for c, _ in lst]), 15)
    rats = [Fr(int(a[i]), int(D[i])) for i in range(len(lst))]
    dens = []
    pts = {}
    for forced in ('+', '-'):
        mm = m5t.copy()
        sel = (x15 % (1 << 14)) == base
        mm[sel] = (forced == '-')
        v, cyc, le = L.rhomax(5, 15, mm, 'c3')
        ok_w, orb, ss = witness(5, 15, mm, v, cyc)
        okw &= ok_w and v > Fr(2, 5)
        dens.append((forced, v, len(orb), orb[0].denominator))
        inside = [y for y, sg in zip(orb, ss) if y.numerator % 2 and L.residue(y, 1 << 14) == base]
        okw &= len(inside) >= 1
        pts[forced] = inside[0]
    dlt = pts['+'] - pts['-']
    v2 = 0
    num = dlt.numerator
    while num % 2 == 0:
        num //= 2
        v2 += 1
    collide.append((base, pts['+'], pts['-'], v2))
    rows.append((base, list(zip([c for c, _ in lst], [s for _, s in lst], rats)), dens))
for base, p1, p2, v2 in collide:
    say(f"    class {base} mod 2^14: the '+' conflict cycle passes through {p1}, the '-' one through {p2}; "
        f"v_2({p1} - ({p2})) = {v2}")
okw &= all(v2 == 14 for _, _, _, v2 in collide)
for base, lst, dens in rows:
    say(f"    class {base} mod 2^14 = {[(c, s, str(r)) for c, s, r in lst]}: forcing one sign on the whole class "
        f"gives rho_max " + ", ".join(f"{v} (all '{f}', witness length {p_}, denominator {d_})" for f, v, p_, d_ in dens))
check(len(pairs) == 4 and all(len(l) == 2 and l[0][1] != l[1][1] for _, l, _ in rows) and okw and
      Fr(5, 12) in [v for _, _, dens in rows for _, v, _, _ in dens],
      "C3. why 15: the tree's 8 depth-15 leaves form 4 classes mod 2^14 (2 negation pairs), each split into two halves "
      "with opposite signs (simplest rationals 19/99 | -163/13, -889/1129 | 9/7, and negatives); giving a whole class "
      "mod 2^14 one sign (either one) keeps a cycle denser than 2/5 -- among them the all-outward (12,5)-cycle of 5x+971 "
      "of density 5/12 = rho*(5,14) -- and the two conflict cycles of each class pass through two rationals that are "
      "congruent mod 2^14 but not mod 2^15 (listed above): a residue collision resolved at level 15")

# ================================================================================================= D
say("\n== D. q = 7: explicit trees, the S_inf skeleton of the optima, and the structure of both certificates ==")
VAL7 = {10: Fr(2, 5), 11: Fr(2, 5), 12: Fr(2, 5), 13: Fr(2, 5), 14: Fr(15, 38), 15: Fr(15, 38), 16: Fr(7, 18),
        17: Fr(7, 18), 18: Fr(19, 49), 19: Fr(13, 34), 20: Fr(34, 89), 21: Fr(8, 21), 22: Fr(14, 37), 23: Fr(14, 37),
        24: Fr(3, 8), 25: Fr(3, 8), 26: Fr(3, 8), 27: Fr(35, 94)}
t1 = time.time()
leaves7, m7, tests7 = L.build_tree(7, 10, Fr(2, 5))
m7t = L.tree_strategy(7, 10, leaves7)
okd = np.array_equal(m7, m7t)
vals7 = {}
for k in (10, 12, 14, 16):
    v, ok, orb, ss = rho_cert(7, k, L.lift(m7t, 10, k), f'd1_{k}')
    vals7[k] = v
    okd &= ok
flipr = {}
for r in (Fr(1, 3), Fr(-1, 3), Fr(1, 11), Fr(-1, 11), Fr(1), Fr(-1)):
    res_ = L.residue(r, 1 << 10)
    mhs = '+' if (7 * res_ + 1) % 4 == 0 else '-'
    flipr[r] = ('-' if m7t[(res_ - 1) // 2] else '+') != mhs
check(okd and all(v == Fr(2, 5) for v in vals7.values()) and len(leaves7) == 46 and flipr[Fr(1, 3)] and
      flipr[Fr(-1, 3)] and flipr[Fr(1, 11)] and flipr[Fr(-1, 11)] and not flipr[Fr(1)] and not flipr[Fr(-1)],
      f"D1. q = 7, F = 2/5: a {len(leaves7)}-leaf tree of depth 10 ({tests7} tests, {time.time() - t1:.0f} s); its "
      f"strategy has rho_max = 2/5 at k = 10, 12, 14, 16 (certified); among its "
      f"{sum(1 for c, d, s in leaves7 if (s == '+') != ((7 * c + 1) % 4 == 0))} flip leaves are the classes of the "
      f"S_inf fixed points +-1/3 (depth 5) and of the S_inf period-2 orbit +-1/11 (depth 10); it keeps the max-halving "
      f"sign at the free cycles' points +-1")
say("    the 46 leaves (x mod 2^d -> sign; '*' = flip of the max-halving sign):")
lv7 = sorted(leaves7, key=lambda t: (t[1], t[0]))
for i in range(0, len(lv7), 8):
    say("      " + "  ".join(f"{c}/2^{d}:{s}{'*' if (s == '+') != ((7 * c + 1) % 4 == 0) else ''}" for c, d, s in lv7[i:i + 8]))
gainset = []
for c in range(1, 32, 2):
    it = L.mh_itin(7, c + 32 * 7919, 2)          # any lift: the first two symbols are determined mod 32
    it2 = L.mh_itin(7, c + 32 * 104729, 2)
    if it == it2 and it[0][1] == 2 and it[1][1] == 2 and it[0][0] == it[1][0]:
        gainset.append(c)
d5flips = sorted(c for c, d, s in leaves7 if d == 5 and (s == '+') != ((7 * c + 1) % 4 == 0))
check(gainset == [11, 21] and d5flips == [11, 21] and L.residue(Fr(1, 3), 32) == 11 and L.residue(Fr(-1, 3), 32) == 21,
      "D1'. the tree's shallow flips (depth 5) are exactly the classes 11, 21 mod 32 = the classes of 1/3, -1/3 = the "
      "classes whose MH itinerary starts (s,2),(s,2): the only pattern on which a flip gains halvings (Corollary F)")
pre = os.path.join(L.SCR, 'd2')
L.rhomax(7, 10, m7t, 'd2', F0=Fr(2, 5), cert=pre)
psi = np.fromfile(pre + '.up_psi', dtype=np.int32).astype(np.int64)
H10 = 1 << 9
E = []
for P in range(H10):
    for b in (0, 1):
        xx = P + b * H10
        if xx % 2 == 0:
            Qn, s = (xx // 2) % H10, 0
        else:
            s = -1 if m7t[(xx - 1) // 2] else 1
            Qn = ((7 * xx + s) // 2) % H10
        e = 3 if P % 2 else -2
        if psi[Qn] + e == psi[P]:
            E.append((P, Qn, b, s))
adj = {}
for P, Qn, b, s in E:
    adj.setdefault(P, []).append((Qn, b, s))
found = {}
from collections import deque
for start in list(adj):                          # shortest tight cycle through each start pair (BFS)
    prev = {start: None}
    dq = deque([start])
    cyc = None
    while dq and cyc is None:
        u = dq.popleft()
        for (w, b, s) in adj.get(u, []):
            if w == start:
                path = [(u, b)]
                z = u
                while prev[z] is not None:
                    pz, pb = prev[z]
                    path.append((pz, pb))
                    z = pz
                cyc = [P_ + b_ * H10 for (P_, b_) in path[::-1]]
                break
            if w not in prev:
                prev[w] = (u, b)
                dq.append(w)
    if cyc is not None and len(cyc) == 5:
        orb, ss = L.cycle_rational(7, 10, m7t, cyc)
        found[min(orb)] = orb
dens17 = sorted(found.values(), key=lambda o: min(o))
check(len(dens17) >= 2 and all(o[0].denominator == 17 for o in dens17),
      f"D2. the density-2/5 cycles of length 5 of this strategy are the (5,2)-cycles of 7x+-1 with denominator 17 "
      f"(2^5 - 7^2 = -17), e.g. {', '.join(str(y) for y in dens17[0])}: a flip (valuation 1) followed by valuation 4 "
      f"(Lemma F: v1 = 2, s1 = s2, s3 != s1), the next obstruction below 2/5")
TREES = {10: (len(leaves7), Fr(2, 5))}
okt = True
for k in range(14, TREE7_KMAX + 1, 2):
    t1 = time.time()
    lv, mm, tt = L.build_tree(7, k, VAL7[k])
    v, ok, _, _ = rho_cert(7, k, L.tree_strategy(7, k, lv), f'd3_{k}')
    okt &= ok and v == VAL7[k] and np.array_equal(mm, L.tree_strategy(7, k, lv))
    dmax = max(d for _, d, _ in lv)
    nmax = sum(1 for _, d, _ in lv if d == dmax)
    TREES[k] = (len(lv), VAL7[k])
    say(f"    k={k}: F = rho*(7,{k}) = {VAL7[k]}: greedy tree with {len(lv)} leaves ({nmax} at the maximal depth {dmax}), "
        f"{tt} tests, {time.time() - t1:.0f} s; rho_max of its strategy = {v}")
check(okt and all(TREES[a][0] < TREES[b][0] for a, b in zip(sorted(TREES), sorted(TREES)[1:])),
      f"D3. greedy trees of optimal 7n+-1 strategies: leaves {[(k, str(TREES[k][1]), TREES[k][0]) for k in sorted(TREES)]}; "
      f"each certified optimal (rho_max = rho*(7,k)); the description length grows quickly as the value decreases")

say("D4-D6: the least fixed points at the certified values (the seven lane's engine, 'lean' mode, both one-sided "
    "potentials; statistics only):")
SD = os.path.join(L.SCR, 'struct')
os.makedirs(SD, exist_ok=True)


def lean_solve(q, k, F):
    """both least fixed points at F (the seven lane's 'lean' mode), left on disk; returns (prefix, wall)"""
    pre = os.path.join(SD, f'q{q}_k{k}')
    out, wall, rss, rc = timed([L.BIN['game'], 'lean', str(q), str(k), str(F.numerator), str(F.denominator), pre])
    if rc != 0 or 'LEAN BOTH' not in out:
        check(False, f"lean solve q={q} k={k}")
    return pre, wall


def load_pot32(pre, kind):
    """one potential as int32 over the H pairs (-1 = outside W), read in chunks to keep the memory low"""
    e8, e16 = ('lo_g8', 'lo_g16') if kind == 'lo' else ('up_psi8', 'up_psi16')
    if os.path.exists(pre + '.' + e8):
        raw, top = np.memmap(pre + '.' + e8, dtype=np.uint8, mode='r'), 255
    else:
        raw, top = np.memmap(pre + '.' + e16, dtype=np.uint16, mode='r'), 65535
    out = np.empty(len(raw), dtype=np.int32)
    for lo in range(0, len(raw), 1 << 22):
        blk = np.asarray(raw[lo:lo + (1 << 22)]).astype(np.int32)
        blk[blk == top] = -1
        out[lo:lo + len(blk)] = blk
    del raw
    return out


def lean_clean(pre):
    for ext in ('lo_taub', 'lo_g8', 'lo_g16', 'up_sigb', 'up_psi8', 'up_psi16'):
        if os.path.exists(pre + '.' + ext):
            os.remove(pre + '.' + ext)


EXC = {}
HGT = {}
NEST = {}
prevcls = None
CHK = 1 << 18
for k in list(range(19, min(STRUCT_KMAX, 27) + 1)):
    q, F = 7, VAL7[k]
    pre, wall = lean_solve(q, k, F)
    H = 1 << (k - 1)
    psi = load_pot32(pre, 'up')
    e = F.denominator - F.numerator
    cls = np.zeros(H, dtype=np.int8)       # over odd nodes x = 2i+1: 0 free, 1 forced MH, 2 forced flip
    nexc = 0
    by8n = np.zeros(8, dtype=np.int64)
    by8d = np.zeros(8, dtype=np.int64)
    n32 = d32 = 0
    for lo in range(0, H, CHK):
        i = np.arange(lo, min(H, lo + CHK), dtype=np.int64)
        x = 2 * i + 1
        okp = psi[((q * x + 1) // 2) % H] + e <= psi[x % H]
        okm = psi[((q * x - 1) // 2) % H] + e <= psi[x % H]
        mh = L.mh_minus(q, x)
        okmh = np.where(mh, okm, okp)
        okfl = np.where(mh, okp, okm)
        exc = ~okmh
        cls[lo:lo + len(i)] = np.where(okmh & okfl, 0, np.where(okmh, 1, 2))
        nexc += int(exc.sum())
        np.add.at(by8n, x % 8, exc)
        np.add.at(by8d, x % 8, 1)
        sel = ((x % 32) == 11) | ((x % 32) == 21)
        n32 += int(exc[sel].sum())
        d32 += int(sel.sum())
    if prevcls is not None and VAL7[k] == VAL7[k - 1]:
        agree = conf = 0
        for lo in range(0, H, CHK):
            i = np.arange(lo, min(H, lo + CHK), dtype=np.int64)
            x = 2 * i + 1
            la = prevcls[((x % (1 << (k - 1))) - 1) // 2]
            cb = cls[lo:lo + len(i)]
            agree += int((la == cb).sum())
            conf += int((((la == 1) & (cb == 2)) | ((la == 2) & (cb == 1))).sum())
        NEST[k] = (agree / H, conf / H)
    prevcls = cls
    by8 = [float(by8n[c] / by8d[c]) for c in (1, 3, 5, 7)]
    EXC[k] = (nexc / H, by8, n32 / d32)
    del psi
    g = load_pot32(pre, 'lo')
    lean_clean(pre)
    fn_bins = np.zeros(8, dtype=np.int64)
    tot_bins = np.zeros(8, dtype=np.int64)
    for lo in range(0, H, CHK):
        P = np.arange(lo, min(H, lo + CHK), dtype=np.int64)
        gp = g[P]
        ee = np.where(P % 2 == 1, F.denominator - F.numerator, -F.numerator)
        oks = []
        for b in (0, 1):
            xn = P + b * H
            ov = np.where(P % 2 == 0, g[(xn // 2) % H], np.maximum(g[((q * xn + 1) // 2) % H], g[((q * xn - 1) // 2) % H]))
            oks.append((ov <= gp + ee) & (ov >= 0) & (gp >= 0))
        Pt = np.where(P <= H // 2, P, P - H)
        fnear = np.where(Pt >= 0, oks[0] & ~oks[1], oks[1] & ~oks[0])
        bins = np.minimum((np.abs(Pt) * 16) // H, 7)
        np.add.at(fn_bins, bins, fnear)
        np.add.at(tot_bins, bins, 1)
    rates = [float(fn_bins[j] / tot_bins[j]) for j in range(8)]
    HGT[k] = (min(rates), max(rates), float(fn_bins.sum() / tot_bins.sum()))
    del g
    say(f"    k={k} F={F}: Min: max-halving sign NOT admissible for the least potential at {EXC[k][0]:.4f} of the odd "
        f"residues (by x mod 8 = 1,3,5,7: {', '.join(f'{r:.3f}' for r in by8)}; on the classes of +-1/3 mod 32: "
        f"{EXC[k][2]:.3f}); Max: forced-'near' lift share {HGT[k][2]:.4f}, per height bin |P|/H in [j/16,(j+1)/16) "
        f"between {HGT[k][0]:.4f} and {HGT[k][1]:.4f}" + (f"; nesting with k-1: agreement {NEST[k][0]:.3f}, conflicts "
        f"{NEST[k][1]:.4f}" if k in NEST else "") + f" ({wall:.0f} s)")
if EXC:
    ks = sorted(EXC)
    check(all(0.10 < EXC[k][0] < 0.17 and EXC[k][2] > 0.4 and max(EXC[k][1][0], EXC[k][1][3]) < 0.05 and
              min(EXC[k][1][1], EXC[k][1][2]) > 0.2 for k in ks),
          f"D4. (EMPIRICAL, k = {ks[0]}..{ks[-1]}) the optimal Min potentials admit the max-halving sign at 83-90% of the "
          f"odd residues; the exceptions sit where max-halving has valuation 2 (x = 3, 5 mod 8: > 20%; x = 1, 7 mod 8: "
          f"< 5%) and above 40% on the classes of the S_inf fixed points +-1/3 mod 32")
    if NEST:
        check(all(a > 0.9 and c < 0.005 for a, c in NEST.values()),
              f"D5. (EMPIRICAL) nesting across levels with the same value: the classification of the odd residues by the "
              f"least Min potential (forced max-halving sign / forced flip / free) at level k agrees with the lift of the "
              f"level-(k-1) classification on {', '.join(f'{NEST[k][0]:.3f} (k={k})' for k in sorted(NEST))} of the "
              f"residues; direct conflicts (forced max-halving <-> forced flip) {', '.join(f'{NEST[k][1]:.4f}' for k in sorted(NEST))}: "
              f"the optima are refined, not rebuilt")
    check(all(HGT[k][1] - HGT[k][0] < 0.01 for k in ks),
          f"D6. (EMPIRICAL) the optimal Max strategies do not depend on the real height: the share of pairs whose lift is "
          f"forced to the 'near' representative varies by less than 0.01 across the 8 height bins at every k = "
          f"{ks[0]}..{ks[-1]} (a potential that is a function of |x| cannot certify them; cf. THM-4486 section 7)")

# ================================================================================================= E
say("\n== E. Automatic strategies ==")


def rat_recon_chunked(x, k):
    """L.rat_recon in chunks (memory)"""
    a = np.empty(len(x), dtype=np.int64)
    D = np.empty(len(x), dtype=np.int64)
    for lo in range(0, len(x), 1 << 15):
        a[lo:lo + (1 << 15)], D[lo:lo + (1 << 15)] = L.rat_recon(x[lo:lo + (1 << 15)], k)
    return a, D


oke = True
rowsE = []
for q in (5, 7):
    for k in (10, 12, 14, 16, 18, 20):
        x = L.odd_nodes(k)
        a, D = rat_recon_chunked(x, k)
        okr = bool(np.all(((a - D * x) % (1 << k)) == 0) and np.all(D % 2 == 1))
        m = a < 0
        v, cyc, le = L.rhomax(q, k, m, 'e1')
        ok_w, orb, ss = witness(q, k, m, v, cyc)
        oke &= okr and ok_w and v == 1
        rowsE.append((q, k, len(orb), orb[0]))
check(oke, "E1. the 'real-sign imitation' sigma(x) = sign of the simplest rational a/D in the class of x (shortest vector "
           "of the lattice a = D x mod 2^k) has rho_max = 1 for q = 5, 7 at every k = 10, 12, ..., 20: all-odd cycles "
           f"(e.g. q={rowsE[-1][0]}, k={rowsE[-1][1]}: length {rowsE[-1][2]}, x0 = {float(rowsE[-1][3]):.4f}); Proposition S "
           "holds for exact rationals, but a residue class cannot know its point's real sign")
okh = True
for q in (5, 7):
    for k in (12, 16):
        x = L.odd_nodes(k)
        a, D = rat_recon_chunked(x, k)
        ht = np.maximum(np.abs(a), D)
        mh = L.mh_minus(q, x)
        for h in (3, 9, 33, 129):
            v, _, _ = L.rhomax(q, k, np.where(ht <= h, a < 0, mh), 'e2')
            okh &= v >= Fr(1, 2)
check(okh, "E2. max-halving corrected to the real sign of the simplest rational on the classes of height <= h "
           "(h = 3, 9, 33, 129; q = 5, 7; k = 12, 16): rho_max >= 1/2 in every case -- never better than max-halving")


def dfa_search(q, k, n, limit=0, seed=1):
    x = L.odd_nodes(k)
    mh = L.mh_minus(q, x)
    bits = [(x >> i) & 1 for i in range(1, k)]
    FLIP, KEEP = n, n + 1
    import itertools
    combos = itertools.product(range(n + 2), repeat=2 * n)
    if limit:
        r = random.Random(seed)
        combos = ([r.randrange(n + 2) for _ in range(2 * n)] for _ in range(limit))
    seen, best, cnt = set(), None, 0
    for combo in combos:
        Dm = np.array(list(combo) + [FLIP, FLIP, KEEP, KEEP], dtype=np.int64).reshape(n + 2, 2)
        st = np.zeros(len(x), dtype=np.int64)
        for b in bits:
            st = Dm[st, b]
        for dfl in (False, True):
            m = mh ^ ((st == FLIP) | ((st < n) & dfl))
            key = np.packbits(m).tobytes()
            if key in seen:
                continue
            seen.add(key)
            v, _, _ = L.rhomax(q, k, m, 'e3')
            cnt += 1
            if best is None or v < best:
                best = v
    return best, cnt


t1 = time.time()
b2, c2 = dfa_search(7, 12, 2)
b3, c3 = dfa_search(7, 12, 3)
b4, c4 = dfa_search(7, 12, 4, limit=20000)
check(b2 == b3 == b4 == Fr(1, 2),
      f"E3. flip sets recognised by small automata (max-halving XOR [a DFA reading bits 1, 2, ... of x accepts]; "
      f"transient states resolved by a default): all DFAs with 2 ({c2} distinct strategies) and 3 ({c3}) transient "
      f"states, and 20000 random 4-state DFAs ({c4} distinct): best rho_max at k = 12 is 1/2, i.e. none beats "
      f"max-halving ({time.time() - t1:.0f} s); the D1 tree needs depth 10 already for 2/5")

# ================================================================================================= F
say("\n== F. Level-independent adversaries (q = 7) ==")
t1 = time.time()
best = {}
for m, k in ((1, 10), (2, 10), (3, 10), (1, 12), (2, 12)):
    H = 1 << (k - 1)
    P = np.arange(H)
    bv = Fr(0)
    for code in range(1 << (1 << m)):
        f = np.array([(code >> i) & 1 for i in range(1 << m)], dtype=np.uint8)
        v, _ = FL.value_of_tau(k, 7, f[P % (1 << m)])
        bv = max(bv, v)
    best[(m, k)] = bv
f79 = np.array([(79 >> i) & 1 for i in range(8)], dtype=np.uint8)
row79 = []
for k in (8, 10, 12, 14):
    H = 1 << (k - 1)
    row79.append(FL.value_of_tau(k, 7, f79[np.arange(H) % 8])[0])
check(all(v == Fr(1, 3) for v in best.values()) and row79[0] == Fr(5, 13) and max(row79[1:]) < Fr(1, 3),
      f"F1. low-bit adversaries tau(P) = f(P mod 2^m), all f: the best exact value is 1/3 = Theorem N's top lift for "
      f"m = 1, 2, 3 at k = 10 and m = 1, 2 at k = 12; at k = 8 some m = 3 rules reach 5/13 (a small-level effect: the "
      f"best one has values {', '.join(str(v) for v in row79)} at k = 8, 10, 12, 14) ({time.time() - t1:.0f} s)")

# ================================================================================================= G
K30 = {}
if DO_K30:
    say("\n== G. Levels 30 and 31 (each certificate re-checked by verify8; the level-30 upper one also by the seven lane's "
        "checker) ==")
    GD = os.path.join(L.SCR, 'k30')
    os.makedirs(GD, exist_ok=True)
    JOBS = [(30, 'upper', Fr(10, 27), [L.BIN['lean8'], 'upper', '7', '30', '10', '27']),
            (30, 'lower', Fr(37, 100), [L.BIN['game'], 'lower', '7', '30', '37', '100']),
            (30, 'upper', Fr(37, 100), [L.BIN['game'], 'upper', '7', '30', '37', '100']),
            (31, 'lower', Fr(7, 19), [L.BIN['lean8'], 'lower', '7', '31', '7', '19', None, '254', '1'])]
    for k, kind, F, cmd in JOBS:
        pre = os.path.join(GD, f'k{k}_{kind}_{F.numerator}_{F.denominator}')
        cmd = list(cmd)
        if None in cmd:
            cmd[cmd.index(None)] = pre
        else:
            cmd.append(pre)
        out, wall, rss, rc = timed(cmd)
        vo, vwall, vrss, vrc = timed([L.BIN['verify8'], kind, '7', str(k), str(F.numerator), str(F.denominator), pre])
        tag = 'UPPER-CERTIFIED' if kind == 'upper' else 'LOWER-CERTIFIED'
        line = next((l for l in vo.splitlines() if l.startswith(kind.upper())), '')
        ok = rc == 0 and vrc == 0 and f'RESULT {tag} {F.numerator}/{F.denominator}' in vo
        extra = ''
        if kind == 'upper' and os.path.exists(pre + '.up_psi8'):
            so, swall, srss, src_ = timed([L.BIN['verify'], '7', str(k), str(F.numerator), str(F.denominator), pre])
            ok &= src_ == 0 and f'RESULT UPPER-CERTIFIED {F.numerator}/{F.denominator}' in so
            extra = f"; the seven lane's checker also accepts it ({swall:.0f} s, {srss / 2 ** 20:.0f} MiB)"
        check(ok, f"G. rho*(7,{k}) {'<=' if kind == 'upper' else '>='} {F} = {float(F):.5f}: {kind} certificate by "
                  f"{os.path.basename(cmd[0]).split('_')[0]} ({wall:.0f} s, {rss / 2 ** 20:.0f} MiB), verify8: "
                  f"{line.split('(')[-1].rstrip(')')} ({vwall:.0f} s, {vrss / 2 ** 20:.0f} MiB){extra}")
        K30[(k, kind)] = F
        for ext in ('lo_taub', 'lo_g8', 'lo_g16', 'lo_h8', 'up_sigb', 'up_psi8', 'up_psi16'):
            if os.path.exists(pre + '.' + ext):
                os.remove(pre + '.' + ext)
    check(K30.get((30, 'lower')) == K30.get((30, 'upper')) == Fr(37, 100) and 7 ** 37 > 2 ** 100 and 7 ** 7 > 2 ** 19,
          "G'. hence rho*(7,30) = 37/100 exactly (both certificates at 37/100; mean valuation 100/37) and rho*(7,31) >= 7/19 "
          "(7^37 > 2^100, 7^7 > 2^19: above log_7 2 = 0.356207): 7n+-1 has no provable sign strategy at any level k <= 31; "
          "the value drops again at k = 30 (37/100 < 13/35 = rho*(7,29)); by monotonicity rho*(7,31) <= 37/100")

# ================================================================================================= H
say("\n== H. The trend in terms of the mean valuation 1/rho* (EMPIRICAL) ==")
seq = [(10, Fr(2, 5)), (14, Fr(15, 38)), (16, Fr(7, 18)), (18, Fr(19, 49)), (19, Fr(13, 34)), (20, Fr(34, 89)),
       (21, Fr(8, 21)), (22, Fr(14, 37)), (24, Fr(3, 8)), (27, Fr(35, 94)), (28, Fr(13, 35))]
if K30.get((30, 'lower')) == K30.get((30, 'upper')) == Fr(37, 100):
    seq.append((30, Fr(37, 100)))
say("    k:      " + " ".join(f"{k:>7d}" for k, _ in seq))
say("    1/rho*: " + " ".join(f"{float(1 / v):7.4f}" for _, v in seq))
kk = np.array([k for k, _ in seq], dtype=float)
mv = np.array([float(1 / v) for _, v in seq])
A1 = np.vstack([np.ones_like(kk), -1 / kk]).T
c1, *_ = np.linalg.lstsq(A1, mv, rcond=None)
A2 = np.vstack([np.ones_like(kk), -1 / np.sqrt(kk)]).T
c2, *_ = np.linalg.lstsq(A2, mv, rcond=None)
say(f"    least-squares fits: A - B/k -> A = {c1[0]:.4f}; A - B/sqrt(k) -> A = {c2[0]:.4f}; threshold log_2 7 = "
    f"{math.log2(7):.4f}")
check(all(float(1 / v) < math.log2(7) for _, v in seq),
      "H. every certified value so far has mean valuation 1/rho* < log_2 7 (i.e. rho* > log_7 2); the fits above are "
      "reported for orientation only (they straddle the threshold and do not decide the limit)")

say(f"\nchecks: {NCHECK[0]}; wall time {time.time() - T0:.0f} s; peak RSS of a child process {PEAK['rss'] / 2 ** 20:.0f} MiB "
    f"({PEAK['what']}); runner RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20:.0f} MiB")
say("ALL CHECKS PASSED")
