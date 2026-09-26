#!/usr/bin/env python3
"""
procgen_seven_20260926_run.py -- runner of lane "seven" (session collatz-procgen-20260922, 2026-09-26):
does 7n +- 1 ever have a bounded-lookahead-provable sign strategy, i.e. is rho*(7,k) < log_7 2 for some k?

Every printed claim is a check(...) that raises on failure; the output ends with ALL CHECKS PASSED.

Engines (compiled into scratch/procgen_seven/, not committed):
  procgen_seven_20260926_game.c    exact min-max solver (pair form of the THM-4486 game, negation quotient,
                                   compact uint16 sweeps); writes both certificates at the value
  procgen_seven_20260926_verify.c  independent checker of the two certificates (all 2^k nodes, exact)
  procgen_seven_20260926_tight.c   tight cores of certificates (structure only)
Python: procgen_seven_20260926_potential.py (exact corrected potentials, q = 9, 11), and the floor lane's
procgen_floor_20260926_lib.py (read-only reuse: its engine, its certificate verifiers, value_of_tau).

Environment: SEVEN_KMAX7 (default 27), SEVEN_KBR7 (1: k = 28, 29, 30), SEVEN_KMAXQ (default 24; q = 9, 11),
SEVEN_KSTRUCT (24).  Peak memory: about 605 MB (the checker of the k = 30 lower certificate).
HINTS: Farey brackets from exploratory runs for the two largest levels of q = 7; they only shorten the search (a wrong
hint makes the search fail, it cannot produce a value without two checked certificates).
"""
import os
import sys
import time
import math
import hashlib
import platform
import resource
import subprocess
from fractions import Fraction as Fr

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
SCR = os.path.join(ROOT, 'scratch', 'procgen_seven', 'run')
sys.path.insert(0, HERE)
import procgen_floor_20260926_lib as FL            # noqa: E402  (floor lane, read-only)
import procgen_seven_20260926_potential as PT      # noqa: E402
import procgen_seven_20260926_cycles as CY         # noqa: E402

KMAX7 = int(os.environ.get('SEVEN_KMAX7', 27))       # exact values of rho*(7,k) up to this level
KBR7 = int(os.environ.get('SEVEN_KBR7', 1))          # 1: add the one-sided certificates at k = 28, 29, 30
KMAXQ = int(os.environ.get('SEVEN_KMAXQ', 24))
KSTRUCT = int(os.environ.get('SEVEN_KSTRUCT', 24))
T0 = time.time()
NCHECK = [0]
LOG2 = math.log(2)


def say(*a):
    print(*a, flush=True)


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)
    NCHECK[0] += 1
    say('  [ok] ' + msg)


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        h.update(f.read())
    return h.hexdigest()


def below_log(q, F):
    """F = a/p < log_q 2  <=>  q^a < 2^p (exact)"""
    return q ** F.numerator < 2 ** F.denominator


# ------------------------------------------------------------------------------------------------ build
os.makedirs(SCR, exist_ok=True)
BIN = {}
for name in ('game', 'verify', 'tight'):
    src = os.path.join(HERE, f'procgen_seven_20260926_{name}.c')
    tag = sha(src)[:12]
    exe = os.path.join(SCR, f'{name}_{tag}')
    if not os.path.exists(exe):
        subprocess.run(['cc', '-O2', '-Wall', '-o', exe + '.tmp', src], check=True)
        os.replace(exe + '.tmp', exe)
    BIN[name] = exe


def timed(cmd):
    """run cmd under /usr/bin/time -l; returns (stdout, wall seconds, max RSS bytes)"""
    t = time.time()
    p = subprocess.run(['/usr/bin/time', '-l'] + cmd, capture_output=True, text=True)
    rss = 0
    for line in p.stderr.splitlines():
        if 'maximum resident set size' in line:
            rss = int(line.split()[0])
    return p.stdout, time.time() - t, rss, p.returncode


PEAK = {'rss': 0, 'what': ''}


def note_rss(rss, what):
    if rss > PEAK['rss']:
        PEAK['rss'], PEAK['what'] = rss, what


CERT_EXT = ('lo_taub', 'lo_g8', 'lo_g16', 'up_sigb', 'up_psi8', 'up_psi16', 'lo_tau', 'lo_g', 'up_sig', 'up_psi')


def clean(prefix):
    """remove certificate files of a prefix (so that no stale file of another format can be read)"""
    for ext in CERT_EXT:
        if os.path.exists(prefix + '.' + ext):
            os.remove(prefix + '.' + ext)


def solve(q, k, lo=(0, 1), hi=(1, 1), ku=(0, 0), keep=False):
    """exact rho*(q,k) by the lean search; both certificates re-checked by the independent checker.
    Returns (value, info)."""
    prefix = os.path.join(SCR, f'c_q{q}_k{k}')
    clean(prefix)
    out, wall, rss, rc = timed([BIN['game'], 'lsearch2', str(q), str(k), str(lo[0]), str(lo[1]), str(hi[0]),
                                str(hi[1]), str(ku[0]), str(ku[1]), prefix, '16'])
    note_rss(rss, f'solver q={q} k={k}')
    val = None
    for line in out.splitlines():
        if line.startswith('VALUE'):
            _, a, b = line.split()
            val = Fr(int(a), int(b))
    check(rc == 0 and val is not None, f"solver found a value for q={q}, k={k}")
    vout, vwall, vrss, vrc = timed([BIN['verify'], str(q), str(k), str(val.numerator), str(val.denominator), prefix])
    note_rss(vrss, f'checker q={q} k={k}')
    ok = vrc == 0 and f"RESULT CERTIFIED {val.numerator}/{val.denominator}" in vout
    info = {'wall': wall, 'rss': rss, 'vwall': vwall, 'vrss': vrss, 'tests': out.count('\n  '),
            'lower': next(l for l in vout.splitlines() if l.startswith('LOWER')),
            'upper': next(l for l in vout.splitlines() if l.startswith('UPPER')),
            'digest': next(l for l in vout.splitlines() if l.startswith('DIGEST')).split()[1], 'prefix': prefix}
    check(ok, f"q={q} k={k}: rho* = {val} certified by the independent checker (lower and upper certificates, all "
              f"{1 << k} nodes); {info['lower'].split('(')[1].split(',')[0]}")
    if not keep:
        for ext in ('lo_taub', 'lo_g8', 'lo_g16', 'up_sigb', 'up_psi8', 'up_psi16'):
            p = prefix + '.' + ext
            if os.path.exists(p):
                os.remove(p)
    return val, info


say("procgen_seven_20260926_run.py -- lane seven, collatz-procgen-20260922, 2026-09-26")
say(f"python {platform.python_version()}, numpy {np.__version__}; KMAX7={KMAX7} KMAXQ={KMAXQ} KSTRUCT={KSTRUCT}")
for fn in ['procgen_seven_20260926_game.c', 'procgen_seven_20260926_verify.c', 'procgen_seven_20260926_tight.c',
           'procgen_seven_20260926_potential.py', 'procgen_seven_20260926_cycles.py', 'procgen_seven_20260926_run.py',
           'procgen_floor_20260926_lib.py']:
    say(f"  sha256 {fn} = {sha(os.path.join(HERE, fn))}")

# ================================================================================================= A
say("\n== A. The lean solver and the independent checker ==")
say("A1. rho*(q,k) from the lean solver equals the floor lane's certified values (THM-4486), every value re-checked by "
    "the independent checker:")
FLOOR = {7: {8: Fr(3, 7), 9: Fr(3, 7), 10: Fr(2, 5), 12: Fr(2, 5), 14: Fr(15, 38), 15: Fr(15, 38), 16: Fr(7, 18),
             17: Fr(7, 18), 18: Fr(19, 49)},
         5: {10: Fr(3, 7), 12: Fr(5, 12), 14: Fr(5, 12), 15: Fr(2, 5), 18: Fr(2, 5)},
         9: {12: Fr(5, 12), 14: Fr(2, 5), 16: Fr(2, 5)},
         11: {12: Fr(9, 22), 14: Fr(13, 33), 16: Fr(21, 55)}}
n = 0
for q, row in FLOOR.items():
    for k, v in row.items():
        prefix = os.path.join(SCR, f'a1_q{q}_k{k}')
        clean(prefix)
        out, wall, rss, rc = timed([BIN['game'], 'lsearch2', str(q), str(k), '0', '1', '1', '1', '0', '0', prefix, '16'])
        got = next(Fr(int(l.split()[1]), int(l.split()[2])) for l in out.splitlines() if l.startswith('VALUE'))
        vout, _, _, vrc = timed([BIN['verify'], str(q), str(k), str(got.numerator), str(got.denominator), prefix])
        if not (got == v and vrc == 0 and 'RESULT CERTIFIED' in vout):
            check(False, f"A1 mismatch q={q} k={k}: {got} vs {v}")
        n += 1
check(True, f"{n} (q,k) pairs (q = 5, 7, 9, 11; k = 8..18): lean-solver value = THM-4486 value, both certificates "
            f"accepted by the independent checker")

say("A2. the independent checker agrees with the floor lane's Python verifiers (a third implementation) on the same "
    "certificates, converted to node form (x = P + tau(P) H carries g(P); psi(x) = Psi(x mod H)):")


def load_lean(prefix, k):
    H = 1 << (k - 1)
    out = {}
    for kind, bits, pot in (('lo', 'lo_taub', ('lo_g8', 'lo_g16')), ('up', 'up_sigb', ('up_psi8', 'up_psi16'))):
        b = np.unpackbits(np.fromfile(prefix + '.' + bits, dtype=np.uint8), bitorder='little')[:H]
        if os.path.exists(prefix + '.' + pot[0]):
            v = np.fromfile(prefix + '.' + pot[0], dtype=np.uint8).astype(np.int64)
            v[v == 255] = -1
        else:
            v = np.fromfile(prefix + '.' + pot[1], dtype=np.uint16).astype(np.int64)
            v[v == 65535] = -1
        out[kind] = (b.astype(np.uint8), v)
    return out


n = 0
for q, k in [(7, 12), (7, 14), (7, 16), (5, 15), (11, 14)]:
    prefix = os.path.join(SCR, f'a1_q{q}_k{k}')
    F = FLOOR[q][k]
    c = load_lean(prefix, k)
    N, H = 1 << k, 1 << (k - 1)
    tau, g = c['lo']
    W_nodes = np.zeros(N, dtype=bool)
    f_nodes = np.zeros(N, dtype=np.int64)
    P = np.arange(H)
    inW = g >= 0
    xs = P + tau.astype(np.int64) * H
    W_nodes[xs[inW]] = True
    f_nodes[xs[inW]] = g[inW]
    sig, psi = c['up']
    flip = np.zeros(N, dtype=np.uint8)
    flip[1::2] = sig
    psi_nodes = psi[np.arange(N) % H]
    ok = FL.verify_lower(k, q, W_nodes, tau, f_nodes, F) and FL.verify_upper(k, q, flip, psi_nodes, F)
    if not ok:
        check(False, f"A2 floor-lane verifier rejects q={q} k={k}")
    n += 1
check(True, f"{n} certificates (q = 5, 7, 11; k = 12..16) are accepted by the floor lane's verify_lower/verify_upper")

say("A3. negative tests: the checker rejects corrupted certificates:")
q, k = 7, 16
prefix = os.path.join(SCR, f'a1_q{q}_k{k}')
F = FLOOR[q][k]
raw_g = open(prefix + '.lo_g8', 'rb').read() if os.path.exists(prefix + '.lo_g8') else None
raw_p = open(prefix + '.up_psi8', 'rb').read() if os.path.exists(prefix + '.up_psi8') else None
check(raw_g is not None and raw_p is not None, "the q=7, k=16 certificate uses the uint8 format")
c = load_lean(prefix, k)
H = 1 << (k - 1)
tau, g = c['lo']
# a tight lower edge: P -> Q with g(Q) = g(P) + e; raising g(Q) by 1 must break it
e_odd, e_even = F.denominator - F.numerator, -F.numerator
found = None
for Pp in range(H):
    if g[Pp] < 0:
        continue
    x = Pp + int(tau[Pp]) * H
    e = e_odd if x % 2 else e_even
    tg = [(x // 2) % H] if x % 2 == 0 else [((q * x + 1) // 2) % H, ((q * x - 1) // 2) % H]
    for Qq in tg:
        if g[Qq] == g[Pp] + e:
            found = Qq
            break
    if found is not None:
        break
bad = bytearray(raw_g)
bad[found] = bad[found] + 1
with open(prefix + '.lo_g8', 'wb') as f:
    f.write(bytes(bad))
vout, _, _, vrc = timed([BIN['verify'], str(q), str(k), str(F.numerator), str(F.denominator), prefix])
with open(prefix + '.lo_g8', 'wb') as f:
    f.write(raw_g)
r1 = vrc != 0 and 'LOWER FAIL' in vout and 'UPPER OK' in vout
sig, psi = c['up']
tightQ = None
N = 1 << k
for x in range(N):
    e = e_odd if x % 2 else e_even
    Qq = (x // 2) % H if x % 2 == 0 else ((q * x + (-1 if sig[(x - 1) // 2] else 1)) // 2) % H
    if psi[Qq] + e == psi[x % H] and psi[Qq] < 254:
        tightQ = Qq
        break
bad = bytearray(raw_p)
bad[tightQ] += 1
with open(prefix + '.up_psi8', 'wb') as f:
    f.write(bytes(bad))
vout2, _, _, vrc2 = timed([BIN['verify'], str(q), str(k), str(F.numerator), str(F.denominator), prefix])
with open(prefix + '.up_psi8', 'wb') as f:
    f.write(raw_p)
r2 = vrc2 != 0 and 'UPPER FAIL' in vout2 and 'LOWER OK' in vout2
vout3, _, _, vrc3 = timed([BIN['verify'], str(q), str(k), '5', '13', prefix])
vout4, _, _, vrc4 = timed([BIN['verify'], str(q), str(k), str(F.numerator), str(F.denominator), prefix])
check(r1 and r2 and vrc3 != 0 and vrc4 == 0,
      "raising g by 1 at the head of a tight lower edge -> LOWER FAIL; raising psi by 1 at the head of a tight upper edge -> "
      "UPPER FAIL; the q=7, k=16 certificates do not certify 5/13; restored files pass again")

say("A5. the certified strategies re-evaluated by Karp's algorithm (the floor lane's pure-Python rho_max_exact for "
    "sigma; its C Karp on the reachable part of G^tau for tau):")
n = 0
for q, k, F in [(7, 8, Fr(3, 7)), (7, 10, Fr(2, 5)), (5, 10, Fr(3, 7)), (9, 10, Fr(5, 12)), (11, 9, None)]:
    prefix = os.path.join(SCR, f'a5_q{q}_k{k}')
    clean(prefix)
    out, _, _, _ = timed([BIN['game'], 'lsearch2', str(q), str(k), '0', '1', '1', '1', '0', '0', prefix, '16'])
    got = next(Fr(int(l.split()[1]), int(l.split()[2])) for l in out.splitlines() if l.startswith('VALUE'))
    c = load_lean(prefix, k)
    Nk, Hk = 1 << k, 1 << (k - 1)
    sig, psi = c['up']
    flip = np.zeros(Nk, dtype=np.uint8)
    flip[1::2] = sig
    r_sigma = FL.rho_max_exact(k, q, flip)
    tau, g = c['lo']
    G = FL.Game(k, q)
    P0 = int(np.nonzero(g >= 0)[0][0])
    m_tau = G.min_mean_tau(tau, P0 + int(tau[P0]) * Hk)
    if not (r_sigma == got and m_tau is not None and m_tau >= got and (F is None or got == F)):
        check(False, f"A5 Karp mismatch q={q} k={k}: {got}, {r_sigma}, {m_tau}")
    n += 1
check(True, f"{n} instances (k = 8..10): rho_max of the certified sign strategy (exhaustive Karp) equals the value, and "
            f"the least cycle density Min can reach against the certified lift strategy (Karp) is >= the value")

say("A4. agreement with the floor lane's engine (a different solver) on fresh instances:")
n = 0
for q, k in [(7, 11), (7, 13), (13, 12), (15, 12), (17, 11), (21, 12)]:
    res = FL.rho_star(k, q)
    FL.certify(k, q, res)
    prefix = os.path.join(SCR, f'a4_q{q}_k{k}')
    clean(prefix)
    out, _, _, _ = timed([BIN['game'], 'lsearch2', str(q), str(k), '0', '1', '1', '1', '0', '0', prefix, '16'])
    got = next(Fr(int(l.split()[1]), int(l.split()[2])) for l in out.splitlines() if l.startswith('VALUE'))
    vout, _, _, vrc = timed([BIN['verify'], str(q), str(k), str(got.numerator), str(got.denominator), prefix])
    if not (got == res['rho'] and vrc == 0):
        check(False, f"A4 mismatch q={q} k={k}")
    n += 1
check(True, f"{n} instances (q = 7, 13, 15, 17, 21; k = 11..13): identical values from both engines, both certified")

# ================================================================================================= B
say("\n== B. Exact min-max densities beyond THM-4486 (every entry: both certificates re-checked on all 2^k nodes) ==")
VAL = {}


HINTS = {(7, 27): ((16, 43), (3, 8))}


def chain(q, k1, k2, ku0, known=None):
    ku = ku0
    for k in range(k1, k2 + 1):
        lo, hi = HINTS.get((q, k), ((0, 1), (1, 1)))
        if (q, k) in HINTS:
            say(f"    (hint for q={q}, k={k}: Farey bracket [{lo[0]}/{lo[1]}, {hi[0]}/{hi[1]}])")
        val, info = solve(q, k, lo=lo, hi=hi, ku=(ku.numerator, ku.denominator),
                          keep=(q == 7 and k in (KSTRUCT, KMAX7)))
        VAL[(q, k)] = val
        say(f"    q={q:2d} k={k:2d}: rho* = {str(val):>7s} = {float(val):.5f}   (lean search {info['wall']:.1f} s, "
            f"RSS {info['rss'] / 2 ** 20:.0f} MB; checker {info['vwall']:.1f} s, RSS {info['vrss'] / 2 ** 20:.0f} MB; "
            f"digest {info['digest']}; {info['lower'].split('(')[1].split(',')[2].strip()}, "
            f"{info['upper'].split('(')[1].split(',')[1].strip().rstrip(')')})")
        if known and k in known:
            check(val == known[k], f"q={q} k={k}: equals the THM-4486 value {known[k]}")
        check(val <= ku, f"q={q} k={k}: non-increasing ({val} <= {ku})")
        ku = val


say("B1. q = 7:")
chain(7, 19, KMAX7, Fr(19, 49), known={19: Fr(13, 34), 20: Fr(34, 89), 21: Fr(8, 21), 22: Fr(14, 37)})
ks7 = sorted(k for (q, k) in VAL if q == 7)
check(all(not below_log(7, VAL[(7, k)]) for k in ks7),
      f"7^a > 2^p for rho*(7,k) = a/p at every k = 19..{KMAX7}: class (i) is EMPTY for 7n+-1 at every level k <= {KMAX7} "
      f"(k <= 18 by THM-4486 and monotonicity)")
if KMAX7 >= 27:
    check(VAL[(7, 23)] == Fr(14, 37) and all(VAL[(7, k)] == Fr(3, 8) for k in (24, 25, 26)) and VAL[(7, 27)] < Fr(3, 8),
          f"q = 7: 14/37 at k = 22, 23; a plateau at exactly 3/8 for k = 24, 25, 26; below 3/8 at k = 27 "
          f"({VAL[(7, 27)]} = {float(VAL[(7, 27)]):.5f})")
    gap = float(VAL[(7, KMAX7)]) - LOG2 / math.log(7)
    check(VAL[(7, 27)] == Fr(35, 94) and gap > 0.01,
          f"q = 7: rho*(7,27) = 35/94; the value at k = {KMAX7} ({VAL[(7, KMAX7)]}) still exceeds log_7 2 = "
          f"{LOG2 / math.log(7):.6f} by {gap:.4f}")


def one_sided(q, k, F, side):
    """a lower-only (side = 'lower') or upper-only ('upper') certificate at threshold F, re-checked"""
    prefix = os.path.join(SCR, f'{side}_q{q}_k{k}')
    clean(prefix)
    out, wall, rss, rc = timed([BIN['game'], side, str(q), str(k), str(F.numerator), str(F.denominator), prefix])
    note_rss(rss, f'{side} q={q} k={k}')
    check(rc == 0 and f"{side.upper()}-WRITTEN" in out, f"{side} energy fixed point at F = {F} finite for q={q}, k={k} "
                                                         f"({wall:.1f} s, RSS {rss / 2 ** 20:.0f} MB)")
    vout, vwall, vrss, vrc = timed([BIN['verify'], str(q), str(k), str(F.numerator), str(F.denominator), prefix])
    note_rss(vrss, f'checker {side} q={q} k={k}')
    tag = 'LOWER-CERTIFIED' if side == 'lower' else 'UPPER-CERTIFIED'
    line = next(l for l in vout.splitlines() if l.startswith('LOWER' if side == 'lower' else 'UPPER'))
    check(vrc == 0 and f"RESULT {tag} {F.numerator}/{F.denominator}" in vout,
          f"q={q} k={k}: rho* {'>=' if side == 'lower' else '<='} {F} certified by the independent checker "
          f"({vwall:.1f} s, RSS {vrss / 2 ** 20:.0f} MB; {line.split('(')[1].rstrip(')')})")
    for ext in ('lo_taub', 'lo_g8', 'lo_g16', 'up_sigb', 'up_psi8', 'up_psi16'):
        if os.path.exists(prefix + '.' + ext):
            os.remove(prefix + '.' + ext)


if KBR7 and KMAX7 >= 27:
    say("B1'. q = 7 at k = 28, 29, 30 by one-sided certificates (compact engine; a lower and an upper certificate at "
        "the same threshold certify the value):")
    one_sided(7, 28, Fr(13, 35), 'lower')
    one_sided(7, 28, Fr(13, 35), 'upper')
    VAL[(7, 28)] = Fr(13, 35)
    one_sided(7, 29, Fr(13, 35), 'lower')
    VAL[(7, 29)] = Fr(13, 35)
    one_sided(7, 30, Fr(7, 19), 'lower')
    check(Fr(13, 35) < VAL[(7, 27)] and not below_log(7, Fr(13, 35)) and not below_log(7, Fr(7, 19)),
          "rho*(7,28) = 13/35 (both certificates at 13/35) < 35/94 = rho*(7,27): the value drops again; rho*(7,29) = 13/35 "
          "(lower certificate at k = 29, upper bound by monotonicity); rho*(7,30) >= 7/19 = 0.36842. Since 7^13 > 2^35 and "
          "7^7 > 2^19, class (i) is EMPTY for 7n+-1 at every level k <= 30 (monotonicity below)")
say("B2. q = 9 and q = 11 (THM-4486 computed k <= 18):")
chain(9, 17, KMAXQ, Fr(2, 5), known={17: Fr(12, 31), 18: Fr(5, 13)})
chain(11, 17, KMAXQ, Fr(21, 55), known={17: Fr(29, 77), 18: Fr(41, 110)})
for q in (9, 11):
    ks = sorted(k for (qq, k) in VAL if qq == q)
    check(all(not below_log(q, VAL[(q, k)]) for k in ks),
          f"q={q}: rho*(q,k) > log_q 2 = {LOG2 / math.log(q):.5f} exactly at every k <= {KMAXQ} (class (i) empty); "
          f"values {[str(VAL[(q, k)]) for k in ks]}")

# ================================================================================================= C
say(f"\n== C. Structure: the 3/8 plateau (k = {KSTRUCT}) and the level k = {KMAX7} ==")


def core_edges(kind, prefix, q, k, F):
    out, _, rss, rc = timed([BIN['tight'], str(q), str(k), str(F.numerator), str(F.denominator), prefix, kind])
    note_rss(rss, f'tight core q={q} k={k} {kind}')
    E = [tuple(int(t) for t in l.split()[1:]) for l in out.splitlines() if l.startswith('E ')]
    size = int(out.splitlines()[0].split()[1])
    return size, E


def cycles_in_core(E, H, q, maxlen=128, tries=4000):
    """shortest cycles through sampled core nodes (BFS in the tight core); returns [(p, a, x0, orbit)]"""
    from collections import deque
    adj = {}
    for P, Q, b, s_ in E:
        adj.setdefault(P, []).append((Q, b, s_))
    found = {}
    for start in list(adj)[:tries]:
        prev = {start: None}
        dq = deque([start])
        cyc = None
        while dq and cyc is None:
            u = dq.popleft()
            for (v, b, s_) in adj.get(u, []):
                if v == start:
                    path = [(u, b, s_)]
                    w = u
                    while prev[w] is not None:
                        pw, pb, ps = prev[w]
                        path.append((pw, pb, ps))
                        w = pw
                    cyc = path[::-1]
                    break
                if v not in prev and v in adj:
                    prev[v] = (u, b, s_)
                    dq.append(v)
        if cyc is None or len(cyc) > maxlen:
            continue
        xs = [P + b * H for (P, b, s_) in cyc]
        ss = [s_ for (P, b, s_) in cyc]
        p_, a_ = len(xs), sum(1 for x in xs if x % 2)
        c_, oa = 0, a_
        for j, x in enumerate(xs):
            if x % 2:
                oa -= 1
                c_ += ss[j] * q ** oa * 2 ** j
        x0 = Fr(c_, 2 ** p_ - q ** a_)
        orb, y = [x0], x0
        for j in range(p_ - 1):
            y = y / 2 if xs[j] % 2 == 0 else (q * y + ss[j]) / 2
            orb.append(y)
        found[min(orb)] = (p_, a_, x0, orb)
    return sorted(found.values(), key=lambda t: t[0])


def sig_of(prefix, k):
    H = 1 << (k - 1)
    b = np.unpackbits(np.fromfile(prefix + '.up_sigb', dtype=np.uint8), bitorder='little')[:H]
    return lambda i: int(b[i])


if KMAX7 >= KSTRUCT >= 24:
    q, k, F = 7, KSTRUCT, VAL[(7, KSTRUCT)]
    N, H = 1 << k, 1 << (k - 1)
    prefix = os.path.join(SCR, f'c_q7_k{k}')
    check(F == Fr(3, 8), f"k = {k} lies on the plateau (rho* = 3/8)")
    c = load_lean(prefix, k)
    sig, psi = c['up']
    x0 = Fr(-169, 87)
    orb = [x0]
    for _ in range(8):
        y = orb[-1]
        orb.append(y / 2 if y.numerator % 2 == 0 else (q * y + 1) / 2)
    rr = [CY.residue(y, N) for y in orb]
    signs_ok = all(sig[(r - 1) // 2] == 0 for r in rr[:8] if r % 2 == 1)
    e_of = lambda r: F.denominator - F.numerator if r % 2 else -F.numerator
    tight_ok = all(psi[rr[i + 1] % H] + e_of(rr[i]) == psi[rr[i] % H] for i in range(8))
    a8 = sum(1 for y in orb[:8] if y.numerator % 2)
    check(orb[8] == x0 and a8 == 3 and len(set(rr[:8])) == 8 and signs_ok and tight_ok and Fr(-109, 87) in orb,
          f"C1. the 7x+1 cycle (-169, -548, -274, -137, -436, -218, -109, -338)/87 (shape (8,3), 2^8 - 7^3 = -87, density "
          f"3/8) is a cycle of the certified optimal sign strategy at k = {k}, and every edge of it is tight for the "
          f"certified potential psi")
    size, E = core_edges('lo', prefix, q, k, F)
    cyc = cycles_in_core(E, H, q, maxlen=64, tries=3000)
    shapes = sorted(set((p_, a_) for p_, a_, _, _ in cyc))
    dens87 = [x for p_, a_, x, _ in cyc if (p_, a_) == (8, 3)]
    check(size > 0 and len(cyc) > 0 and all(Fr(a_, p_) == F for p_, a_ in shapes) and (8, 3) in shapes and
          all(87 % x.denominator == 0 for x in dens87),
          f"C2. lower certificate at k = {k}: tight core of {size} pairs; Min's best responses against the certified Max "
          f"strategy (tight cycles) all have density exactly 3/8, with shapes {shapes[:8]}; the (8,3) ones are "
          f"rational cycles with denominators dividing 87; all are real-bounded (|x| <= "
          f"{float(max(max(abs(y) for y in o) for _, _, _, o in cyc)):.0f})")
    say("C3. census of the simple rational cycles of small shape kept by the certified optimal strategies (a cycle is "
        "kept iff sigma agrees with its sign at the residue of each of its odd points):")
    SHAPES = [(2, 1), (5, 2), (7, 3), (8, 3), (13, 5), (11, 4), (14, 5)]
    CYC = {sh: CY.shape_cycles(7, *sh) for sh in SHAPES}
    CEN = {}
    for kk in sorted(set((KSTRUCT, KMAX7))):
        sb = sig_of(os.path.join(SCR, f'c_q7_k{kk}'), kk)
        CEN[kk] = {sh: sum(CY.present(kk, sb, cc) for cc in CYC[sh]) for sh in SHAPES}
        say(f"    k={kk} (rho* = {VAL[(7, kk)]}): " + ", ".join(
            f"({p_},{a_}) {'exp' if 7 ** a_ > 2 ** p_ else 'con'} {CEN[kk][(p_, a_)]}/{len(CYC[(p_, a_)])}" for p_, a_ in SHAPES))
    ok = all(CEN[kk][sh] == 0 for kk in CEN for sh in SHAPES if Fr(sh[1], sh[0]) > VAL[(7, kk)])
    check(ok and CEN[KSTRUCT][(8, 3)] > 0 and len(CYC[(8, 3)]) == 54,
          f"every simple cycle of density > rho* among these shapes is broken (as it must be), the plateau strategy keeps "
          f"{CEN[KSTRUCT][(8, 3)]} of the 54 simple (8,3)-cycles"
          + (f", and at k = {KMAX7} none of them" if KMAX7 > KSTRUCT and CEN[KMAX7][(8, 3)] == 0 else "")
          + f"; hundreds of expanding (14,5)-cycles (7^5 > 2^14, density 5/14 = 0.35714) survive in these optimal "
            f"strategies, while a class-(i) strategy would have to break all {len(CYC[(14, 5)])} of them")
    say("C4. the rational relaxation: the real-sign assignment sigma(x) = sgn(x) keeps no expanding cycle (Proposition "
        "S of the note, proved by hand); exhaustive check over all simple cycles of 7x+-1 with p <= 11:")
    nexp = ncon = kexp = kcon = 0
    for p_ in range(1, 12):
        for a_ in range(1, p_):
            for cc in CY.shape_cycles(7, p_, a_):
                kp = CY.sgn_keeps(cc)
                if 7 ** a_ > 2 ** p_:
                    nexp += 1
                    kexp += kp
                else:
                    ncon += 1
                    kcon += kp
    check(kexp == 0 and kcon > 0,
          f"{nexp} expanding simple cycles (p <= 11): sgn keeps none; {ncon} contracting ones: sgn keeps {kcon} (e.g. the "
          f"free cycle (1, 4, 2)); so no finite set of rational cycles obstructs provability at every level")
    for kk in sorted(set((KSTRUCT, KMAX7))):
        pre = os.path.join(SCR, f'c_q7_k{kk}')
        for ext in ('lo_taub', 'lo_g8', 'lo_g16', 'up_sigb', 'up_psi8', 'up_psi16'):
            if os.path.exists(pre + '.' + ext):
                os.remove(pre + '.' + ext)

# ================================================================================================= D
say("\n== D. All-level floors for q = 9 and q = 11 from the negative-integer adversary (corrected potentials) ==")
POT = [(9, 5, 16, 5, 16384, [1, 5, 22, 11, 50, 25, 112, 56, 28, 14, 7, 32, 16, 8, 4, 2]),
       (11, 2, 7, 2, 1024, [1, 6, 3, 16, 8, 4, 2])]
for q, a0, p0, beta, U0, cyc in POT:
    t1 = time.time()
    Phi, A, B, passes = PT.least_potential(q, a0, p0, beta, U0)
    r = PT.check_potential(q, a0, p0, beta, U0, Phi, A, B)
    a, p = PT.tight_cycle_check(q, cyc)
    kmin = next(k for k in range(2, 40) if 2 ** (k - 1) >= r['Hmin'])
    check(r['C1'] and r['C2'] and r['C3'] and r['C4'] and r['C5'],
          f"q={q}: Phi(u) = u^{beta} for u >= {U0}, least solution below (exact rationals, {passes} passes): C1 (E),(O) on "
          f"u < {U0}; C2 Phi(w) <= w^{beta} on [{U0 // 2}, {U0}); C3 ({q}*{U0}+1)^{beta} <= A 2^{beta} {U0}^{beta} with "
          f"A = 2^{int(math.log2(A))}; C4 max Phi <= A {U0}^{beta}; C5 Phi > 0; max Phi(u)/u^{beta} = "
          f"{float(r['cmax']):.4f} ({time.time() - t1:.1f} s)")
    check(Fr(a, p) == Fr(a0, p0) and max(cyc) <= 2 ** (kmin - 1),
          f"q={q}: hence every cycle of U_k has density >= {a0}/{p0} for every k >= {kmin} (H >= {r['Hmin']}), and so "
          f"rho*({q},k) >= {a0}/{p0} for EVERY k (monotonicity below {kmin}); the cycle {tuple(cyc)} of u -> u/2, "
          f"({q}u -+ 1)/2 lies in U_k (k >= {kmin}) and has density exactly {a0}/{p0}: the negative-integer adversary's "
          f"value is exactly {a0}/{p0}")
    direct = [(kk, PT.check_Uk_edges(q, a0, p0, beta, U0, Phi, A, B, kk)[0]) for kk in (kmin, kmin + 1, kmin + 2)]
    check(all(ok_ for _, ok_ in direct),
          f"q={q}: direct exact check (independent of the case analysis) of Phi(v) <= lambda Phi(u) on every edge of U_k "
          f"for k = {kmin}, {kmin + 1}, {kmin + 2} (wrap-arounds included)")
    if q == 9:
        Phi2, A2, B2, _ = PT.least_potential(q, a0, p0, beta, 4096)
        r2 = PT.check_potential(q, a0, p0, beta, 4096, Phi2, A2, B2)
        check(r2['C1'] and not r2['C2'],
              f"q=9 control: with U0 = 4096 the least solution exists but C2 fails (max Phi(u)/u^5 = "
              f"{float(r2['cmax']):.4f} > 1), so a large U0 is needed")
    check(LOG2 / math.log(q + 1) < a0 / p0 < LOG2 / math.log(q) and below_log(q, Fr(a0, p0)),
          f"q={q}: log_{q + 1} 2 = {LOG2 / math.log(q + 1):.5f} < {a0}/{p0} = {a0 / p0:.5f} < log_{q} 2 = "
          f"{LOG2 / math.log(q):.5f}: Theorem N is improved at every level, but the floor stays below the provability "
          f"threshold ({q}^{a0} < 2^{p0})")

# ================================================================================================= E
say("\n== E. Level-independent adversaries tried for q = 7 (none beats 1/3; exact values of frozen lift strategies) ==")
FAM = {'top lift (Theorem N)': lambda k, H: np.ones(H, dtype=np.uint8)}
for cfr in (Fr(1, 8), Fr(1, 4), Fr(1, 2), Fr(3, 4), Fr(7, 8)):
    FAM[f'window -x in [{1 - cfr}H, {2 - cfr}H)'] = (lambda c: lambda k, H: (np.arange(H) <= int(c * H)).astype(np.uint8))(cfr)


def switched(classes):
    def f(k, H):
        t = np.ones(H, dtype=np.uint8)
        for c in classes(H):
            t[c % H] = 0
        return t
    return f


FAM['top lift switched at the class of -1'] = switched(lambda H: [H - 1])
FAM['top lift switched at -1, -2, -4'] = switched(lambda H: [H - 1, H - 2, H - 4])
rng = np.random.default_rng(20260926)
FAM['uniformly random lifts'] = lambda k, H: rng.integers(0, 2, H).astype(np.uint8)
best = Fr(0)
SW = []
for name, fam in FAM.items():
    row = []
    for k in (8, 10, 12):
        H = 1 << (k - 1)
        v, _ = FL.value_of_tau(k, 7, fam(k, H))
        row.append(v)
        best = max(best, v)
        if name.startswith('top lift switched') and k >= 10:
            SW.append(v)
    say(f"    {name:34s}: value at k = 8, 10, 12: {[str(v) for v in row]}")
check(best == Fr(1, 3) and all(v < Fr(1, 3) for v in SW),
      "the best of these frozen adversaries is exactly 1/3 = Theorem N (EMPIRICAL at k = 8, 10, 12); switching the top "
      "lift at the free cycle's classes (breaking (-1,-4,-2)) makes the adversary WORSE at k = 10, 12 "
      f"({[str(v) for v in SW]}): the switched class is 2-adically next to the free cycle and gives Min a free ascent")

# ================================================================================================= F
say("\n== F. Summary ==")
for q in (7, 9, 11):
    ks = sorted(k for (qq, k) in VAL if qq == q)
    say(f"  q={q:2d} (log_q 2 = {LOG2 / math.log(q):.5f}): " +
        ", ".join(f"k={k}:{VAL[(q, k)]}" for k in ks if k == ks[0] or VAL[(q, k)] != VAL[(q, k - 1)]) +
        (f" (exact to k = {ks[-1]}; rho*(7,30) >= 7/19)" if q == 7 and KBR7 else f" (exact to k = {ks[-1]})"))
say(f"\nchecks: {NCHECK[0]}; wall time {time.time() - T0:.0f} s; peak RSS of a child process "
    f"{PEAK['rss'] / 2 ** 20:.0f} MB ({PEAK['what']}); runner RSS "
    f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20:.0f} MB")
say("ALL CHECKS PASSED")
