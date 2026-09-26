#!/usr/bin/env python3
"""
procgen_seven2_20260926_lib.py -- helpers of lane "seven2" (session collatz-procgen-20260922, 2026-09-26):
structure of the optimal sign strategies of q n +- 1 (q = 5, 7), automatic strategies, adversaries.

Setting (THM-4486).  Level k, odd q, N = 2^k, H = 2^(k-1).  A sign strategy sigma gives the odd residues mod N a
sign; T_sigma(x) = x/2 (x even), (q x + sigma(x))/2 (x odd).  rho_max(sigma) = the largest odd density of a cycle of
the parity graph G_sigma (nodes Z/N, edges x -> both lifts of T_sigma(x) mod H).  A strategy is stored as a boolean
array `minus` over the odd nodes x = 2i+1 (i = 0..H-1), True = sign '-'.

Engines (compiled into scratch/procgen_seven2/run/, not committed):
  procgen_seven2_20260926_rhomax.c    exact rho_max of one strategy (witness cycle + potential certificate)
  procgen_seven2_20260926_restrict.c  Min's energy game with restricted signs (ctypes; search engine only)
  procgen_seven_20260926_game.c       the seven lane's exact game solver (read-only reuse)
  procgen_seven_20260926_verify.c     the seven lane's independent certificate checker (read-only reuse)
Nothing an engine prints is used as a proof step without an exact re-check: witness cycles are re-walked here with
exact integers and identified with rational cycles (fractions.Fraction); potential certificates are re-checked by
the seven lane's independent checker.
"""
import os
import ctypes
import hashlib
import subprocess
from fractions import Fraction as Fr

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
SCR = os.path.join(ROOT, 'scratch', 'procgen_seven2', 'run')
BIN = {}
LIBR = None


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        h.update(f.read())
    return h.hexdigest()


def build():
    """compile the engines into SCR (source hash in the name); returns the dict of executables"""
    global LIBR
    os.makedirs(SCR, exist_ok=True)
    for name, src in (('rhomax', 'procgen_seven2_20260926_rhomax.c'), ('game', 'procgen_seven_20260926_game.c'),
                      ('verify', 'procgen_seven_20260926_verify.c'), ('tight', 'procgen_seven_20260926_tight.c')):
        p = os.path.join(HERE, src)
        exe = os.path.join(SCR, f'{name}_{sha(p)[:12]}')
        if not os.path.exists(exe):
            subprocess.run(['cc', '-O2', '-Wall', '-o', exe + '.tmp', p], check=True)
            os.replace(exe + '.tmp', exe)
        BIN[name] = exe
    p = os.path.join(HERE, 'procgen_seven2_20260926_restrict.c')
    so = os.path.join(SCR, f'restrict_{sha(p)[:12]}.so')
    if not os.path.exists(so):
        subprocess.run(['cc', '-O2', '-shared', '-fPIC', '-o', so + '.tmp', p], check=True)
        os.replace(so + '.tmp', so)
    L = ctypes.CDLL(so)
    ll = ctypes.c_longlong
    L.r_init.argtypes = [ctypes.c_int, ll]
    L.r_solve.argtypes = [ll, ll, ll]
    L.r_solve.restype = ll
    L.r_mask_class.argtypes = [ll, ctypes.c_int, ctypes.c_int]
    L.r_get_mask.argtypes = [ctypes.c_void_p]
    L.r_set_mask.argtypes = [ctypes.c_void_p]
    L.r_strategy.argtypes = [ctypes.c_void_p]
    LIBR = L
    return BIN


# ------------------------------------------------------------------------------------------------ strategies
def odd_nodes(k):
    return np.arange(1, 1 << k, 2, dtype=np.int64)


def mh_minus(q, x):
    """max-halving (MH): the sign s with q x + s = 0 mod 4; True where s = -1"""
    return (q * x + 1) % 4 != 0


def lift(minus, k, k2):
    """the level-k strategy viewed at level k2 >= k (same map T_sigma)"""
    x = odd_nodes(k2)
    return minus[((x % (1 << k)) - 1) // 2]


def write_sig(minus, path):
    np.packbits(np.asarray(minus, dtype=np.uint8), bitorder='little').tofile(path)


def rhomax(q, k, minus, tag='s', F0=None, cert=None):
    """exact rho_max(sigma) by the C engine.  Returns (value, cycle, le) where value is a Fraction, cycle the witness
    (list of nodes of G_sigma) and le = True if the engine only reports rho_max <= F0 (no witness)."""
    path = os.path.join(SCR, f'{tag}_q{q}_k{k}.sig')
    write_sig(minus, path)
    cmd = [BIN['rhomax'], str(q), str(k), path]
    if F0 is not None or cert is not None:
        F0 = F0 if F0 is not None else Fr(0)
        cmd += [str(F0.numerator), str(F0.denominator)]
        if cert:
            cmd += [cert]
    out = subprocess.run(cmd, capture_output=True, text=True).stdout
    val, cyc, le = None, None, False
    for line in out.splitlines():
        if line.startswith('RHOMAX '):
            a, b = line.split()[1:]
            val = Fr(int(a), int(b))
        elif line.startswith('RHOMAX_LE'):
            a, b = line.split()[1:]
            val, le = Fr(int(a), int(b)), True
        elif line.startswith('CYCLE'):
            cyc = [int(t) for t in line.split()[1:]]
    if val is None:
        raise RuntimeError('rhomax failed: ' + out)
    return val, cyc, le


def check_cycle(q, k, minus, cyc):
    """exact check that cyc is a closed walk of G_sigma; returns (odd count, length)"""
    H = 1 << (k - 1)
    for i, x in enumerate(cyc):
        y = cyc[(i + 1) % len(cyc)]
        if x % 2 == 0:
            t = (x // 2) % H
        else:
            s = -1 if minus[(x - 1) // 2] else 1
            t = ((q * x + s) // 2) % H
        if y % H != t:
            raise AssertionError('not a closed walk of G_sigma')
    return sum(1 for x in cyc if x % 2), len(cyc)


def residue(r, N):
    return (r.numerator * pow(r.denominator, -1, N)) % N


def cycle_rational(q, k, minus, cyc):
    """the rational periodic point of T_sigma coded by a closed walk (Lemma C of the note), with an exact check that
    its orbit points have the residues of the walk; returns (orbit, signs)"""
    N = 1 << k
    ss = [(-1 if minus[(x - 1) // 2] else 1) if x % 2 else 0 for x in cyc]
    p, a = len(cyc), sum(1 for x in cyc if x % 2)
    c, oa = 0, a
    for j, x in enumerate(cyc):
        if x % 2:
            oa -= 1
            c += ss[j] * q ** oa * 2 ** j
    x0 = Fr(c, 2 ** p - q ** a)
    orb, y = [], x0
    for j in range(p):
        orb.append(y)
        if residue(y, N) != cyc[j] % N:
            raise AssertionError('residue mismatch')
        y = y / 2 if cyc[j] % 2 == 0 else (q * y + ss[j]) / 2
    if y != x0:
        raise AssertionError('not periodic')
    return orb, ss


def rat_recon(x, k):
    """for each odd residue x mod 2^k: a short vector (a, D), D odd > 0, of the lattice {(a, D): a = D x mod 2^k}
    (Lagrange reduction; then the shortest small combination with D odd) -- 'the simplest rational in the class'"""
    N = 1 << k
    x = np.asarray(x, dtype=np.int64)
    u0 = np.full(len(x), N, dtype=np.int64)
    u1 = np.zeros(len(x), dtype=np.int64)
    v0 = x.copy()
    v1 = np.ones(len(x), dtype=np.int64)
    for _ in range(300):
        nu = u0 * u0 + u1 * u1
        nv = v0 * v0 + v1 * v1
        sw = nv > nu
        u0, v0 = np.where(sw, v0, u0), np.where(sw, u0, v0)
        u1, v1 = np.where(sw, v1, u1), np.where(sw, u1, v1)
        nv = np.where(sw, nu, nv)
        m = np.rint((u0 * v0 + u1 * v1) / nv).astype(np.int64)
        if not np.any(m != 0):
            break
        u0 = u0 - m * v0
        u1 = u1 - m * v1
    ba = np.zeros(len(x), dtype=np.int64)
    bD = np.zeros(len(x), dtype=np.int64)
    bn = np.full(len(x), np.iinfo(np.int64).max)
    for cu, cv in ((0, 1), (1, 0), (1, 1), (1, -1), (1, 2), (1, -2), (2, 1), (2, -1)):
        a = cu * u0 + cv * v0
        D = cu * u1 + cv * v1
        n = a * a + D * D
        ok = (D % 2 != 0) & (n < bn)
        ba, bD, bn = np.where(ok, a, ba), np.where(ok, D, bD), np.where(ok, n, bn)
    sg = np.sign(bD)
    return ba * sg, bD * sg


# ------------------------------------------------------------------------------------------------ MH itineraries
def mh_itin(q, x, steps):
    """exact MH accelerated itinerary of the 2-adic integer x (a Python int, known exactly): [(s, v), ...]"""
    out = []
    for _ in range(steps):
        s = 1 if (q * x + 1) % 4 == 0 else -1
        t = q * x + s
        v = (t & -t).bit_length() - 1
        out.append((s, v))
        x = t >> v
    return out


def sinf_periodic_points(q, n):
    """the periodic points of period dividing n of the v=2 set S_inf of MH: for each sign word (s_1..s_n) the fixed
    point of x -> (q x + s)/4 composed, kept iff its MH itinerary is (s_i, 2) at every step (exact rationals)"""
    pts = []
    for w in range(1 << n):
        signs = [1 if (w >> i) & 1 else -1 for i in range(n)]
        # x = (q^n x + sum_i s_i q^(n-1-i) 4^i) / 4^n  ==>  x (4^n - q^n) = sum ...
        c = sum(signs[i] * q ** (n - 1 - i) * 4 ** i for i in range(n))
        x0 = Fr(c, 4 ** n - q ** n)
        y, ok = x0, True
        for i in range(n):
            if y.denominator % 2 == 0 or y.numerator % 2 == 0:
                ok = False
                break
            ys = y.numerator * pow(y.denominator, -1, 8) % 8
            s = 1 if (q * ys + 1) % 4 == 0 else -1
            if s != signs[i] or (q * ys + s) % 8 != 4:
                ok = False
                break
            y = (q * y + s) / 4
        if ok and y == x0:
            pts.append((x0, tuple(signs)))
    return pts


# ------------------------------------------------------------------------------------------------ restricted game
def build_tree(q, k, F, cap=None):
    """greedy top-down decision tree on the low bits: classes c mod 2^d (c odd) are processed breadth first; each
    class is forced to one sign (MH sign first) if Min's restricted energy game at F stays finite everywhere,
    otherwise split by the next bit.  Returns (leaves [(c, d, sign)], strategy `minus`, number of tests).  The final
    strategy is an exploratory object: its rho_max is re-established by rhomax (witness + certificate)."""
    cap = cap or 20 * F.denominator + 200
    L = LIBR
    L.r_init(k, q)
    L.r_reset_values()
    H = 1 << (k - 1)
    if L.r_solve(F.numerator, F.denominator, cap) != 0:
        raise RuntimeError('infeasible at the start')
    mask = np.full(H, 3, dtype=np.uint8)
    leaves, queue, tests, head = [], [(1, 1)], 0, 0
    while head < len(queue):
        c, d = queue[head]
        head += 1
        if d >= 2:
            m1 = 1 if (q * c + 1) % 4 == 0 else 2
            tries = [m1, 3 - m1]
        else:
            tries = [1, 2]
        done = False
        for m in tries:
            L.r_save()
            L.r_get_mask(mask.ctypes.data)
            L.r_mask_class(c, d, m)
            tests += 1
            if L.r_solve(F.numerator, F.denominator, cap) == 0:
                leaves.append((c, d, '+' if m == 1 else '-'))
                done = True
                break
            L.r_restore()
            L.r_set_mask(mask.ctypes.data)
        if not done:
            if d >= k:
                raise RuntimeError('a single residue admits no sign')
            queue.append((c, d + 1))
            queue.append((c + (1 << d), d + 1))
    minus = np.zeros(H, dtype=np.uint8)
    L.r_strategy(minus.ctypes.data)
    return leaves, minus.astype(bool), tests


def tree_strategy(q, k, leaves):
    """the strategy of a list of leaves (c, d, sign) at level k >= all depths (the leaves must partition the odd
    residues); returns `minus`"""
    x = odd_nodes(k)
    out = np.zeros(len(x), dtype=np.int8) - 1
    for c, d, s in leaves:
        sel = (x % (1 << d)) == c
        if np.any(out[sel] >= 0):
            raise AssertionError('overlapping leaves')
        out[sel] = 1 if s == '-' else 0
    if np.any(out < 0):
        raise AssertionError('leaves do not cover')
    return out.astype(bool)
