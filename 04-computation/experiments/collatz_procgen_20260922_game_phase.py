#!/usr/bin/env python3
"""Mantissa-phase law for Althofer's 3n+-1 game.

Exact identity (any play): if n = n_0 -> n_1 -> ... -> n_T = 1 with n_{t+1} = (3 n_t + e_t)/2^{j_t},
e_t = +-1, then
    log2 n + T log2 3 = S - sum_t log2(1 + e_t/(3 n_t)),     S = sum_t j_t  (an integer),
so the phase  phi = frac(log2 n + T log2 3)  equals  frac(-sum_t log2(1 + e_t/(3 n_t))).
Every move advances frac(log2 n) by log2 3 (mod 1) up to a correction of size about 1/(3 n_t).

This script reads the remoteness table written by collatz_procgen_20260922_game_remote.c
(rem_<CAP>.bin: uint64 CAP, then uint16 per odd non-multiple of 3, 0 = unresolved, else rem+1),
follows the principal variation (winner: fastest win, ties -> descending move; loser: slowest loss,
ties -> descending move), and reports per octave:
  * min/max of phi(n) = frac(log2 n + T*(n) log2 3) (signed to (-1/2,1/2]);
  * the most frequent endgames (last three positions of the principal variation);
  * how well the value (parity of T*) is predicted by theta = frac(log2 n) alone
    (majority vote in 64/256/1024 bins).
usage: game_phase.py rem_file [kmin kmax]
"""
import sys, math
from collections import Counter
import numpy as np

L3 = math.log2(3)


def load(fn):
    raw = open(fn, 'rb').read()
    cap = int.from_bytes(raw[:8], 'little')
    R = np.frombuffer(raw[8:], dtype=np.uint16)
    return cap, R


def main():
    fn = sys.argv[1]
    cap, R = load(fn)
    kmin = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    kmax = int(sys.argv[3]) if len(sys.argv) > 3 else 40

    def oddpart(x):
        while x % 2 == 0:
            x //= 2
        return x

    def dch(n):
        return oddpart(3 * n + 1) if n % 4 == 1 else oddpart(3 * n - 1)

    def uch(n):
        return (3 * n - 1) // 2 if n % 4 == 1 else (3 * n + 1) // 2

    def rem(y):
        if y == 1:
            return 0
        if y > cap:
            return -1
        r = int(R[y // 3])
        return r - 1 if r else -1

    def rem_start(n):
        if n == 1:
            return 1
        if n % 3:
            return rem(n)
        ra, rb = rem(dch(n)), rem(uch(n))
        ev = [r for r in (ra, rb) if r >= 0 and r % 2 == 0]
        if ev:
            return 1 + min(ev)
        if ra >= 0 and rb >= 0:
            return 1 + max(ra, rb)
        return -1

    def pv(n):
        path = [n]
        x = n
        while x != 1:
            a, b = dch(x), uch(x)
            ra, rb = rem(a), rem(b)
            rx = rem_start(x)
            if rx % 2 == 1:
                pa = (a == 1) or (ra >= 0 and ra % 2 == 0)
                pb = rb >= 0 and rb % 2 == 0
                if a == 1:
                    nx = a
                elif pb and (not pa or rb < ra):
                    nx = b
                else:
                    nx = a
            else:
                nx = b if rb > ra else a
            path.append(nx)
            x = nx
            if len(path) > 2000:
                break
        return path

    # fully resolved prefix
    B = 3
    while rem_start(B) >= 0:
        B += 2
    print(f"# rem table cap={cap}; every odd start < {B} has a capped remoteness")
    print("# octave k: [min phi, max phi] over odd n in [2^k,2^(k+1)); top endgames (last 3 PV positions);"
          " accuracy of theta-bin majority predictor (64/256/1024 bins)")
    allphi_min, allphi_max = 1, -1
    for k in range(kmin, kmax):
        lo, hi = 1 << k, 1 << (k + 1)
        if hi > B:
            break
        ns = np.arange(lo + 1, hi, 2, dtype=np.int64)
        T = np.array([rem_start(int(n)) for n in ns])
        th = np.log2(ns.astype(np.float64)) % 1.0
        phi = (np.log2(ns.astype(np.float64)) + T * L3) % 1.0
        phi = np.where(phi > 0.5, phi - 1, phi)
        allphi_min = min(allphi_min, phi.min())
        allphi_max = max(allphi_max, phi.max())
        P = (T % 2 == 0)
        accs = []
        for bins in (64, 256, 1024):
            b = np.minimum((th * bins).astype(int), bins - 1)
            cntP = np.bincount(b, weights=P, minlength=bins)
            cnt = np.bincount(b, minlength=bins)
            maj = cntP * 2 >= cnt
            correct = np.where(maj[b], P, ~P).mean()
            accs.append(correct)
        # endgames on a sample
        ends = Counter()
        step = max(1, len(ns) // 4000)
        for n in ns[::step]:
            p = pv(int(n))
            ends[tuple(p[-4:-1])] += 1
        top = ', '.join(f"{'>'.join(map(str, e))}:{c}" for e, c in ends.most_common(4))
        print(f"k={k:2d}: phi in [{phi.min():+.4f},{phi.max():+.4f}]  P={P.mean():.4f}  "
              f"theta-predictor acc {accs[0]:.4f}/{accs[1]:.4f}/{accs[2]:.4f}  endgames: {top}")
    print(f"# over all octaves checked: phi in [{allphi_min:+.4f},{allphi_max:+.4f}]")
    # endgame corrections: exact values of -sum log2(1+e/(3 n_t)) for common endings
    print("# exact phase contributions log2(3n/(3n+e)) of small moves (descending and ascending):")
    for n in (3, 5, 7, 9, 11, 13, 17, 19, 21, 23, 43, 85):
        for kind, e in (("descending", 1 if n % 4 == 1 else -1), ("ascending", -1 if n % 4 == 1 else 1)):
            num = 3 * n + e
            print(f"   {n} -> {oddpart(num)} ({kind}, 3n{'+' if e > 0 else '-'}1={num}): "
                  f"{math.log2(3 * n / num):+.5f}")
    print("#   e.g. endgame 19->7->5->1 contributes %+.5f" % (math.log2(57 / 56) + math.log2(21 / 20) + math.log2(15 / 16)))


if __name__ == '__main__':
    main()
