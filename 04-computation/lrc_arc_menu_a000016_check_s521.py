#!/usr/bin/env python3
"""
lrc_arc_menu_a000016_check_s521.py   claudebox-2026-06-01-S521

Self-contained extension + closed-form check for HYP-1993.

(1) Recompute the geometric LRC arc menu for m movers (m=4..10) using the
    machinery in lrc_arc_menu_confirm_s521.py (exact difference-constraint
    feasibility + refinement-canonical iso count), and
(2) check it against the closed form
        A000016(m) = (1/2m) * Σ_{d|m, d odd} φ(d) 2^{m/d}.

A000016 = number of binary necklaces of length m under rotation + complement
(distinct output sequences of an m-stage shift register feeding back the
complement of the last stage).  Conjecture (HYP-1993): menu(m)=A000016(m), m>=4.
(m=3 is the degenerate L=1/2 boundary: menu(3)=1, A000016(3)=2.)

Run m=10 takes ~1 min; m>=11 is slow (skipped by default, pass an arg to include).
"""
from __future__ import annotations
import sys
from fractions import Fraction
from importlib import util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = util.spec_from_file_location("conf", os.path.join(HERE, "lrc_arc_menu_confirm_s521.py"))
C = util.module_from_spec(spec)
sys.argv = [sys.argv[0]]            # stop conf.main() from parsing our args
spec.loader.exec_module(C)

def phi(n):
    r, nn, p = n, n, 2
    while p * p <= nn:
        if nn % p == 0:
            while nn % p == 0: nn //= p
            r -= r // p
        p += 1
    if nn > 1: r -= r // nn
    return r

def A000016(m):
    s = sum(phi(d) * 2 ** (m // d) for d in range(1, m + 1) if m % d == 0 and d % 2 == 1)
    assert s % (2 * m) == 0, (m, s)
    return s // (2 * m)

def menu(m):
    L = Fraction(m - 1, m + 1)       # n=m+1, L=(n-2)/n=(m-1)/(m+1)
    feas = [S for S in C.upsets(m) if C.realizable(m, S, L)]
    raw = set(C.build_adj(m, S) for S in feas)
    classes = set(C.canon_fast(adj) for adj in raw)
    return len(feas), len(classes)

def main():
    mmax = 11 if (len(sys.argv) > 1 and sys.argv[1] == "--slow") else 10
    print("HYP-1993 closed-form check: menu(m) vs A000016(m)\n")
    print(f" {'m':>2} {'n':>2} {'#feasS':>7} {'2^(m-1)':>8} {'menu':>5} {'A000016':>8} {'match':>6}")
    for m in range(4, mmax + 1):
        nf, mn = menu(m)
        a = A000016(m)
        print(f" {m:>2} {m+1:>2} {nf:>7} {2**(m-1):>8} {mn:>5} {a:>8} {str(mn==a):>6}")
        sys.stdout.flush()
    print("\nA000016(m), m=3..13:", [A000016(m) for m in range(3, 14)])
    print("(menu(3)=1 is the L=1/2 boundary exception; A000016(3)=2.)")

if __name__ == "__main__":
    main()
