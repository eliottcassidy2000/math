#!/usr/bin/env python3
"""procgen_petersen_20261001_census26.py -- census of all-odd anti-circulant tournaments on Z/26 (m = 13).

An anti-circulant tournament T_s on Z/2m (m odd) has x -> y iff s(y - x) = (-1)^x, where s(-d) = -(-1)^d s(d).
All 2^13 sign patterns are reduced modulo multipliers u in (Z/26)^* and global negation (both give isomorphic
tournaments); every orbit representative is tested with the one-array parity engine procgen_petersen_20261001_anti.c
(256 MB, about 3 s each), the all-odd ones are canonicalised with nauty, and the prime-derived classes QR_p[mu26]
(p = 3 mod 4, p < 3000) and QR_27 - 0 are located among them.
Prints to stdout only; ends with 'ALL CHECKS PASSED'. Runtime about 20 minutes, RSS about 270 MB.
"""
import itertools
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import procgen_petersen_20261001_lib as L  # noqa: E402

T0 = time.time()
WORK = os.path.join(os.path.dirname(os.path.dirname(HERE)), "scratch", "procgen_petersen", "build")
os.makedirs(WORK, exist_ok=True)
ANTI = os.path.join(WORK, "procgen_petersen_anti")
subprocess.run(["gcc", "-O2", "-o", ANTI, os.path.join(HERE, "procgen_petersen_20261001_anti.c")], check=True)
NCHECK = 0


def check(cond, msg):
    global NCHECK
    NCHECK += 1
    if not cond:
        print("CHECK FAILED:", msg)
        sys.exit(1)


m, n = 13, 26
reps = {}
for free in itertools.product((1, -1), repeat=m):
    s = L.anti_full_s(m, free)
    reps.setdefault(L.anti_pattern_key(m, s), s)
print("N = 26: %d sign patterns, %d orbit representatives under multipliers and negation" % (2 ** m, len(reps)))
sys.stdout.flush()
found = []
for i, (key, s) in enumerate(sorted(reps.items())):
    A = L.anti_tournament(m, s)
    H2, par = L.anti_engine(ANTI, A)
    check(H2 == 1, "Redei")
    check(par[m][2] == 1, "Theorem 6.1: antipodal arcs odd")
    if all(p3 == 1 for (_, _, p3) in par.values()):
        found.append((key, A))
        print("   all-odd representative #%d: s(1..13) = %s" % (i, "".join("+" if s[d] == 1 else "-" for d in range(1, 14))))
        sys.stdout.flush()
cans = L.labelg([L.digraph6(n, L.tour_arcs(A)) for _, A in found])
classes = sorted(set(cans))
print("all-odd orbit representatives: %d; distinct isomorphism classes: %d" % (len(found), len(classes)))
A27 = L.f27_tournament()
c27 = L.labelg([L.digraph6(n, L.tour_arcs(A27))])[0]
check(c27 in classes, "QR27 - 0 is among the all-odd classes")
prime_classes = {}
for p in range(7, 3000):
    if p % 4 != 3 or not L.is_prime(p) or ((p - 1) // 2) % m:
        continue
    c = L.labelg([L.digraph6(n, L.tour_arcs(L.mu_tournament(p, m)))])[0]
    prime_classes.setdefault(c, []).append(p)
for c in classes:
    tag = []
    if c == c27:
        tag.append("QR27 - 0")
    if c in prime_classes:
        tag.append("QR_p[mu26] for p = %s" % prime_classes[c])
    print("   class %s  %s" % (c, "; ".join(tag) if tag else "(not QR27 - 0, not QR_p[mu26] for p < 3000)"))
print("prime-derived classes (p < 3000): %d, of which all-odd: %d" % (
    len(prime_classes), sum(1 for c in prime_classes if c in classes)))
print()
print("%d checks, %.0f s" % (NCHECK, time.time() - T0))
print("ALL CHECKS PASSED")
