"""aud_c_selftest.py - the C counter aud_pairstate (both of its methods) against plain permutation enumeration
(aud_core.brute_state) on random pairs of orders 1..8, and against known transitive values at orders 12..20."""
import random, subprocess, os, sys
from aud_core import *

HERE = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.join(HERE, "aud_pairstate.exe")
rng = random.Random(4242)


def rec(mode, A, O):
    n = len(A)
    return " ".join([str(mode), str(n)] + [str(sum(1 << j for j in range(n) if A[i][j])) for i in range(n)]
                    + [str(sum(1 << o for o in O))])


cases = []
for n in range(1, 9):
    for _ in range(60 if n <= 7 else 15):
        A = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                if rng.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
        O = [i for i in range(n) if rng.random() < 0.5]
        cases.append((A, O))
out = subprocess.run([EXE], input="\n".join(rec(1, A, O) for A, O in cases) + "\n", capture_output=True, text=True).stdout.split("\n")
bad = 0
for (A, O), o in zip(cases, out):
    a, b = o.split("|")
    a = tuple(map(int, a.split()))
    b = tuple(map(int, b.split()))
    bs = brute_state(A, O)
    if not (a[1:] == bs == b):
        bad += 1
        print("MISMATCH", A, O, o, bs)
print(("PASS  " if bad == 0 else "FAIL  ") + "C counter (recursion and literal generation) = permutation enumeration on %d random pairs of orders 1..8" % len(cases))

# transitive controls at larger orders: s(TT_n, {p}) = (1, 2^(n-p), 1, [p = n]); s(TT_a, {i, i+delta}) = (r_X, (2^delta+1) 2^m, 1, [m = 0])
lines = []
exp = []
for n in (12, 14, 16, 18, 20):
    for p in (1, 3, n):
        lines.append(rec(0, transitive(n), [p - 1]))
        exp.append((n, 1, 2 ** (n - p), 1, 1 if p == n else 0))
    for (i, delta) in ((1, 1), (2, 2), (3, 3), (n - 2, 2), (n - 1, 1)):
        m = n - i - delta
        lines.append(rec(0, transitive(n), [i - 1, i - 1 + delta]))
        exp.append((n, -1 if delta == 1 else 1 - 2 ** (delta - 1), (2 ** delta + 1) * 2 ** m, 1, 1 if m == 0 else 0))
out = subprocess.run([EXE], input="\n".join(lines) + "\n", capture_output=True, text=True).stdout.split("\n")
bad = sum(1 for e, o in zip(exp, out) if tuple(map(int, o.split())) != e)
print(("PASS  " if bad == 0 else "FAIL  ") + "C counter on transitive pairs of orders 12..20 (%d cases)" % len(exp))
