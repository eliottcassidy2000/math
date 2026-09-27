#!/usr/bin/env python3
"""Orchestrator audit of lane `cnpos` (positive Hamiltonicity of the Collatz alphabet C_n). The lane's
construction is used only as a black box to PRODUCE candidate paths (certificates); every path is then verified by
this script's own checker, written from the definition of C_n (x ~ y iff x + y is a power of 2 or of 3):
a permutation of 1..n with every consecutive sum a power of 2 or 3.

  1. Clean n = T1 - 1 at the levels a = 2..11 (n = 7 .. 131071): paths verified.
  2. Window right ends (Conjecture B4): C_(3^a - 1) at a = 6, 7, 9 (two-power levels) and C_(3^a) at a = 8, 10
     (one-power levels), where the construction gives them.
  3. Theorem F's reduction is exercised by the construction; in addition the modulus m equals |3^a - 2^k|.
"""
import os, sys, time
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)
    print("  ok:", msg, flush=True)


def is_pow(x, b):
    while x % b == 0:
        x //= b
    return x == 1


def target(s):
    return s >= 2 and (is_pow(s, 2) or is_pow(s, 3))


def verify(n, seq):
    if len(seq) != n or sorted(seq) != list(range(1, n + 1)):
        return False
    return all(target(seq[i] + seq[i + 1]) for i in range(n - 1))


sys.argv = [sys.argv[0]]
src = open(os.path.join(here, "procgen_cnpos_20260926_run.py")).read()
cut = src.index("# B.1 the clean right ends at every level a <= AMAX")
ns = {"__name__": "cnpos_run_defs", "__file__": os.path.join(here, "procgen_cnpos_20260926_run.py")}
exec(compile(src[:cut], "cnpos_run_defs", "exec"), ns)
construct = ns["construct"]
C = ns["C"]

done = []
for a in range(2, 12):
    q = ns["clean_levels"](a)
    T1 = min(3 ** a, 2 ** q)
    n = T1 - 1
    t = time.time()
    seq = construct(n, want_seq=True)
    seq = seq[2] if isinstance(seq, tuple) and seq[0] == 'PATH' else None
    assert seq is not None and verify(n, list(seq)), (a, n)
    m = C.reduced_data(n)["m"]
    assert m == abs(3 ** a - 2 ** q) or m in (abs(3 ** a - 2 ** (q - 1)), abs(3 ** a - 2 ** (q + 1))), (a, n, m)
    done.append((a, n, round(time.time() - t, 1)))
check(True, "clean n = T1 - 1 verified Hamiltonian by the orchestrator's own checker: " + ", ".join(f"a={a}: n={n}" for a, n, _ in done))
ends = []
for a, n in ((6, 3 ** 6 - 1), (7, 3 ** 7 - 1), (9, 3 ** 9 - 1), (8, 3 ** 8), (10, 3 ** 10)):
    if n == 3 ** a:
        # one-power level: the leaf 3^a hangs on m = 2^q - 3^a; a path of C_(3^a - 1) ending at m extends
        n1 = n - 1
        mm = C.reduced_data(n1)["m"]
        res = construct(n1, end_at=mm, want_seq=True)
        seq = list(res[2]) + [n] if res[0] == "PATH" else None
    else:
        res = construct(n, want_seq=True)
        seq = list(res[2]) if res[0] == "PATH" else None
    if seq is not None and verify(n, seq):
        ends.append(n)
check(len(ends) == 5, f"window right ends verified Hamiltonian by the orchestrator's own checker (Conjecture B4 instances): {ends}")
