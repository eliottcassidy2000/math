#!/usr/bin/env python3
"""Adversarial LRC(14) gap-engineering engine.

Exact M(V) = max_t min_i ||v_i t|| for 13-speed families via the
breakpoint-moduli method: the maximizer is at t = k/q with
q in {v_i+v_j} u {|v_i-v_j|} u {2 v_i}.  All arithmetic exact (integers,
Fraction only at the comparison layer).

Diagnosis: every escaping witness (q, k) with min_i dist >= 3q/41 among
breakpoint moduli, plus the pinning profile for design moduli.
"""
import numpy as np
from fractions import Fraction
from math import gcd
from itertools import combinations

THR = Fraction(3, 41)   # gap ceiling
FLOOR = Fraction(1, 14)

def breakpoint_moduli(V):
    Q = set()
    for a, b in combinations(V, 2):
        Q.add(a + b)
        if a != b:
            Q.add(abs(a - b))
    for v in V:
        Q.add(2 * v)
    Q.discard(0)
    return sorted(Q)

def W_of_q(V, q):
    """max over k in [1,q-1] of min_i dist_q(k*v_i); returns (W, argk).
    Exact: integer arithmetic."""
    if q < 2:
        return 0, 0
    v = np.array([x % q for x in V], dtype=np.int64)
    ks = np.arange(1, q, dtype=np.int64)
    R = (v[:, None] * ks[None, :]) % q
    D = np.minimum(R, q - R)
    m = D.min(axis=0)
    j = int(m.argmax())
    return int(m[j]), int(ks[j])

def exact_M(V, extra_moduli=()):
    """Exact M as Fraction, plus witness (q,k). Includes optional extra moduli
    (harmless: every t=k/q lower-bounds M; breakpoints already attain it)."""
    best = Fraction(0, 1); wit = (1, 0)
    for q in sorted(set(breakpoint_moduli(V)) | set(extra_moduli)):
        W, k = W_of_q(V, q)
        f = Fraction(W, q)
        if f > best:
            best, wit = f, (q, k)
    return best, wit

def escapes(V, thr=THR, qmax=None):
    """All breakpoint moduli q carrying a witness with min dist/q >= thr.
    Returns list of (q, k, W, Fraction)."""
    out = []
    for q in breakpoint_moduli(V):
        if qmax and q > qmax:
            continue
        W, k = W_of_q(V, q)
        if Fraction(W, q) >= thr:
            out.append((q, k, W, Fraction(W, q)))
    return out

def covering_check(V):
    missing = [q for q in range(2, 14) if all(v % q for v in V)]
    return missing

def pin_profile(V, q):
    """For modulus q: (covered?, worst depth over units, worst unit)."""
    if any(v % q == 0 for v in V):
        return True, 0, None
    worst, wa = -1, None
    for a in range(1, q):
        if gcd(a, q) != 1:
            continue
        m = min(min((a * v) % q, q - (a * v) % q) for v in V)
        if m > worst:
            worst, wa = m, a
    return False, worst, wa

def d_gap(q):
    """Pinning depth ceiling for M < 3/41: need min dist <= ceil(3q/41)-1."""
    return -(-3 * q // 41) - 1

def stack_report(V, qlo=14, qhi=48):
    """Design-stack report over [qlo,qhi]: which moduli covered, which pinned,
    which VIOLATED (would push M >= 3/41)."""
    bad = []
    for q in range(qlo, qhi + 1):
        cov, worst, wa = pin_profile(V, q)
        if not cov and worst > d_gap(q):
            bad.append((q, worst, wa))
    return bad

def report(name, V):
    V = sorted(V)
    assert len(set(V)) == 13, f"{name}: need 13 distinct"
    miss = covering_check(V)
    M, (q, k) = exact_M(V)
    fars = sum(1 for v in V if v > 14)
    pos = ("IN GAP (1/14,3/41) *** COUNTEREXAMPLE ***" if FLOOR < M < THR
           else ("BELOW FLOOR *** LRC(14) COUNTEREXAMPLE ***" if M < FLOOR
                 else ("== 1/14 (tight)" if M == FLOOR else f">= 3/41 (fails, M-3/41 = {float(M-THR):+.5f})")))
    print(f"[{name}] V={V}")
    print(f"    M = {M} = {float(M):.6f}  witness t = {k}/{q}   {pos}")
    print(f"    fars(>14) = {fars}, v_max = {V[-1]}, covering-missing = {miss}")
    return M, (q, k)

if __name__ == "__main__":
    # THM-1290 referee gates
    gates = [
        ("ladder", list(range(1, 14)), Fraction(1, 14)),
        ("1..12,14", list(range(1, 13)) + [14], Fraction(1, 13)),
        ("1..12,26", list(range(1, 13)) + [26], Fraction(2, 27)),
        ("1..11,13,24", list(range(1, 12)) + [13, 24], Fraction(1, 14)),
        ("1..11,13,36", list(range(1, 12)) + [13, 36], Fraction(3, 41)),
        ("8/105-family", [1, 2, 3, 5, 7, 8, 9, 10, 11, 12, 17, 19, 104], Fraction(8, 105)),
        ("deep-well 1..12,182", list(range(1, 13)) + [182], Fraction(14, 183)),
    ]
    ok = True
    for name, V, expect in gates:
        M, wit = exact_M(V)
        s = "OK " if M == expect else "FAIL"
        if M != expect:
            ok = False
        print(f"gate {s} {name}: M = {M} (expect {expect}) witness {wit[1]}/{wit[0]}")
    print("ALL GATES PASS" if ok else "GATE FAILURE — engine unsound")
