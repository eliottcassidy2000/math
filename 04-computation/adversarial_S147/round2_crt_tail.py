#!/usr/bin/env python3
"""Round 2: CRT-tuned tail on 12-element cores (the survivor architecture).

core C (12 distinct values, heights <= 60) + tail h >= 65.
Constraints assembled exactly:
  A. covering 2..13 for C u {h}   (h may carry leftover duties)
  B. [14,54] stack for C u {h}: for each q, either some element ≡ 0 (mod q)
     or depth-pinning at d_gap(q).  For fixed C this defines the set
     R(q) = { r in [0,q) : C u {r} passes at q }  -> h mod q must be in R(q).
  C. tail-band safety: for every v in C, the breakpoint moduli h+v (and h-v
     when >= 55) must satisfy W_{C u {h}}(q) <= ceil(3q/41)-1; likewise 2h
     (and h itself as modulus? h appears as modulus only if h = pair sum of
     others - impossible since h > max pair? no: pair sums of C reach 110+;
     handled by full exact-M at the end anyway).
  D. exact M of the full 13-set for every candidate that clears A-C.
Cores: all 12-subsets of each atlas survivor (59 x 13, deduped), plus the
known near-floor cores {1..11,13} etc.
"""
import sys
from math import gcd
from fractions import Fraction
from itertools import combinations
sys.path.insert(0, ".")
from lrc_engine import (exact_M, covering_check, d_gap, W_of_q, THR, FLOOR)

def pass_at_q(vals, q):
    """cover-or-pin test at modulus q for the value multiset."""
    res = [v % q for v in vals]
    if 0 in res:
        return True
    d = d_gap(q)
    if d <= 0:
        return False
    for a in range(1, q):
        if gcd(a, q) != 1:
            continue
        m = min(min((a * r) % q, q - (a * r) % q) for r in res)
        if m > d:
            return False
    return True

def residue_options(C, q):
    """R(q): residues r mod q such that C u {r} passes at q."""
    ok = []
    for r in range(q):
        if pass_at_q(list(C) + [r], q):
            ok.append(r)
    return ok

def leftover_cover(C):
    return [q for q in range(2, 14) if all(v % q for v in C)]

def core_tail_candidates(C, hmin=65, hmax=2600):
    """Yield h satisfying A (covering) and B ([14,54] via residue sets)."""
    need = leftover_cover(C)
    # build per-q admissible residue sets for q in [14,54]
    Rs = {}
    for q in range(14, 55):
        R = residue_options(C, q)
        if not R:
            return  # core dead: no tail residue can fix q
        Rs[q] = set(R)
    for h in range(hmin, hmax + 1):
        if h in C:
            continue
        if any(h % q for q in need):
            continue
        ok = True
        for q, R in Rs.items():
            if h % q not in R:
                ok = False
                break
        if ok:
            yield h

def tailband_ok(C, h):
    """Check C u {h} at moduli h+v, |h-v| (>=55), 2h exactly."""
    V = sorted(set(C) | {h})
    for v in C:
        for q in (h + v, h - v):
            if q >= 55:
                W, k = W_of_q(V, q)
                if W > d_gap(q):
                    return False, (q, k, W)
    q = 2 * h
    W, k = W_of_q(V, q)
    if W > d_gap(q):
        return False, (q, k, W)
    return True, None

def main():
    SURV = [tuple(int(x) for x in l.split(":")[1].split())
            for l in open("survivors.txt") if ":" in l]
    cores = set()
    for S in SURV:
        for i in range(13):
            cores.add(tuple(sorted(S[:i] + S[i+1:])))
    # known near-floor 12-cores as extra seeds
    cores.add(tuple(list(range(1, 12)) + [13]))
    cores.add(tuple(range(1, 13)))
    cores.add(tuple(list(range(1, 11)) + [12, 13]))
    cores = sorted(cores)
    print(f"{len(cores)} candidate 12-cores")

    n_dead = n_hB = n_hC = 0
    results = []
    for ci, C in enumerate(cores):
        got = False
        for h in core_tail_candidates(C):
            got = True
            n_hB += 1
            ok, kill = tailband_ok(C, h)
            if not ok:
                continue
            n_hC += 1
            V = sorted(set(C) | {h})
            if len(V) != 13:
                continue
            g = 0
            for v in V:
                g = gcd(g, v)
            if g != 1 or covering_check(V):
                continue
            M, (q, k) = exact_M(V)
            results.append((M, tuple(V), q, k, h))
            if FLOOR < M < THR:
                print(f"*** IN GAP *** {V} M={M}")
        if not got:
            n_dead += 1
        if (ci + 1) % 100 == 0:
            print(f"  ..{ci+1} cores, {n_hB} AB-tails, {n_hC} ABC-tails, "
                  f"{len(results)} exact-M evals", flush=True)
    print(f"cores: {len(cores)}, dead(A/B infeasible): {n_dead}, "
          f"tails passing A+B: {n_hB}, passing A+B+C: {n_hC}")
    results.sort()
    seen = set()
    print("\nbest exact-M results:")
    for M, V, q, k, h in results[:40]:
        if V in seen:
            continue
        seen.add(V)
        tag = ("*** IN GAP ***" if FLOOR < M < THR else
               ("TIGHT 1/14" if M == FLOOR else ""))
        print(f"M={str(M):>8} ={float(M):.6f} wit {k}/{q} tail={h} V={list(V)} {tag}")

if __name__ == "__main__":
    main()
