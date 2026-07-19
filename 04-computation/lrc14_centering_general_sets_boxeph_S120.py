#!/usr/bin/env python3
"""
boxeph-2026-07-18-S120 — the centering witness on GENERAL 12-sets: probing uniqueness of
the loneliness minimizer {1,...,12} (Tao's n=12 optimistic conjecture / INVcov).

BACKGROUND. LRC(13) (proved for the project): every 12 distinct positive speeds have
M = max_t min_i ||v_i t|| >= 1/13, with equality (conjecturally, HYP-4382 / Tao n=12)
ONLY for dilates of {1,...,12}.  So {1,...,12} is the UNIQUE minimizer of loneliness.
Uniqueness <=> every 12-set C != c*{1,...,12} has some t with min-dist > 1/13 -- a WITNESS.

S118's centering witness produced such a t for APs via q = v_min + v_max (a pairwise SUM) and
t = d^{-1}/q.  THIS SCRIPT tests whether structured moduli -- pairwise SUMS q = v_i+v_j --
supply a beat-1/13 witness for NON-AP sets, i.e. whether uniqueness reduces to a finite
per-set check over pairwise sums.

CONJECTURE UNDER TEST (centering rigidity):
  For every 12-set C that is not a dilate of {1,...,12}, some pairwise-sum modulus q=v_i+v_j
  yields a witness t=p/q with min_k ||v_k t|| > 1/13.  {1,...,12} alone has NO such witness.

We compute, for each C:
  - true M (exact, tent-pole denominators),
  - bestSum = max over q in {v_i+v_j} and p of min_k ||v_k (p/q)||  (the pairwise-sum witness),
  - which q, and whether bestSum > 1/13,
and test the conjecture on: {1,...,12}; all single-element perturbations; reflective
non-AP sets; and random 12-sets.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations


def frac_dist(x: F) -> F:
    r = x - (x.numerator // x.denominator)
    return min(r, 1 - r)


def fmin(speeds, t: F) -> F:
    return min(frac_dist(F(v) * t) for v in speeds)


def true_M(speeds):
    Ds = set()
    n = len(speeds)
    for i in range(n):
        for j in range(i + 1, n):
            Ds.add(speeds[i] + speeds[j])
            Ds.add(abs(speeds[i] - speeds[j]))
    best, bestt = F(0), F(0)
    for D in Ds:
        if D <= 0:
            continue
        for m in range(1, D):
            val = fmin(speeds, F(m, D))
            if val > best:
                best, bestt = val, F(m, D)
    return best, bestt


def best_sum_witness(speeds):
    """max over q = v_i+v_j (pairwise sums) and residues p of min_k ||v_k p/q||."""
    qs = sorted({speeds[i] + speeds[j] for i in range(len(speeds))
                 for j in range(i + 1, len(speeds))})
    best, bq, bp = F(0), None, None
    for q in qs:
        for p in range(1, q):
            if gcd(p, q) != 1:
                continue
            val = fmin(speeds, F(p, q))
            if val > best:
                best, bq, bp = val, q, p
    return best, bq, bp


def best_any_witness_q(speeds, qmax):
    """max over ALL q<=qmax and residues p (to see if non-sum moduli ever beat sums)."""
    best, bq = F(0), None
    for q in range(2, qmax + 1):
        for p in range(1, q):
            if gcd(p, q) != 1:
                continue
            val = fmin(speeds, F(p, q))
            if val > best:
                best, bq = val, q
    return best, bq


THIRT = F(1, 13)
print("=" * 78, flush=True)
print("(A) BASELINE {1,...,12}: is it uniquely witness-less at 1/13?", flush=True)
print("=" * 78, flush=True)
base = list(range(1, 13))
M, t = true_M(base)
bs, bq, bp = best_sum_witness(base)
print(f"  {base}", flush=True)
print(f"  true M = {M} at t={t};  best pairwise-sum witness = {bs} (q={bq},p={bp})", flush=True)
print(f"  M == 1/13: {M == THIRT};  best-sum-witness > 1/13: {bs > THIRT}", flush=True)

print("", flush=True)
print("=" * 78, flush=True)
print("(B) SINGLE-ELEMENT PERTURBATIONS of {1,...,12}: does a pairwise-sum q beat 1/13?", flush=True)
print("=" * 78, flush=True)
print(f"{'replace':>10} {'->':>3} {'newset (short)':>22} {'trueM':>8} {'bestSum':>8} "
      f"{'q*':>5} {'M>1/13':>7} {'wit>1/13':>8}", flush=True)
fails_conj = []
tested = 0
for m in range(1, 13):
    for mp in list(range(1, 30)):
        if mp in base and mp != m:
            continue
        if mp == m:
            continue
        C = sorted(set(base) - {m} | {mp})
        if len(C) != 12:
            continue
        g = 0
        for v in C:
            g = gcd(g, v)
        # skip dilates check: C is a dilate of {1..12} only if C == c*{1..12}; single-move never is
        M, _ = true_M(C)
        bs, bq, bp = best_sum_witness(C)
        tested += 1
        witok = bs > THIRT
        Mok = M > THIRT
        if not witok:
            fails_conj.append((m, mp, C, M, bs, bq))
        if mp <= 16 and (mp > 12 or m in (1, 6, 12)):
            short = "{" + ",".join(str(x) for x in C[:2]) + ",..," + str(C[-1]) + "}"
            print(f"{m:>10} {'->':>3} {short:>22} {str(M):>8} {str(bs):>8} "
                  f"{str(bq):>5} {str(Mok):>7} {str(witok):>8}", flush=True)
print(f"\n  tested {tested} single-element perturbations", flush=True)
print(f"  perturbations where NO pairwise-sum witness beats 1/13: {len(fails_conj)}", flush=True)
for f in fails_conj[:15]:
    print(f"    replace {f[0]}->{f[1]}: set={f[2]} trueM={f[3]} bestSum={f[4]} (all q tried)", flush=True)

print("", flush=True)
print("=" * 78, flush=True)
print("(C) For the CONJECTURE FAILURES (if any), does a NON-sum modulus witness beat 1/13?", flush=True)
print("=" * 78, flush=True)
if fails_conj:
    for f in fails_conj[:20]:
        C = f[2]
        ba, bqa = best_any_witness_q(C, 3 * max(C))
        print(f"    replace {f[0]}->{f[1]} {C}: trueM={f[3]}  best-ANY-q witness={ba} at q={bqa}  "
              f">1/13={ba > THIRT}", flush=True)
else:
    print("    (none) — every single-element perturbation has a pairwise-sum witness > 1/13.", flush=True)

print("", flush=True)
print("=" * 78, flush=True)
print("(D) REFLECTIVE non-AP sets (symmetric v_i+v_{13-i}=const but NOT an AP)", flush=True)
print("=" * 78, flush=True)
# build symmetric sets: pick 6 values, mirror them around a center S
refl_tests = [
    ([1, 2, 3, 4, 5, 6], 13),     # -> {1..12} (AP, control)
    ([1, 2, 3, 4, 5, 7], 15),     # 7 mirrors to 8; {1,2,3,4,5,7,8,10,11,12,13,14}
    ([1, 2, 3, 4, 6, 7], 15),
    ([1, 2, 4, 5, 7, 8], 17),
    ([1, 3, 4, 6, 7, 9], 19),
]
print(f"{'set (short)':>28} {'reflective?':>11} {'AP?':>5} {'trueM':>8} {'bestSum':>8} {'q*':>5} {'wit>1/13':>8}", flush=True)
for half, S in refl_tests:
    C = sorted(set(half) | {S - x for x in half})
    if len(C) != 12:
        continue
    reflective = all((C[i] + C[11 - i]) == C[0] + C[11] for i in range(6))
    diffs = {C[i + 1] - C[i] for i in range(11)}
    isap = (len(diffs) == 1)
    M, _ = true_M(C)
    bs, bq, bp = best_sum_witness(C)
    short = "{" + ",".join(str(x) for x in C) + "}"
    if len(short) > 27:
        short = short[:24] + "..}"
    print(f"{short:>28} {str(reflective):>11} {str(isap):>5} {str(M):>8} {str(bs):>8} "
          f"{str(bq):>5} {str(bs > THIRT):>8}", flush=True)

print("", flush=True)
print("=" * 78, flush=True)
print("(E) RANDOM-ish 12-sets (structured spread): pairwise-sum witness vs true M", flush=True)
print("=" * 78, flush=True)
rnd = [
    [1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233],   # Fibonacci-ish
    [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048],  # powers of 2
    [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37],    # primes
    [1, 4, 9, 16, 25, 36, 49, 64, 81, 100, 121, 144],  # squares
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13],         # near-AP (one off)
]
print(f"{'family':>18} {'trueM':>10} {'bestSum':>10} {'q*':>6} {'sum=trueM?':>10} {'wit>1/13':>8}", flush=True)
names = ["fib", "pow2", "primes", "squares", "near-AP"]
for nm, C in zip(names, rnd):
    M, _ = true_M(C)
    bs, bq, bp = best_sum_witness(C)
    print(f"{nm:>18} {str(M):>10} {str(bs):>10} {str(bq):>6} {str(bs == M):>10} "
          f"{str(bs > THIRT):>8}", flush=True)

print("\nDONE.", flush=True)
