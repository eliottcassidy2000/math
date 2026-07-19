#!/usr/bin/env python3
"""
kind-pasteur-2026-07-19-S128c87 -- HYP-7915: the FINITE-CHECK FEASIBILITY LEDGER
for LRC(14) via the Rosenfeld / Sungkawichai-Trakulthongchai architecture.

Method (from the primary sources, fetched this session):
  Rosenfeld arXiv:2509.14111 (k=7 speeds, 8 runners):
    - MSS Cor.3 product bound: any minimal counterexample has
      prod v_i <= B_k := (C(k+1,2)^(k-1)/k)^k     [k=7: < 7.4e54]
    - Lemma 4: lcm(2..k+1) divides the counterexample product.
    - Lemma 6/7: for prime p, IF no improper k-tuple exists mod (k+1)p
      (a finite covering check), THEN p | prod v_i of any counterexample.
    - Kill: verify enough primes that lcm(2..k+1)*prod(verified p) > B_k.
    - k=7: 27 primes 31..163; p=163 took ~32 h on one core.
  S-T arXiv:2604.23906 (k=10,11,12):
    - Same architecture + intermediate sieves/backward projection + mod-p
      symmetry (cost ~ p^k/(2^k (k-1)!) representatives) + polynomial
      shortcut (Prop 4.4) WHEN k+1 IS PRIME (k=10,12) avoiding the most
      expensive lifts; k=11 (k+1=12 composite) needed the full lift chain
      x2,x2,x2,x2,x3,x3 (cost c^k * |S| per lift).
    - k=12: 73 primes in [167,733]; ln(prod p) > 547 > 546 > ln(B_12);
      est. ~40 days on a 10-core Apple M4; heuristic per-prime cost
      p^((k+1)/2)/(k 2^k).
    - k=13 bottleneck (their Section 7): efficient computation of
      I(13,p,1) (the initial improper-tuple sieve); k+1=14=2*7 is
      COMPOSITE so the polynomial shortcut fails and the lift chain
      needs x7 steps (7^13 per surviving tuple).

This script computes the exact ledger arithmetic:
  (1) B_k for k=7,12,13 (verify k=7,12 against the papers' numbers);
  (2) the k=13 prime budget (greedy smallest usable primes until the
      product beats B_13, with and without the lcm(2..14) credit);
  (3) the cost extrapolation via S-T's own heuristic (ratio to their
      k=12 set = the calibrated 40-day baseline);
  (4) the x7-lift and I(13,p,1) wall arithmetic.
All exact integer arithmetic where possible; ln/log10 for display.
"""
from math import log, log10, comb, gcd
from functools import reduce

def primes_upto(n):
    sieve = bytearray([1])*(n+1)
    sieve[0:2] = b'\x00\x00'
    for i in range(2, int(n**0.5)+1):
        if sieve[i]:
            sieve[i*i::i] = b'\x00'*len(sieve[i*i::i])
    return [i for i in range(2, n+1) if sieve[i]]

PR = primes_upto(3000)

def lcm_range(a, b):
    v = 1
    for x in range(a, b+1):
        v = v*x//gcd(v, x)
    return v

def B(k):
    # MSS Cor.3: (C(k+1,2)^(k-1)/k)^k  -- keep exact as a fraction power
    num = comb(k+1, 2)**(k-1)
    return num, k  # represent as (num/k)^k

def lnB(k):
    num, den = B(k)
    return k*(log(num) - log(den))

print("== (1) MSS product bounds ==")
for k in (7, 12, 13):
    lb = lnB(k)
    print(f"  k={k:2d}: B_k = (C({k+1},2)^{k-1}/{k})^{k}"
          f"  ln B = {lb:9.2f}  log10 B = {lb/log(10):8.2f}")
print(f"  [verify] k=7: log10 = {lnB(7)/log(10):.2f} (paper: < 7.4e54 -> 54.87)  "
      f"k=12: ln = {lnB(12):.1f} (paper: < 546)")
L14 = lcm_range(2, 14)
print(f"  lcm(2..14) = {L14}  (ln = {log(L14):.3f})")

print("\n== (2) the k=13 prime budget ==")
# k=12 reconstruction: greedy smallest primes >= 167 until sum ln p > ln B_12
target12 = lnB(12)
ps12, s = [], 0.0
for p in PR:
    if p < 167: continue
    ps12.append(p); s += log(p)
    if s > target12: break
print(f"  k=12 reconstruction: greedy primes from 167: count={len(ps12)}, "
      f"largest={ps12[-1]}, sum ln p = {s:.1f} (paper: 73 primes to 733, >547)")

# k=13: polynomial shortcut needs p > k^2+k; but k+1=14 composite so the
# shortcut is unavailable anyway -- keep the same lower-cutoff convention
# (first usable prime after the k^2+k analogue) for comparability, and show
# sensitivity to the cutoff.
target13 = lnB(13)
for cutoff, note in [(191, "first prime > k^2+k = 182"),
                     (167, "same cutoff as k=12 (optimistic)"),
                     (400, "pessimistic: small primes fail the improper check")]:
    ps13, s = [], 0.0
    for p in PR:
        if p < cutoff: continue
        ps13.append(p); s += log(p)
        if s > target13 - log(L14):  # credit lcm(2..14) as Rosenfeld does
            break
    print(f"  k=13 cutoff {cutoff} ({note}): count={len(ps13)}, largest={ps13[-1]}, "
          f"sum ln p = {s:.1f} (target {target13 - log(L14):.1f} after lcm credit)")

# canonical set for the cost extrapolation: cutoff 191
ps13, s = [], 0.0
for p in PR:
    if p < 191: continue
    ps13.append(p); s += log(p)
    if s > target13 - log(L14): break

print("\n== (3) cost extrapolation (S-T's own heuristic p^((k+1)/2)/(k 2^k)) ==")
def cost(p, k):
    return p**((k+1)/2) / (k * 2**k)
c12 = sum(cost(p, 12) for p in ps12)
c13 = sum(cost(p, 13) for p in ps13)
R = c13/c12
days13 = 40*R
print(f"  k=12 baseline: {len(ps12)} primes, heuristic mass {c12:.3e}  == ~40 days on 10-core M4")
print(f"  k=13 target:   {len(ps13)} primes to {ps13[-1]}, heuristic mass {c13:.3e}")
print(f"  RATIO = {R:.1f}x  ->  ~{days13:,.0f} days on the same 10-core M4")
print(f"                     = {days13/365:.1f} machine-years (10-core)")
print(f"                     = {days13*10/365:.0f} core-years")
print(f"  at 1,000 cores:  ~{days13*10/1000:.0f} days;  at 10,000 cores: ~{days13*10/10000:.1f} days")

print("\n== (4) the two structural walls the heuristic does NOT price ==")
sevenlift = 7**13
print(f"  (a) x7 lift (k+1=14=2*7 composite; the k=10/12 polynomial shortcut needs")
print(f"      k+1 prime and FAILS): one x7 lift costs 7^13 = {sevenlift:,} (~{sevenlift:.2e})")
print(f"      tuple-checks per surviving tuple, vs 2^13 = {2**13:,} for a x2 lift")
print(f"      ratio (7/2)^13 = {(7/2)**13:.2e}. Feasible ONLY if the surviving set |S|")
print(f"      before the x7 step is tiny: |S| * 7^13 <= 1e15 requires |S| <= {1e15/sevenlift:,.0f}.")
p12max, p13max = 733, ps13[-1]
i12 = 12*log10(p12max) - log10(2**12) - log10(reduce(lambda a,b: a*b, range(1,12), 1))
i13 = 13*log10(p13max) - log10(2**13) - log10(reduce(lambda a,b: a*b, range(1,13), 1))
print(f"  (b) I(k,p,1) representative space (mod-p symmetry-reduced ~ p^k/(2^k (k-1)!)):")
print(f"      k=12 @ p={p12max}: 10^{i12:.1f}   k=13 @ p={p13max}: 10^{i13:.1f}")
print(f"      growth factor 10^{i13-i12:.1f} -- and S-T name THIS (not the lifting) as the")
print(f"      k=13 bottleneck; it is a covering-configuration classification mod p.")

print("\n== (5) sanity: Rosenfeld k=7 replication ==")
target7 = lnB(7)
L8 = lcm_range(2, 8)
ps7, s = [], 0.0
for p in PR:
    if p < 31: continue
    ps7.append(p); s += log(p)
    if s > target7 - log(L8): break
print(f"  greedy primes from 31: count={len(ps7)}, largest={ps7[-1]} "
      f"(paper: 27 primes, 31..163) -> {'MATCH' if len(ps7)==27 and ps7[-1]==163 else 'CHECK'}")
