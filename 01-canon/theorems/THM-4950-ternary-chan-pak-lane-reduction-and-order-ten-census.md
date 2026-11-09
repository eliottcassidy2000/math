---
id: THM-4950
title: "Ternary Chan-Pak: lane reduction and the order-10 strong census. Let M3 be the multiplicative monoid of the integers v = 1 (mod 3), with 10 and 31 removed, and let H_3(T) = I(Omega(T), 3). (1) Base-3 Chan-Pak (HYP-9530 (A): every element of M3 is a value H_3(T)) holds iff every irreducible of M3 is H_3 of a strong tournament. (2) The irreducibles of M3 are exactly: the primes p = 1 (mod 3) other than 31 (lane P); the products pq of two primes = 2 (mod 3), other than 10 (lane Q); 10p (p = 1 mod 3 prime, p != 31); 31p (p = 1 mod 3 prime, including 961); 31pq (pq as in lane Q); and five sporadic values 40, 250, 310, 9610, 29791. (3) Every v = 1 (mod 3) with v <= 29,035, other than 10 and 31, is a value H_3(T), already as an order join of strong tournaments of order <= 10 (previous range 5,398); the 10,070 strong order-10 values each have a witness re-checked by independent code. (5) Strong spectrum: below 28,981 the values 1 < v = 1 (mod 3) that are H_3 of no strong tournament are exactly 10, 31 and 64 (exact below 337 by H_3 - 1 >= (3/2)(H - 1) and Busch's floor; constructive above), confirming the strong-spectrum sentence of HYP-9541 (4) on that range."
status: >
  PROVED: (1) and (2), elementary, for all v.
  FINITE-EXACT:
  - (3), by explicit witnesses: each value is an order join of strong tournaments of order <= 10. This is the same
    kind of statement as THM-4940 (9). Pre-registered Q1-Q3 (PREREGISTRATION.md) all CONFIRMED: m* = 9679 lies in
    the predicted [5000, 20000], and every irreducible below 1 + 3m* is a piece value.
  - (4): the H_3 spectrum of the strong tournaments of order 10, by complete one-vertex extension of every order-9
    class. Every order-10 value was re-derived by independent code from a stored witness (10,070 of 10,070).
  - (5): exact below 337, constructive on [337, 28,980].
  - The shape theorem (2) was also checked by direct factorisation below 200,000.
  KNOWN and used:
  - THM-4940 (1): H_3 is multiplicative over order joins.
  - THM-4940 (3), (4): 10 and 31 are never values.
  - THM-4094: the binary lane reduction, which this mirrors in base 3. S86's README (section 8) proposed the ternary
    lane reduction as the natural route.
  - kind-pasteur's classgen.c (order-9 classes) and canonical form.
  - S86's pieces: the alpha_1 <= 80 inventory, one-apex words, and Session 16's intercyclic rays.
  - Busch's formula for the strong floor: f(3q) = 3 * 5^(q-1), f(3q + 1) = 5^q, f(3q + 2) = 9 * 5^(q-1), so
    f(11) = 225 (EJC 13 (2006) N3, Cor. 1; THM-4744, THM-4870). THM-1370 (3): f is nondecreasing (Moon).
  - Concurrent lane mac-mini-2026-11-08-ternary, THM-4952 (E) (commit 844706ab1c, 14 minutes before (5) was added
    here). It introduced the method of (5): H = I(Omega, 2) <= I(Omega, 3) = H_3, so Busch's floor bounds the order of
    a strong realiser. With it, it proved that 64 has no strong realiser (orders <= 8 suffice). HYP-9541 (4)
    conjectures that the values 1 (mod 3) attained but never by a strong tournament of order >= 3 are exactly 1 and
    64. Equivalently, 64 is the only such value v > 1. (5) checks that sentence; it does not address the
    ribbon-completeness sentence of HYP-9541 (4).
  NOT CLAIMED: ternary Chan-Pak itself, which stays OPEN as HYP-9530 (A). The reduction does not make it easier than
  binary Chan-Pak: lanes P and Q are infinite.
  AUDIT: independent audit PASS WITH CORRECTIONS
  (05-knowledge/results/ternary_metagraph_20261109/audit/AUDIT.md). Every number was reproduced with independent
  code, the order-10 census by a different method without canonical forms. The proof labels of (1) are fixed (m1),
  (3) is retyped FINITE-EXACT (m2), and credit and scope are added to (5) (m3). The binary comparison is corrected
  (m4). MISTAKE-742.
session: opus-2026-11-09-S87
source: 05-knowledge/results/ternary_metagraph_20261109/README.md
related:
  - THM-4940
  - HYP-9530
  - THM-4094
  - THM-1370 (h-spectrum-omits-7-21-all-n)
  - THM-4942
  - THM-4952
  - HYP-9541
  - THM-4744
  - THM-4870
---

# THM-4950 - ternary Chan-Pak: lane reduction and the order-10 census

## Setting

- **Base-3 count.** For a tournament T, H_3(T) = I(Omega(T), 3) = sum_k alpha_k(T) 3^k is the number of 3-coloured
  odd-cycle packings (THM-4940). H_3 = 1 (mod 3), and H_3 is multiplicative over order joins (THM-4940 (1)).
- **Spec_3.** Spec_3 is the set of values H_3(T) over all tournaments of all orders. 10 and 31 are not in Spec_3
  (THM-4940 (3), (4)).
- **The monoids.** M = {v >= 1 : v = 1 (mod 3)} under multiplication, and M3 = M minus {10, 31}.
  - M3 is a submonoid: 10 = 2 * 5 and 31 are M-atoms, so a product of two elements of M3 other than 1 has at least
    two M-atoms and is neither 10 nor 31.
  - Spec_3 is contained in M3 and is closed under multiplication.

## Statement

**(1) Reduction.** Spec_3 = M3 if and only if every irreducible element of M3 is H_3 of a strong tournament.

**(2) Shape theorem.** The irreducibles of M3 are exactly the elements of these lanes:

| lane | elements |
|---|---|
| P | p prime, p = 1 (mod 3), p != 31 |
| Q | pq with p <= q primes = 2 (mod 3), pq != 10 |
| 10P | 10p, p prime = 1 (mod 3), p != 31 |
| 31P | 31p, p prime = 1 (mod 3), including 961 = 31^2 |
| 31Q | 31pq with pq as in lane Q |
| sporadic | 40 = 2^3 5, 250 = 2 5^3, 310 = 10 * 31, 9610 = 10 * 31^2, 29791 = 31^3 |

Below 200,000 there are 23,925 irreducibles: P 8987, Q 12882, 10P 1123, 31P 413, 31Q 515, sporadic 5.

**(3) Coverage to 29,035.** Every v = 1 (mod 3) with v <= 29,035, other than 10 and 31, lies in Spec_3.
- Each such v is an order join of strong tournaments of order <= 10.
- In the notation of THM-4940, Gap_3 contains no alpha-sum m <= 9678 other than 3 and 10. The previous range was
  m <= 1799.
- S86's other pieces add nothing below 29,038 (audit, closure_check.out): the alpha_1 <= 80 inventory, the strong
  one-apex tournaments and Session 16's intercyclic rays. The threshold m* = 9679 is the same with or without them.

**(4) Order-10 census.**
- The strong tournaments of order 10 have exactly 10,070 distinct H_3 values. The least is 199 and the largest is
  31,549.
- 4,063 of these values were not given by S86's pieces.
- The construction threshold is m* = 9679, value 29,038 = 2 * 14,519, the least lane-Q irreducible not attained.
- The first missing irreducible in each lane is:

  | lane | Q | P | 31P | 31Q | 10P |
  |---|---|---|---|---|---|
  | first missing | 29,038 | 29,191 | 30,721 | 31,093 | 31,210 |

- All five sporadic irreducibles are attained; 29,791 = 31^3 at order 10.

**(5) Strong spectrum.** Let v = 1 (mod 3) with v < 28,981, v > 1. Then v is H_3 of no strong tournament iff
v is in {10, 31, 64}.
- Every other such v is H_3 of a strong tournament of order <= 10.
- The first value with no strong realiser of order <= 10 is 28,981.
- Below 337 this is exact: no larger order can contribute.
- HYP-9541's convention leaves order 1 out of the strong spectrum (H_3 = 1). With that convention, (5) confirms the
  strong-spectrum sentence of HYP-9541 (4) on [1, 28,980]: the values attained but never by a strong tournament of
  order >= 3 are exactly 1 and 64.

**Calibration of the threshold.**

| strong pieces up to order | m* | value 1 + 3m* |
|---|---|---|
| 8 | 370 | 1,111 |
| 9 | 1,800 | 5,401 |
| 10 | 9,679 | 29,038 |

The threshold sits at about 92% of the largest order-10 value.

## Proof

**(1).**
- **Necessity.** If Spec_3 = M3, an irreducible v of M3 equals H_3(T) for some T.
  - v is the product of H_3 over the strong components of T (THM-4940 (1)).
  - Each factor lies in Spec_3, which is contained in M3.
  - Since v is irreducible, all factors but one equal 1. So v is H_3 of a strong component.
- **Sufficiency.** Every element of M3 is a product of irreducibles, by induction on size. Spec_3 is closed under
  products.

**(2).** Every v in M is a product of M-atoms. The M-atoms are the primes = 1 (mod 3) and the products of two primes
= 2 (mod 3). The number of atoms in any M-factorisation is the same, l(v) = w1(v) + w2(v)/2, where w1 and w2 count
the prime factors = 1 and = 2 (mod 3) with multiplicity. A factorisation v = xy in M3 with x, y > 1 is a split of an
M-factorisation into two nonempty parts, neither part equal to 10 or 31.

- **l(v) >= 4.** Every v is reducible: split off two atoms. Both parts have length >= 2, so neither is 10 or 31.
- **l(v) = 3.** v is reducible iff some M-factorisation of v contains an atom outside {10, 31}: split off that atom;
  the rest has length 2. So v is irreducible iff every M-factorisation uses only the atoms 10 and 31.
  - 10^3 = (4)(25)(10) and 10^2 * 31 = (4)(25)(31) have other factorisations, so they are reducible.
  - 10 * 31^2 = 9610 and 31^3 = 29791 have no other factorisation, so they are irreducible.
- **l(v) = 2.** v is irreducible iff every M-factorisation into two atoms uses 10 or 31.
  - v = p1 p2 (both = 1 mod 3): the factorisation is unique, so v is irreducible iff 31 divides v. This gives lane
    31P.
  - v = p q1 q2: the factorisation {p, q1 q2} is unique. v is irreducible iff p = 31 or q1 q2 = 10. This gives lanes
    31Q and 10P, plus 310 when both hold.
  - v = q1 q2 q3 q4: v is irreducible iff each of the three pairings of the multiset contains the pair {2, 5}. Write
    the multiset {2, 5, x, y}.
    - The pairing {2, x}{5, y} contains {2, 5} iff x = 5 or y = 2.
    - The pairing {2, y}{5, x} contains {2, 5} iff y = 5 or x = 2.
    - The solutions are {2, 5, 5, 5} = 250 and {2, 2, 2, 5} = 40. 100 = (4)(25) is reducible.
- **l(v) = 1.** v is an M-atom other than 10 and 31. It is irreducible. This gives lanes P and Q.

`ternary_lanes.py` confirms that the irreducibles below 200,000 are exactly these lane elements. The audit's
independent classifier agrees.

**(3), (4).**
- `ternary_order10.c` extends each of the 191,536 order-9 classes by a vertex in all 2^9 ways and keeps the strong
  extensions. Every class of order 10 is such an extension. For each, after canonical relabelling, it computes:
  - the odd-cycle counts per vertex subset (subset DP);
  - the packing vector alpha_k (subset DP);
  - H_3 = sum alpha_k 3^k.
- Check: on every 256th class, sum alpha_k 2^k = H by an independent Hamiltonian-path DP. There were 277,185 checks
  and 0 failures.
- The slices deduplicate separately. Their 70,960,990 processed classes include repeats across slices. Only the set
  of values is used, so repeats are harmless.
- Witness mode re-runs each slice and stores one canonical code per value. `verify_witness10.py` (independent Python)
  checks all 10,070 witnesses: decode, strong, H_3 recomputed from scratch, OCF. There were 0 failures.
- The audit recomputed the order-10 value set by incremental extension, with no canonical forms, and got the same
  set.
- `ternary_coverage10.py` forms the multiplicative closure of all piece values, using S86's closure routine. Every
  product of piece values is realised by the order join of the pieces.

**(5).** This is THM-4952 (E)'s method, run on the complete census of orders <= 10.
- For every tournament, H_3 - 1 = sum_k alpha_k 3^k >= (3/2) sum_k alpha_k 2^k = (3/2)(H - 1), since
  3^k >= (3/2) 2^k for k >= 1.
- A strong tournament of order n >= 11 has H >= f(n) >= f(11) = 225 (Busch; THM-1370 (3)). So its H_3 is at least
  1 + (3/2) * 224 = 337. A strong realiser of a value below 337 therefore has order <= 10.
- The strong classes of order <= 10 are complete: S86's census for orders <= 9 and (4) for order 10. Among their H_3
  values below 337, exactly 10, 31 and 64 are missing.
- Above 336, every v = 1 (mod 3) up to 28,978 occurs among them.
- The check is `strong_spectrum10.py` (output `strong_spectrum10.out`); the audit's `compare_census.out` agrees.

## Remarks

- **Why order 10 helps.** The lane-P and lane-Q irreducibles have no smaller factorisation, so each must be a single
  strong value.
  - H_3(T) <= (3/2)^K H(T), where K is the largest k with alpha_k(T) > 0. At order 10, K <= 3 and max H = 15,745, so
    every piece value is below 53,200.
  - The threshold m* therefore tracks the top of the newest order's spectrum, moving about 5 times further per order.
- **Comparison with base 2.**
  - In the odd numbers without 7 and 21 (THM-4094; S86 README section 5) the irreducibles are:
    - the primes p != 7;
    - 7p with p not 3 or 7;
    - the carry atoms 49, 63 = 9 * 7 and 343 = 7^3.
  - Base 3 has the counterparts 40 = 4 * 10 and 29791 = 31^3 among its sporadic values.
  - Base 3 has more infinite lanes because M is not factorial: products of two primes = 2 (mod 3) are atoms. The
    excluded value 10 is itself such a product, which creates the lane 10P and the sporadic 40, 250.
- **Status of HYP-9530 (A).** It is FINITE-EXACT for alpha-sums <= 9678 (values <= 29,035) and OPEN beyond.
