# The Hamiltonian-path spectrum of the 3-adic chart class

kind-pasteur-2026-11-13-S27, THM-5152. **Chan–Pak remains open.**

## 0. The class

Let G be the grammar generated from one vertex by ordinal joins X => Y and cyclic joins C3[X,Y,Z] (X -> Y -> Z -> X).
Session 312 ([CLUSTER_CHART](../session312_cluster_pythagorean_symmetry_20261113/CLUSTER_CHART.md)) proves the
following. A tournament is the orientation chi(det(v_i, v_j)) of projectively distinct points v_i in Q_3^2 (chi =
leading 3-adic unit) iff some vertex switch of it lies in G. Call these tournaments the **chart class** C. They are
the tournaments whose arcs are read off from the digit-reversed address chart (valuation = common ancestry, leading
digit = orientation).

The spectrum of G itself, without switching, is KNOWN: it is Session 268's composable-tournament census
(session268_composable_tournaments_20261111/modular, n <= 16), which names 13 as outside it. Only the switching
closure C is studied here.

## 1. A finite decision procedure (proved)

**Lemma 1 (heredity).** C is closed under induced subtournaments. This is Session 312's remark, since a subset of a
point configuration is a point configuration.

**Lemma 2 (reduction).** Let T be in C with strong components S_1 => ... => S_k. Deleting the singleton components
keeps T in C (Lemma 1) and keeps H(T), because H is multiplicative over ordinal sums. Each nontrivial S_i is a strong
member of C, and H(T) is the product of the H(S_i).

**Lemma 3 (Busch).** A strong tournament of order n >= 3 has H >= f(n) = min{3^a 5^b : 2a + 3b = n − 1}. The values
are 3, 5, 9, 15, 25, 45, 75, 125, 225, 375 for n = 3..12, and f increases. This is Busch's theorem, attained by
interval chains (THM-4744).

**Proposition 4.** Fix B and let n_0 be largest with f(n_0) < B. The set of odd values below B attained by C is the
multiplicative closure, inside [1, B), of the values attained by strong members of C of order at most n_0.

*Proof.* Combine Lemmas 2 and 3. A strong member with H < B has f(n) <= H < B, so n <= n_0. QED.

So whether C attains v < B is a finite computation. In particular, a prime or 7p value missed by every strong member
of order <= n_0 is missed at **every** order, since 7 is never a count (THM-1370).

## 2. The computation (FINITE-EXACT)

**Generation.** G is enumerated up to isomorphism by canonical modular-decomposition trees: linear nodes flattened
into ordinal chains, C3 nodes up to rotation, both cyclic orders allowed. Gallai's decomposition is unique, so these
trees are exactly the isomorphism classes. The counts are 1, 1, 2, 4, 9, 24, 65, 182, 528, 1,558 for orders 1..10;
4,674 at order 11 (7,048 classes in all).

The read-only sweep of this session regenerated the counts through order 10 independently and matched Session 268's
classes and H sets. It also brute-forced all labelled tournaments of orders 5 and 6 as a control.

**Switching.** For each class and each of its 2^(n−1) switch sets (a set and its complement give the same
tournament), compute H by the subset dynamic programme (scripts/swspec.c).

**Result.** Every class and every switch through order 11 was computed (48 minutes for order <= 11). The run was
pre-registered as P1 at 961f30e87b. A first launch was stopped unread, and the class count was corrected before
the result was read.

- **Values below 3,000.** 125 odd values are not attained through order 11
  (scripts/miss11_order_le_11_below_3000.txt). Values above f(12) = 375 are not yet decided.
- **Below f(12) = 375 the omissions are decided**, by Proposition 4 together with the multiplicative closure
  (scripts/closure.py). Besides 7, the values omitted at every order are:
  - 10 primes: 31, 53, 113, 127, 149, 173, 181, 251, 263, 277;
  - 5 values 7p: 21, 35, 203, 217, 287;
  - 2 other values: 341 = 11 · 31 and 343 = 7^3.

  The omitted primes are 11 of the 73 odd primes below 375, counting 7.
- **Pre-registered P1: CONFIRMED.** The ten omissions below 225 persist, and all six candidates in (225, 375)
  (251, 263, 277, 287, 341, 343) are omitted at every order. 31 and 127 remain omitted.
- **Independent check.** The read-only sweep of this session reproduced the order-10 spectrum with its own code
  (sweep27c/sw.c). It certified the omissions 31, 35, 53, 113, 127, 149, 173, 181 and 203 by the same reduction.

## 3. Consequences

- **A chart-only selector cannot be complete.** It misses the primes and the 7p targets listed in section 2 at every
  order. The 7p omissions (21, 35, 203, 217, ...) are exactly the kind of target THM-4094 §6 leaves to strong cores.
- **Most omissions are primes or 7p.** Below 225 every omission is prime or 7p. Below 375 two more appear, 341 =
  11 · 31 and 343 = 7^3. In both, every factorization uses an omitted value (31 or 7), as Proposition 4 requires.
- **Mersenne curiosity.** 7, 31 and 127, the Mersenne primes 2^p − 1 for p = 3, 5, 7, are all omitted, while 3, 63,
  511 and 1023 are attained. No mechanism is known. The next test is 8191, which needs strong members of order up to
  n with f(n) <= 8191, that is n <= 18.
- **Relation to Session 312's boundary.** That boundary shows the neutral partner D of the order-8 core B is not
  chartable. The omissions here show that even values attained by small ordinary tournaments, such as 31 at order 6,
  need non-chart structure. Module interiors, or Session 313's tournament-specific port data, are unavoidable for
  selection.
- **Not closed under adding a sink.** The sweep observed that 13 is a chart value at order 5 but not at order 6.
  This is a property of the class: an ordinal join with a sink need not be a switch of a grammar member. Proposition 4
  is unaffected, because it only deletes singleton components.
