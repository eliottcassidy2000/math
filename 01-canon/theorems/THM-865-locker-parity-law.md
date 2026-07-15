---
id: THM-865
title: THE LOCKER PARITY LAW IS FALSE — H(D_11) = 4027 ≡ 3 (mod 4): refutation of THM-853(III)'s conjecture, with the exact localization (α₁ stratifies over top vertices; t_m odd first at m = 11) and the post-mortem of the divisor-pairing mechanism (pair-matching is real at small composite m, dies at m = 20; primes are lawless)
status: (i) stratification PROVED (all n); (ii) refutation MACHINE-EXACT, doubly verified (independent DFS + subset-DP for α₁; Held-Karp for H; THM-466 digit identity cross-checked at every n ≤ 15); (iii)-(iv) exact censuses m ≤ 20
source: mac-mini-2026-07-15-S109; owner directive 2026-07-15 ("prove the locker parity law via the divisor-pairing involution") — the honest outcome is the refutation
depends_on:
  - THM-002 (OCF), THM-466 (m=2 digit corollary: H ≡ 1 + 2α₁ mod 4)
  - THM-853(III) (D_n definition + the conjecture, kind-pasteur S128 c20)
related:
  - opus-S316 probe (ii) (the 2-adic tower above Rédei decays at n = 7 — this file is its arithmetic-family witness)
scripts: 04-computation/locker_parity_law_census_macmini_S109.py, locker_parity_law_verify_macmini_S109.py -> 05-knowledge/results/locker_parity_law_{census,verify}_macmini_S109.out (+ locker_parity_t17_19_macmini_S109.out)
---

# THM-865 — the locker parity law is false

D_n = the locker (divisibility) tournament of THM-853(III): vertices 1..n; for u < v,
v → u iff u | v or v = u+1, else u → v. Conjectured there: H(D_n) ≡ 1 (mod 4) for all n,
via "the divisor-pairing involution d ↔ x/d on odd cycles". By THM-466 (m = 2),
H ≡ 1 + 2α₁ (mod 4), so the conjecture says: α₁(D_n) = #{directed odd cycles} is even.

## (i) Stratification (PROVED, all n)

D_m is the induced subtournament of D_n on {1..m} (the arc rule depends only on the pair),
so every cycle has a well-defined top vertex and

> α₁(D_n) = Σ_{m ≤ n} t_m,  t_m := #{directed odd cycles through vertex m in D_m}.

Vertex 1 has out-degree 0 (loses every arc), so no cycle uses it. The out-set of the top
vertex m is {proper divisors of m} ∪ {m−1}, and exits m → 1 continue nowhere; so every
cycle through m enters from a non-divisor u ≤ m−2 and exits to a proper divisor ≥ 2 of m
or to m−1. The conjecture ⟺ t_m even for every m. ∎

## (ii) Refutation (machine-exact, two independent methods)

**t_11 = 451 is odd.** Hence α₁(D_11) = 909 is odd and

> **H(D_11) = 4027 ≡ 3 (mod 4).**

The law holds for n = 3..10 and first fails at n = 11. Full rows (Held-Karp H, exact):

| n | 3..10 | 11 | 12 | 13 | 14 | 15 |
|---|-------|----|----|----|----|----|
| H(D_n) | 1, 1, 5, 9, 33, 109, 469, 1721 | 4027 | 28851 | 83817 | 400569 | 3141317 |
| H mod 4 | all 1 | **3** | **3** | 1 | 1 | 1 |

H ≡ 1 + 2α₁ (mod 4) verified at every n ≤ 15 (THM-466 cross-check, 0 failures). The
n = 5..9 row 5, 9, 33, 109, 469 reproduces THM-853(III) exactly (same convention).

## (iii) The t_m parity table (exact, m ≤ 20)

| m | 3..10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 |
|---|-------|----|----|----|----|----|----|----|----|----|----|
| t_m | 0,0,2,2,8,28,96,322 | **451** | 5272 | **9999** | 53840 | 423430 | **4101613** | 5037556 | 88568440 | 148162470 | **2818298633** |
| parity | even | **odd** | even | **odd** | even | even | **odd** | even | even | even | **odd** |

Odd exactly at m ∈ {11, 13, 16, 20}: two primes, a square, and a number that is neither.
**No prime/composite/square law.** Induced H(D_n) mod 4, n = 3..20:
1,1,1,1,1,1,1,1,3,3,1,1,1,3,3,3,3,1.

## (iv) Post-mortem of the divisor-pairing mechanism

Stratify t_m by exit cell (the divisor b with m → b, or b = m−1):

- **Pair-matching is real at small composite m**: the cells of a divisor pair d ↔ m/d have
  EQUAL parity at m = 8 (2↔4), 12 (2↔6 and 3↔4), 14 (2↔7), 16 (2↔8), 18 — so paired cells
  cancel mod 2 — and at the squares the fixed cell d = √m carries the parity: even at
  m = 9, **odd at m = 16** (the locker signature: t_16's oddness sits exactly on the
  self-paired divisor). Through m ≤ 18 the composite data looked exactly like the
  conjectured mechanism, with mismatched pairs (m = 6: 2↔3; 10: 2↔5; 15: 3↔5) always
  compensated by an odd m−1 cell.
- **The mechanism dies at m = 20**: pair 2↔10 mismatches (odd/even), pair 4↔5 matches
  (odd/odd), and the m−1 = 19 cell is EVEN — no compensation. t_20 is odd at a
  non-square composite.
- **Primes never had the mechanism**: a prime top vertex has the single live exit m−1
  (divisor exits: only the dead vertex 1), so nothing pairs; t_p = 2, 8, 451, 9999,
  5037556, 148162470 at p = 5, 7, 11, 13, 17, 19 — parities 0,0,1,1,0,0, lawless. The
  global law was first killed by the prime channel (m = 11), and would have been killed
  by the composite channel at m = 20 anyway.

## (v) Reading

The τ-parity of the integer locker does NOT transfer to the H-parity of the locker
tournament. This is the expected face of THM-466(iv) and opus-S316 probe (ii): digit 0 of
H is Rédei-constant, but digit 1 (α₁ mod 2) is not tamed — not by scores, and now not by
the divisibility lattice either, even on a single canonical arithmetic family. The 2-adic
tower's decay has an arithmetic witness: the locker tournament inherits the locker's
DEFINITION but not its PARITY LAW. (What survives of the intuition: the exit-cell
pair-matching at m = 8..18 and the square-cell signature at m = 16 — recorded above as
data, unproved, and already known to be non-universal.)
