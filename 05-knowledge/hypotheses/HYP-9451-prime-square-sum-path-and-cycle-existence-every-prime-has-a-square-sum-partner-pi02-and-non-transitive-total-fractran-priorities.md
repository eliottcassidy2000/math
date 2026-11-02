---
id: HYP-9451
title: "After THM-4860: (1) [existence PROVED by THM-4860 (f): the first 50,396 primes admit a square-sum arrangement] minimality also PROVED (THM-4860 (e′): 50,396 is the least count); the remaining question: do infinitely many bounds admit an arrangement (every one needs Chebyshev's bias mod 4 in {0,1,2} and six mod-36 pair conditions, which below 5·10^6 hold only on two windows); a Hamiltonian cycle needs bias exactly 1 and no leaf prime and has no candidate below 700,000; (2) every prime p has a prime q with p + q a square (the infinite prime square-sum graph has no isolated vertex; Pi^0_2, Bunyakovsky type for the family s^2 - p); (3) a FRACTRAN program whose conflict rule is a non-transitive tournament on its fractions can be total only when its denominators share factors, and PRIMEGAME's do; does any total non-transitive rule exist for PRIMEGAME's fourteen denominators?"
status: >
  OPEN (klein-2026-11-01) apart from the existence and minimality in (1), PROVED (THM-4860 (e′), (f)). (1) The necessary
  conditions are exact (THM-4860 (d), (d′), (d″), (e), (e′)); the flow-plus-merging construction decided the first
  survivor in seconds. Chebyshev's bias is negative only in rare windows (Leech 1957: first crossing 26,861; second window 616,841 to
  633,798, CITED), so candidates are rare and the question is genuinely about those windows. (2) is weaker than
  Bunyakovsky for s^2 - p (one prime value instead of infinitely many) but I know no proof; the data (one isolated prime
  below 10^5, namely 57,649, first partner 118,751; every prime below 100 has a partner below 500) is consistent with
  the heuristic that s^2 - p is prime for some s <= C sqrt(p) log p. (3) is a definition-level question (README §1
  (D1)): with pairwise coprime denominators totality forces transitivity; PRIMEGAME's realisable applicable sets are
  the up-sets of divisibility by its denominators and may admit a non-transitive total rule.
session: klein-2026-11-01 (codex/session-ramsey-20261030)
source: 05-knowledge/results/machines_automata_fractran_tournaments_dictionary_20261101/README.md
related:
  - 01-canon/theorems/THM-4860-the-prime-square-sum-graph-is-bipartite-by-residue-mod-4-two-primes-1-mod-3-are-never-adjacent-chebyshev-bias-obstruction-to-prime-square-sum-paths-none-below-615731-and-the-25-primes-below-100-by-ogg-layer.md
  - 01-canon/theorems/THM-4856-goldbachgame-fermatgame-legendregame-halts-iff-counterexample-fractran-programs-built-on-primegames-divisor-loop.md
  - 01-canon/theorems/THM-4852-primegame-phase-law-both-published-lists-agree-selector-cellular-automaton-in-base-29-primorial.md
---

# HYP-9451 — after THM-4860: the prime square-sum path, the partner question, and tournament-priority FRACTRAN

Terms are those of THM-4860: `G_N` is the square-sum graph on the primes below `N`; `A = {2} ∪ {p ≡ 1 (4)}`, `B = {p ≡ 3 (4)}`.

## 1. The prime square-sum path: existence PROVED, minimality and infinitude OPEN

**Existence (PROVED, THM-4860 (f)).** The first 50,396 primes (all primes below 617,059) admit a Hamiltonian path of `G_N`; the `Σ⁰₁` search of the dictionary note halted. Every bound below 617,053 is excluded by THM-4860 (d), (d′), (e).

**Minimality (PROVED, THM-4860 (e′)).** The bound 617,053 (the first 50,395 primes) passes the 2-matching test but violates the segment condition `c_25 − c_11 ≤ e_25 = 0`; so 50,396 is the exact least count. Nothing is open here.

**Infinitude (OPEN).** Do infinitely many bounds admit an arrangement? Each needs Chebyshev's mod-4 bias in `{0,1,2}` and the six mod-36 pair conditions of THM-4860 (d″); Littlewood's theorem gives infinitely many sign changes of the mod-4 bias, but the joint condition is a question about the simultaneous fluctuations of eight residue classes, and below `5·10^6` only the windows `(100, 27,031]` and `[617,039, 634,199]` pass the counts.

**Cycle version.** A Hamiltonian cycle of `G_N` needs `|A| = |B|` (bias exactly 1) and minimum degree 2; no bound below 700,000 qualifies, because a leaf prime is present at every bias-1 bound in the window (at 700,000 the only leaf is 582,451, whose single partner is 362,333). The next Chebyshev windows are the only places to look.

## 2. Every prime has a square-sum prime partner (OPEN; `Π⁰₂`)

**Claim.** For every prime `p` there is `s` with `s² − p` prime. Equivalently the infinite prime square-sum graph has no isolated vertex. Data: isolated primes below `10^5` are reduced to 57,649 (partner 118,751); 111,031 is the only isolated prime below 120,000 (partner 136,973). The statement is `∀p ∃s`, so, unlike the `Π⁰₁` prime conjectures of THM-4856, it is the totality of a machine and not the halting of one.

## 3. Non-transitive total priorities for FRACTRAN (OPEN; definition-level)

With "apply the applicable fraction that beats all other applicable ones" under a tournament `R`, the step is total iff every realisable applicable set `App(n) = {i : d_i | n}` induces a subtournament with a source. For pairwise coprime denominators this forces `R` transitive (README §1 (D1)). PRIMEGAME's denominators `91, 85, 51, 38, 33, 29, 23, 19, 17, 13, 11, 14, 2, 1` share factors; **is there a non-transitive `R` that is total on PRIMEGAME?** If yes, PRIMEGAME's "first applicable" order is a choice, not a necessity, and the resulting programs are new objects to compare with THM-4852's phase law.
