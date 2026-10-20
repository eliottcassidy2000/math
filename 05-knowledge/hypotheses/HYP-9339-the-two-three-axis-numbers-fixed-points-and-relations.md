---
id: HYP-9339
title: "The 2-3 axis: numbers, fixed points and relations. Thesis: the arithmetic of the repository's three conjectures is organised by the two smallest primes acting in three roles - multiplicatively (the lattice 2^a 3^b: 3-smooth ladders, Pythagorean scales, the critical slope log_3 2, Furstenberg's x2x3 rigidity), additively (the only coincidences |2^a - 3^b| = 1 of Gersonides/Stormer and Catalan's 9 = 8 + 1, which are the Collatz free cycles, the Catalan carry and the first anchor collision), and relationally (binary information = tournaments = Z/2, ternary information = oriented graphs = Z/3-tournaments and 3-hypertournaments, the tie state existing iff 2 is invertible mod k) - with fixed points as the common currency: Burnside counts isomorphism classes by fixed points of relabellings, Collatz anchors are fixed points of letters, the Farey map's two cusps 0 and 1 are the even and odd anchors, the golden ratio is the fixed point phi^2 = phi + 1 where the two operations agree, and the 2-adic order of 3 (period 2^(K-2)) and the 3-adic order of 2 (clock tower 2 x 3^(n-1)) are the two clocks every Collatz word obeys; conjectures: (1) every integer 3x+1 cycle beyond Gersonides' four has an unbalanced word, (2) the Stormer pairs are the only slopes at which an integer cycle can be balanced, (3) the k-tournament counts T_k(n) for even k satisfy a coloured 'tournaments = even graphs' identity, (4) Furstenberg's x2x3 rigidity is the measure-theoretic shadow of the Collatz critical slope: the only ergodic x2x3-invariant measures are Lebesgue or atomic, and the only Collatz-invariant 'slopes' with integer cycles are the convergents of log_3 2"
status: >
  Synthesis hypothesis; every row typed. Multiplicative/additive rows PROVED or KNOWN in the cited theorems (THM-4484,
  THM-4640, THM-4716, THM-4717, THM-4719, THM-4591, THM-4633; Stormer 1897 and Gersonides 1343 CITED; Furstenberg 1967 and
  Rudolph 1990 CITED). Relational row PROVED (THM-4720). Fixed-point row: THM-4600 (anchors), THM-4714 (Farey cusps), the
  golden self-duality (HYP-9337). Conjecture (1) is half proved (THM-4716: balanced => free => Gersonides), (2) is a
  restatement of THM-4716 plus the Stormer classification (the balanced cycles exist only at |2^q - 3^p| = 1, i.e. the
  consecutive 3-smooth pairs 1|2, 2|3, 3|4, 8|9), hence PROVED; (3) OPEN; (4) a reading, OPEN. Pre-audit.
session: mac-mini-2026-10-19-twothree
source: 05-knowledge/results/two_three_axis_20261019_numbers_fixed_points_relations.md
related:
  - THM-4720, THM-4719, THM-4717, THM-4716, THM-4714, THM-4640, THM-4591, THM-4556, THM-4554, THM-4522 (x2x3 lonely spectrum), THM-4521 (dictionary of sixes), THM-4581 6(b) (grand orbits = Z[1/6] x <2,3>), THM-1980, THM-4600, HYP-9337
---

# HYP-9339 — the 2–3 axis: numbers, fixed points and relations

## The three roles of 2 and 3

| role | objects | repository theorems | status |
|---|---|---|---|
| multiplicative | the lattice `2^a 3^b`; 3-smooth ladders; Pythagorean scales; windows; the critical slope `log_3 2 = 0.6309`; `x2x3` on the circle | THM-4716/4717/4719; THM-4591; THM-4633; THM-4522 (lonely spectrum of `2^j 3^k`); Furstenberg 1967, Rudolph 1990 (CITED) | PROVED / CITED |
| additive | `|2^a − 3^b| = 1` only for `(1,0), (1,1), (2,1), (3,2)` (Gersonides; Stormer's theorem for `S = {2,3}`); `9 = 8 + 1` is the Catalan carry and the first anchor collision; the Gersonides shapes are the four free `3x+1` cycles | THM-4484, THM-4640 (5), THM-4713 | KNOWN / PROVED |
| p-adic clocks | the 2-adic order of 3 (`2^(K−2)`: the Mersenne word period) and the 3-adic order of 2 (`2·3^(n−1)`: the clock tower); the level operator on units mod `3^n` with spectrum on the half-circle | THM-4556 (iv), THM-4554, THM-4520 | PROVED |
| relational | binary information on pairs = tournaments (`Z/2`); ternary = oriented graphs (`Z/3`, tie `= 1/2 mod 3`) and 3-hypertournaments (cyclic orientations of triples); the tie state exists iff `2` is invertible mod `k` | THM-4720; THM-001 (Redei), THM-1980 (2-adic edge of the Redei count) | PROVED |
| fixed points | Burnside = fixed points of relabellings; anchors `c_a = 1/(2^a − 3)` = fixed points of letters; the Farey map's cusps `0, 1` = the even/odd anchors; `φ² = φ + 1`; Minkowski `?` fixes `0, 1/2, 1` | THM-4600, THM-4714, THM-4717, HYP-9337 | PROVED |
| the conductor 6 | Collatz lives on `Z[1/6]`: grand orbits are the `Z[1/6] ⋊ ⟨2,3⟩` orbits up to null sets; the mod-6 Collatz graph; the dictionary of sixes; A073631 excludes 561 through the prime 3 | THM-4581 6(b), THM-4521, THM-4713 note | PROVED / KNOWN |

## Theory-building statements

1. **Coincidence principle.** Additive coincidences of the multiplicative atoms are finite and tiny (Stormer), and each one is a structural object elsewhere: `2 − 1 = 1` (the fixed point 0 / the even anchor), `3 − 2 = 1` (the odd anchor `−1`), `4 − 3 = 1` (the trivial cycle `{1, 2}`, the unit carry), `9 − 8 = 1` (the `{−5, −7, −10}` cycle, the Catalan carry, the `(4,1,1)` re-anchoring block). Nothing larger: the fifth coincidence does not exist, and that absence is why `3x+1` has no further balanced cycle (THM-4716).
2. **Critical-slope principle.** The one slope at which the two multiplicative actions balance, `a = log_3 2`, is the slope of every integer cycle and the slope at which the Sturmian circles fail to exist (THM-4717). Subcritical slopes (the golden one included) carry positive circles, supercritical ones negative circles. Furstenberg's rigidity says the `x2x3` semigroup has no invariant structure between "everything" and "finite"; the Collatz reading is that cycles can only hide exactly at the balance point.
3. **Parity principle for relations.** Going from binary to `k`-valued relations, the structure splits by the parity of `k`: even `k` behaves like tournaments (no symmetric state, only odd-cycle permutations fix anything, rigidity defect `k^(−2n)`), odd `k` like graphs (a tie state, every permutation contributes, defect `k^(−n)`). The arity direction (triples) keeps the two-state sign structure (even/odd return permutation).
4. **Fixed points as the common currency.** Counting (Burnside), dynamics (anchors, cusps, the golden point) and relations (the tie state `x = 1 − x`) all reduce to fixed points of an involution or a letter; the involution `x ↦ 1 − x` on `Z/k`, the involution `x ↦ 1 − x` on `[0,1]` under which `?` is symmetric, and the letter `a` fixing `c_a` are the same kind of object.

## Conjectures

* **(1)** Every integer cycle of `3x+1` other than `{0}, {−1}, {1,2}, {−5,−7,−10}` has an unbalanced parity word (half proved: THM-4716).
* **(2)** PROVED form: balanced integer cycles of `3x+1` exist exactly at the Stormer pairs `(2^q, 3^p) ∈ {(2,1), (2,3), (4,3), (8,9)}`.
* **(3)** For even `k`, `T_k(n)` equals the number of isomorphism classes of `(k/2)`-coloured even graphs (a coloured version of "tournaments and even graphs are equinumerous"); OPEN, testable by Burnside.
* **(4)** The measure-theoretic rigidity of `x2x3` (Furstenberg; Rudolph for positive entropy) and the arithmetic rigidity of Collatz cycles at the critical slope are two faces of one statement about the semigroup `⟨2, 3⟩`; OPEN as mathematics, offered as the organising thesis.
