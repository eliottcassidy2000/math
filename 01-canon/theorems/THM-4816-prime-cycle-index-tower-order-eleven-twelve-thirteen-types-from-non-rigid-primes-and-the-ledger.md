---
id: THM-4816
title: "The prime cycle-index tower: at every order the rigid primes are counted by the ledger and only the non-rigid primes need enumerating, by invariant enumeration over the prime-order cycle types the cycle lemma allows; the order-11 types (12 rows over 796,876,429 primes; 12,155 non-rigid), the order-12 types (14 rows over 142,980,435,791 primes; 285,705 non-rigid) and the order-13 types (22 rows over 46,453,867,277,461 primes; 11,256,145 non-rigid), each confirmed by the labelled count"
status: >
  PROVED (A): the tower principle (cycle lemma, invariant parametrisation, rigid self-converse anti-automorphism type,
  ledger closure, the labelled-count identity). FINITE-EXACT (B)-(E): the non-rigid census at orders 3-13
  (nonrigid_primes.c), the tower tables at orders 3-13 (tower.py), with A(n) and SC(n) recomputed by Burnside's lemma
  (burnside_counts.py; equal to OEIS A000568 and A002785 through n = 15). Two independent confirmations: at orders
  3-10 the tower reproduces klein's gentourng census (THM-4814 (A)) row for row with no gentourng input, and at every
  order 3-13 the labelled prime count from 2^C(n,2) equals the sum of n!/|Aut| over the tower's rows. Not
  independently audited at filing.
session: mac-mini-2026-10-27-dipole
source: 05-knowledge/results/dipole_principle_and_prime_cycle_index_tower_20261027.md (section 6)
scripts:
  - 04-computation/experiments/prime_cycle_index_tower_20261027/nonrigid_primes.c (non-rigid census; out/nr{n}.txt, out/nonrigid{12,13}.log)
  - 04-computation/experiments/prime_cycle_index_tower_20261027/tower.py (ledger + types; out/t{n}.txt, tower_n14.out)
  - 04-computation/experiments/prime_cycle_index_tower_20261027/burnside_counts.py (A(n), SC(n) by Burnside)
  - 04-computation/experiments/prime_cycle_index_tower_20261027/format_types.py (types_tables_n11_n13.md)
related:
  - 01-canon/theorems/THM-4814-order-ten-prime-cycle-index-types-order-eleven-ledger-skew-hadamard-nuclear-norm-extremality-and-black-core-distance-of-regular-tournaments.md
  - 01-canon/theorems/THM-4803-intransitivity-is-the-relation-modulo-its-modular-tree-and-a000568-from-the-prime-cycle-index-census.md
  - 01-canon/theorems/THM-4817-order-twelve-thirteen-fourteen-ledgers-prime-self-converse-intransitivity-and-labelled-counts.md
  - 05-knowledge/hypotheses/HYP-9422-regular-tournaments-exactly-m-black-moves-from-a-singleton-plus-strong-core-nuclear-norm-maximizers-off-the-hadamard-orders-and-the-order-11-types.md
  - 05-knowledge/hypotheses/HYP-9423-dipole-endpoint-conjectures-two-connected-descent-and-the-prime-tower-predictions.md
---

# THM-4816 — The prime cycle-index tower; the order-11, 12 and 13 types

A **type** of a prime tournament `Q` is the multiset of cycle types of `Aut(Q)` together with, when `Q` is
self-converse, the multiset of cycle types of its anti-automorphisms (klein, THM-4814 (A)). The species identity of
session56 / THM-4803 needs only the types of all primes of orders `< n` to count everything of order `n` except the
order-`n` primes themselves. THM-4814 obtained the types through order 10 from a gentourng census of all 9,733,056
classes; the same route at order 11 needs all 903,753,248 classes (HYP-9422 (4): "hours on eight cores"). This
theorem replaces the census by a tower that enumerates only the non-rigid primes, which are few.

## A. The tower principle (PROVED)

**A1 (cycle lemma; klein, THM-4814 (A)).** An automorphism of prime order `p` of a prime tournament on `n > p`
vertices has at least two `p`-cycles. Automorphisms of tournaments have odd order, so a non-rigid prime tournament has
an automorphism `g` of odd prime order whose cycle type is `1^f p^c` with `n = cp + f` and either `c >= 2`, or `c = 1`
and `f = 0` (`n = p`, a circulant).

Consequences at small orders: at orders 8-13 no prime tournament has an element of order 7 (a single 7-cycle needs
`n = 7`, two need `n >= 14`); no element of order 9 at order 11 (a 9-cycle with two fixed points is excluded by the same module argument), but elements of order 9
with cycle type `3^1 9^1` (order 12) and `1^1 3^1 9^1` (order 13) do occur.

**A2 (invariant parametrisation).** Fix `g` of type `1^f p^c` with fixed set `F` and cycles `C_1, ..., C_c`. A
tournament is `g`-invariant iff it is constant on the `<g>`-orbits of pairs; since `|<g>|` is odd no orbit contains a
pair together with its reverse, so each orbit carries one free bit. The orbits are: the pairs inside `F` (the
tournament `T[F]`); for each fixed vertex and cycle, one orbit (`v -> C_i` or `C_i -> v`); inside each cycle, one
orbit per difference `1, ..., (p-1)/2` (a circulant on `Z_p`); between two cycles, one orbit per offset `e` in `Z_p`
(`c_{i,s} -> c_{j,s+e}` for all `s`, or all reversed). Hence `cf + c(p-1)/2 + C(c,2)p` bits besides `T[F]`. Since
`Sym(F)` centralises `g`, `T[F]` may be taken up to isomorphism. Every non-rigid prime class has a labelled member
invariant under the representative of the conjugacy class of one of its prime-order automorphisms, so enumerating
these invariant tournaments, keeping the prime ones and deduplicating by canonical form lists every non-rigid prime
class exactly once.

**A3 (rigid self-converse primes).** If `Q` is rigid and self-converse, its unique anti-automorphism `psi` satisfies
`psi^2 in Aut(Q) = 1`, and an anti-automorphism has at most one fixed point (two fixed points `u, v` would give
`u -> v` iff `v -> u`). So `psi` has cycle type `1^1 2^((n-1)/2)` (n odd) or `2^(n/2)` (n even).

**A4 (ledger closure).** Given the types of orders `< n`, the species identity determines the number of strong
non-prime classes and of non-strong classes of order `n`; with `A(n)` (all classes) it gives the number `P(n)` of
prime classes. The self-converse decomposition identity of THM-4814 (A) with `SC(n)` (self-converse classes) gives
the number `SCP(n)` of self-converse primes. The rigid rows of the order-`n` table are then forced:
`rigid self-converse = SCP(n) - (non-rigid self-converse)`,
`rigid not self-converse = P(n) - SCP(n) - (non-rigid not self-converse)`.
`A(n)` and `SC(n)` are recomputed here by Burnside's lemma over cycle types (`s(T) = T` counted with
`2^#pair-orbits` when no pair orbit is reversed, Davis 1954; `s(T) = T^op` when every orbit has even length and no
reversal or odd length and a reversal, which counts classes with `T ≅ T^op`). They agree with OEIS A000568 and
A002785 through `n = 15`.

**A5 (labelled-count identity).** `P_lab(n) = sum over prime classes Q of n!/|Aut(Q)|`, and `P_lab(n)` is also
`s_lab(n)` (labelled strong tournaments, from `2^C(n,2)` by `1 + T = 1/(1 - S)`) minus the labelled strong non-primes,
`n! sum_k P_lab(k)/k! [x^n] (T(x) - 1)^k` over `k < n` (THM-4814 (B)). The second expression uses no unlabelled data;
the first uses the tower's table, whose rigid rows came from `A(n)`. Their equality, the Burnside defect
`n! P(n) - P_lab(n) = n! sum over non-rigid classes of (1 - 1/|Aut|)`, is therefore an independent check of the
non-rigid census, of `A(n)` and of the ledger at once.

So the whole type table of every order follows from the non-rigid primes alone, and those have a labelled member
among `sum_{allowed (p,c)} t(f) 2^(cf + c(p-1)/2 + C(c,2) p)` invariant tournaments, where `t(f)` is the number of
tournament classes on `f` vertices: 657,440 at order 11, 15,736,832 at order 12, 373,358,656 at order 13.

## B. Validation (FINITE-EXACT)

Concurrently, klein (THM-4818 (B), pushed the same morning) enumerated the order-11 *symmetric* primes, non-rigid or
self-converse, by the same invariant method. It also included the anti-automorphism types, which this tower does not
enumerate. Its 270,999 classes split into 12,155 non-rigid and 259,307 self-converse, with the same twelve types. This
confirms the order-11 table and, by direct enumeration, the ledger's self-converse prime count `SCP(11) = 259,307`.
The theorem-level table was left to this entry (THM-4818's status line).


At orders 3-10 the non-rigid census gives 1, 0, 1, 1, 5, 12, 102, 895 classes (orders 3..10), identical type lines
to the non-rigid rows of klein's census, and the tower's full tables (rigid rows from A4) are identical to klein's
type files `prime_types_n{3..10}.txt` row for row; neither program uses the other's output. The identity A5 holds at
every order 3-13. The run times are one second at order 11, 35 seconds at order 12 and 26 minutes at order 13.

## C. The order-11 types (FINITE-EXACT)

796,876,429 prime classes (= `P(11)` of THM-4814 (B), now confirmed by A5); 12,155 non-rigid (463 self-converse),
with `|Aut|` in {1, 3, 5, 11, 55}: 3 for 12,128 classes, 5 for 23, 11 for the three non-Paley circulants, 55 for
`P_11`. HYP-9422 (4) predicted fewer than two hundred thousand non-rigid primes and about 12,200 from the Burnside
defect; the count is 12,155.

| \|Aut\| | Aut cycle types (elements) | self-converse | anti-automorphism cycle types (elements) | classes |
|---|---|---|---|---|
| 1 | `1^11` (1) | no | – | 796,605,430 |
| 1 | `1^11` (1) | yes | `1^1 2^5` (1) | 258,844 |
| 3 | `1^11` (1), `1^2 3^3` (2) | no | – | 1,224 |
| 3 | `1^11` (1), `1^2 3^3` (2) | yes | `1^1 2^5` (3) | 96 |
| 3 | `1^11` (1), `1^5 3^2` (2) | no | – | 10,456 |
| 3 | `1^11` (1), `1^5 3^2` (2) | yes | `1^1 2^2 6^1` (2), `1^1 2^5` (1) | 176 |
| 3 | `1^11` (1), `1^5 3^2` (2) | yes | `1^1 2^5` (3) | 176 |
| 5 | `1^11` (1), `1^1 5^2` (4) | no | – | 12 |
| 5 | `1^11` (1), `1^1 5^2` (4) | yes | `1^1 2^5` (1), `1^1 10^1` (4) | 5 |
| 5 | `1^11` (1), `1^1 5^2` (4) | yes | `1^1 2^5` (5) | 6 |
| 11 | `1^11` (1), `11^1` (10) | yes | `1^1 2^5` (11) | 3 |
| 55 | `1^11` (1), `1^1 5^2` (44), `11^1` (10) | yes | `1^1 2^5` (11), `1^1 10^1` (44) | 1 |

The circulant count is Burnside's: the multiplier group `Z_11^*` acts on the 32 connection sets, only the squares fix
a set (the residues), so there are `(32 + 4·2)/10 = 4` classes. Every row obeys A1 and THM-4814's anti-automorphism
lemma (no odd cycle of length `>= 3`, at most one fixed point).

## D. The order-12 types (FINITE-EXACT)

142,980,435,791 prime classes; 285,705 non-rigid (773 self-converse), `|Aut|` in {1, 3, 5, 9}: 3 for 285,601
classes, 5 for 72, 9 for 32 (8 with `Aut ≅ Z_9`, generated by a `3^1 9^1` element; 24 with `Aut ≅ Z_3^2`).

| \|Aut\| | Aut cycle types (elements) | self-converse | anti-automorphism cycle types (elements) | classes |
|---|---|---|---|---|
| 1 | `1^12` (1) | no | – | 142,978,979,910 |
| 1 | `1^12` (1) | yes | `2^6` (1) | 1,170,176 |
| 3 | `1^12` (1), `1^3 3^3` (2) | no | – | 14,572 |
| 3 | `1^12` (1), `1^6 3^2` (2) | no | – | 267,488 |
| 3 | `1^12` (1), `1^6 3^2` (2) | yes | `2^3 6^1` (2), `2^6` (1) | 304 |
| 3 | `1^12` (1), `1^6 3^2` (2) | yes | `2^6` (3) | 304 |
| 3 | `1^12` (1), `3^4` (2) | no | – | 2,774 |
| 3 | `1^12` (1), `3^4` (2) | yes | `2^6` (1), `6^2` (2) | 71 |
| 3 | `1^12` (1), `3^4` (2) | yes | `2^6` (3) | 88 |
| 5 | `1^12` (1), `1^2 5^2` (4) | no | – | 72 |
| 9 | `1^12` (1), `1^3 3^3` (2), `3^1 9^1` (6) | no | – | 8 |
| 9 | `1^12` (1), `1^3 3^3` (2), `3^4` (6) | no | – | 4 |
| 9 | `1^12` (1), `1^6 3^2` (4), `3^4` (4) | no | – | 14 |
| 9 | `1^12` (1), `1^6 3^2` (4), `3^4` (4) | yes | `2^6` (3), `6^2` (6) | 6 |

## E. The order-13 types (FINITE-EXACT)

46,453,867,277,461 prime classes; 11,256,145 non-rigid (20,925 self-converse), `|Aut|` in {1, 3, 5, 9, 13, 39}: 3
for 11,255,339 classes, 5 for 432, 9 for 368 (16 with `Aut ≅ Z_9`), 13 for five circulants and 39 for the circulant
whose connection set is a union of cosets of the cubes `{1, 3, 9}` (Burnside: `(64 + 4 + 4)/12 = 6` circulant classes).
The full 22-row table is in `types_tables_n11_n13.md`; its rigid rows are 46,453,763,046,048 (not self-converse) and
92,975,268 (self-converse, anti-automorphism `1^1 2^6`).

## F. What the tower gives next

The order-13 table closes the order-14 ledger (THM-4817): `P(14) = 27,687,212,017,612,715`, `SCP(14) = 767,626,075`,
and the defect `sum over non-rigid order-14 primes of (1 - 1/|Aut|) = 3,574,801,896/7`, whose denominator 7 says that
order-14 primes with an automorphism of order 7 exist (the first order where A1 allows `7^2`). The non-rigid counts
1, 0, 1, 1, 5, 12, 102, 895, 12,155, 285,705, 11,256,145 (orders 3-13) are dominated by the type `1^(n-6) 3^2`, whose
invariant count `2^(C(n-6,2) + 2(n-6) + 5)/(12 (n-6)!)` times the prime fraction of the invariant tournaments matches
the census to about 1% at orders 11-13 (HYP-9423 (5)). An order-14 census by the present parametrisation would need
about `1.7 · 10^10` invariant tournaments (`6880 · 2^21` for the type `1^8 3^2` alone), so the next order wants a centraliser reduction or canonical augmentation.
