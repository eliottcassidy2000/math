---
id: THM-4644
title: "Fermat primes are exactly the primes all of whose circulant tournaments are rigid; Wieferich primes are exactly the primes at which the circulant tournaments on Z_(p^2) invariant under the odd part of doubling explode from 2^((p-1)/d') to 2^((p+1)(p-1)/(2d'))"
status: PROVED (elementary; part 1 forward direction = THM-871 via Burnside's theorem on transitive groups of prime degree; exhaustive checks in the script)
source: opus-2026-10-14-S33
depends_on: [THM-871 (Fermat-rung rigidity), Burnside 1901 (a transitive group of prime degree is doubly transitive or contained in AGL(1,p))]
verification: 04-computation/experiments/number_tournament_bridges_20261014.py (sections 1-2; output 05-knowledge/results/number_tournament_bridges_20261014.out)
note: 05-knowledge/results/number_tournament_bridges_20261014.md
---

# THM-4644 - Fermat and Wieferich primes as circulant-tournament phenomena

A circulant tournament on `Z_n` is given by a connection set `S ⊂ Z_n \ {0}` with `S ∩ (−S) = ∅` and `S ∪ (−S) = Z_n \ {0}`. The arcs are `i → j` iff `j − i ∈ S`.

## Statement

1. **Fermat primes.** For an odd prime `p`, the following are equivalent:
   * every circulant tournament on `Z_p` has `Aut = Z_p` (translations only);
   * `p − 1` is a power of 2, i.e. `p` is a Fermat prime.

   If `p − 1` has an odd prime factor `r`, some circulant tournament on `Z_p` has `|Aut| >= p r`.
2. **Wieferich primes.** Let `p` be an odd prime and set:
   * `d = ord_p(2)` and `d'` the odd part of `d`;
   * `g = 2^(d/d')`.

   Then:
   * `<g>` has odd order in `(Z/p^2)^*`.
   * The circulant tournaments on `Z_(p^2)` whose connection set satisfies `gS = S` number exactly `2^E`, where

     `E = (p − 1)/d'` if `2^(p−1) ≢ 1 (mod p^2)`, and `E = (p + 1)(p − 1)/(2 d')` if `p` is a Wieferich prime.

   * So `p` is Wieferich iff the count exceeds `2^((p−1)/d')`. At the two known Wieferich primes the exponent jumps:
     * `1093`: `d = 364`, `d' = 91`; `E = 6564` instead of `12`;
     * `3511`: `d = d' = 1755`; `E = 3512` instead of `2`.
3. **Paley primes.** If `ord_p(2) = (p − 1)/2` and `p ≡ 3 (mod 4)` (so `d = d'`), the doubling-invariant circulant tournaments on `Z_p` are exactly the Paley tournament and its reverse. Examples:
   * `p = 7`, where `<2> = {1, 2, 4}`;
   * `p = 3511`, a Wieferich prime with `ord(2) = 1755`.

## Proof

**Part 1, rigidity of Fermat primes.**
- Let `T` be a circulant tournament on `Z_p`. `Aut(T)` is transitive of prime degree `p`. It is not doubly transitive, since no automorphism maps an arc `(a, b)` to `(b, a)`.
- By Burnside's theorem, `Aut(T) <= AGL(1, p)`. So `Aut(T) = Z_p ⋊ M` with `M = {a ∈ (Z/p)^* : aS = S}`.
- `−1 ∉ M`, because `−S` is the complement of `S`.
- If `p − 1 = 2^k`, the cyclic 2-group `(Z/p)^*` has `−1` as its unique involution. Every subgroup `M ≠ 1` contains `−1`, so `M = 1` (THM-871).

**Part 1, the converse.**
- Let `r` be an odd prime dividing `p − 1`, and `H` the subgroup of order `r`. Then `−1 ∉ H`.
- Negation permutes the cosets of `H` without fixed points, since `−aH = aH` would force `−1 ∈ H`.
- Choose one coset from each pair `{aH, −aH}` and let `S` be their union. Then `S ∩ (−S) = ∅`, `S ∪ (−S) = (Z/p)^*` and `HS = S`. So `Aut ⊇ Z_p ⋊ H`.

**Part 2, the orbit count.**
- **Orders.** `ord_(p^2)(2) ∈ {d, p d}`, and it equals `d` iff `2^(p−1) ≡ 1 (mod p^2)`. Hence `ord_(p^2)(g) = d' p^ε`, with `ε = 0` exactly for Wieferich primes; this is odd.
- **Unit orbits.** `<g>` acts freely on the units `(Z/p^2)^*` (order `p(p−1)`), with `p(p − 1)/(d' p^ε)` orbits.
- **Non-unit orbits.** On the non-zero non-units `pZ/p^2 ≅ (Z/p)^*` it acts with orbit size `d'`, giving `(p − 1)/d'` orbits.
- **Pairing.** Negation maps orbits to orbits and fixes none, since `<g>` has odd order and so does not contain `−1` modulo `p^2` or modulo `p`.
- **Count.** A `g`-invariant tournament is a choice of one orbit from each negation pair. So the count is `2^(N/2)`, where `N = p(p−1)/(d' p^ε) + (p−1)/d'`, which gives the two values of `E`.

**Part 3.** If `<2>` has index 2 and odd order `(p−1)/2`, it is the subgroup of squares. Its two cosets are `QR_p` and `−QR_p`, i.e. the Paley tournament and its reverse. ∎

## Checks

`.out` §1–2:
- the brute-force automorphism groups of all circulant tournaments for `p = 3, 5, 7, 11`;
- the multiplier stabilizers of all circulant tournaments for every prime `p <= 31` (rigid exactly at 3, 5, 17);
- the exponent `E` for every prime `p < 4000`, which jumps exactly at 1093 and 3511.

## Relation to other results

- **Part 1, forward direction.** This is THM-871 (Fermat-rung rigidity; there via Alspach's theorem). The converse and the equivalence are what this theorem adds.
- **Part 2.** This is the tournament count of the same doubling-orbit jump as:
  - THM-4631 (`c_2 = p c_1` exactly at Wieferich primes: 1093 has 3 and 3279 fast modes, 3511 has 2 and 7022);
  - codex's same-day split/grow tower (`05-knowledge/results/fermat_wieferich_order_20261014e.md`).

  The new ingredient is the non-unit orbits, which make `E` an exact tournament count.
- **A different tournament criterion.** Codex's same-day path-defect packet (`05-knowledge/results/tournament_path_defects_20261014e.md`) shows that `p` is Wieferich iff `p` divides the number `(2^p − 2)/p` of free marked-path orbits of `E_(p+2)`. That action is not by automorphisms.
- **Fermat numbers as path counts.** The Fermat numbers themselves appear as `H(E_n) = 1 + 2^(n−2)` (`05-knowledge/results/arithmetic_seams_20260921_tournaments.md` §5).

## Readings (not part of the theorem)

- **Fermat primes.** The five known ones are the five known primes at which the circulant tournaments are rigid. "Are there infinitely many Fermat primes?" is the same question as "are there infinitely many such primes?"
- **Wieferich primes and tournament counting.** The Euler quotient of 2 vanishes mod `p` exactly at Wieferich primes, and it enters the Burnside count `T(p) mod p` of all tournaments (THM-307; the mac-mini Wilson/{2,4,8} dictionary, THM-4641). Part 2 makes the same primes visible as a jump in a count of structured tournaments.
- **The `{7, 21}` link is NUMEROLOGY.** `7` and `21 = |Aut(P_7)|` are the vertex count and symmetry order of the first Paley prime, and 3511 is a Paley prime too. No mechanism connects Wieferich primes to the forbidden Hamiltonian-path counts `{7, 21}`.
