---
id: THM-4644
title: "Fermat primes are exactly the primes all of whose circulant tournaments are rigid; Wieferich primes are exactly the primes at which the circulant tournaments on Z_(p^2) invariant under the odd part of doubling explode from 2^((p-1)/d') to 2^((p+1)(p-1)/(2d'))"
status: PROVED (elementary assembly, audited - MISTAKE-608. Part 1 at the multiplier level is Akin 2018, arXiv:1806.11241, Thms 3.13-3.14 (both directions); the full automorphism group follows from Alspach's theorem (as in THM-871) or from Burnside's, which answers Akin's open remark about 17 and the other Fermat primes. Part 2 is Akin's invariant-subset count for <g> on Z_(p^2) plus the Wieferich order dichotomy. Exhaustive checks in the script.)
source: opus-2026-10-14-S33
depends_on: [THM-871 (Fermat-rung rigidity, via Alspach), Akin 2018 (arXiv:1806.11241, Thms 3.12-3.14), Burnside 1901 (a transitive group of prime degree is doubly transitive or solvable with a normal Sylow p-subgroup), Alspach 1970]
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
- By Burnside's theorem, `Aut(T)` is solvable with a normal Sylow `p`-subgroup of order `p` (as `p ∥ p!`). The translation group `Z_p` has order `p`, so it is that normal subgroup. Hence `Aut(T) <= N_(S_p)(Z_p) = AGL(1, p)` in the given coordinates. (Alspach's theorem gives the same directly.)
- So `Aut(T) = Z_p ⋊ M` with `M = {a ∈ (Z/p)^* : aS = S}`.
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

- **Part 1.**
  - At the level of multipliers, both directions are Akin's Theorems 3.13–3.14 (arXiv:1806.11241, 2018), stated for every odd modulus: `Z^*_(2n+1)` acts freely on connection sets iff `2n+1` is a squarefree product of Fermat primes, and every odd-order multiplier subgroup fixes some connection set.
  - Akin then writes that he does not know whether `Aut = Z_p` for every rotational tournament at 17 or the other known Fermat primes.
  - The forward direction for the full automorphism group is THM-871 (via Alspach's theorem); the Sylow argument above gives it via Burnside. Either answers Akin's question for every Fermat prime.
  - This theorem assembles these; nothing in part 1 beyond the assembly is new.
- **Part 2.** This is the tournament count of the same doubling-orbit jump as:
  - THM-4631 (`c_2 = p c_1` exactly at Wieferich primes). THM-4631 counts `<2>`-orbits: 1093 has 3 and 3279 fast modes, 3511 has 2 and 7022. Part 2 counts `<g>`-orbits: 12 and 13116 unit orbits at 1093, coinciding with THM-4631 at 3511, where `g = 2`. Both counts jump by the factor `p`; the exponent `E` jumps by `(p+1)/2`.
  - codex's same-day split/grow tower (`05-knowledge/results/fermat_wieferich_order_20261014e.md`).

  The count `2^m` of invariant connection sets is Akin's Theorem 3.14 construction. The new ingredients are the Wieferich dichotomy for `<g>` and the non-unit orbits, which together make `E` an exact tournament count.
- **A different tournament criterion.** Codex's same-day path-defect packet (`05-knowledge/results/tournament_path_defects_20261014e.md`) shows that `p` is Wieferich iff `p` divides the number `(2^p − 2)/p` of free marked-path orbits of `E_(p+2)`. That action is not by automorphisms.
- **Fermat numbers as path counts.** The Fermat numbers themselves appear as `H(E_n) = 1 + 2^(n−2)` (`05-knowledge/results/arithmetic_seams_20260921_tournaments.md` §5).

## Readings (not part of the theorem)

- **Fermat primes.** The five known ones are the five known primes at which the circulant tournaments are rigid (Akin asked this for 17 and the other known Fermat primes). "Are there infinitely many Fermat primes?" is the same question as "are there infinitely many such primes?"
- **Not doubling-invariant.** At 1093, `−1 = 2^182 ∈ <2>` modulo `p^2`, so no circulant tournament on `Z_(p^2)` is invariant under doubling itself. That is why the theorem uses the odd part `g`.
- **Wieferich primes and tournament counting.** The Euler quotient of 2 vanishes mod `p` exactly at Wieferich primes, and it enters the Burnside count `T(p) mod p` of all tournaments (THM-307; the mac-mini Wilson/{2,4,8} dictionary, THM-4641). Part 2 makes the same primes visible as a jump in a count of structured tournaments.
- **The `{7, 21}` link is NUMEROLOGY.** `7` and `21 = |Aut(P_7)|` are the vertex count and symmetry order of the first Paley prime, and 3511 is a Paley prime too. No mechanism connects Wieferich primes to the forbidden Hamiltonian-path counts `{7, 21}`.
