---
id: HYP-9414
title: "After THM-4804: injective fibres (distinct tilings of a carousel class lead to distinct partner classes, for every odd m), the pair-module partners and the limit proportion of prime partners, connectivity of the fibre's line neighbourhood, the orbit-size law tc(T_w) = H |Orb_G(w)|/(3^r |G|) for every decorated fibre of Theorem G, and which classes of the line metagraph are injective in general"
status: "NUMERICAL ((1), (2), (3) exhaustive at m = 3 / n = 12); (4) PROVED for the carousel case by THM-4804 (2), conjectured in general; OPEN"
session: opus-2026-10-26-S59
source: 05-knowledge/results/module_shattering_lines_carousel_finite_field_20261026.md
related:
  - 01-canon/theorems/THM-4804-complement-lines-shatter-every-module-and-the-carousel-fibre-as-frobenius-orbits-of-a-normal-basis.md
  - 05-knowledge/results/cyclictriple_composition_20261026/README.md
  - 05-knowledge/results/tournament_composition_classification_20261026/CLASSIFICATION.md
  - 05-knowledge/hypotheses/HYP-9352-presentation-count-laws-strong-spectrum-forbidden-fibre-sizes-and-the-three-move-law-at-n-11.md
---

# HYP-9414: after module shattering

1. **Injective fibres.** For every odd `m >= 3` and every necklace class `T_w` of the carousel fibre, distinct tilings of `T_w` lead to non-isomorphic line partners, so every black line of `T_w` has class multiplicity one in THM-4731's black multigraph and the line neighbourhood of `T_w` has exactly `tc(T_w) = H_m d/(m 3^m)` classes. (Exhaustive at `m = 3`.) General form: a class whose partners are prime-rich in this sense is injective; at `n <= 9` the census has black multiplicities above one, so injectivity is not universal; conjecture: it holds for every class all of whose line partners are prime. Hostile to look for: a class with two tilings whose partners are isomorphic through an isomorphism that is not induced by an automorphism of the class.

2. **Pair-module partners.** The partners of the carousel fibre that are not prime have exactly one proper non-singleton module, of size 2 (exhaustive at `m = 3`: `6480` or `6453` of `462,915` presentations per class). Conjecture: for every odd `m` the non-prime partners of `T_w` have only a pair module, formed by two consecutive path vertices of the same block whose remaining in- and out-neighbourhoods coincide after the reversal, and their proportion among the presentations tends to zero as `m` grows.

3. **Shared partners and connectivity.** At `m = 3` the four fibre nodes share `296` partner classes (`148` between the members of each converse pair, `148` between an `|Aut| = 81` and an `|Aut| = 27` class), so the fibre is connected to itself through length-two metagraph paths. Conjecture: for every odd `m` the carousel fibre lies in one connected component of the line metagraph through such shared partners, and the shared partner classes between converse-paired fibre nodes are exactly the self-converse partner classes.

4. **Orbit-size law.** For every prime quotient `Q` of order `r` with `G = Aut(Q)` and the decorated fibre `Q[A_(w_i)]`, `tc(T_w) = H |Orb_G(w)| / (3^r |G|)`: the presentation count is proportional to the size of the decoration's `G`-orbit (THM-4804 (2) proves the cyclic case, where the orbit size is the field degree; the general case follows from Theorem G's group formula once `H` is constant on the fibre, which F4 gives). The same holds for child types of any common cyclic-triple type, with `3^r` replaced by the product of the child automorphism orders.

5. **F_1 reading.** For `q` child types with `q` not a prime power the necklace structure has no field; conjecture (definitional): the class counts `N_q(m)` are polynomials in `q` whose value at `q = 1` is 1 and whose coefficients count rotation types, so the one-type fibre is the `F_1`-point of the family of decorated fibres in the sense of counting polynomials.

Non-consequences: none of these bears on black-core connectivity (HYP-9346 (2)); (1)-(3) are statements about the carousel fibre only.
