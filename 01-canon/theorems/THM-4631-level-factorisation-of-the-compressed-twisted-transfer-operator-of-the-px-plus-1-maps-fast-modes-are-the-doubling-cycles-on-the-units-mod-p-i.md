---
id: THM-4631
title: "Level factorisation of the compressed twisted transfer operator of the maps x/2, (px+1)/2 (p an odd prime; p = 3 is THM-4629): on the class means modulo 2^j p^k, det(I - M_(j,k)(s)) = [1 - u/2 - v/2] prod_(i=1)^k [1 - (u/2)^(L_i)]^(phi(p^i)/L_i) with u = zeta^alpha e^(-s), v = zeta^(alpha+beta) e^(s(theta_p - 1)), L_i = ord_(p^i)(2): the graded pieces of the p-adic depth filtration are the Cayley graphs of multiplication by 2^(-1) on (Z/p^i)^x (disjoint cycles of length L_i, phi(p^i)/L_i = [(Z/p^i)^x : <2>] of them), the odd branch lowers depth, and the 2-adic (de Bruijn) part contributes only the even density 1/2; hence the fast modes (decay 1/2 per octave) have multiplicity [(Z/p^i)^x : <2>] at level i, which is 1 at every level iff 2 is a primitive root modulo p^2 (Artin primes that are not Wieferich), jumps by the factor p at level 2 for the Wieferich primes 1093 and 3511, and is 2 at level 1 for p = 7 along the quadratic-residue split {1,2,4} | {3,5,6} of the Paley tournament"
status: >
  PROVED: statements 1-3 (the filtration argument of THM-4629, proof 2(b), with p in place of 3; the only new input is
  the cycle decomposition of multiplication by 2^(-1) on (Z/p^i)^x into cosets of <2>). PROVED: the corollaries of
  statement 4 (elementary group theory: 2 primitive root mod p^2 iff primitive root mod every p^i; ord_(p^2)(2) = ord_p(2)
  iff p is Wieferich; <2> = squares mod 7). FINITE-EXACT: the identity verified symbolically (sympy, symbols u, v) in
  12 cases (p = 3, 5, 7, 11, 13, 17; sizes to 18) and by exact rational determinants at three random rational points in
  25 cases (sizes to 343, including p = 7 at k = 3 with three doubled levels and j up to 3). NUMERICAL: complex
  determinants at random (u, v) in 8 further cases including the Wieferich primes 1093 and 3511 at level 1 (relative
  error at double precision). CITED (from memory, not re-fetched): Artin's primitive-root conjecture and its constant
  0.3739558 (Hooley 1967 conditional on GRH); the Wieferich primes 1093 and 3511 (Meissner 1913, Beeger 1922); the
  Paley tournament. For p = 3 the dynamical meaning is THM-4629's (tau, omega are finite for every n under the Collatz
  conjecture's hypothesis of convergence to 1); for p >= 5 the operator is the formal compression of the finite-horizon
  step-count identities (orbits may diverge), and the theorem is a statement about the operator. OPEN: unchanged from
  THM-4629 (the 3-adic closure); whether the fast modes of the px+1 maps are detectable in integer data (for p = 3 they
  are not, THM-4629 statement 7).
session: opus-2026-10-13-S30
source: 05-knowledge/results/graph_theory_as_common_carrier_20261013.md
scripts:
  - 04-computation/experiments/graph_carriers_level_factorisation_20261013.py (+ .out: Part 1 the determinants; Part 2 ord_(p^i)(2), cycle counts, Artin and Wieferich; Part 5 the de Bruijn spectrum)
related:
  - THM-4629 (the case p = 3; the operator, the filtration proof, the scalar law and its resonances)
  - HYP-9261 (the graph carriers of the repository's threads), HYP-9260 (the Riemann-sphere reading of the scalar law)
  - 05-knowledge/results/ramanujan_heegner7_trivial_cycle_20261006.md and collatz_paley_bridge_20261001.md (the {1,2,4} | {3,5,6} split), THM-4628 (residue colourings modulo conjugation-fixed primes)
  - 05-knowledge/results/graph_theory_as_common_carrier_20261013.md
---

# THM-4631 — the fast modes of the twisted transfer law are the doubling cycles on the units

## Setting

`p` an odd prime, `T_p(n) = n/2` (`n` even), `(pn+1)/2` (`n` odd); `theta_p = log2 p`; a character `(alpha, beta)` modulo `m`, `zeta = e^(2 pi i/m)`; `u = zeta^alpha e^(-s)`, `v = zeta^(alpha+beta) e^(s(theta_p - 1))`. `W_(j,k)` is the space of functions of the class `(r mod 2^j, rho mod p^k)`. The compressed transfer operator (THM-4629 statement 2 with `p` for `3`):

`(M v)(r, rho) = (u/2) [v(r/2, rho 2^(-1)) + v(r/2 + 2^(j-1), rho 2^(-1))]` for `r` even, `= (v/2) [v(s_r, c(rho)) + v(s_r + 2^(j-1), c(rho))]` for `r` odd, with `s_r = (pr+1)/2 mod 2^(j-1)` and `c(rho) = (p rho' + 1)/2 mod p^k`, `rho'` the odd representative of `rho` modulo `2 p^k` (for `j = 0` both branches act on the `p`-adic classes with weight `1/2`).

`L_i = ord_(p^i)(2)`, `phi(p^i) = (p-1) p^(i-1)`, `c_i = phi(p^i)/L_i = [(Z/p^i)^x : <2>]`.

## Statements

1. **(Exact class-mean identities; PROVED.)** THM-4629 statement 1 holds verbatim with `p` in place of `3`: halving maps the class `(r, rho)` onto `(r/2 mod 2^(j-1), rho 2^(-1) mod p^k)` and the odd step maps `(r, rho)` onto `((pr+1)/2 mod 2^(j-1), c~(rho) mod p^(k+1))` with `c~(rho) = (p rho' + 1)/2`, because `2` is invertible modulo `p^k` and `n = rho' (mod 2p^k)` determines `(pn+1)/2` modulo `p^(k+1)`. `M_(j,k)` is the compression of these identities under the `p`-adic closure `H_k` of THM-4629.
2. **(Level factorisation; PROVED.)** For every `j >= 0`, `k >= 0` and every character,
   `det(I - M_(j,k)(s)) = [1 - u/2 - v/2] prod_(i=1)^(k) [1 - (u/2)^(L_i)]^(c_i)`.
3. **(Graph reading; PROVED.)** Let `V_i` be the functions of `p`-adic depth `<= i`. On the graded piece `V_i/V_(i-1)` (`i >= 1`), spanned by the additive characters `e(a rho/p^i)` with `a` a unit, the odd branch is zero (its image has depth `<= i-1`) and the halving branch acts as `(u/2)` times the permutation `a -> a 2^(-1)` of `(Z/p^i)^x` tensored with the de Bruijn averaging `2E_j` (spectrum `{1, 0, ..., 0}`). The permutation is the Cayley graph `Cay((Z/p^i)^x, {2^(-1)})`: `c_i` disjoint directed cycles of length `L_i`, the cosets of `<2>`. Its characteristic polynomial is the displayed level factor. On `V_0` (functions of `r` alone) the determinant is the scalar law.
4. **(Corollaries; PROVED.)**
   - (a) The fast modes of level `i` sit at `|u/2| = 1`, i.e. `Re s = -ln 2` (decay `1/2` per octave), at the frequencies `f = alpha/m - l/L_i (mod 1)`, each with multiplicity `c_i = [(Z/p^i)^x : <2>]`.
   - (b) `c_i = 1` for every `i` iff `2` is a primitive root modulo `p^2` iff `2` is a primitive root modulo `p` and `p` is not a Wieferich prime. Among the primes below `100` this holds for `p = 3, 5, 11, 13, 19, 29, 37, 53, 59, 61, 67, 83` (12 of 24; Artin's constant predicts `9.0` on average).
   - (c) For a Wieferich prime (`2^(p-1) = 1 mod p^2`: `1093`, `3511`) `L_2 = L_1`, so `c_2 = p c_1`: `1093` has `3` fast modes at level 1 and `3279` at level 2; `3511` has `2` and `7022`. For every other prime `L_i = L_1 p^(i-1)` and `c_i = c_1` for all `i`.
   - (d) `p = 7`: `L_1 = 3`, `<2> = {1, 2, 4}` is the subgroup of squares, the two cycles are the quadratic residues and the non-residues, the out-neighbourhood and in-neighbourhood of `0` in the Paley tournament on `7` vertices and the two halves of the code of `ramanujan_heegner7_trivial_cycle_20261006.md`; the level-1 factor is `(1 - (u/2)^3)^2`. In general `c_1 = 2` with the two cycles equal to the residues and non-residues iff `ord_p(2) = (p-1)/2` (`2` is then a square): `p = 7, 17, 23, 41, 47, 71, 79, 97` below `100`, and the Paley tournament exists for those `= 3 mod 4`: `7, 23, 47, 71, 79` (Part 2; the other indices below `100` are `3` at `43`, `6` at `31`, `8` at `73` and `89`).
5. **(Verification; FINITE-EXACT, NUMERICAL.)** `graph_carriers_level_factorisation_20261013.out`, Part 1: the identity of statement 2 as polynomials in `(u, v)` in the 12 symbolic cases `(p, k, j) = (3,1,0..2), (3,2,0..1), (5,1,0..1), (7,1,0..1), (11,1,0), (13,1,0), (17,1,0)`; exact equality at three random rational `(u, v)` in 25 cases up to size `343` (`p = 3` to `k = 4`, `p = 5` to `k = 3`, `p = 7` to `k = 3`, `p = 11, 13, 17, 19, 23` to `k = 2`, `j` to `3`); double-precision agreement in 8 numeric cases including `p = 1093` and `3511` at level 1. Part 5: `E_j` has the single nonzero eigenvalue `1/2` and `rank(E_j^q)` falls to `1` at `q = j` for `j = 1..7`.

## Proof

**1.** As THM-4629, proof 1, with `p` for `3`: `n = n_0 + 2^j p^k t` gives `(pn+1)/2 = (pn_0+1)/2 + p^(k+1) 2^(j-1) t`. ∎

**2, 3.** THM-4629, proof 2(b), with `p` for `3`. The filtration `V_i` is preserved by both branches: halving because `rho 2^(-1) mod p^i` depends on `rho mod p^i`; the odd branch because `c(rho) mod p^i` depends only on `rho mod p^(i-1)` (the `p` lifts of a class modulo `p^(i-1)` have odd representatives differing by multiples of `p^(i-1)`, so `p n + 1` differs by multiples of `p^i`, and `2` is invertible modulo `p^i`). Hence on `V_i/V_(i-1)` the odd branch vanishes and `M` acts as `(u/2) (2E_j) (x) Pbar_i`, where `E_j v(r) = [r even] (1/2)(v(r/2) + v(r/2 + 2^(j-1)))` and `Pbar_i` is the pullback of `rho -> rho 2^(-1)` on the characters of exact conductor `p^i`, i.e. the permutation `a -> a 2^(-1)` of the units modulo `p^i`. The orbits of multiplication by `2^(-1)` on the group `(Z/p^i)^x` are the cosets of the subgroup `<2>`, each of size `L_i`; so `Pbar_i` is a permutation matrix of cycle type `(L_i^(c_i))`, with characteristic polynomial `prod_(cycles) (1 - x^(L_i))`. `2E_j` has the eigenvalue `1` once (its rank-one part, the even density) and is otherwise nilpotent (Part 5; for `j = 0` it is the scalar `1`). The determinant over the graded piece is therefore `prod_(eta) (1 - (u/2) eta)` over the `L_i`-th roots of unity of each cycle, `= (1 - (u/2)^(L_i))^(c_i)`. On `V_0` the operator is THM-4629's level-`j` 2-adic compression with the odd-branch target `s_r = (pr+1)/2`; its only nonzero eigenvalue is `(u + v)/2` (the two lifts of every image class enter with equal weight, so `M|V_0` has rank one on the constants' complement... precisely: `M|V_0` is the weighted de Bruijn operator whose non-constant part is nilpotent, THM-4629 proof 2(b) for `j >= 1` and the scalar law for `j = 0`), giving the factor `1 - u/2 - v/2`. The determinant of a block-triangular operator is the product over the graded pieces. ∎

**4.** (a) `|u/2| = 1` at `Re s = -ln 2`; the arguments of the `L_i`-th roots of unity give the frequencies. (b) If `2` is a primitive root modulo `p^2` it is one modulo every `p^i` (the kernel of `(Z/p^i)^x -> (Z/p^2)^x` is a `p`-group generated by `1 + p^2`, and `2^(L_2)` generates it: standard lifting), so `c_i = 1`; conversely `c_1 = 1` needs `2` primitive modulo `p`, and `c_2 = 1` needs `L_2 = p(p-1)`, i.e. `2^(p-1) != 1 (mod p^2)`. The list below `100` is Part 2 of the output. (c) Wieferich means `2^(p-1) = 1 (mod p^2)`, so `L_2 | p - 1`, hence `L_2 = L_1` (always `L_1 | L_2 | L_1 p`); then `c_2 = p(p-1)/L_1 = p c_1`. Otherwise `L_2 = p L_1` and inductively `L_i = p^(i-1) L_1` (Part 2 checks `i = 2` for `p < 100`). (d) `2^3 = 8 = 1 (mod 7)`; `{1, 2, 4}` are the squares modulo `7`; the Paley tournament on `F_7` joins `x -> x + d` for `d` a nonzero square. The general statement is the definition of `c_1 = 2` with `<2>` of index `2`, i.e. `<2>` = the squares. ∎

**5.** Computation (script above). ∎

## Scope and boundary

- The theorem concerns the compressed operators; the step from the finite compressions to the integers is THM-4629's open 3-adic closure (for `p = 3`) and is not touched. For `p >= 5` the maps are expanding on average (`theta_p > 2`) and `tau` is not defined on divergent orbits; the class-mean identities are then statements about finite-horizon step counts and the operator is the object.
- The multiplicities `c_i` are the only place where the arithmetic of `p` enters the spectrum; the scalar law knows `p` only through `theta_p`. In particular the Artin and Wieferich dichotomies are properties of the fast, undetectable modes, not of the slow lines; the theorem explains why the `3`-adic hierarchy of THM-4629 is "silent" (one cycle per level, every eigenvalue at modulus `1/2` per octave) and predicts the same silence with multiplicity for every `p`.
- Related readings: `graph_theory_as_common_carrier_20261013.md`, sections 2 and 6 (the Cayley/de Bruijn carrier; the Mersenne clocks `p | 2^n - 1 iff L_1 | n` as the same cycles read as periods).
