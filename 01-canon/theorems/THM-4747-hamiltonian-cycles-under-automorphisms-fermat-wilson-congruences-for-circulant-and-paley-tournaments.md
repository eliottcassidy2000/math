---
id: THM-4747
title: "Hamiltonian cycles and paths under (anti-)automorphisms: Fermat and Wilson congruences. (1) A nontrivial automorphism that maps a directed Hamiltonian cycle to itself acts semiregularly; so an automorphism of prime order q with a fixed vertex forces q | C(T), where C(T) is the number of directed Hamiltonian cycles. (2) For a circulant tournament Cay(Z_p, S), p prime, with multiplier stabilizer M_S = {a : aS = S} of order mu: C(T) = (p-1)/2 (mod p*mu); for the Paley tournament, C(P_p) = (p-1)/2 (mod p(p-1)/2). (3) The number rho_0 of Hamiltonian paths reversed by the reflection x -> -x satisfies rho_0 = mu (mod 2 mu); in particular it is odd, and every tournament with an involutory anti-automorphism has an odd number of Hamiltonian paths reversed by it. (4) Summing over all labelled tournaments on Z_p recovers Wilson's theorem. (5) FINITE-EXACT: all 264 circulant classes of odd order 3 <= n <= 23 satisfy (2)-(3); at the composite orders 9, 15, 21 no circulant class satisfies the naive prime-order congruence C = #(S cap Z_n^*) (mod n) (0/4, 0/16, 0/88). Exact Paley data: (H, C, H/|Aut|, rho_0) = (3, 1, 1, 1), (189, 24, 9, 9), (95095, 5505, 1729, 185), (1172695746915, 34358763933, 6857869865, 573057), (15760206976379349, 374127973957716, 62293308207033, 63871533) for q = 3, 7, 11, 19, 23."
status: >
  PROVED: statements 1-4 (elementary group actions; |Aut| | H is the classical free action, LEM-003). FINITE-EXACT:
  statement 5, from 04-computation/experiments/tournament_pseudoprimes_20261021/circ_counts.c (subset DP from vertex 0,
  64-bit counts; H(P_11) = 95095 and H(P_7) = 189 reproduce A038375, H(C_3[C_3]) = 3159 and the n <= 9 values reproduce
  the gentourng census tour_census.c). NOT CHECKED: novelty (statements 1-4 are routine orbit counting and are recorded
  because the repo had not stated them; the Paley cycle counts beyond q = 11 are not checked against the literature).
session: mac-mini-2026-10-21-fanopsp
source: 05-knowledge/results/tournament_pseudoprimes_camion_fermat_20261021.md
scripts:
  - 04-computation/experiments/tournament_pseudoprimes_20261021/circ_counts.c (+ circ_small.out, circ_13_17.out, circ_19.out, circ_21.out, circ_23.out, paley19.out, paley23.out)
  - 04-computation/experiments/tournament_pseudoprimes_20261021/check_congruences.py (+ .out)
depends_on:
  - LEM-003 (Aut(T) acts freely on Hamiltonian paths) — CITED in-repo
related:
  - 01-canon/theorems/THM-4746-camion-pseudoprimes-strong-tournament-counts-are-quasi-periodic-mod-p-korselt-criterion-carmichael-numbers-and-the-triangle.md
  - 01-canon/theorems/THM-4741-presentation-counts-single-view-classes-are-k1-c3-stacks-rotation-inequality-and-forbidden-fibre-sizes-7-and-21.md (tc = H/|Aut| and the rotation inequality)
  - 01-canon/theorems/THM-4730-strong-pure-blue-alphabet-complete-points-of-symmetry-are-k1-c3-t5.md (paths reversed by anti-automorphisms)
  - 05-knowledge/hypotheses/HYP-9357-paley-hamiltonian-path-maximality-and-the-kantor-galois-exceptional-primes-3-7-11.md
---

# THM-4747 — Hamiltonian cycles under automorphisms

## 1. The fixed-vertex lemma (PROVED)

**Lemma.** If `g ∈ Aut(T)`, `g ≠ 1`, maps a directed Hamiltonian cycle `Z = (v_0 … v_(n−1))` to itself, then `g` acts on `Z` as a rotation `v_i ↦ v_(i+k)` with `k ≢ 0`. So `g` is semiregular: all its cycles on `V(T)` have length `n/gcd(n,k)`, and it fixes no vertex.

* *Proof.* An orientation-preserving automorphism of a directed `n`-cycle is a rotation. ∎
* **Corollary.** If `T` has an automorphism of prime order `q` that fixes a vertex, then `q | C(T)`.
* **Corollary.** In general, `C(T) ≡ #{Z : g(Z) = Z} (mod q)` for every automorphism `g` of prime order `q`.

## 2. Circulant and Paley tournaments (PROVED)

**Setting.** `T = Cay(Z_p, S)`, `p` prime. `A = {x ↦ ax + b : aS = S} ≤ Aut(T)` has order `pμ`, where `μ = |M_S|`.

**Which Hamiltonian cycles have nontrivial stabilizers.**

* A non-translation in `A` fixes the vertex `b/(1−a)`, so by §1 it fixes no Hamiltonian cycle.
* A translation `x ↦ x+b` fixes `Z` iff `v_(i+k) = v_i + b` for all `i`. Since `k` is invertible mod `p`, this gives `v_i = v_0 + i·k^(−1)b`: `Z` is the *arithmetic cycle* with step `s = k^(−1)b ∈ S`.

**Counting orbits.**

* There are `(p−1)/2` arithmetic cycles. `x ↦ ax+b` maps step `s` to step `as`, and `M_S` acts freely on `S`. So they form `(p−1)/(2μ)` orbits of size `μ`.
* Every other Hamiltonian cycle has trivial stabilizer, so its orbit has size `pμ`.

Hence

    C(T) ≡ (p−1)/2   (mod pμ).

* For the Paley tournament (`μ = (p−1)/2`): `C(P_p) = ((p−1)/2)(1 + p·t)` with `t ∈ Z_(≥0)`.
* `t = 0, 1, 100, 200928444, 1478766695485` for `q = 3, 7, 11, 19, 23`.

**Paths.** `|A|` divides `H` because `A` acts freely on Hamiltonian paths (LEM-003). For the Paley tournament, `p(p−1)/2 | H(P_p)`.

## 3. Reflection-symmetric paths (PROVED)

**The signed group.** `φ(x) = −ax + b` with `a ∈ M_S` is an anti-automorphism. Together with `A` these maps form a group `A^±` of order `2pμ`. It acts on Hamiltonian paths:

* automorphisms by `P ↦ g(P)`;
* anti-automorphisms by `P ↦ φ(P)^rev`.

**Which elements fix paths.**

* An anti-automorphism that fixes a path squares to an automorphism fixing it, so it is an involution.
* `−1 ∉ M_S` because `−S ≠ S`. So the involutions are exactly the reflections `x ↦ 2c − x`.
* The reflections are conjugate under translations, so each fixes the same number `ρ_0` of paths.

**Burnside.** `(H + pρ_0)/(2pμ)` is an integer. Write `H = pμ·t` with `t` odd (`H` is odd by Rédei). Then

    ρ_0 ≡ μ (mod 2μ).

In general, for any tournament with an involutory anti-automorphism `σ`, `P ↦ σ(P)^rev` is an involution on the odd set of Hamiltonian paths. So `σ` reverses an odd number of them.

## 4. Wilson's theorem from tournaments (PROVED; a remark)

**Count pairs.** Count pairs `(T, Z)`: a labelled tournament on `Z_p` and one of its directed Hamiltonian cycles.

* There are `(p−1)!·2^(C(p,2)−p)` such pairs.
* The rotations fix exactly the pairs where `Z` is arithmetic with step `s` (`p−1` choices) and `T` is a circulant with `s ∈ S` (`2^((p−3)/2)` choices).

**Conclude.** Hence `(p−1)!·2^(C(p,2)−p) ≡ (p−1)·2^((p−3)/2) (mod p)`. Fermat gives `2^(C(p,2)−p) ≡ 2^((p−3)/2)`, and so `(p−1)! ≡ −1`.

* This is the classical necklace proof, weighted by tournaments.
* For composite `n > 4`, `(n−1)! ≡ 0 (mod n)`, so the summed cycle count has no pseudoprimes. This contrasts with the strong-tournament count of THM-4746.

## 5. Computations (FINITE-EXACT)

**Census.** `circ_counts.c` enumerates every circulant class (connection sets up to multipliers) of odd order `3 ≤ n ≤ 23`. For each it computes:

* `H`, from paths starting at vertex 0, times `n`;
* `C`;
* `ρ_0`, by a subset DP over the `±` pairs;
* `A = #(S ∩ Z_n^*)`.

`check_congruences.out`:

* Every one of the 264 classes has `H ≡ 0 (mod n)` and `ρ_0` odd.
* At prime `n`, every class satisfies `H ≡ 0 (mod nμ)`, `C ≡ (n−1)/2 (mod nμ)` and `ρ_0 ≡ μ (mod 2μ)`.
* At `n = 9, 15, 21`, no class satisfies `C ≡ A (mod n)`. The exact congruence there is the necklace one, `Σ_(e|n) μ(e) F(e) ≡ 0 (mod n)`, where `F(e)` counts the cycles fixed by the rotation subgroup of order `e`.

**Paley data.**

| `q` | `H(P_q)` | `C(P_q)` | `tc = H/(q(q−1)/2)` | `ρ_0` |
|---|---|---|---|---|
| 3 | 3 | 1 | 1 | 1 |
| 7 | 189 = 3³·7 | 24 | 9 | 9 |
| 11 | 95095 = 5·7·11·13·19 | 5505 | 1729 = 7·13·19 | 185 |
| 19 | 1172695746915 | 34358763933 | 6857869865 = 5·7·11·23·774463 | 573057 |
| 23 | 15760206976379349 | 374127973957716 | 62293308207033 = 3·167·4567·27225299 | 63871533 |

* `tc(P_q) mod q = 1, 2, 2, 4, 14`: no pattern is claimed.
* That `tc(P_11) = 1729` is a Carmichael number is NUMEROLOGY: `tc(P_19)` fails Korselt at 7.
* The circulant H-maxima are recorded in HYP-9357.
