---
id: THM-2477
title: "Rigid strong tournaments are C3 through m = 62; modular reduction and the prime minimal counterexample"
status: >
  PROVED (unconditional, 2 <= m <= 62: the only strongly connected
  tournament with H(T) = |Aut(T)| is C3) + PROVED (all m: Theorem R,
  modular reduction -- H = |Aut| descends factor-by-factor through
  the Gallai decomposition via the concatenation/Aut-sequence
  squeeze; Corollary R-prime -- every minimal counterexample is a
  PRIME strong tournament with |Aut| >= f(m) and imprimitive-or-
  intransitive Aut) + PROVED (all m: any counterexample with
  PRIMITIVE Aut is excluded -- affine bound beats the Busch floor
  for m >= 25, finite check below; corollary tc(T_p) > 1 for ALL
  primes p >= 5) + PROVED (LEM-003+: Aut acts freely on Hamiltonian
  paths for EVERY tournament, adversarially re-verified on all
  2,396,745 labeled tournaments m <= 7) + OPEN historical auxiliary
  statement (lemma ML: for m >= 63, an odd non-primitive group of order
  >= f(m) admits no prime strong invariant tournament; crossovers
  are co-finite). SUPERSEDED rigidity frontier (2026-10-20): the
  only nontrivial strong path-rigid tournament is C3 at EVERY
  order, by the Hamiltonian-bypass proof below; ML is not needed.
  Consequently rigid-SC(n)=pal13(n) holds for every n. The
  rotational non-rigidity corollary does NOT settle purity:
  regular cyclic T5 is pure-blue with tc=3.
  Cited: Busch 2006 (via THM-1370/MISTAKE-055), Dixon 1967,
  Burnside prime-degree, Feit-Thompson (primitive branch only).
source: kind-pasteur-2026-07-27-S135
depends_on:
  - THM-2453-palindromic-narayana-law-for-the-rigid-stratum
  - THM-1370-h-spectrum-omits-7-21-all-n
related:
  - THM-2454-pure-blue-alphabet-and-center-law
  - THM-2450-rigid-self-converse-classes-are-cyclic-ternary-towers
scripts: >
  04-computation/goal_thm_crossover_dp_kps_S134.py,
  aut_max_vs_busch_theory_kps.py, aut_max_witness_verify_kps.py,
  aut_max_window27_kps.py, aut_max_exhaustive_m7_kps.cpp,
  goal_adv_verify_kps_S135.py (adversarial re-implementation)
hash_basis: working-tree bytes (LF); produced by a four-agent
  workflow (bounds / structure / adversarial verify / synthesis)
  with the two-independent-implementation rule met at every
  computational layer (three at the DP layer)
---

# THM-2477 -- rigidity forces the 3-cycle, sixty-two vertices deep

## Current correction and supersession — 2026-10-20

The historical through-62 result is superseded by an independent,
**all-order** classification: the only nontrivial strongly connected
path-rigid tournament (`H=|Aut|`) is `C3`; the singleton is the other
rigid strong atom. See the [complete proof](../../05-knowledge/results/merged_metagraph_finite_pure_atoms_20261020/FINITE_ATOMS_AND_RIGID_CLOSURE.md#3-rigid-strong-tournaments-are-completely-classified).
Its mechanism is the one-orbit closing property of Hamiltonian paths,
contradicted by Benhocine–Wojda's (1983) Hamiltonian-bypass existence theorem
at orders at least seven, with direct small-order arguments. It does not
use the historical modular/group estimates below or lemma ML.

Therefore the full rigid class is the ordered `{1,C3}` stack class, and
rigid self-converse classes are its palindromes at **every order**.
The old `n<=62` frontier, the assertion that ML is needed to finish this
rigidity goal, and the bounded consequences below are **SUPERSEDED**.
The separately formulated auxiliary ML statement is not proved or
re-audited here; only its claimed necessity for this goal is removed.

There is also a logical correction: `tc>1` means **nonrigid**, not
**non-pure-blue**. Regular cyclic `T5` has `H=15`, `|Aut|=5`, `tc=3`,
and all three fibre tilings are blue. Thus the old claim that a
non-rigidity corollary settles THM-2454's rotational-purity question is
withdrawn. The [cycle reduction](../../05-knowledge/results/merged_metagraph_purity_and_strong_core_20261020/PURE_BLUE_STRONG_CYCLE_REDUCTION.md)
and the all-order rigid theorem are different statements; neither alone
classifies every odd nonrigid pure-blue strong atom. The separate
[uniform prefix-rerouting proof](../../05-knowledge/results/pure_blue_atom_classification_20261020/UNIFORM_PREFIX_REROUTING_CLASSIFICATION.md)
now supplies that classification: the pure-blue strong atoms are exactly
`1,C3,regular C5`. This independently resolves the old rotational-purity
question while leaving the nonrigidity-to-impurity inference invalid.

The original proof skeleton, computations, and historical consequences
below are retained without being relabelled as newly audited arguments.


**PROVED + GAP** as itemized. Proof skeleton (each step labeled in
the workflow synthesis, adversarially audited):

1. `|Aut|` divides `H` with free action (LEM-003+, all tournaments);
   `|Aut|` odd; strong forces `H >= f(m)` (Busch floor, cited +
   exhaustive m <= 7). So a rigid strong T needs `|Aut| >= f(m)` --
   an exponentially large odd automorphism group.
2. **Theorem R (all m):** for strong non-prime `T = Q[M_1..M_k]`
   (Gallai), `H(T) = H(Q) prod H(M_i)` squeezes between the
   concatenation injection and the Aut exact sequence, so
   `H = |Aut|` descends to `Q` and every module. **R-prime:** a
   minimal counterexample is prime (if `Q = C3` by minimality, the
   run-transfer 2-run walk gives strict surplus).
3. Dixon's `o(m)` (max odd subgroup of S_m) vs `f(m)`: crossovers
   on `[3, 62]` are exactly `{3, 9, 27, 54}` (three independent
   exact DPs). `m = 9`: the Lagrange window sieves to the full
   Sylow-3; all 4 invariant tournaments are `T9 = C3[C3,C3,C3]`,
   non-prime, killed by R-prime. `m = 27`: kernel pinch forces
   `Syl_3(S_27)`; all 8 invariant tournaments are `T27`, non-prime.
   `m = 54`: forced product group makes cross arcs uniform -- a
   stack, not strong.
4. Primitive Aut, any m: Burnside/affine bounds beat `f(m)` for
   `m >= 25` (checked to k = 5000 with monotone argument beyond);
   below 25 the DP covers. Corollary: `tc(T_p) > 1` for every
   prime `p >= 5` -- the rotational family is never rigid beyond
   C3, settling THM-2454's family question at the tc level.
5. Remaining case = **(ML)**: `m >= 63`, odd non-primitive
   `|G| >= f(m)`: show no G-invariant prime strong tournament.
   The frontier is `m = 63` (ratio 1.679, group shapes with an F21
   top factor -- the first non-Sylow stress test); the adopted
   mechanism is the subdirect-deficiency dichotomy (uniform block
   arcs => module or stack; proper subdirect linking costs order
   `3^t`), with the m = 162 index-3 subdirect subgroups as the
   known danger shape. Crossovers are CO-FINITE (every m >= 243
   achievable, >= 189 by upper bounds; largest upper-bound
   non-crossover found: 188), so counting can never finish the
   job: ML is the theorem's true content at scale.

## Consequences

- **THM-2453 upgraded:** rigid-SC(n) = pal13(n) (Narayana
  palindromes) is PROVED for all `n <= 62` -- every strong
  component of a rigid tournament on <= 62 vertices is C3 or a
  point by this theorem plus Theorem R's descent.
- **THM-2454:** the rigid stratum of the pure-blue alphabet law is
  now mechanism-complete through n = 62; only the nonrigid-atom
  completeness (T5-uniqueness among centers) remains census-bound.
- The proof pattern is the HYP-9030 atom-floor program executed in
  the H-world: floor (Busch) vs group order (Dixon), crossovers
  closed by classification -- the template for the Keller
  degree-semigroup atom question.

## Reproduction

The six ported scripts recompute: the o(m)-vs-f(m) DP (m <= 70),
the m = 3..7 exhaustive census (two engines), the m = 9/27/54
window closures, and the affine inequality sweep. The full
synthesis transcript with the step-by-step audit is in the S135
workflow journal (four agents, zero unresolved disagreements;
m = 27 kernel pinch flagged as thinnest step and independently
reconstructed).
