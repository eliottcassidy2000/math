# The rank-11 AP-core is the achiral vertex where codex's rank-or-Euler frontier meets

> **CORRECTED BY MISTAKE-227.** The adjacent rows `d_k` span an index-`11!`
> path/Jacobi sublattice, not a basis of `L(AP)`. THM-2052's private anchor
> rows are a different thirteen-coordinate bounded code; no map makes the two
> objects equal. The 60 minimal vectors are signed off-diagonal Schur triples,
> not four-term additive energy. Rank/Euler convergence at the AP and the claim
> that rank twelve is the sole terminal are withdrawn; THM-2053 supplies a
> rank-eleven tangent-disk terminal. Read the remainder as idea provenance.

*boxeph-2026-07-21-S214. Owner: work incoming LRC progress (pull often); explore rank-11 / "11 private-
coordinate relations", through the lens "relations between a set of things is a tournament." Builds on
incoming codex THM-2052 (relation-rank descent, PROVED) + HYP-8841 (rank-or-Euler frontier, OPEN); mac-mini
S25 (rank-11 relation lattice L(AP)); my THM-1750 (reify ladder), S208 (braid arrangement), S211 (additive
energy = CT-moment), S212/S213 (mirror-parity / chirality = the Euler branch, HYP-8845). Verified in
`04-computation/rank11_relation_lattice_is_the_transitive_tournament_boxeph_S214.py`.*

## Incoming: the rank-or-Euler frontier (pulled this session)

codex just proved **THM-2052** ("finite LRC height forces a high-rank bounded relation code"): a
pigeonhole harvests **`dim W_{91⁶,3}(v) ≥ 11`** independent support-≤3 bounded relations on any
counterexample `v`, and a **12th** independent relation pins `⟨v⟩` by maximal minors (a finite Hadamard
box). So the open step is a single **rank 11 → 12 descent**, packaged as the **rank-or-Euler frontier**
(HYP-8841, OPEN): after strict-volume tests, each peel forces *either*

- **rank branch:** a bounded relation outside the rank-11 code (rank `11→12`, the finite terminal), *or*
- **Euler branch:** an owner-labelled **Euler survivor** — a component of `G_{1/14}` with `χ>0`, using
  **HYP-8845** (my mirror-parity: `χ≥2` on covering rows).

My recent chirality/Euler work *is* the Euler branch. This reflection places the rank-11 object it descends
toward on the same map, through the "relations = tournament" lens the owner asked for — with the honest
caveat codex's theorem forces.

## Three rank-11/12 objects — do not conflate (from the pull)

The number `11` names three different things; the whole program is the gap between the first two:

| | object | rank | status |
|---|---|---|---|
| **A** | harvested **bounded-relation code** `W_{Q,3}(v)` (THM-2052) | `n−s+1 = 11` (→12 terminal) | PROVED |
| **B** | ambient **resonance lattice** `{k∈ℤ¹³:k·v=0}` (single character) | `n−1 = 12` | structural |
| **C** | **AP resonance lattice** `L(AP)={a∈ℤ¹²:Σi·aᵢ=0}` (mac-mini S25) | `\|S\|−1 = 11` | exact lattice; later synthesis corrected |

The AP-core `C` is what the extremal `{1,…,12}` *is*; the descent lives in `A` inside `B`. My contribution
below is a structural reading of `C` — the target vertex — not of the descent carrier `A`.

## The rank-11 AP-core through the tournament lens (verified, object C)

Through "relations = tournament" glasses the 12 AP speeds `{1,…,12}` are the **transitive tournament
`T₁₂`**, and its rank-11 relation lattice `L(AP)` records that structure exactly:

- **Finite-index chain frame.** The rows
  `d_k=(k+1)eₖ−k e_{k+1}` (`k=1..11`) have adjacent support and a
  tridiagonal Gram matrix, but span a sublattice of index `11!`. Saturation is
  load-bearing and restores relations absent from the path frame. They are not
  THM-2052's private-coordinate anchor rows.
- **The long-range relations = additive energy.** The **minimal** vectors of `L(AP)` (norm 3, kissing
  number **60**) are exactly the `±`(additive triples) `vᵢ+vⱼ=v_{i+j}` (30 of them); kissing `= 2·(#triples)
  = additive energy` (mac-mini S25; my S211 CT-moment), which the AP uniquely maximizes.
- **The reify-ladder vertex.** `T₁₂` has score sequence `0,1,…,11` (the AP again), `char_A = x¹²` (fully
  nilpotent — the THM-1750 **nullcone vertex**), and its transitivity Vandermonde `∏(j−i)≠0` is the **braid
  arrangement `A₁₁`** (rank 11, S208). So `rank 11 = 12−1 = rank(A₁₁) = rank(L(AP))` — one rank-11 vertex.
- **It is achiral.** `{1,…,12}` is **palindromic**: `vᵢ+v_{13−i}=13` for all `i` (verified). So the reversal
  `i↦13−i` fixes the AP — the AP is **self-converse = achiral = the fixed point** of the chirality
  involution (S213), and the constant pair-sum `13` feeds codex's pair-sum wall `q∣vᵢ+vⱼ` (THM-2047 §2).

## The synthesis: the AP is the achiral fixed point where both branches meet

The two branches of the frontier meet at one object. The rank branch **descends toward** the rank-11
AP-core (the transitive/nullcone vertex, where the harvestable code saturates). The Euler branch is
**maximal at** the AP: the deep well's good set has `χ=24` = twelve mirror pairs (S212), and the AP is the
`ι`-**fixed** (self-converse) configuration. So:

> **The extremal `{1,…,12}` is the rank-11, achiral, self-converse transitive vertex — simultaneously the
> bottom of the relation-rank descent and the top of the Euler/mirror-parity count.** Wall A ("the AP is the
> unique extremal 12-core") is exactly "the rank-11 achiral transitive vertex is the unique optimizer" —
> the reify-ladder (THM-1750) and chirality (S213) statement.

This is why *both* branches point at the same place: the AP is the fixed point of the reversal `ι`, and both
"rank saturates" and "`χ` is maximal" are properties of that fixed point.

## The honest caveat (forced by THM-2052 and MISTAKE-224)

The tournament lens is **structural/diagnostic, not the proof carrier.** THM-2052 states the descent is
**signed coding theory, not Tournament Analysis**: "any binary orientation discards the signed coefficients
and relation height, so tournament scores, cycles, SCCs, and Hamiltonian paths do not encode" the code. And
MISTAKE-224 forbids the antisymmetric tournament as an LRC *positivity* engine. So:

- The "relations = tournament" reading illuminates the *structure* of the rank-11 AP-core `C` (its Hasse
  spine, its additive-energy shells, its achirality) and correctly identifies it with the reify-ladder
  vertex — a genuine organizing lens (S210/S386/S537).
- But the **rank descent itself** (object `A`) must keep the **signed integer coefficients and heights** a
  tournament orientation throws away. The only antisymmetric content that survives *into the proof* is the
  **pair-sum law** `q∣vᵢ+vⱼ` (THM-2047 §2) and the **`t↦−t` mirror-pairing** (my S212/S213, HYP-8845) — not
  a tournament orientation.

So: tournament lens for intuition and for locating the extremal vertex; signed coding for the descent;
mirror-parity (`χ`) for the Euler branch. They are complementary halves of the one frontier.

## Scope

Verified structural facts about the rank-11 AP-core (object `C`): tridiagonal/Jacobi Gram, 11 adjacent-pair
covering relations = transitive spine, minimal vectors = additive triples = additive energy, `char_A=x¹²`
nullcone vertex, palindromic/self-converse achirality. These are re-verifications and a synthesis, not a new
proof step; the value is **placing the rank-11 extremal on codex's incoming rank-or-Euler frontier** (its
Euler branch is my HYP-8845), identifying it as the achiral fixed point both branches target, and flagging
the THM-2052 caveat that the tournament orientation is a lens, not the signed-coding carrier of the descent.

Links: HYP-8855, THM-2052 (codex), HYP-8841 (rank-or-Euler), HYP-8845, THM-1750,
[[chirality-toothpicks-and-why-tournament-counts-are-even-one-lefschetz-parity-boxeph-S213]],
[[the-good-sets-reversal-symmetry-an-equivariant-mirror-parity-sharpening-of-the-chi-criterion-boxeph-S212]].
