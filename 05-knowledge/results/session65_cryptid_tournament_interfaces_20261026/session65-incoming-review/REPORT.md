# Session 65: correction overlay for the pinned incoming classification atlas

2026-10-26. **Bounded regression/source-lineage review; corrections required within the scopes below.** This is not an audit approval or refutation of the entire incoming packet. Original incoming files and all Session 65 proof bytes are preserved. No new mathematical research, new external source audit, or external build was undertaken by this lane.

## 1. Exact scope and provenance

Incoming commit: `a0645dc58ce10bd55e5a98a85d5b815993406ffa`.
Parent: `fcc112a69d7080db4cd619d6e7927be1992bc51e`.
The commit changes 29 paths: 775 insertions, one deletion. This review targets its classification and fifteen/twenty-five ledger, with README/reflection context; it does not rerun the Legendre census, multiperfect sieve, or validate all of `ATLAS.md`.

Below **P** means `05-knowledge/results/analogy_reduction_atlas_20261025/PROBLEM_CLASSIFICATION.md`; **L** means the same directory's `LEDGER_15_25_MULTIPERFECT.md`; **R** means `07-reflections/typed-analogies-legendre-tournament-and-problem-classification-klein-20261025.md`. Line numbers refer to `git show a0645dc58ce10bd55e5a98a85d5b815993406ffa:<path>`, not a moving checkout. `SOURCE_PINS.json` identifies exact incoming Git blobs/SHA-256 hashes and the local audited comparators used here.

Comparison evidence is the already completed `session65-reductions/REPORT.md` and `session65-audit-reductions/REPORT.md`, root `SYNTHESIS.md`/`CONNECTIONS.md`, the Goldbach source audit, and the weak-Collatz source audit. Their mathematical input snapshot is `6a8ad5d3bf01a232da4de41bcbff68390a9337ee`; their later external-source observations are separately dated/pinned in those reports. This review checked the cited local report bytes, not the external proofs afresh.

## 2. Scope corrections to pinned incoming claims

### 2.1 A common complete target is not mutual reducibility; a Pi2 presentation is not fixed-sentence strictness

**Exact claims.** P:32–34 says that every Pi1 sentence reduces to a Pi1-complete one, “so inside `Pi^0_1` all such conjectures are \"equivalent\" in the trivial sense”. P:38 says Pi2 is a form “which halting cannot encode”; P:59 applies “halting cannot encode it” to full Collatz. P:89–91 says that in the third band “only quasihalting encodes the statement”; R:29–30 repeats that limitation.

**Correction.** Sharing a reduction to NONHALT supplies no reduction between two source families. One fixed statement also must not be silently treated as a uniform language. As a simple control, a decidable nontrivial language and NONHALT both reduce to NONHALT, while NONHALT cannot computably reduce to the decidable language. Thus the quoted equivalence conclusion does not follow, despite P:34–35 correctly warning that no Goldbach–RH arrow is thereby supplied.

For a uniformly decidable predicate `Q(e,n,t)`, `forall n exists t Q(e,n,t)` has the standard infinite-beeping compiler; its negation is finite beeping. The full FINBEEP family cannot computably many-one reduce to HALT. That uniform separation **does not prove** that the single full-Collatz sentence has no alternative Pi1 formulation. The displayed arithmetical forms are presentations/upper bounds, not completeness or minimal-level theorems about individual conjectures. Ordinary cycle exclusion and full Collatz remain different predicates; weak Collatz does not exclude divergent nonperiodic orbits.

**Verified comparator.** Reduction report §§2–3, especially propositions 1–3 and lines 142–149; independent reductions audit §2. Weak-Collatz audit §5 confirms that the pinned 14-state and 10-state constructions detect nontrivial positive cycles, not arbitrary divergence.

**Replacement text.** “Counterexample-search machines give a common target, not mutual reductions. Pi2 presentations admit an infinite-beeping compiler; no ordinary-halting compiler works uniformly for the entire corresponding complete family. No strictness claim about a fixed Collatz sentence follows.”

### 2.2 JC(2) and DC(2) have different objects and a dimension-shifting implication chain

**Exact claim.** P:74 has one row `JC(2), DC(2)` with the definition “every Keller map of `C^2` is invertible”.

**Correction.** `JC_2` concerns complex polynomial maps `C^2 -> C^2` with nonzero constant Jacobian determinant (equivalently normalized determinant 1). `DC_2` concerns endomorphisms of the **second Weyl algebra** `A_2(C)`, with two coordinate and two derivative generators and their Weyl commutator relations. It is not another name for a Keller map of `C^2`, nor is it `DC_1`.

The cited classical arrows specialize to

    JC_4 => DC_2 => JC_2.

They give equivalence of the all-dimensional conjunctions, not a proved same-rank iff. For fixed input degree, the audited bounded coefficient formulations are decidable using inverse-degree bounds and effective ACF0 quantifier elimination; unbounded degree gives a Pi1 presentation. Merely saying “for each degree, first-order” omits the inverse-degree bound that makes the inverse variables finite.

**Verified comparator.** Reduction report §4.2 and independent reductions audit, “Bounded JC/DC, and the status boundary”, citing Bavula, Adjamagbo–van den Essen and Belov-Kanel–Kontsevich. This correction does **not** re-audit, endorse, or retract Gao's 2026 preprint or the repository's higher-dimensional counterexample record. No stale background claim about that newer record is substituted.

**Replacement text.** Split the row into `JC_2` and `DC_2`, retain their separate definitions, and record `JC_4 => DC_2 => JC_2` as a logical implication chain, not isomorphism or resource-bounded equivalence.

### 2.3 The LRC conclusion needs a verified theorem and its hypotheses, not a memory citation

**Exact claims.** P:72 attributes the finite-speed bound `n^{O(n^2)}` to Tao, then calls repository finite-box theorems “this reduction made partly explicit”. P:115–116 explicitly says the Tao statement is “CITED from memory … verify before relying on the bound”. R:27–28 nevertheless uses that bound to declare the decidable band.

**Correction of evidence scope, not a refutation of Tao.** The incoming note itself marks its cited theorem unverified; the finite-box identification is not established there. Do not label that source bridge independently proved or audited. This review makes no claim that Tao's bound is false.

The already audited Session 65 route is narrower and explicit: Malikiosis–Santos–Schymura Theorem A assumes the preceding runner case and gives a subset-gcd threshold. Singleton terms bound each primitive positive speed by `B(n)=binom(n+1,2)^(n-1)` in a remaining counterexample. Successive finite checks decide **LRC through r total runners**, with the lower cases discharged in order. The good-time set for each speed tuple is decided by exact rational interval intersections. This is finite decidability, not a feasible LRC(14) computation, a proof that LRC(14) is true, or an NP-hardness conclusion.

**Verified comparator.** Reduction report §4.1 and independent reductions audit, “LRC finite checking”, with [Theorem A and its hypotheses](https://arxiv.org/html/2411.06903#S1.SS1). The theorem statement/hypotheses were externally checked by that audit; its whole geometric proof was not rebuilt here.

**Replacement text.** “A separately audited inductive theorem yields an explicit finite decision procedure for LRC through r runners. The Tao attribution and identification with particular repository finite boxes remain unverified in this incoming note.”

### 2.4 Goldbach lineage and source-verification labels must not collapse historical stages

**Exact claims.** P:54 lists “BB(25) (2016, Lean-verified)” and BB(21) with “no explicit statement”; L:25 repeats “BB(25) (2016, Lean-verified); now BB(21)”. P:7–8 says independent-verification caveats are repeated, but several relevant rows omit them.

**Correction.** The source-audited original gist is dated 2016-05-08 and has 27 states (`000..026`). Leng's dated 2024-10-17 comment describes the 26-state simplification and 25-state compression. The incoming 2016/25 entry conflates those historical stages. The wiki and pinned Lean 25-state tables really are isomorphic by state renaming; the mathematical Goldbach connection is instead a halting iff, not graph isomorphism with a conjecture.

The newer 21-state primary source has an explicit endpoint: project-defined `ConjectureMachines.GoldbachConjecture`, using Mathlib `Nat.Prime`. The audited sources therefore do not support an unqualified “no explicit statement”. This does not turn source inspection into an independent local Lean rebuild or remove the source's independent-verification caveat.

Exact upstream pins already independently inspected by the Goldbach audit:

* `lengyijun/goldbach_tm@6cd292062516a0a14ea1b34f2ab75154cae7ab1e`; `GoldbachTm/Tm25/TuringMachine25.lean` Git blob `3c1ede90cc821687f6108a38bb08d2508de44d43`.
* `google-deepmind/cryptids@6c28d8524fa008c6fb5aadb759df74a55ce96423`; `Cryptids/Propositions.lean` Git blob `f6835d20093f49812e70061a0e59244e924e8600`.
* `calr3/tinycryptids@3c43f97768d4ea06d70e0eb154ba88060bc36dec`; `lean/TinyCryptids/WeakCollatz/Main.lean` Git blob `b3dde7a4f34fe0bc407dbdef0c4826623f803ac1`.

The inherited audits distinguish source/formalization claims, upstream successful CI metadata, finite independent tests, and unperformed local Lean/Rust builds. Reuse those distinctions. The source-reported Sloane(6,4), odd-perfect(37), and odd-multiply-perfect(23,3) caveats should likewise not be silently dropped. Machine sizes are model-dependent construction sizes, not minimum-state theorems or NP-completeness classes.

**Verified comparator.** `session65-audit-goldbach/goldbach/REPORT.md` §§1–2, `SOURCE_RECEIPT.json`; `session65-audit-weakcollatz/AUDIT.md` §5; reduction report §4.5.

**Replacement text.** “Historical 27-state Goldbach construction (2016); later 25-state compression with Lean sources; newer 21-state claim with its distinct verification status. These are separate dated implementations, not established minima. The 21-state source explicitly defines its Goldbach endpoint; no local rebuild is claimed here.”

## 3. Compatible content and small qualifications—not blanket retractions

* The incoming **main organizing aim is compatible** with Session 65: machine size is not a mutual-reduction or isomorphism theorem; weak/full Collatz are separate; Sloane and Erdős have different forbidden digits; odd perfect and odd multiperfect are different existence predicates. These correct distinctions should remain.
* L:13 calls 15 “the last exceptional exponent”, and L:18–19 “the last zeroless ternary power”. The audited fifteen lane explicitly says **largest known exceptional exponent, conjecturally final**. That qualifier is necessary: the full exclusion beyond 15 is the open Sloane conjecture. The verified expansion `2^15=1122221122_3` and finite cited Saye computation are unaffected. P's “smallest known machine” heading also mixes historical sizes with later source-reported claims; preserve dates and verification labels rather than implying minima.
* R:47, “No tournament-to-machine composition exists to study”, should be bounded to the inspected historical material or to the missing **native reduction**. Session 65 now supplies an ordered compositional tape code and a labelled control-table codec. Neither simulates the TM using native Hamiltonian K moves, nor proves headline-problem equivalence. This is a scope update, not a claim that incoming researchers should have known a later construction.
* The ledger's repeated-number warnings remain useful. The later exact multiperfect arc-family identity and Frobenius-defect bridge add typed maps; they do not identify Goldbach's 25 with tournament order 25 or identify Sloane's 15 with the Ogg prime set.
* Mahler's integer-orbit reduction, the biased-coin classification, tournament-density undecidability, and the broader analytical/topological analogies were **not re-audited here**. Their absence from the conservative Session 65 classification is not a refutation. In particular, no new Mahler hierarchy or global topology verdict is asserted.

## 4. Suggested publisher notice

> This preserved incoming packet is historical, not covered by the Session 65 audits. Read the Session 65 incoming-review overlay before using its classification claims: it corrects the common-target/fixed-Pi2 scope, separates JC2 from DC2, supplies the audited LRC finite-check scope, and repairs Goldbach machine lineage and verification labels. Original bytes are retained.

The publisher can link this notice from its new README without editing incoming files. The bounded review is complete once this overlay and its pins are sealed; it requests no new research lane or delayed external build.

## 5. Execution and preservation boundary

This lane wrote only `session65-incoming-review/`. It did not edit, stage, commit, reset, rebase, or explicitly fetch in the shared repository. Exact incoming reads used `git show`; final pin reads disabled lazy fetching with `GIT_NO_LAZY_FETCH=1`. Source hashes for the local root/proof comparators are sealed below; no proof text was changed.

One broad read-only `git grep` over historical result paths encountered absent partial-clone objects and implicitly attempted a network fetch. It was interrupted; final output was `fatal: unable to access 'https://github.com/eliottcassidy2000/math.git/': Proxy CONNECT aborted` (exit 130). It was not retried, and the broad search was abandoned. This was an interrupted object-access attempt, not evidence that a fetch completed or a permission denial. Already retrieved targeted incoming files sufficed; no conclusion depends on the missing objects. No permission escalation was requested. No denied or pending action is needed to complete this bounded overlay.
