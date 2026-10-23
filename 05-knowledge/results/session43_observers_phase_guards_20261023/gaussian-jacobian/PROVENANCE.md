# Provenance and recovered frontier

Read-only repository discovery on 2026-10-23. Checkout: session35-publication/clean-repo. HEAD observed at final provenance capture: 8b71b746f7a9d90378c702e01608207227b1f83e. The checkout can advance concurrently; source snapshots included below fix exactly what was inspected. No Git writes by this worker.

## What the definitions mean

GMC(n): for complex-coefficient P in n independent standard real Gaussian variables, E(P^m)=0 for every positive integer m implies, for each polynomial Q, E(QP^m)=0 eventually in m. Real-coefficient P makes the premise trivial by E(P^2)=0. For two real variables, normalized circular coordinates Z,W satisfy E(Z^a W^b)=delta_ab a!. This is NOT the signed order moment mu(F) used for the tournament coupling operator, and not a finite sum of path-cover radial moments.

Primary source: Derksen, van den Essen, Zhao, https://arxiv.org/html/1506.05192, Definitions 1.1/2.2 and Theorem 1.6. Their all-dimensional implication to JC does not imply GMC(2) => JC(2). Fixed-dimensional reductions must retain their dimension change.

Primary update: Christopher D. Long, https://arxiv.org/html/2607.18186v1, Theorem 5.1 and Section 7. The explicit three-variable witness has pure moments zero and mixed moments m! at every positive order; padding gives failure in higher dimensions. Section 7 does not settle GMC(2). We use this as a cited primary preprint, not as a fresh theorem or an audit of its entire historical account.

## Repository recovery

- prior-2022.md: THM-2022 claims NC2 and GMC(2), identifies the nullcone as the strict one-sided charge loci. Its updated formalization front door is 04-computation/lean/TournamentH7/TournamentH7/GMC2Main.lean. Source inspection confirms that gmc2 takes only the pure-moment premise and delegates to GMC2DvdKOmegaWiring.gmc2_unconditional. I also inspected the latter file's small-root/splitting-field assembly. I did NOT rebuild Lean, check the complete dependency graph or independently establish its claimed axiom report. Therefore this session records an inherited proof/formalization claim, not an independently certified GMC2 resolution. The public preprint and repo claim have different status and dates; neither should silently overwrite the other.
- prior-2019.md: affine charge-height supports and common radial multipliers have an exact angular/radial factorization. This is the established narrow family motivating preservation of integer charge.
- prior-1795.md: explicitly corrected from a detection-depth theorem to a four-cell recurrence-order conjecture. Leading-coefficient singularities block promotion of recurrence order to a uniform detection bound. This warning is directly relevant to replacing arbitrary powers by finite frame tests.
- prior-3290.md: repository's generalized spherical/Archimedean GMC3 witness family; claims a proof with exact checks and flags immutable independent audit pending. This worker has not redone it.
- prior-4727.md: mixed-weight torus quotient Jacobian identity, special equivariant planar rigidity, and a clear NOT PROVED JC2 boundary. Section 4 of THEOREMS.md independently rederives just its elementary planar weight-(1,-w) argument.
- THM-4626, inspected in source: dual-field surjectivity/local-nilpotence and a dimension-three statement equivalent to planar JC; an equivalence is not a proof of its common truth.
- prior-1810.md: identifier/title collision discovered during broad grep. The actual retrieved document concerns Paley skew determinants and binary-form analogies, despite a gmc-themed filename. It is not evidence for GMC. It contains broad structural claims not audited here. Included as a provenance warning against relying on filenames or shared numbers.
- session42-gaussian/RESULTS.md and the Session42 common-frame summary motivate simultaneous versus independent orbit compression. Our mixed-charge counterexample proves a separate angular loss independently.
- session42-sign-moments/THEOREM.md and session42-coupling/COUPLING_OPERATOR.md concern finite Walsh/observer operators. The new finite-frame theorem explains why their moment terminology alone provides no GMC transfer.

## Geometry and moduli boundaries

Spherical Gaussian polar decomposition is an actual integration identity. Hyperbolic weights in the planar torus theorem mean opposite-sign C* weights, not a metric of constant negative curvature. A finite-angle guard modulo 8 or 9 is not an infinite torus grading. Combining such angular guards gives lcm 72, while the tournament coupling gap 8/9 is an eigenvalue calculation in a different space. These distinctions preserve usable maps rather than dismissing the analogies.

## Root's concurrent bridge (not independently checked here)

Root reports the path-cover expansion F_T(x)=sum_m m! u_m(T)(x-1)^(n-m), hence a finite radial Gaussian realization replacing m! by E|Z|^(2m). This is compatible with the finite-bandwidth compiler, but does not satisfy the fixed-polynomial all-powers hypothesis of GMC by construction. See root's separate packet for proof and cyclotomic guards.

## Deliverable status

THEOREMS.md: fresh elementary scoped deductions, with inherited planar special case marked. check.py: 330 exact checks, normal and -O identical. No general Gaussian theorem, JC2 solution, native metagraph wiring rule, or literature-priority claim is asserted by this session.

## Snapshot link normalization

The copied prior-4727.md had one relative link that would break when relocated. Only that link was normalized to its immutable GitHub URL at commit 8b71b746f7a9d90378c702e01608207227b1f83e; the target was verified with git cat-file. The original repository theorem was not changed. Original unnormalized snapshot SHA256: eebf30e950564595a72f48b62304184930f7a48b3b6503473daebef96c9288df. All other copied prior-*.md Markdown links were checked; no other relative links were present. No mathematical text changed.
