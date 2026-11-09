# Session 232 final synopsis consistency review

**Final verdict: PASS.** The corrected canonical synopses THM-4987, THM-4988, THM-4989 and THM-4990, final README and SYNTHESIS are consistent with the independently audited proofs and retain the required hypotheses, statistics and inheritance boundaries. No mathematical blocker remains.

This is an independent editorial and mathematical consistency review, not a rerun of the heavy calculations or historical THM-4812 certificates. The reviewer made no repository edits. The exact reviewed working-tree snapshots are pinned below; these hashes, not a later mutable branch name, identify this approval.

## Corrections completed during review

- Reversed README divisibility was repaired: the cut-stabilizer order divides each coordinate, not the reverse.
- THM-4990, README and SYNTHESIS retain the specified singleton-padded star closures sharing the same marked core. No unrestricted closure-difference claim remains.
- THM-4989 credits LEM-003 specifically for the two rotation-fixed R5 Hamilton cycles and path quotient 3, rather than potentially attributing the five-triangle calculation to it.
- The copied subgroup/transfer audit's printed commands now use `python3`, not `python=3`.
- README takes k>=1 when calling C3^k nontrivial.
- The transfer synopsis and README explicitly use nonempty digraphs; THM-4988 likewise specifies a nonempty tournament.
- SYNTHESIS says that the full unlabelled action-graph automorphism group *need not* be the specified acting group, rather than incorrectly asserting that they never agree.

The audited mathematical arguments did not require repair; the first item was a real synopsis direction error and the remaining items prevent scope or attribution ambiguity.

## Accepted mathematical scope

1. **THM-4987.** The exact subgroup index, full-group vertex-orbit decomposition, strict q_A<23 threshold, prime valuations, threefold alternative and equality iff structure agree with the audited main proof and SHARP_STRATA. The dangerous subgroup indices are excluded by the prime-cost argument, not by an invalid closure-of-values claim. The 5/7/25 strata remain within their respective strict windows. No dependence on unaudited THM-4812(g) is introduced.
2. **THM-4988.** Endpoint quotients use vertex stabilizers rather than the full group. Minimal orders 6 and 7, witness masks/counts, generator-labelled action graphs and faithful graph encoding agree with the audited proof. Permutation/edge-reversal sign-evenness is not even-degree. The transported counted objects are admissible traces, never ordinary graph Hamilton paths. The encoding's leaves prohibit the latter.
3. **THM-4989.** The free path quotient, coarse packing quotient and weighted action-groupoid mass remain distinct. R5 and its joins support the exact 7/21 countercontrols. The tower uses actual restriction maps, with singleton path fibres and three-point packing fibres. Equal mass is not equivalence; stable images are not eventual constancy. The inherited all-subgroup mass exclusion uses the now-proved THM-4987, not OCF cardinality alone.
4. **THM-4990.** The ratio N*k_a/k_b, pointwise rather than setwise stabilizer, separate mass/coarse initial conditions and sharp tournament terminal bound are correctly retained. The digraph generalization does not import the tournament-only motion bound. Cut-state rigidity is necessary, not sufficient, and does not expand Session231 baseline placement.
5. **README and SYNTHESIS.** They preserve count/action distinctions, numerical universes, finite/all-order boundaries and inherited/new results. The statements do not identify an abstract action graph's full automorphism group with its chosen acting group or confuse graph-automorphism orders with path counts.
6. **Source scope.** Shao/Zhang are credited for methodological inspiration only. No geometric functor, p-adic explanation of 7/21, or imported geometric theorem is asserted. This review inspected SOURCE_NOTES but did not independently re-audit the external papers.
7. **Inheritance.** No synopsis contradicts THM-4812. Its full-group gap and strict vertex-transitive window are inherited with the later audit scope and repairs. All completeness claims remain OPEN; no external priority is claimed.

## Input fidelity and execution boundary

The published subgroup PROOF and SHARP_STRATA, boundary graph PROOF, and orbit-transfer PROOF and SOURCE_NOTES are byte-identical to the corresponding inspected producer originals. The packing proof differs from its audited candidate only by the status promotion recorded in PROMOTION.json. Replacing that one status line reconstructs SHA-256 `13ef71bafa1c035c18776082817fe9633020562f0a9479f44bfdb8a6146b19ec`, exactly the independent audit's candidate pin. All three independent audit reports were inspected. The inherited THM-4812 statement, strict-window lemma, later audit scope and LEM-003 cycle caveat were inspected directly.

REPRODUCE was checked for consistency with the publication adaptations. The inspected final REPLAY.json records PASS for 16 isolated-copy runs of eight programs, each normal/optimized pair equal and equal to its stored receipt. This reviewer inspected that receipt but did not independently execute those 16 runs. Git publication and final package checksum generation remain the parent's separate close-out responsibility.

## Exact final inspected SHA-256 pins

- `01-canon/theorems/THM-4987-subgroup-path-orbit-gaps-and-small-symmetry-cost.md`
  `a6ee9dffdbf42b3c467cbf21bd5a129534341556a38d7f789ca28e013ce4e9f7`
- `01-canon/theorems/THM-4988-endpoint-721-boundaries-and-faithful-automorphism-graphs.md`
  `4e3bf5a495e91f85a0c238e0488effdea052b197effbf9b9018c3f5cfcea38a9`
- `01-canon/theorems/THM-4989-colored-cycle-packing-coarse-orbits-and-stabilizer-mass.md`
  `514fa6058c239b511c0635e090215d46cb5e527a32e0544490f0cb2f6c03906f`
- `01-canon/theorems/THM-4990-stabilizer-corrected-hamilton-prefix-orbit-transfer.md`
  `fdb960c69e2c3c4c83b1ec5c5e155b30682d0c75a670ed39b4e3b850a53e43c2`
- `01-canon/theorems/THM-4812-presentation-fibres-7-and-21-are-absent-at-every-order-orbit-product-lemma-vertex-transitive-lemma-closed-window.md`
  `7c40f6586d9f39dbdddfc9899b4dde9cd8b00e4cea0cb7e987042d52adf466f9`
- `01-canon/theorems/LEM-003-aut-acts-freely-on-ham-paths.md`
  `09289c2ad4c12bcadf5d29a54a198e56b0ef8223853ac3f83a9caa931c65b9b8`
- `05-knowledge/results/session232_automorphisms_20261109/README.md`
  `e5520c3d5ff6c10d3563fc5493134164130f84da6eb48984ec3506d659a197a3`
- `05-knowledge/results/session232_automorphisms_20261109/SYNTHESIS.md`
  `76d4936b05473227d2adf861cc540ae9a53575fc0b08a94785decd04617f8dd1`
- `05-knowledge/results/session232_automorphisms_20261109/REPRODUCE.md`
  `806ebdc03f61c01c240256b67d2468ffe100fca7691a25412c43fa7c2e53811f`
- `05-knowledge/results/session232_automorphisms_20261109/subgroups/PROOF.md`
  `204d1d2060ade9e9f957f11f975a4035bc06731994c354077c12db4e69d8593f`
- `05-knowledge/results/session232_automorphisms_20261109/subgroups/SHARP_STRATA.md`
  `49ddbac6b8f1ca82f68777ab9c6f47ab77c976abb1d99f6939cb58aaa891556a`
- `05-knowledge/results/session232_automorphisms_20261109/boundary_graphs/PROOF.md`
  `5cc45fc818759076f10adbd304996fb63e7dbdb10dfbc3ef7e86c27813c6b4c5`
- `05-knowledge/results/session232_automorphisms_20261109/packing_mass/PROOF.md`
  `8b2688665bc645d178d67592a573c752e165668a2ac74034b484138f6e1a2fef`
- `05-knowledge/results/session232_automorphisms_20261109/packing_mass/PROMOTION.json`
  `9895f9b37caa43bd958ebc99688f97fbd0d90eb506cc30375f6c5de8e07b702c`
- `05-knowledge/results/session232_automorphisms_20261109/orbit_transfer/PROOF.md`
  `22ac3eef230caceecf3bf635bee19e82528db9756ee89e7826b8586075eb4955`
- `05-knowledge/results/session232_automorphisms_20261109/orbit_transfer/SOURCE_NOTES.md`
  `5bc9aa428c01497eb9f828c816b9a65189e596a173c56d0fc08d8656cffe3f70`
- `05-knowledge/results/session232_automorphisms_20261109/subgroup_transfer_audit/AUDIT.md`
  `abb722a109c087af630acd6ad16e04b209505394691cb7ef6d36820f41d8a789`
- `05-knowledge/results/session232_automorphisms_20261109/graph_audit/REPORT.md`
  `e2d40c7d37f4be6257644dedabfabd4be922e0b2d4a304752462dfff4af2f29b`
- `05-knowledge/results/session232_automorphisms_20261109/packing_audit/AUDIT.md`
  `3f72353530837b51fe2e4bee9175f9dce8b03e98ad984e87f8c421fa0b9f8765`
- `05-knowledge/results/session232_automorphisms_20261109/REPLAY.json`
  `fc72bb0aba9d2eb9c8a4e5d73b985bab087182b601313face42bd8abdfac26ee`
