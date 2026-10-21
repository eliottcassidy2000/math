# Session24 signed-exception independent review

Status: signed packet APPROVED mathematically; structural all-order proof and finite verifier APPROVED after written review and both-mode replay. Repository was read only throughout the independent review. The historical audit stages below are preserved; the final approval supersedes their pending language.

Publication integration: source-scope qualifications were made explicitly historical and STRUCTURE.md was promoted to PROVED relative to the cited classical group theorems after approval. No new mathematical assertion was inferred from the finite controls.

## Signed packet

Approved input hashes:

- REPORT.md: 71151bc45770990336a62f62650ccc2341778e263dbe4886f258a72fb97a2c0b
- verify_small_actions.py: ae8982b0a2a8677124f1c9cf198b411c5c840fdc94002512c20dd716819e0146
- receipt.json and receipt-O.json: f4a8120b8b90e9220def8c185f24d409d4bfb233fb6ce956924fa4c11b15a9bc

The report's all-order results are correct deductions relative to Gluck's theorem, Feit–Thompson, and the standard affine structure theorem for primitive solvable groups. Ordinary and signed primitivity coincide for order at least three; the order-two caveat is necessary and present. The Paley proof covers every finite-field prime power congruent to three modulo four, without assuming the full Paley automorphism classification. The vertex-transitivity reduction by orbit complement is correct. Its strongness and equal canonical-module consequences follow.

The source-scope criticism was verified against the inherited pre-session checkout: THM-4743 left full classification open, HYP-9353 had OPEN only-if status, and its source section 4 nevertheless stated completeness. The original audit also read THM-2450's finite census through ten, but the root then identified its superseding current canon: [THM-2477](../../../../01-canon/theorems/THM-2477-rigid-strong-classification-through-sixty-two.md) and the [all-order bypass proof](../../merged_metagraph_finite_pure_atoms_20261020/FINITE_ATOMS_AND_RIGID_CLOSURE.md#3-rigid-strong-tournaments-are-completely-classified) prove strong path-rigid atoms exactly `{1,C3}`, arbitrary rigid stacks, and self-converse palindromes at every order. Therefore THM-2450's historical OPEN language is not the current frontier. “Rigid” means one ordinary Hamilton-path orbit, not trivial automorphism group. Current THM-2454 separately classifies strong pure-blue atoms as `{1,C3,T5}`. Exceptional triangle-substitution towers are different: `C3[C3]` has `h=39`, not one. These repaired scope statements do not retroactively validate the original unsupported signed-completeness sentence.

Integration condition satisfied: the signed report's remaining-OPEN statements are explicitly scoped to the inherited pre-session frontier and superseded by the new structural theorem. They are not objections to the new proof.

## Primary source check

Read the supplied primary Gluck PDF and independently located its Cambridge publisher version:
https://www.cambridge.org/core/services/aop-cambridge-core/content/view/9667B85D9DDD0FE911C8382287B86167/S0008414X00022550a.pdf/trivial-set-stabilizers-in-finite-permutation-groups.pdf

Printed page 65 Theorem 1 has precisely the degree/stabilizer table used in the report. The section's initial hypothesis is primitive solvable, with a regular elementary abelian normal subgroup and irreducible point stabilizer. Applying the twice-odd group-order filter leaves orders 6,10,42 only. The citation is verified; the derived tournament statement is not attributed as a theorem explicitly stated in Gluck.

## FINITE-EXACT replay

Ran the producer script normally and under python -O, writing only review-owned outputs. Both are byte-identical with the input receipt. The scope is all 1098 labelled tournaments of orders two through five, the named order-seven controls, and C3[C3]. Inspected the direct permutation test, subset setwise stabilizers, primitivity block criterion, and dual Hamilton enumeration/DP implementations. No removable assertions support correctness.

The Hamilton carrier controls are sound: P7 has (H,ordinary order,h,fH,signed path orbits)=(189,21,9,3,6), whereas C3[C3] has (3159,81,39,11,25). Fixed binary-palette orbits do not imply fixed Hamilton orbits. The packet does not certify the unavailable incoming order-nine or circulant-fifteen census.

## Independent hostile controls

hostile_controls.py is independent of producer imports. It enumerates every involution type on palettes q=1,...,9 and explicitly acts on nonconstant length-three cyclic necklaces. All 29 controls verify q'=(q^3-q)/3 and f'=f(q-1).

It additionally exhausts cyclic-invariant F3-linear subspaces in dimensions 3 and 5, and every reflection in each cyclic dihedral action. Every proper invariant kernel fails to contain image(I+P_tau), as required by the proposed full-classification kernel argument. These finite controls are not substitutes for the all-order proof.

Reproduce:

    python3 session24-signed/verify_small_actions.py --output session24-signed-review/replay.json
    python3 -O session24-signed/verify_small_actions.py --output session24-signed-review/replay-O.json
    python3 session24-signed-review/hostile_controls.py

## Historical structural candidate audit (superseded by final approval)

The communicated argument has passed an independent line-by-line conceptual audit so far. Key necessary details were sent to its author: minimal nontrivial signed blocks induce primitive block-stabilizer actions; kernel sign is positive because a negative map fixing two odd blocks would fix two vertices; use a sign-matched wreath supergroup before establishing modularity; identify a doubleton marker with its omitted vertex; propagate a kernel basis vector by ordinary transitivity. Awaiting written packet and frozen hashes before final structural approval.

## Written structural proof approval

Reviewed every section of session24-structure/STRUCTURE.md, including the added general palette bound and non-D-finite family generating function. The mathematical proof is APPROVED. It closes HYP-9353 relative to the cited classical group theorems and the independently verified primitive seed classification; there is no circular reliance on the bounded census or the desired classification. The finite verifier remained a separate evidence gate at this stage; it is passed in the final section below.

The crucial section 6 uses a regular coloring for the faithful ordinary quotient Q, not an assumed full quotient group. Every marker subset is ordinary-regular. Exceptionality produces an anti-stabilizer; all quotient anti-stabilizers preserving the block coloring coincide. Their entire lift set is exactly one K-coset. Universal marker coverage implies image(I+P_tau) lies in K. A fixed coordinate of the odd-degree quotient involution and transitivity of K's conjugation action force K to contain every coordinate vector. Independent triangle rotations then force module status. Every quantifier in this step is essential and present.

The majority quotient is used only as a sign-preserving envelope before module status is proved. Larger-group regularity implies smaller-group regularity in exactly the direction used. Minimal non-singleton blocks induce primitive block-stabilizer groups; a primitive subgroup implies the full induced signed tournament group is primitive. Nonexceptional block palettes and q>=3 exceptional palettes each contradict global exceptionality, so only the C3 primitive seed remains. This removes the old nonmodule-imprimitivity frontier rather than presuming it closed.

The C3 recursion, its counts for C3[R5] and C3[P7], and the sparse-family non-D-finiteness argument are correct. Sufficiency already follows from inherited tower substitution; no necessity result is recycled into that argument.

Minor presentation recommendation only: replace 'its moved ordinary orbits are exactly the E(T)-regular subsets' with 'the union of its moved ordinary orbits is exactly the set of E(T)-regular subsets' to distinguish the orbit and subset carriers. After audit, use PROVED relative to cited classical dependencies, rather than suggesting the verified primitive seeds remain an open mathematical condition.

Written structural version reviewed: SHA256 044500da260a52b87e49ce980293bedfca75686caa169ae2fb45fc1006a909b3. Any status-only or the specified wording correction should receive a new frozen hash before publication.

Review-generated finite evidence hashes:
- hostile_controls.py: 6fb845dca52c0701b87feb9bd1f11d8e91cdce36aba83861dfbb8cabf1341806
- hostile_receipt.json: 7465ac6ed1558e1c47e4411e7ac1db8ff6302631dbd259ca5a5e94accb9b1fc4

The extra linear hostile without transitivity has a proper nine-element kernel containing image(I+P_tau) in F3^3. This confirms that transitive propagation is indispensable; image containment alone does not force the full kernel.

Do not publish gluck.txt, which is a local primary-article extraction for source checking only.

## Structural finite-verifier review

The author applied the requested orbit/subset wording correction. Revised STRUCTURE.md SHA256 9a6d9b06ed2e6724c106c5cb101fa8ea95aeeefffc4c0cbb37ed877b0b6b9a1f is APPROVED mathematically. The proof is unchanged.

Reviewed verify_structure.py line by line. It independently enumerates all labelled tournaments of orders two through five, tests constructive orbit toggling and the palette lower bound, inspects named examples through nine, enumerates palette involutions through q=10, counts actual local subset triples for the C3/R5/P7 substitution controls, and enumerates the four cyclic-invariant ternary kernels for three blocks. The latter are signed permutation-subgroup controls, not a claim that each proper subgroup is the full automorphism group of a new tournament.

Normal replay exactly matches the producer receipt. Approved script SHA256 a2a323f78713c5a13ca0a3ed09ca3f3048ec2b0a59ecdfc667ebb29938d05d41; receipt SHA256 8f04dab6a2a0928690ab3f47d35cde03894806ec9cbbf61fc730da1e61c1c0d4. Actual C3[R5] and C3[P7] signed-regular subset counts are 15000 and 222264. Zero, diagonal, sum-zero and full ternary kernels give 420,288,216,0 signed-regular subsets respectively. These checks confirm the proof's finite controls without claiming to computationally establish the all-order result.

Optimized structural replay also matches byte-for-byte. All requested finite replay gates passed. Final audit verdict: APPROVED, with source-scope historical qualifications and no unresolved mathematical correction.

## Frozen publication revisions

The root independently checked the final scope clarification requiring a nontrivial block system. Its approved pre-promotion STRUCTURE.md hash was 675e8943bbf39317af11dce4ab3f64f97dc07375871bdb0ad5e54aca504432bd. The publication-only status and historical-scope revisions have these hashes:

- session24-structure/STRUCTURE.md: 62d039ad51079d719d9a8fb01a2bec17e47788103f24d3a0f7f7971a69f492c7
- session24-signed/REPORT.md: 6863f778c8ec66fef5db9e4c8d6be325c8b81a528c88c164b8b3ce342ec93a84

The earlier input hashes above identify the audited lineage, not the final publication bytes. Verifier and receipt bytes are unchanged.
