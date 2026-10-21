# Independent arithmetic/tournament packet audit

Date: 2026-10-21. Verdict: APPROVED at the exact hashes in APPROVED-SHA256SUMS. No required mathematical or reproducibility correction found. This approves a research packet, not a literature-novelty claim or automatic canon promotion.

Scope: README.md, COPRIME-EXTENSION.md, the three final Python programs, six frozen receipts, and SHA256SUMS. The exploratory probe_joint.py is explicitly excluded. The requested `receipts` location is interpreted as the six JSON files actually present; there is no receipts directory. No Git commands or repository edits were performed. All audit writes are confined to session24-arithmetic-review.

## Proof findings

1. Exact clock maps: L(n)=3n+1 is a bijection modulo341, with inverse114(x-1), and L(4n+1)=4L(n). The golden quotient sends phi to81, not4. Substitution in the quadratic and period relations gives exactly the gcd341 quotient. Sharing a carrier does not identify the two clocks or integer trajectories.
2. Signed-section criterion: for odd m an even order modulo some prime divisor creates a flipped additive-order-prime stratum. Odd prime orders lift through odd prime-power factors, so the converse holds. The orbit count is correct. At341 the unit-only sections for2 and3 coexist with the distinct zero-stratum obstructions 2^5*31=-31 and3^15*11=-11. The subgroup need not contain global -1 for a nonunit to flip.
3. Joint classification: 81/4 has CRT coordinates(1,28), and28 has order15. Together with4 this gives QR11×QR31. Its eight nonzero-difference orbits pair into four antipodal choices. Every antipodal Boolean function on the four nonzero sign quadrants is ±chi11 or ±chi31, and each axis gives one more independent sign. Thus the four normal forms and four gauges each exhaust16 sets. All four common-out profiles were independently reconstructed and agree; arbitrary isomorphisms cannot merge the four classes.
4. Maximal-module criterion: the union of intersecting modules is a module; for distinct overlapping maximal proper modules, union=V. Their three nonempty regions are transitively ordered, contradicting strength. Hence maximal proper modules partition a strong tournament canonically. Transitivity gives uniform block size; translation of the block through0 by any member intersects and equals that block, making it a subgroup. Decomposability therefore forces a nontrivial proper subgroup module. Conversely any such subgroup is a nontrivial module. A subgroup is a module exactly when every other coset has a constant connection sign. The argument is not an inference from an imprimitive exhibited subgroup.
5. Coprime extension: every skew cyclic tournament is strong, since ±1 yields a directed Hamiltonian cycle. Coprime projections arise by integer multiplication, so every subgroup of C_m×C_n is rectangular. The three cases in COPRIME-EXTENSION.md cover every nontrivial proper rectangle: H_A nonzero/H_B proper; H_B full/H_A proper; H_A zero/H_B nonzero proper. Each displayed outsider really lies outside H, and the displayed vertices lie in H and have opposite signs. No primality or factor indecomposability is used. Oddness ensures a and-a are distinct and permits skew signs. The noncoprime C3×C3 diagonal hostile is valid and does not identify that carrier with C9.
6. Full automorphisms: maximum common-out count intrinsically recovers the dominant-coordinate blocks in all four forms. In a reversed form q(q-1)>q makes the block-majority direction unambiguous. The minority arcs are exactly the equality matching. Preservation forces every within-block coordinate map tau_a to be the same tau, and the internal arcs force tau in Aut(P_q); the induced map is in Aut(P_p). Conversely each product pair preserves all arcs. The factor backtracking search has complete domains, injectivity and all assigned pair constraints, fixes0, and finds exactly the square multipliers. This is a structural proof with finite-exact factor-group identification, not an unproved all-prime automorphism classification.
7. Signed subset: the acting maps are ordinary permutations of subsets, never subset complementation. Negation is anti, so the full signed group is G union(-G). The three corners have row and column multiplicities2,1, forcing both distinguished coordinates fixed individually. Affine coordinate maps fixing0 and1 are identities; an anti-map cannot fix the endpoints of a directed edge. Enumeration confirms the trivial stabilizer. Every2-set is fixed by its anti-reflection x→u+v-x, every singleton has its conjugate75-element ordinary stabilizer, and the empty set has the full stabilizer. Thus minimum size is exactly3, with signed orbit51,150.
8. Free action gives25,575 dividing Hamiltonian-path and two-leaf spanning-out-tree totals, but does not compute either total. The packet preserves this scope.
9. Fano: rank-four incidence image, affine base96, relative syndrome and one-line repair are consistent. The checker independently uses DP and every7! ordering, and identifies all16 zero-syndrome orientations with Paley7. This does not inherit the bracket escape predicate or distinguish the prime triple from another independent prime triple.

## Inherited-source status audit

Read the applicable AGENTS.md, START-HERE, tournament frontier, Research Protocol and active guardrails. The parent task's explicit read-only/no-Git boundary overrides mutating startup/closeout routines.

- procgen_numerology_20260923_odd_square_brackets.md §§1–2 proves exactly the doubling escape set; procgen_brackets_20260924_pairings_transitions.md §1 gives the multiplier interval and six integer pairs. The Sep24 correction in MISTAKES distinguishes the real cutoff from the integer cutoff; the packet preserves it.
- collatz_mod6_20260917_wild_typing.md §3 and scaffolding_audit.md §6 support the trunk, primitive-prime refutation, exact order criteria and85 base4/base2 separation. The packet does not reuse the source's stronger literature claims as new dependencies.
- golden_signed_relation_geometry_20261019.md explicitly labels its elementary proofs as a research note, not audited canon. Its §3 proves the specified quotient; §4 separates tripling order5 and golden-root37 order55 modulo121. The packet states those statuses correctly.
- THM-4592-the-eleven-squares-fibonacci-torus-is-the-period-five-collatz-cycle-set.md has scoped elementary equivariance/n=5 results and finite-exact all-odd bijectivity evidence only through22. MISTAKE-584 explicitly corrects the old all-odd statement, P/Q naming and factor2 projection twist. The packet correctly follows the correction rather than the overbroad historical sentence.
- THM-4591-integer-cycles-of-the-terras-map-to-period-301993-are-the-five-pythagorean-cycles.md labels its census FINITE-EXACT and separates period11 from Golay/torus11; seam_mills_20260925.md §5 proves the distinct parity mechanism. Both are represented faithfully.
- THM-4725-no-uniform-dilation-rule-for-prime-circulants-signed-set-burnside-and-the-p-mod-4-dichotomy-with-the-paley-rule.md is PROVED(pre-audit). Only its elementary fixed-section observation is needed; the packet does not elevate its other arithmetic or literature assertions.
- arithmetic_braids2_20260917_fano_code.md explicitly separates the proved incidence-code mechanism from FINITE-EXACT128-orientation Hamiltonian counts, including base96 and171/189. These are recovered, not newly claimed.
- THM-1745-leaf-graded-arborescence-filtration-721-shadow.md §1 contains the scoped free-action argument. The exact slug is essential because another unrelated file has the same numeric ID. The packet uses the correct one.

These checks verify the claims actually inherited here; they are not independent audits of every other theorem in those source files or external literature.

## Reproduction and independent hostiles

Ran check.py, paley_automorphisms.py and coprime_extension.py both normally and with python3 -O. Every output matches its frozen original byte-for-byte, including745 explicit main checks, factor search nodes52/572 and the64/8/4 extension controls. All11 original SHA256SUMS entries pass.

New hostiles.py uses no packet imports. Both modes agree over58,019 exception-based checks:
- every subset of every cyclic skew tournament at orders3,5,7,9: actual nontrivial module existence agrees with subgroup-module existence (decomposable counts0,0,0,4);
- transitive T3 has overlapping maximal proper modules, showing why the strong hypothesis matters;
- all341 conjugacy and inverse identities; explicit81≠4 and differing orders;
- every57,970 pair is stabilized by its ordinary-subset anti-reflection;
- direct CRT adjacency reconstruction gives all four published common-out profiles and fibre equivalences;
- a translation inside just one block is an automorphism in the aligned forms but fails in the reversed forms, isolating the matching constraint;
- both factors decomposable at orders9 and25: all seven nontrivial proper rectangular subgroups fail the module test.

The generalized theorem was checked logically above; finite samples are only hostile controls. No universal noncoprime theorem, universal product full-Aut theorem,341 Hamiltonian count, or Collatz/LRC conclusion is approved.

## Approval boundary

APPROVED-SHA256SUMS names exactly the12 approved packet files, relative to session24-arithmetic. Only those exact bytes are approved. Any changed byte needs reconciliation; probe_joint.py is not approved by this audit. AUDIT-SHA256SUMS separately inventories the audit artifacts.
