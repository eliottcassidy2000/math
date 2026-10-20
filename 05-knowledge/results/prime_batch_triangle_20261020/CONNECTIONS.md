# Retained connections and their limits

These are object-level transfers or explicitly labelled analogies. They are
not implications between unsolved problems. Current proved packets outrank
the historical drafts reviewed for provenance.

## 1. What was recovered rather than re-proved

* [Arithmetic braids, divisor note](../arithmetic_braids_20260917_divisors.md),
  DB1--DB3, already proves the exact F=S+U classification. F counts proper
  nontrivial divisors, S counts those squarefree, U counts those prime;
  U is a subset of S, so this is not a disjoint partition. The whole exponent
  box, including removal of its apex, is essential.
* [Finite divisor-profile controller](../seam_prime_balance_20260925.md)
  already has eight live exponent profiles and an absorbing multiplication
  reject state. It also proves the reject state is not preserved by Collatz
  halving. Rebuilding that classifier would give a weaker result, not progress.
* [Divisor balance family](../collatz_mod6_20260917_divisor_balance_family.md)
  preserves unreached nonsquarefree divisors, distinguishing exponent-box
  occupancy from support size. Our smooth-partition bound uses the analogous
  extra coordinate of weighted occupancy, with a different counting predicate.
* [Source orientation and boundary contexts](../source_orientation_boundary_contexts_20261018k.md)
  and [carry transport](../carry_descent_transport_branches_20261018l.md)
  retain the original source row and complete six-colour boundary constraints.
  The finite 375-branch survival result remains unchanged. Neither a full
  six-colour transfer nor a positive formal affine slope proves entry.
* [Defect deficiency](../defect_transitions_20261020/strata.md) concerns
  D=2k+1-h, not delta=3f-s. Its bounded exceptional graph is useful as a
  presentation method. Its variable must not be substituted for our slack.

## 2. Incoming work read during this session

Baseline was 7d21b9fdcf; the session fetched 2e8ff23dc3, 71f9f0bcba and
b3d59a8495 while the first tests ran.

1. [Cycle-changing descent](../black_descent_20261020c/cycles.md), especially
   its prime-orbit blocks and additive boundary lemma: it can change the
   Hamilton cycle to recover a descending endpoint pair. The same flexibility
   is not available for a fixed Collatz source. Our replacement changes an
   *upper-bound partition*, leaving the physical valuation word untouched.
2. [105/17 maps](../black_descent_20261020c/105_17.md): composing the matching
   incidence and junction observations loses more than either intended
   interface alone. This motivated testing the square of the triangle transfer
   before attempting a geometric chain as a forcing gadget. Its complete
   boundary relation becomes full at four colours, a decisive stopping reason.
3. [Odd blue spectrum](../odd_blue_spectrum_hollow_cubes_20261020/ODD_BLUE_SPECTRUM_AND_HOLLOW_CUBES.md):
   the hollow anti-involution frame retains matching completion bits.
   The packet's preserved invariant is blue degree, not the original blue
   adjacency. Our two synchronized components give another precise test:
   forgetting relative cyclic orientation preserves f,s,|G|,H but changes B.
4. [Paired-black leaves](../paired_black_leaf_reduction_20261020/BLACK_QUOTIENT_LEAF_CLASSIFICATION.md)
   reduces mixed leaves to rigid outer words and a strong deficit-two center.
   Its historical header leaves that center open; the later
   [completed classification](../deficit_two_centers_20261020/README.md) and
   [audit completion](../deficit_two_centers_20261020/session15-deficit-centers/AUDIT_COMPLETION.md)
   supersede that boundary and give the unique center S4. We read the cutoff,
   full residual census, exact leaf series and approval record. The resulting
   leaf/slack intersection and synchronized degree series are in README §3.

The newest result actually narrows our grammar: non-self-converse leaves are
zero-slack words modulo reversal; mixed leaves form the mirrored S4 subfamily
inside slack three. Neither synchronized slack-six center is a leaf. This is
a genuine subset relationship between independently defined predicates.

## 3. Connection contracts

### Named valuations -> smooth continuation kernels

**Source:** processed cycle-length multiset. **Map:** strip prime factors
greater than the current p, but retain every remaining exponent, multiplicity,
old/new colour and total committed degree. **Preserved:** every future
gcd interaction and recurrence weight after its already accumulated part.
**Lost:** which stripped prime produced a kernel. **Sidecar:** accumulated
rational coefficient and actual degree. **Test:** the exact smooth-multiset
DP from the inherited packet; the new two-colour partition product and common
degree prime-bank Gram matrix. Keeping only prime support fails at lengths
3 and 9 when a future query sees gcd with 9.

### Two oriented triangles -> a two-state matching orientation

**Source:** two directed 3-cycles and a bijection between their vertices.
**Map:** compare the cyclic orders through that bijection, giving epsilon=+/-.
**Preserved:** whether the matching preserves cyclic orientation.
**Lost:** the shared C3 phase, vertex identities and all exterior contacts.
**Sidecar:** matched vertex addresses and full cross-orientation matrix.
**Test:** S_+ and S_- have identical H=45 but blue-orbit counts 9 and 3.

The same quotient S3/C3 governs the relative orientation of two ordered
three-colour palettes. This is a common small moduli object, not a graph
isomorphism between the tournament and the unit-distance strip. At four
colours one spare colour supplies intermediate triples in both orientation
classes. Thus the proposed palette obstruction disappears unless a separate
geometric context really forces the common three-colour palette.

### Group order -> divisor balance, with a retained clock

**Source:** full permutation automorphism group G. **Map:** its odd order M,
then the proper divisor counts F(M),S(M),U(M). **Preserved:** the exponent-box
predicate F=S+U. **Lost:** the action, element orders, ports and cycles.
**Sidecar:** the cycle index or explicit marked action. **Test:** C3^2 and
C9 have the same order nine but minimum possible tournament slacks zero
and 21. C3^3 and C27 both give the balanced order 27 but slacks zero and 246.
The mixed balanced shape first becomes possible at slack 15, witnessed by
order 315. This is a new use of the old balance theorem, not its rediscovery.

### Fixed Collatz word -> ordered multiset-block envelopes

**Source:** original positive odd n and its exact exponent word. **Map:**
retain consecutive blocks, their exponent multisets, original cylinder and
integer clock; evaluate each decreasing rearrangement only as an affine upper
bound. **Preserved:** a sufficient strict descent certificate. **Lost:** within-
block positions in the envelope, but not in the authenticated source word.
**Sidecar:** original word/residue and block cuts. **Test:** one-block failure
at source27 is repaired by four blocks; source63 needs five at its first exit.
Refinement is monotone because it restricts allowable rearrangements.

This is an actual source-preserving paid-rule mechanism. A theorem that every
source enters a finite paid prefix is still the missing universal statement.
Prime factor shape does not determine the tested budget: balanced 27 needs
four first-exit blocks, unbalanced 63 needs five, and unbalanced 75 needs one.
No conclusion from F=S+U alone is available.

## 4. What should be tried next

* **Narrow continuation width:** retain smooth-part mass and use a sharper
  restricted-partition estimate or an adaptive cut. A polynomial-width claim
  for the unchanged prime-batch interface is now excluded by Proposition 2.
* **Deduplicate slack-six presentations:** canonicalize the four endpoints
  and gap words while keeping the actual homogeneous triangles and their ports.
  The current local recognizer is exact but counts presentations with overlap.
* **Force a palette before using orientation:** the strip gives a clean
  unit-distance host on which to test extra cross-interior contacts. Reject
  candidate quotients whose complete six-port relation has already become full.
* **Bound refinement demand:** use the exact DP (14) to search for long
  source-specific sequences of needed cuts, then ask for an arithmetic bound
  on these cuts. First-exit four-block universality is refuted by63; any later-
  prefix or growing-budget conjecture needs a separate statement and test.

Existing META-PATTERNS cards used: compose witness fibres before scalar
projection; retain operation and carrier; test the residual independently.
No generic new card is promoted from one session's analogies.
