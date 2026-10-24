# Independent review of the Session 46 synthesis

24 October 2026. **PASS after one pre-publication normalization correction.**
The result is a cohesive conditional theory, not an unconditional all-size
native wiring theorem. The reviewed synthesis and all antecedents are pinned
in `SOURCE_PINS.json`. No Git operation or checkout mutation was performed.

## Important correction found and resolved

The first draft of decomposition-theory Section 6 called F the complete actual
Laplacian perturbation while testing it against the already multiplied map
ell_b. That notation could apply the multiplicity twice. The final proof now
explicitly writes actual Delta L = 2^b F and tests the divided integer F against
ell_b = 2^b pi_2, or tests actual Delta L against the unmultiplied ell_0. If no
common positive power of two divides every column, it uses b=0.

This matters mathematically. In the weighted two-vertex graph L=4B, where
B=[[1,-1],[-1,1]], take actual Delta L=2B and h=(0,1) modulo 4. Then Lh=0 but
Delta L h=(2,2) is nonzero. Applying both doublings wrongly reports zero;
either corrected convention detects the failure. The new controls verify both
correct conventions on 225 weighted two-vertex cases and find 36 false
positives from the wrong convention. The native image equations already used
the correct single multiplicity and are unchanged.

The producer corrected this wording before first publication, refreshed its
manifest, and the earlier theory audit was repinned. This audit checked that
new freeze rather than relying on the superseded proof hash.

## Proof review, independent of the numerical outputs

- The common-target reduction needs disjoint actual source/target roles,
  equal marked target charge, complete current rows, exact integer margins,
  joint pins, disjoint available normalized resources, and no omitted changing
  factor. Those assumptions now appear explicitly in the synthesis. Minimal
  balanced-subtrade extraction gives conformal endpoint-union paths, but not
  a promise that atoms are small.
- Every pairing of a primitive balanced trade gives a minimal zero-sum
  sequence, hence support at most twice the Davenport constant. The reverse
  implication for a single pairing is false. The augmentation-ideal proof
  correctly gives D((C2)^4+C128)=132. The bound 264 has neither native sharpness
  nor arbitrary-target or order-independent scope. This is the classical
  p-group formula, not a new general Davenport theorem.
- The partial-sum buffer criterion cancels the precise margin and charge at
  each event and returns its buffer. Native supply, joint compatibility and
  availability remain hypotheses; group labels do not create resources.
- The full incidence theorem retains physical rows after combining overlapping
  roles. For its exact equality/congruence lattice, Dickson's lemma proves
  finiteness of conformally minimal nonzero elements in finitely many orthants.
  Subtracting an atom lowers the integer l1 norm. Sign-compatible partial sums
  stay between the endpoints, so all equalities, congruences and coordinate
  bounds remain valid.
- The proposed one-hot encoding is lossless when it retains every actual
  allowed record, owner value, resource slack and exact factor tuple. Integral
  occurrence consistency remains exact on cyclic factor graphs. It makes no
  assertion about a fractional local relaxation. Slot one-hot coordinates
  ensure each changed slot switches directly once, while an atom may still
  involve all slots. Finiteness does not imply an efficient compiler.
- The full28 move-lattice proof has the correct ambient zero-margin lattice:
  certified basis B is contained in the move lattice M, which is contained in
  the charge kernel K, and both outer indices are 2048. Thus B=M=K. This does
  not make its generators feasible in every unit-capacity box. The hostile
  charge fiber permits all 28 slots to vary, not only its eight changed slots;
  its two disjoint four-subsets require a direct eight-slot event.
- In the same labelled ground set, distinct group labels forbid a one-for-one
  charge-preserving exchange. Every non-singleton fixed-cardinality charge
  fiber therefore fails the ordinary matroid basis-exchange axiom. This narrow
  statement excludes no recoded or generalized matroid representation.

## Native claims, scope and provenance

The complete-domain, larger-circuit, compensation and move-lattice proofs and
independent audits were reviewed. The 43-clock and 37-clock eighteen-source
packets remain distinct. Minimum 37 and uniqueness up to sign apply only to
eighteen-source primitive relations in the supplied 28-option inventory. The
cost-36 census contains 1,129 relations but no primitive one; the synthesis
now says this precisely rather than apparently excluding all relations.

The literal three-cycle at source 2835 fails on its second application to the
initial path. The synthesis's misleading “third iterate” wording was corrected.
The reciprocal polynomial and half-odd-cost statements retain their precise
univariate/bivariate scopes.

Final compensation claims were reconciled with both producer and independent
auditor after their freezes. The six-slot lower bound applies throughout the
28-source supply. Event/work/peak minima instead belong to the specified
selected complete fibers: for the 37-clock circuit with buffer 2071, four events
and work 49 are attained together at peak 13; peak 12 requires five events and
work 57 in its 38-state menu. These minima are not silently transferred to
larger menus. All positive restored-buffer routes use reversible capacity-one
borrowing. Strict temporal nonreturn retains the original primitive atom.

Nodewise all-dyadic current preservation is not every ordinary per-port phase,
odd-primary preservation, kernel equality, or a fully occupied atlas rewiring.
The immutable remainder and complete shared-coordinate records remain real
conditions. The synthesis's opening now distinguishes the new fixed-order
native certificates from the conditional all-size theory and inherited
higher-order examples.

The inherited Gray-cycle and four-layer proofs were read in the repository.
Their shared resources violate the pairwise-disjoint-alternatives hypothesis.
The repair enlarges the source, changes resource inventory, releases its stated
auxiliary pins, and assumes injection is the only changing cross-slot condition.
This is a proof-level comparison with already inherited work, not a new incoming
commit claim or an independent replay of its complete native certificates.

## Executed checks versus reviewed evidence

This worker executed `check.py` in normal and optimized Python: **31,667 explicit
checks**, with byte-identical receipts. It imports no producer checker. Besides
the multiplicity controls, it exercises:

- 6,722 ordered endpoint pairs across 1,712 small exact lattice fibers, including
  explicit one-hot capacity and residual-factor encodings;
- a cyclic factor system whose fractional local constraints are feasible but
  whose integral exact assignment set is empty;
- the Gray resource skeleton at r=2,3,4, exhaustively enumerating 16, 256 and
  65,536 binary assignments and checking the enlarged vacancy routes;
- 3,956 failed basis-exchange instances in 98 nontrivial distinct-label fibers.

These are abstract finite controls, not new native realizations or substitutes
for the all-size proofs.

`verify_inputs.py` passes **643 checks** in both modes, with identical receipts:
33 reviewed source hashes, 192 manifest-listed artifacts from 11 packets,
32 byte-identical normal/optimized receipt pairs, and explicit synthesis claim
cross-checks. These are provenance/receipt checks, not producer replays.

The native enumeration, edge-level multiplication, complete dyadic filtration,
sharp-boundary census and larger lattice enumeration were executed by their
named independent auditors and inspected here through their proofs, scripts
where relevant, receipts and hashes. This worker does not relabel those runs
as newly executed. Likewise the earlier theory audit's 40 controls and
965,633/18,321-check replays are separate antecedents, not part of this worker's
31,667 checks.

## Reproduction

From the published audit-synthesis directory:

    python3 -B check.py --out normal.json
    python3 -O -B check.py --out optimized.json
    cmp normal.json optimized.json
    python3 -B verify_inputs.py --result-root .. --repo ../../../.. --out inputs.normal.json
    python3 -O -B verify_inputs.py --result-root .. --repo ../../../.. --out inputs.optimized.json
    cmp inputs.normal.json inputs.optimized.json

Before assembly, use `--workspace WORKSPACE` instead of `--result-root`; the
workspace layout is the sibling session46 packet directories. The explicit
`MANIFEST.json` lists this audit's deliverables only and excludes itself.
