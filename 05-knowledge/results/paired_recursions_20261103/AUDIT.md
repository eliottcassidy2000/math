# Verification and scope audit

**PASS: self-audit with independent algorithms.** No external reviewer or
formal proof assistant is claimed. The all-order arguments are in PAIRINGS.md
and FOUR_BITS.md; bounded checks provide controls, not their quantifiers.

## Computation universe and paths

- Exact Fraction Burnside counts for every n=0,...,40, plus separate
  full sums for the valuation hostiles40,48,64.
- New four-bit formulas at every n=3,...,200 against the inherited
  bounded-cycle modular evaluator; against full sums through40.
- Independent signed unordered-pair cycles of a representative of
  **every integer partition** through12, with the actual anti-action
  endpoint-swap condition. These agree with the separate channel formula.
- Every labelled Boolean orientation through6 assigned to a canonical
  ordinary isomorphism class. Each class has actual subset-DP H,
  converse class, and reachability strongness. Through5, literal
  permutation path counts provide an independent H check.
- Explicit commuting free actions at their maximum predicted ranks and
  H-preserving involutions with exactly the predicted fixed-point defect.
  The checker verifies generator involutions, pairwise commutativity,
  commutation with converse, and all free orbit sizes.
- Both H-weighted SCC recurrences checked against the complete class census.
- All256 order-five matching frames, all15 distinct full blocks, all455
  block triples, and a concrete four-block cover. This is exact optimality,
  not a greedy upper bound. The old reflection and MISTAKES entry are repaired.
- Weighted CRT counts checked for505 even moduli (all499 even values4..1000
  plus six extras), at five binary precisions, against direct gcd scans.
  The full even four-bit formula is independently checked on the same505.
- Actual81/27 tournaments reconstructed; common67 cyclic triples,
  H=462915 for both, automorphism orders81 and27, port8 after a sink gives
  2390013 versus1754919. This inherited hostile is freshly replayed.

Normal and python -O runs give byte-identical JSON receipts. All checks
raise explicit exceptions; Python assertion removal cannot disable them.
The mathematical verification uses the standard library. Pillow is used
only for the static diagram, which was rendered and visually inspected.
Timing output is separate from deterministic receipts.

## Proof audit

P1 necessity uses invariance of the R-fixed locus under commuting actions;
sufficiency explicitly numbers free translation blocks. It imposes neither
metagraph adjacency nor H preservation. The c=0 and d=0 boundaries are
stated. P2 keeps F_2 separate. P3 is fibrewise and constructs an attaining
involution, so its minimum is an iff/exact optimum rather than a lower
bound presented as a classification.

F1 and F2 classify all multiplicity patterns that can survive modulo16,
including n3,4,6,10. Gcds are odd; raising a gcd from1 to3 changes the
valuation by2, not1. Odd factors in multiplicity factorials are retained.
F4 divides by odd d only and retains the full modulus in the interval
boundary; radical reduction is used solely for the prime-exclusion set.

The weighted strong-component formulas are inherited. No independent
classification of all strong H values, no proof of the7/21 completeness
conjecture, and no generic low-dimensional source-preserving carrier follow.
The finite matching-frame obstruction is not equated with an LRC theorem.

## Attribution and literature

The April valuation theorem and October precision/channel theorems were
recovered before deriving the extensions. Incoming Session149 and the
tournament_joint_parity packet were read at a28303347b and are credited.
The new packet makes no literature-priority claim. OEIS A000568 was checked
for the standard initial counts, not used as an independent proof.
Applegate--Pol--Sloane's primary arXiv HTML and Royle et al.'s primary
abstract were opened; RECOVERY.md links them beside the limited cited claims.

The recovered June greedy-cover claim failed and is explicitly corrected.
The stronger totient valuation formula was a local scout guess, not an
accepted theorem; its first failure in the declared even range is recorded.
No other canonical claim is promoted from a computational prefix.

See manifest.json for hashes of artifacts and the actual local dependencies.

