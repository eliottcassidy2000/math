# Scoped independent review: Session274

**Result: PASS WITH SCOPE CORRECTION, applied.** The first review found the
then-stated proofs sound. After the all-size paired-root compiler was added,
the complement-defect section required an explicit n≥4 scope; this is now
present in the proof and synthesis. At parent order2 an internal birth has
a one-edge child chord graph, so its defect is a cut, not a nonzero sector. This review does not certify novelty against the entire repository,
re-audit inherited theorems, prove Chan–Pak, or validate unrelated research.
Reviewed source hashes are in `REVIEWED_SHA256SUMS`.

## Proof checks

1. **Whole-profile rooted fan.** The path-block deletion argument is bijective.
   A new vertex can exit the new chain only to the root; each nonroot colour
   therefore has one terminal new block, while the root colour has a possible
   pre-root block and terminal block. Allocations between the blocks are legal
   even when the new labels interleave: old nonroot vertices beat every new
   vertex. The singleton-root case forces exclusion of its terminal block.
   This validates all three modes, including k=0 and t=1 (0^0 convention).
2. **Observer quotient.** With U=Q1−tQ0 and
   V=Q2−2tQ1+t²Q0, U=A−C and V=A+C. Thus the two division-by-two
   formulas are exact polynomial identities; the sequence is completely
   determined by lengths0,1,2 for this one operation language. No converse
   classification of rooted tournaments is claimed.
3. **Cut-space geometry.** The cut/cycle structure is explicitly inherited.
   For n≥4, connected old G forces the child gauge to agree with the old one
   up to a global constant. Internal birth imposes one common value on
   {predecessor,root} and one on its complement. First-root birth imposes one
   common value off the root. Modulo constants each kernel has dimension1.
   Rank-nullity therefore proves the exact two-to-one sector fan. The n≥4
   scope is load-bearing and correctly stated.
4. **Complement transport defect.** At a first-root birth the new chord mask
   is a newborn star, hence a cut. Internally the released predecessor-root
   chord has defect1 while all old chords have defect0; old connectedness
   precludes this being a cut. The proposed closed-path detector is valid.
5. **Forest means and moments.** A single path cover's chord graph is a forest;
   its edge equations are independent. A tuple of paths can instead create
   cycles; incidence rank and cycle parity give the exact probability in
   the all-moment formula. Conflicting uses of a repeated chord and backward
   compulsory path arcs are correctly removed.
6. **Explicit mean and automaton.** In the permutation expansion, selected
   runs along the numerical path must be uniformly directed. Their two signs
   sum to zero for odd runs and to two for even runs. Contracting selected
   edges leaves (n−|S|)! possibilities. The three-state weighted-run update
   correctly permits only completed even runs at acceptance.

7. **Paired-root compiler and synthesis.** Fixing the old root and newborn
   gauge bits to zero leaves exactly n−1 unknowns. Newborn chords read every
   unknown except the internal predecessor, whose bit is read from the newly
   freed predecessor-root chord. These independent coordinates give the
   explicit inverse at every n≥1. This is a bijection of framed tilings,
   not ordinary isomorphism classes or Hamiltonian counts. The synthesis
   preserves those distinctions and the remaining arithmetic selector.

## Independent computational route

Run from the packet directory:

    python3 review/independent.py

The script does not import the production scripts.

* Fan profiles are computed by enumerating Hamiltonian permutations on subsets
  and convolving them over named colours. E and L are computed independently
  by choosing the root-colour subset and explicitly enumerating its terminal
  or split paths. Full formula, cubic recurrence and integral recovery of A,C
  pass on all labelled n≤3 tournaments, every root, t=1..3, k=0..3:
  **87 input/root/colour cases, four fan lengths each**.
* The four order5 second moments are independently evaluated by summing over
  ordered path pairs, solving their binary edge equations and weighting by
  2^(−rank), rather than using subset path-count DP on orbit members.
  Results by representative0,1,2,3 are221/2,227/2,217/2,109.
* Actual born adjacency matrices, reframed by inserting the newborn before
  its root, confirm the two-per-child-sector fan for every root in every
  parent orbit at n=4,5: **24 orbit/root cases**.

* The paired-root compiler inverse is independently tested on EVERY child
  tiling at every root for parent n=1,...,5: **5,405 cases**. The test decodes
  gauge bits, extracts the parent, re-encodes, and recovers the child exactly.

`RESULTS.json` records this scope. Tests support the new proofs, not their
all-order conclusion by extrapolation. No larger numerical coverage census
was used. The unresolved exact-H selector is preserved throughout.
