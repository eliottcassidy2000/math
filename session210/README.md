# Session 210: exact carriers for the two-gap spectrum

**Status: independently audited structural results and finite checks.
Chan–Pak spectrum completeness remains OPEN.** No literature-priority claim
is made. The classical oddness theorem, exclusions `7,21`, odd-cycle formula,
and one-reversed-arc count are inherited, not new results of this session.

## Four packets

1. **[Prior-work and structure atlas](structure-atlas/STRUCTURE_ATLAS.md)**
   recovers the native meaning of the indices `3,10`, distinguishes actual
   tournament realizability from arbitrary counting objects, and reconciles
   old scope errors. The latest inherited coverage is independently audited
   through `2^32`; the larger `2^35` fixed-bank result is source-reported from
   one largest-scale run, not independently rerun there. This session does
   not extend that coverage.

2. **[Projective cycle carrier](projective-cycle-carrier/PROOF.md)** realizes
   `(H(T)-1)/2` as the total number of cells of a natural regular CW complex:
   nonempty disjoint odd-cycle packings with signs, modulo simultaneous sign
   reversal. Its cell grading retains the complete packing polynomial; its
   canonical mod-two cohomology class, for a nonempty carrier, has largest
   nonzero power at exponent `nu-1`, where `nu` is the maximum disjoint-
   triangle packing size. Triangle SCC words give real
   projective spaces, and adding one `S4` component gives the stated
   projective-space–sphere wedge. Equal counts do not determine topology:
   at order six, `H=9` can give four points or a circle.

3. **[Local gates and closure responses](graded-spectrum-gates/PROOF.md)**
   gives a genuinely graded path/cycle correspondence. For
   `Gamma(T): a -> b -> T -> a`, the arc `a -> b` is recognized by
   `d+(a)=d-(b)=1`; deleting the gate bijects Hamiltonian cycles of the
   closure with paths of the core. The gate is unique from order four.
   More strongly, the ordered path-cover polynomial obeys

       P_Gamma(T)(q) = q[P_T(q+1)+qP_T(q)-P_T(q-1)].

   For known core order `n`, the first `floor((n-1)/2)` responses
   `H(T), H(Gamma(T)), ...` recover `P_T` exactly. This is reconstruction of
   an invariant, not of the tournament or a minimum-query claim.

4. **[Coherent permutohedra and subspace cells](coherent-permutohedron/PROOF.md)**
   bijects tournaments with orientations of the permutohedron agreeing on
   parallel edge classes; Hamilton paths are exactly sinks. Face sink counts
   factor over their blocks. Dropping coherence admits every sink count
   through `n!/2` for `n>=2`, whereas generic linear objectives admit only one.
   For the classical extreme-arc-reversal family, feedback energy gives a
   Gaussian-binomial refinement and an explicit bijection from
   **field-decorated** nontrivial paths to all finite-field subspaces.
   Bare paths index the cells, not all their field-valued points.

## The cross-packet consequence

Let `m_f` count permutation words with `f` forward adjacent pairs, and let
`F_T(z)` sum `z^dim(face) * sinks(face)` over all nonempty permutohedron faces,
including the whole polytope. The gate and geometric packets prove

    P_T(q) = sum_f m_f binom(q+f,n),
    F_T(z) = sum_f m_f (1+z)^f.

Both transforms are invertible at known order. Therefore the finite nested-
closure responses recover the **entire face–sink incidence polynomial**.
This is an immediate synthesis of the audited packet theorems, not a new
independent census. It recovers dimension-aggregated incidence counts, not
the sink count of each individually labeled face.

The [sharp cross-carrier boundary](cross-carrier-boundary/PROOF.md) proves
that even **all** gate responses and the full face-incidence polynomial
cannot recover projective topology. Two explicit order-six tournaments
have `H=17` and the same `P(q)=q^6+8q^4+8q^2`, but their projective carriers
are eight points and a circle plus four points. Order six is the first
possible failure. This is an exact independently checked obstruction, not
a conjecture.

The projective construction supplies a complementary positive cell count
for the odd index. Its join arithmetic is `a+b+2ab`, not addition. None of
these transports removes the central obligation: realizing every allowed
integer by an actual tournament. No independently classified object family
with proven spectrum `N_0 minus {3,10}` was found.

## Verification and reviews

Each packet contains its own verifier, exact outputs, source manifest and
release hashes. Normal and optimized Python replays agree. Independent
reviews are [atlas](structure-atlas/audit/ATLAS_REVIEW.md),
[projective carrier](projective-cycle-carrier/audit/REVIEW.md),
[gates](graded-spectrum-gates/audit/REVIEW.md), and
[permutohedron](coherent-permutohedron/audit/REVIEW.md).

The projective census covers every labeled tournament through order six;
gate, face–sink and cover controls cover every labeled core through order
five, with further separately specified hostile and field tests. These
finite tests support, but do not replace, the all-order proofs.
