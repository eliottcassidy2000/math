# Independent joint path/cycle interface audit

2026-11-05. **PASS for the exact five-coordinate transfer, linear minimality,
and physical equal-profile collision. No correction was needed.**

Frozen proof SHA256:
`d16a02db91d90d37fc0621602029cdc1a1c5e2811394b983e06bdd0c512883f3`.
Verifier SHA256:
`a9ef42072e32c0e91b5be55d8a295a1b399b6b50f4e6f149d387bcc3cd6e7822`.
Certificate SHA256:
`07b9475e64112f437a0f3a9707969a3e80fb2a9447838a5bb0d4dc49c6594cd1`.

## Mathematical scope and proof

The previously independently audited one-seam identities give the stated
integer linear update on (h,e,R,s,C). A path starting at the new root must
immediately use the sole reverse seam, giving s'=d*s. Every cycle must
use that seam, giving C'=E*s. The triangle seed avoids the sole order-two
cycle exception for every continuation. The old C coordinate is needed
as an immediate output, although it does not affect later states.

All five displayed input graphs really are reachable from the fixed
triangle seed by the stated one-step gadgets, and all five are strong.
The five continuation queries are legal for the same grammar. Their
evaluation matrix has determinant -24. Any rational d-dimensional linear
realization with fixed gadget maps and linear H,C outputs factors this
matrix through a d-dimensional space, regardless of how the input states
are encoded. Its rank therefore forces d>=5. The physical transfer gives
the matching upper bound. This is precisely a linear-realization theorem;
it makes no claim about arbitrary nonlinear encodings or algorithms.

The two collision graphs are genuine strong tournaments. Their shared
(H,e,R,deletion H,C) profile is (37,5,71,13,4), but their start counts are
9 and 8. The same singleton seam gives path state (103,32,108) for both,
and cycle counts 9 and 8. Thus the proposed old data do not determine the
next cycle response on all rooted tournaments. The proof does not assert
that these collision graphs belong to the triangle-seeded grammar, and
does not use that assertion in the rank theorem. No minimum dimension
for a model also required to output deletion H is claimed.

## Independent checks

The independent graph constructor puts gadget vertices before old vertices,
unlike the source. A separate endpoint recurrence computes h,e,R,s and
deletion counts; a memoized anchored-cycle recursion computes C. Literal
vertex-permutation path and cycle enumeration independently verifies all
25 rank-witness outputs and eight collision outputs. Rational Gaussian
elimination, rather than the source's signed permutation expansion, gives
determinant -24. The input-state determinant is -12 and the continuation
functional determinant is 2, independently checking the factorization.

Every rooted labelled gadget of orders one through four is physically
applied to each of the five reachable witness inputs: all 1,425 resulting
graphs agree with the full five-coordinate transfer. These include 590
zero-cycle controls and 320 cases with nonunit gadget deletion count.
Thus the checks include non-strong and zero-endpoint gadgets as well as
the strong rank witnesses. The proof, not these finite controls, supplies
the all-order identity.

Normal and optimized source runs are byte-identical and reproduce the
frozen certificate. The independent script is `audit.py`, and its result
is `AUDIT.json`; it uses the sibling independently audited endpoint helper.
Author source files were not edited. No unresolved mathematical,
implementation, or scope issue was found.
