# Independent audit of the Eisenstein–Bernstein tensor proof

Status: **ACCEPTED**. The arithmetic route proves the same tensor lemma as
the frozen elementary induction, without using that induction. It does not
independently establish the tournament identities, an arbitrary-modulus
analogue, a new stability cutoff, or any Riemann-hypothesis conclusion.

Reviewed source: `../stability/EISENSTEIN_BERNSTEIN_TENSOR_PROOF.md`, also
present as Section 2 of `../stability/INDEPENDENT_TENSOR_MATCHING_REVIEW.md`.
The current exact bytes are identified in `CURATED_FILES.json`.

## 1. Polynomial and Bernstein identities

For `P(x)=sum_I c_I x_I`, a one-variable monomial coefficient pair `(a,b)`
contributes `P^2+(P')^2=(a+bx)^2+b^2`. Tensoring this identity gives

    J(P;x)=c^T tensor_i [[1,x_i],[x_i,1+x_i^2]] c.

To reduce to nonnegative coordinates, replace the polynomial by
`P_tilde(y)=P(s_1 y_1,...,s_r y_r)` with `s_i` the signs of the original
coordinates. This preserves coefficient parity and preserves the jet norm
at the corresponding points, since mixed derivatives change only by signs.
It is not a claim that the jet norm of a fixed polynomial is an even function.

With Bernstein weights `((1-x)^2,2x(1-x),x^2)`, the three matrix
coefficients are exactly

    Q0=I,   Q1=[[1,1/2],[1/2,1]],   Q2=[[1,1],[1,2]].

The target `1+x^2` has coefficients `(1,1,2)`. A tensor coefficient with
`k` type-2 coordinates must therefore be at least `2^k`.

## 2. Type-2 transformation and type-0 slices

The integer matrix `M=[[1,1],[0,1]]` satisfies `Q2=M^T M`. Apply `M` on
every type-2 coordinate of the coefficient tensor. Its mod-2 input is the
single tensor basis vector `e1 tensor ... tensor e1`; `M e1=(1,1)`.
Consequently, after transformation its parity is the tensor product of:

- `(1,1)` at each type-2 coordinate;
- `e1` at each type-0 or type-1 coordinate.

The type-0 and transformed type-2 identity factors split the quadratic
form into slices, each carrying `Q1` on its remaining type-1 coordinates.
Fix every type-0 slice index to 1 and allow every type-2 index. There are
exactly `2^k` retained slices. Each has an odd top coefficient and all other
coefficients even. Other slices need not have this parity, but their
quadratic forms are nonnegative: `Q1` has eigenvalues `1/2` and `3/2`,
and its tensor powers are positive definite. Discarding them is valid.
The argument also covers zero type-0 coordinates, zero type-2 coordinates,
and zero remaining type-1 coordinates.

## 3. Phase average and F4 nonvanishing

Let `zeta=exp(i*pi/3)`. Independently choose every `z_i` from
`{zeta,conjugate(zeta)}` with equal probability. The single-coordinate
Gram matrix of `(1,z_i)` is `Q1`: diagonal entries are 1, and the average
of either off-diagonal entry is `Re(zeta)=1/2`. Independence therefore gives

    E |P(z_1,...,z_m)|^2 = c^T Q1^(tensor m) c.

All evaluations lie in `Z[zeta]`, where `zeta^2-zeta+1=0`. Modulo 2,
this ring becomes `F2[z]/(z^2+z+1)`, a field with four elements because
the quadratic has no root in F2. Both possible phase residues are nonzero.
An odd-top/even-rest polynomial reduces to their product, hence has nonzero
residue and cannot be zero in `Z[zeta]`. Its squared complex magnitude is

    N(a+b zeta)=a^2+ab+b^2,

an integer which is positive for `(a,b)!=(0,0)`, so it is at least 1.
This proves each retained slice contributes at least 1. Thus every tensor
Bernstein coefficient of `J(P;x)-product_i(1+x_i^2)` is nonnegative, and
the nonnegative Bernstein weights prove the lemma on the cube.

There is no invalid use of separate formal phase variables: all phases
live in one Eisenstein ring, and the product of their residues is nonzero
in the same F4 field. The norm calculation uses the 60-degree generator;
with this convention the cross term is `+ab`, as in the supplied proof.

## 4. Controls

The producer's `verify_eisenstein_bernstein.py` was inspected and replayed
in both normal and optimized modes; replay receipts are retained here.
It uses exact integers/Fractions and explicit errors, not Python assertions.
Its 762 checks, 18 phase cases, and 484 mixed Bernstein coefficient cases
are algebraic controls, not a substitute for the proof.

`verify_arithmetic_slices.py` is an independent implementation that checks
intermediate tensor-M transformations and the complete identity-slice
partition, rather than only the final coefficient inequality. It checks
all type patterns through dimension 5 on two seeded odd-coset vectors per
dimension, including all retained-slice parities, the exact number of
retained slices, positivity of discarded slices, direct-vs-sliced tensor
quadratic forms, F4 residues, integral norms, and phase-average identities.
Its normal and optimized receipts are required to agree byte for byte.

The curated manifest includes the final primary theorem, both reviewed
arithmetic documents, exact verifiers and receipts, and independent audits.
It deliberately excludes exploratory float probes and unrelated cutoff or
zeta packets. No claim in this audit changes those frozen earlier packets.
