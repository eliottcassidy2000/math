---
id: THM-2449
title: "Switching fourth moment, signed four-cycles, and Gram energy"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For a
  symmetric zero-diagonal sign matrix A and
  Q_A(x)=sum_(i<j) a_ij x_i x_j on the Rademacher cube, the second
  moment is m=binom(n,2), while the fourth moment is
  3m^2-2m+24 times the signed C4 sum. The same C4 sum is exactly
  (||AA^T||_F^2-n(n-1)(2n-3))/8. Hence the fourth-moment relaxation
  of the switching min-max problem is precisely a Gram-energy
  minimization, with an elementary parity floor. A skew tournament
  matrix has x^T K x=0 and is not this object. The identity connects,
  but does not solve, the cited plus-minus extremum or follow from a
  faster algorithm for computing XX^T.
source: mac-mini-2026-07-26-switching-c4-gram
depends_on:
  - THM-2412-delta-exponential-and-central-newton-layer-split
related:
  - THM-2294-anchored-plucker-tournament-and-kakeya-address-bank
script: 04-computation/switching_fourth_moment_signed_c4_gram_thm2449.py
output: 05-knowledge/results/switching_fourth_moment_signed_c4_gram_thm2449.out
script_sha256: 636fa1d11a23553ac6668b2fd071de5e84f1dbc61a247e00aa11bac98c4ab515
output_sha256: 341a32d8738d1f7bad972f2913601015d6995a9956dbbb2847d59193af2c79d3
hash_basis: working-tree bytes (LF)
---

# THM-2449 -- the first nonconstant even switching moment is Gram energy

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Let `A=(a_ij)` be a symmetric `n x n` matrix with zero diagonal and
`a_ij in {+-1}` for `i!=j`.  Put

```text
m=binom(n,2),

Q_A(x)=sum_(i<j) a_ij x_i x_j=(1/2)x^T A x,

x uniform in {+-1}^n.                                           (1)
```

Multiplying row and column `i` of `A` by `x_i` switches the signed
complete graph.  Thus `Q_A(x)` is the upper-triangular edge imbalance
of the switched matrix `D_x A D_x`.

This makes the MathOverflow min-max object a signed-graph switching
problem, not a tournament quadratic form: if `K` is skew-symmetric,
then `x^T K x=0` identically.

## 1. Exact Rademacher moments

For an unoriented four-cycle `C`, let

```text
a(C)=product_(ij in E(C)) a_ij
```

and sum once over the `3 binom(n,4)` four-cycles of `K_n`:

```text
S_4(A)=sum_C a(C).                                               (2)
```

Then

```text
E Q_A(x)^2=m,                                                    (3)

E Q_A(x)^4=3m^2-2m+24 S_4(A).                                  (4)
```

For (3), two edge characters have nonzero expectation only when the
edges agree.  For (4), an ordered four-tuple of edges has even degree
at every vertex in exactly three cases:

- one edge occurs four times: `m` tuples;
- two distinct edges occur twice each: `6 binom(m,2)` tuples; or
- four distinct edges form a simple four-cycle: `4!` orderings of
  each cycle, carrying its sign product.

The first two counts sum to

```text
m+6 binom(m,2)=3m^2-2m,
```

which proves (4).  This classification retains repeated edges; treating
the fourth moment as a count of simple cycles alone is false.

## 2. The identical Gram invariant

Because `A=A^T`,

```text
||AA^T||_F^2=||A^2||_F^2=tr(A^4).                              (5)
```

The diagonal entries of `A^2` all equal `n-1`.  For `i!=j`,

```text
(A^2)_ij=sum_(k notin {i,j}) a_ik a_kj.
```

Squaring and summing gives the exact closed-walk census

```text
||AA^T||_F^2
 =n(n-1)(2n-3)+8S_4(A).                                        (6)
```

Indeed, the constant terms contribute

```text
n(n-1)^2+n(n-1)(n-2)=n(n-1)(2n-3),
```

and every signed four-cycle has four choices of ordered opposite
vertices and coefficient two.

Combining (4) and (6),

```text
E Q_A(x)^4
 =3m^2-2m
  +3(||AA^T||_F^2-n(n-1)(2n-3)).                               (7)
```

Consequently

```text
max_x |Q_A(x)|
 >=(E Q_A(x)^4)^(1/4),                                         (8)
```

and minimizing this fourth-moment lower bound over sign matrices is
exactly the same optimization as minimizing `||AA^T||_F^2`.
This is a relaxation of the min-max problem, not an equality between
the maximum and the fourth moment.

## 3. Switching and parity boundaries

For every `x`,

```text
(D_x A D_x)^2=D_x A^2 D_x.
```

Thus the Gram energy and every signed four-cycle product are switching
invariants.  They retain the first even closed-walk information which
the raw edge sum forgets.

There is also a universal integer-parity floor.  Each off-diagonal
entry of `A^2` is a sum of `n-2` signs.  Therefore

```text
||AA^T||_F^2 >=

  n(n-1)^2                         if n is even,

  n(n-1)^2+n(n-1)=n^2(n-1)         if n is odd.                  (9)
```

Equality in the even case would require all distinct rows to be
orthogonal; equality in the odd case would require every distinct-row
inner product to have magnitude one.  Equation (9) does not assert
that either boundary is attained for every order.

The triangular number `m` enters because the switching cube has one
binary sign per unordered pair.  This is the same edge count behind
`2^m` labelled tournament orientations, but a symmetric sign edge and
a skew tournament arc are different representations.

## 4. Relation to the cited problems

The MathOverflow question asks about the asymptotic minimum, over all
edge signs, of the maximum switching imbalance after normalization by
`n^(3/2)`.  Equations (7)--(9) give an exact switching-invariant moment
coordinate, but their scale alone is only a moment lower bound and does
not establish the requested limit.

The paper *XX^t Can Be Faster* concerns the arithmetic complexity of
computing a matrix times its transpose.  Substituting `X=A` makes its
output the Gram matrix in (5), so the algebraic object is genuinely the
same.  Algorithmic speedup does not classify the sign matrices
minimizing (5), and no extremal conclusion is imported from that paper.

## 5. Exact companion and hostile audit

Run

```text
python3 04-computation/switching_fourth_moment_signed_c4_gram_thm2449.py
python3 -O 04-computation/switching_fourth_moment_signed_c4_gram_thm2449.py
```

The dependency-free referee:

- exhausts all `1,098` signed complete graphs through `n=5`;
- evaluates `41,544` Rademacher switch vectors, with deterministic
  larger controls through `n=9`;
- independently computes moments, signed cycles, and Gram products;
- checks switching invariance and the parity floor; and
- keeps the skew-tournament zero-form as a hostile representation
  control.

Normal and optimized runs reproduce

```text
05-knowledge/results/switching_fourth_moment_signed_c4_gram_thm2449.out
```

byte-for-byte.  A separate read-only brute-force audit reconstructed
both identities before this companion was written and confirmed the
same switching/skew scope boundary.  QED.
