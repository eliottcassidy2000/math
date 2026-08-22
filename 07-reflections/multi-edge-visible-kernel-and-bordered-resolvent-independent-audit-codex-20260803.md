# Independent audit of multi-edge kernel visibility and the bordered compiler

**Status:** **VERIFIED / FINITE-EXACT audit**, with no theorem ID.  The audit
source is
[the independent script](../04-computation/gmc_multi_edge_visible_kernel_bordered_resolvent_independent_audit_20260803.py),
with
[stored output](../05-knowledge/results/gmc_multi_edge_visible_kernel_bordered_resolvent_independent_audit_20260803.out).
It audits the
[multi-edge primary](../04-computation/gmc_multi_edge_visible_kernel_helly_bordered_resolvent_scout_20260803.py),
its stored output, and the accompanying synthesis note by pinned hashes.

The audit concerns only finite static symmetric relation-walk graphs.  It does
not establish a tournament, chronology, response-composition, `GMC`, `FC`, or
`LRC` result.  It also makes no asymptotic or bit-complexity claim.

## Independence boundary

The audited primary is neither imported nor executed.  The audit executes only
the older pinned THM-3287 relation source and reconstructs every decorated
adjacency directly from its row-edge relations.  The older half-transfer
source is hash-pinned as an inheritance target but is not used as executable
mathematical logic.

The two main calculations use different exact routes from the primary:

1. The primary obtains the all-ones kernel projection from a nullspace basis
   and its Gram matrix.  The audit instead computes the cyclic minimal
   polynomial of `T=A^2` on the all-ones vector.
2. The primary expands the local formal inverse `(C-xK(x))^{-1}`.  The audit
   instead solves a local Dyson recursion for the updated boundary vector.

Thus matching results are not just a replay of the same helpers.

## Cyclic certificate for zero visibility

Let `q` be the first exact Krylov relation for `T=A^2` on `1`.  If `q(0)` is
nonzero, the relation explicitly writes `1` as an element of `im(T)`, so it is
orthogonal to `ker(A)`.  If

```text
q(t)=t r(t),
```

then symmetry makes `T` diagonalizable and

```text
P_ker(A) 1 = r(T)1/r(0).
```

The audit verifies both `A P_ker(A)1=0` and the exact projection identity

```text
1^T P_ker(A)1 = ||P_ker(A)1||^2
```

for every visible case.  Kernel dimensions are computed separately as exact
ranks over `QQ`.

This route reproduces all `2^11=2048` lattice records byte for byte.  The
serialized census has `37,894` bytes and SHA-256

```text
53f740b8fd3e8b1baa88a460e210662d51b6aa69ff8a157fcdfc9fd5d8510ad0.
```

It independently confirms:

- all ten visible pairs contain the special edge `(7,18)`, so no new minimal
  circuit appears at arity two;
- there are exactly six inclusion-minimal visible sets, with arity profile
  `(one: 1, two: 0, three: 5)`;
- the six primitive labelled supports and masses agree exactly with the
  primary: one mass `3/5`, four masses `12/35`, and one mass `3/7`.

The support check is stronger than a scalar-mass comparison: it checks every
nonzero coordinate and sign of each primitive direction.

## Persistence versus fragility

For the singleton `{(7,18)}`, the independently recovered primitive direction
has

```text
1^T v=3,  v^T v=15.
```

Every one of the ten remaining row-edge adjacency deltas annihilates `v`.
The same vector therefore stays in the kernel on the entire upper interval.
The exhaustive `1024`-superset census has minimum zero mass exactly `3/5`,
confirming the algebraic lower-bound certificate.

For the first arity-three circuit

```text
{(2,10),(13,18),(17,22)},
```

only `(10,16)` and `(16,21)` annihilate its primitive direction.  The eight
fourth-edge masses independently reproduce

```text
(3,22):12/59   (7,18):14/15   (10,16):12/35   (10,22):0
(11,13):12/37  (13,22):0      (16,21):12/35   (19,22):0.
```

This confirms the primary's important distinction: visible zero mass is not
monotone under edge addition.  The special singleton supplies a persistent
direction, while the triple direction can be preserved, rotated, replaced, or
destroyed by the next edge.

## Independent bordered-resolvent check

In the fixed 26-node universe, write

```text
A=B+U C U^T,
```

where `B` is the selected-tree adjacency padded by isolates and `C` is the
block-swap matrix for the newly added undirected decorated edges.  Put

```text
s(x)=U^T(I-xB)^{-1}w,
K(x)=U^T(I-xB)^{-1}U,
h(x)=U^T(I-xA)^{-1}w.
```

The audit uses the Dyson equations

```text
h=s+xKCh,
G_A-w^T(I-xB)^{-1}w = x s^T C h.
```

Coefficientwise, this gives

```text
h_n=s_n+sum_(i+j=n-1) K_i C h_j,
correction_n=sum_(i+j=n-1) s_i^T C h_j.
```

This is algebraically equivalent to the primary's bordered inverse formula but
is implemented through a different recursion.  After restoring the explicit
degree-zero contribution from newly active observer vertices, the first `24`
coefficients agree with direct powers of `A` for all `55` two-edge updates.
The maximum local dimension is ten, and the update-shape census agrees exactly.
The `55` direct 24-term sequences serialize to `11,560` bytes with SHA-256

```text
d969c893735eca95eb16fec7c574bfdb9e24a68160156d7a8f2eaf12a3050b8d.
```

## Verdict and remaining boundary

The primary's finite claims, mechanism description, and scope warning are
verified.  In particular, the data support the procedural lesson that joint
edge supports can create a visible object absent from every singleton and
pair, and that a useful sequence compiler should retain active births and a
local border before scalar elimination.

What remains deliberately outside the audit is just as important: the
decorated adjacency is symmetric and static; no intrinsic binary orientation
or chronological transition law has been supplied.  The finite local
dimension bound is for these 55 updates, not an asymptotic theorem.  Any use in
another frontier must therefore transfer the *procedure*—minimal-support
closure and sidecar retention—without identifying the mathematical objects.
