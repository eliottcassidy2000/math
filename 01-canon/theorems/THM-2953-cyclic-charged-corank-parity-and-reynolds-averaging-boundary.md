---
id: THM-2953
title: "Cyclic charged corank parity and Reynolds-averaging boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the
  real augmentation representation of an odd
  prime cycle, every invariant kernel has even dimension.  Therefore a
  nonzero codimension-one minor forces injectivity; at p=7 one nonzero
  charged 5-by-5 minor suffices.  Reynolds averaging of T^*T has kernel
  equal to the intersection of the seven rotated kernels, so its
  fifth-compound positivity proves only joint injectivity of the
  rotated views.  A rank-one difference functional has an invertible
  averaged cycle Laplacian on the augmentation space, providing the
  sharp owner/carrier-loss hostile.  No THM-2941 or LRC(14) closure is
  claimed.
source: root-lrc-compound-transfer-2026-07-29
audit: >
  An independent hostile audit rederived the stable-kernel parity,
  codimension-one minor gate, Reynolds kernel intersection, and every
  rank-one hostile constant.  It replayed normal, optimized, and
  stored output byte-for-byte and matched both LF hashes.  The final
  wording repairs the distinction between THM-2608's equivariant
  matching rule and a genuinely equivariant physical response.
depends_on: []
related:
  - THM-2947-conjugate-pair-corank-parity-and-one-minor-resultant-gate
  - THM-2951-fifth-compound-reconstruction-and-v-four-phase-scalarization-boundary
  - THM-2608-alternative-rail-clock-collapse-and-missing-transition-index
script: 04-computation/lrc14_cyclic_charged_corank_reynolds_boundary_thm2953.py
output: 05-knowledge/results/lrc14_cyclic_charged_corank_reynolds_boundary_thm2953.out
script_sha256: a4f0030c8fb4d4ce3cf12c9bb1e92cde609af2f47a269c2fd0fcbb2081881a20
output_sha256: c422cf07c99c2313416fd1b4c3a08ce5c60075ec1af989b18f3dfd65930dc234
hash_basis: LF-normalized bytes
---

# THM-2953 -- cyclic charged corank parity and Reynolds-averaging boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Prime-cyclic charged parity

Let `p` be an odd prime, let `S` be the cyclic shift on `R^p`, and let

```text
V={x in R^p:sum_j x_j=0}.                              (1)
```

Thus `V` is the real augmentation representation of `C_p` and has
dimension `p-1`.  Its complexification is the direct sum of the
nontrivial character lines.  Complex conjugation pairs characters
`k` and `-k`, so over the reals

```text
V=direct_sum_(k=1)^((p-1)/2) V_{+-k},    dim_R V_{+-k}=2. (2)
```

The summands in `(2)` are pairwise nonisomorphic irreducible real
`C_p`-modules.  Hence every `S`-stable real subspace of `V` is a direct
sum of some of them and has even dimension.

Let

```text
T:V->U                                                   (3)
```

be any real linear response into a finite-dimensional real inner
product space.  If `ker(T)` is `S`-stable, then

```text
nullity(T) is even,
rank(T) in {0,2,...,p-1}.                               (4)
```

In particular,

```text
one nonzero (p-2)-minor of T  ==>  T is injective.       (5)
```

For `p=7`, the charged space has dimension six, and `(5)` says that
one nonzero real `5`-minor forces rank six.  Equivalently, the sum of
squares of all `5`-minors is positive exactly when `T` is injective,
provided the stable-kernel hypothesis is retained.

The conclusion is unchanged by a lawful cyclic twist.  If `ker(T)` is
stable, then

```text
ker(T S^a)=S^-a ker(T)=ker(T).                           (6)
```

Left multiplication by an invertible shift on the codomain likewise
does not change the kernel.

## 2. What Reynolds averaging actually proves

No stability hypothesis is needed to form the positive semidefinite
Reynolds average

```text
Abar=(1/p) sum_(j=0)^(p-1) S^j T^*T S^-j : V->V.        (7)
```

It commutes with `S`.  More importantly, for every `v in V`,

```text
<Abar v,v>
 =(1/p) sum_j ||T S^-j v||^2.                           (8)
```

Consequently

```text
ker(Abar)=intersection_j S^j ker(T).                    (9)
```

The kernel in `(9)` is `S`-stable, so the even-nullity and
codimension-one-minor gates apply to `Abar`.  But the conclusion is
only

```text
Abar injective
 iff the seven rotated views of T are jointly injective. (10)
```

It does **not** imply that `T` itself has rank at least `p-2`, has
stable kernel, or retains a named owner/carrier.

If `ker(T)` is already stable, then every term in `(9)` is the same
and Reynolds averaging is faithful at the level of injectivity.
Without that hypothesis, `(9)` is the exact information loss.

## 3. Sharp rank-one hostile

Take `p=7`, put

```text
v=e_0-e_1,
T(x)=<v,x>.                                             (11)
```

Then `T:V->R` has rank one and its kernel is not shift-stable.  Since
`T^*T=vv^*`, formula `(7)` gives the normalized cycle Laplacian

```text
Abar=(2I-S-S^-1)/7.                                    (12)
```

On the six-dimensional augmentation space its eigenvalues are

```text
(2-zeta_7^k-zeta_7^-k)/7,       1<=k<=6,               (13)
```

all positive.  Their product is

```text
det(Abar|V)
 =|product_(k=1)^6(1-zeta_7^k)|^2/7^6
 =49/7^6
 =1/2401.                                              (14)
```

In the integral augmentation basis `b_i=e_i-e_6`,
`0<=i<=5`, the exact `6`-by-`6` matrix has `30` nonzero
`5`-minors among its `36` cofactors.  One diagonal cofactor is

```text
3/2401.                                                (15)
```

Thus the averaged fifth compound is emphatically nonzero although
the original response has rank one.  The example is sharp for the
logical distinction in `(9)--(10)`: positive charged energy after
averaging can be assembled from seven different one-dimensional
views.

The same construction works for every odd prime.  The nonzero
Laplacian eigenvalue product gives

```text
det(Abar|V)=p^(3-p),                                   (16)
```

while `rank(T)=1`.

## 4. LRC boundary

At a balanced seven-cell wall, the centered cell indicators have Gram
matrix proportional to

```text
7I-J.                                                  (17)
```

Its restriction to the label augmentation space is already invertible.
The fifth-compound gate therefore repackages the regular-simplex
geometry; it does not obstruct realization by physical combs.

To use `(5)` on a separate response, its **actual kernel** must be
proved `C_7`-stable.  An endpoint/link matrix, a nonstationary
seven-slot Gram matrix, or a selected-owner current does not acquire
that property merely by averaging.  Applying `(7)` would prove only
joint injectivity of rotated observations and would erase the
owner/carrier distinction exposed by `(11)--(15)`.

This is also compatible with THM-2608, but the typing is important:
that theorem supplies an equivariant **matching rule**, while its
physical clock-block weights vary.  It does not supply a high-rank
equivariant physical response.  Its zero sevenfold Boolean product
still shows why equivariance of the indexing rule is not pointwise
carrier survival.

Accordingly this theorem supplies a cheaper exact rank test **when**
equivariance or stable kernel is independently present.  It proves no
critical-wall exclusion, no THM-2941 statement, and no decrement of
the LRC(14) ledger.

## 5. Exact companion

Run

```text
python 04-computation/lrc14_cyclic_charged_corank_reynolds_boundary_thm2953.py
python -O 04-computation/lrc14_cyclic_charged_corank_reynolds_boundary_thm2953.py
```

The companion uses exact rational arithmetic and explicit exception
gates.  It verifies:

1. the augmentation-basis matrix in `(12)`;
2. shift equivariance and rank six on `V`;
3. determinant `1/2401`;
4. the `30/36` nonzero fifth-minor census and `(15)`;
5. the fifth-minor squared-energy control;
6. rank one and nonstable kernel for `(11)`; and
7. the general determinant formula `(16)` for `p=3,5,7,11`.

The stable-kernel quantifier, the intersection identity `(9)`, the
augmentation-basis normalization, and the owner/carrier scope were
independently hostile-audited.
