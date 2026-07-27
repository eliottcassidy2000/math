---
id: THM-2602
title: "Commutative vertex insertion and ordered-transition curvature no-go"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  Common-base Boolean multiplication, linear marginalization, and terminal
  Fourier readout form a one-way vertex calculus: before the readout every
  root operation is diagonal, all such insertions commute, and every
  chartwise selector contributes only a telescoping coboundary.  A four-atom
  table has mixed Boolean difference one while its insertion commutator is
  zero.  Identity and translation edge correspondences have identical
  source and target marginals but sevenfold monodromies I and S_(7a), so
  vertex data cannot determine curvature.  An ordered nonnegative 13 by 13
  edge bank can cancel THM-2542's class only if its product has a positive
  q to q-7a path, equivalently if K S_(7a) has positive diagonal; full
  deterministic trivialization is K=S_(-7a).  Diagonal vertex insertions
  fail identically, while seven supplied shifts S_(-a) pass and fewer than
  seven leave (7-t)a.  THM-2594's translated primitive response is an exact
  nonzero label-only seam control, not the missing physical correspondence.
  Neither THM-2592 nor THM-2594 currently supplies composable adjacent-edge
  root gluing.  No row is excluded and LRC(14) remains open.
source: wild-holotopy-2026-07-28-ordered-transition-curvature
depends_on:
  - THM-2591-theta-zero-selector-cech-coboundary-and-c91-holonomy-no-go
  - THM-2593-charged-target-section-atlas-and-minimal-c91-holonomy-trivialization
related:
  - THM-2592-fallback-rail-digit-diagonal-pullback-and-primitive-bockstein
  - THM-2594-realized-theta-slaved-contraction-at-the-r5-window
script: 04-computation/lrc14_vertex_insertion_transition_curvature_thm2602.py
output: 05-knowledge/results/lrc14_vertex_insertion_transition_curvature_thm2602.out
script_sha256: bdaea9644d575554682d816ab2ed04d8437ebd5e4b34c6827f8b019b9f0825e8
output_sha256: a487cc5ff3b071f6e94ed4e5f22718b3b86cf53ee0365fda6a8ad8e73eec164d
hash_basis: LF-normalized bytes
---

# THM-2602 -- commutative vertex insertion versus ordered transition curvature

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**

THM-2591 proves that choosing a root separately at each of seven charts is
only a Cech zero-cochain and cannot cancel the nonzero root-deck holonomy.
THM-2592 subsequently built an actual common-`x` pullback, and THM-2594 built
a nonzero mixed contraction on one common ancestry base.  This raises the
next precise question: does a common carrier plus a mixed coefficient already
contain the transition two-cell demanded by THM-2591?

No.  The missing datum is order, not another vertex statistic.  The theorem
below isolates the exact algebraic boundary and gives the cheapest decisive
test for future constructions.

## 1. The one-way vertex-insertion calculus

Let `K` be a field, let `A` be a commutative `K`-algebra of functions on a
common ancestry base, and put

```text
Q=F_13,                    V=A^Q.                         (1)
```

A root-labelled Boolean or weighted gate `g=(g_q)_(q in Q)` acts by the
physical-root diagonal operator

```text
(D_g v)_q=g_q v_q.                                        (2)
```

Define the **one-way vertex calculus** to allow:

1. finite sums and compositions of the operators (2);
2. application of common-base linear marginals `lambda:A -> K`; and
3. finite Fourier transforms in root, clock, or auxiliary labels only as
   terminal linear readouts.

The word "terminal" is load-bearing.  A Fourier output is not fed back as a
new physical-root operator before another gate.  Such a feedback would be a
new non-diagonal transition operation and lies outside this calculus.

For every pair of gates,

```text
D_f D_g=D_(fg)=D_g D_f.                                  (3)
```

Hence the complete insertion algebra is diagonal and commutative in the
physical root basis.  A simultaneous change of basis preserves (3).  Linear
marginals and terminal Fourier maps can detect correlations in the resulting
vector, but they do not turn that vector into an ordered endomorphism
`q -> q'`.

If any possibly nonlinear rule applied to these outputs chooses a chart root

```text
h_ell in F_13,              ell in F_7,                   (4)
```

then its edge contribution is still

```text
(dh)_ell=h_ell-h_(ell-1),
sum_(ell in F_7)(dh)_ell=0.                               (5)
```

This includes changing between positive THM-2586 rails chart by chart.  The
choice can be adaptive and need not be linear: once it assigns only one value
to each vertex, (5) is formal.  The complete 1,312-selector bank of THM-2591,
including the selector that always takes the last available rail, was checked
again exactly.

## 2. Mixed interaction is not transition curvature

The smallest hostile has four atoms indexed by `(i,j) in F_2^2`.  Put

```text
p=1_(i=0),             q=1_(j=0),             f=pq.       (6)
```

Then the mixed Boolean/Mobius difference is

```text
f(0,0)-f(0,1)-f(1,0)+f(1,1)=1,                           (7)
```

while

```text
[D_p,D_q]=0.                                               (8)
```

Thus a nonzero mixed ANOVA coefficient, a double centring, or a primitive
Fourier coefficient can certify genuine dependence without supplying an
ordered square commutator.  This is the exact finite hostile to reading
THM-2594's nonzero `Psi` as transition curvature.

## 3. Same vertex marginals, different monodromy

For `a in F_13^*`, let `S_a` be the `13 x 13` permutation matrix

```text
S_a(q,q')=1_(q'=q+a).                                     (9)
```

Consider the two nonnegative edge measures

```text
E^I(q,q')  =13^(-1) I(q,q'),
E^a(q,q')  =13^(-1) S_a(q,q').                           (10)
```

Both source marginals and both target marginals are the uniform measure on
`F_13`.  Consequently every statistic obtained independently from either
endpoint agrees, including all one-vertex Fourier data.  Nevertheless their
seven-edge conditional monodromies are

```text
I^7=I,                    S_a^7=S_(7a) != I.              (11)
```

This is a sharp information-theoretic hostile: even exact source and target
marginals at every edge do not determine the ordered gluing.  The joint
correspondence, and then its composition order, are indispensable.

## 4. The exact twisted-return criterion

Suppose a future construction supplies seven ordered nonnegative kernels

```text
K_ell(q,q') >= 0,          ell in F_7,                    (12)
```

on one lawfully composable common ancestry carrier.  Use row-vector
conventions and write

```text
K=K_0 K_1 ... K_6.                                         (13)
```

Because every summand in a matrix product is nonnegative,

```text
(K S_(7a))(q,q)>0

iff

there is a positive ordered path q=q_0 -> ... -> q_7=q-7a. (14)
```

Thus a necessary first test for cancelling the THM-2542 base holonomy `7a`
is

```text
diag(K S_(7a)) is not identically zero.                   (15)
```

A deterministic global trivialization has the stronger exact form

```text
K=S_(-7a),              equivalently K S_(7a)=I.          (16)
```

Every product of physical-root vertex insertions is diagonal.  Since
`7a != 0`, its left side in (15) is identically zero.  Therefore no operation
in Section 1 can pass even the necessary path test.

The boundary is sharp.  Supplying the genuine transition `S_(-a)` on every
edge gives

```text
S_(-a)^7=S_(-7a)                                             (17)
```

and passes (16).  Supplying it on only `t` of the seven ordered edges leaves

```text
S_(-ta) S_(7a)=S_((7-t)a),                                (18)
```

which is the identity exactly for `t=7`.  This recovers THM-2593's exact
partial-fill boundary in matrix form.

Equations (12)--(15) deliberately separate algebra from physical typing.
For the matrix product to represent an LRC carrier rather than a formal
coupling, a candidate must additionally provide:

```text
source:     seven edge tables on retained common ancestry data;
gluing:     the outgoing root of edge ell is the incoming root of ell+1;
order:      the seven joint labels survive before marginalization/DFT;
positivity: a path counted in (14) is one positive fibre product;
target:     its total correction is -7a.                   (19)
```

Matching one-edge marginals is insufficient by Section 3.  Independently
choosing seven positive entries is also insufficient: (19) requires one
composable span, not a post hoc coupling.

## 5. A nonzero label-only seam is available

The missing algebra is not mysterious.  Let `B(ell,theta)` be a labelled
array on `F_7 x F_13`, let `xi` and `zeta` have orders seven and thirteen,
and use the Fourier convention

```text
Bhat(beta,alpha)
 =sum_(ell,theta) B(ell,theta) xi^(-beta ell)zeta^(-alpha theta). (20)
```

Define the formal translation

```text
(T_a B)(ell,theta)=B(ell+1,theta+a).                      (21)
```

Then

```text
(T_a B)^hat(beta,alpha)
 =xi^beta zeta^(alpha a) Bhat(beta,alpha),                (22)

((T_a-I)B)^hat
 =(xi^beta zeta^(alpha a)-1)Bhat,                         (23)

((T_a^7-I)B)^hat
 =(zeta^(7 alpha a)-1)Bhat.                               (24)
```

For `a,alpha != 0` and `beta != 0`, both multipliers in (23)--(24)
are nonzero.  The companion checks all `12*12*6=864` primitive cases exactly
in `F_547`, using an element of order 91.

THM-2594 supplies a concrete common-base response bundle with nonzero
primitive modes, so (23)--(24) give a nonvacuous formal seam on its labels.
But the operator `T_a` was inserted by hand.  THM-2594 does not realize it as
a physical map between adjacent ancestry fibres.  This is the desired
positive control: once an ordered translation is supplied, curvature is
visible immediately; the contraction alone does not supply the translation.

## 6. What the current common carriers do and do not contain

THM-2592 is an actual common-`x` fibre product, but its normalized target
section `q` is a vertex label.  It has no second adjacent target root `q'`
and no ordered chart-to-chart kernel of type (12).

THM-2594 retains two root labels `u,q` in one common-base table
`N(u,q,ell,theta)`.  That table is stronger raw material than separate
marginals, but it is a **cospan**, not yet a seven-edge path object.  The
theorem supplies no lawful identification of the `q` output of one clock
with the `u` input of the next, and its factors live at linked stalk nodes
rather than one transition point.  Marginalizing and Fourier transforming
`N` cannot create the missing gluing after it has been forgotten.  Its exact
fixed-absolute-root hostile also shows that nonzero `Psi` is not by itself a
causal certificate for the affine slaving.

Therefore neither result currently defines the product (13).  This is not a
proof that their common ancestry tables cannot be refined into such kernels;
it identifies the additional datum that refinement must retain and the exact
test it must pass.

## 7. Exact companion and scope

Run

```text
python 04-computation/lrc14_vertex_insertion_transition_curvature_thm2602.py
python -O 04-computation/lrc14_vertex_insertion_transition_curvature_thm2602.py
```

The dependency-free companion checks:

1. all 1,312 admissible THM-2591 selectors and all 15,744 nonzero
   marker/selector pairs;
2. all 169 commutators of physical-root diagonal basis insertions and all
   169 commutators after one simultaneous exact Fourier conjugation;
3. the four-atom mixed-difference hostile;
4. equal endpoint marginals and different monodromies for every
   `a in F_13^*`;
5. all 864 primitive label-translation edge and seven-edge seam factors;
6. all 96 partial-fill cases, 156 diagonal twisted-return zeros, and 156
   full twisted returns.

Normal and optimized executions byte-match the stored transcript.  All
decisions are exact; no floating point or random fixture is used.

The theorem proves a stopping boundary and a design criterion.  It does not
construct a physical `K_ell`, identify the THM-2592 and THM-2594 clocks,
transport a semantic owner, recover a relation current, exclude a typed row,
or prove LRC(14).  The next positive object is now sharply specified: one
ordered common-ancestry root correspondence satisfying (19), with a nonzero
twisted return (15).

QED (candidate; independent audit pending).
