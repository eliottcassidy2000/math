---
id: THM-3390
title: "All-iterate commuting completion norm bound"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Let c>=2, let T act on a finite-dimensional
  complex Hilbert space, and suppose an isometric embedding V and a
  contraction Q give uniformly bounded defects
  E_n=c V^*Q^{*n}V-T^{*n}, n>=1, each commuting with T.  Then ||T||<=c.
  Boundedness of the defects first makes T power bounded, which is the
  missing tail justification in the singular-vector telescope.  The
  constant is sharp for T=cJ_2 using an explicit bilateral-shift unitary
  dilation.  Q=0,T=MJ_2 gives the minimal noncommuting hostile for every
  M>c; a scalar M gives the corresponding unbounded-defect hostile.  If T
  is the raw compression V^*QV, the first commutator already forces T
  normal, so a nontrivial application needs a distinct positive completion
  layer.  The c=2 mechanism is inspired by Lorist--Schwenninger
  arXiv:2608.03841v1, a very recent cited preprint claiming Crouzeix's
  conjecture; the present theorem is self-contained and does not promote
  that external claim to established literature.
source: root-2608-crouzeix-puzzle-2026-08-14
audit: independent recurrence, singular-vector, square-completion, tail, sharp-dilation, omission-hostile, and raw-compression line audit
related:
  - THM-2222-scalar-transfer-parity-tower-and-four-checkpoint-survivor-reduction
  - THM-3315-tournament-cut-switching-centered-coronal-walk-compiler
script: 04-computation/all_iterate_commuting_completion_norm_bound_thm3390.py
output: 05-knowledge/results/all_iterate_commuting_completion_norm_bound_thm3390.out
script_sha256: 8d5b315015a6a2cfdfabd156c79a684c811b6db468cf9c3ad70c4e4ca96f7d20
output_sha256: 08f1f377a07b894cd77809624ef3971ae8265e77f2a25bd223ecb247bff3b148
semantic_sha256: ad59ba4dea8c755a57c1947ff114593e383caf8597ec185497b4c952f4969c88
hash_basis: LF-normalized bytes
---

# THM-3390 -- a bounded commuting completion of every adjoint power controls the norm

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. The theorem

Let `H` be a finite-dimensional complex Hilbert space, let `T in L(H)`,
and let `c>=2`.  Suppose there are a Hilbert space `K`, an isometry

```text
V:H -> K,
```

and a contraction `Q in L(K)` such that, for every integer `n>=1`,

```text
E_n=c V^* Q^{*n} V-T^{*n}                              (1)
```

satisfies

```text
sup_(n>=1)||E_n||<infinity,             E_n T=T E_n.    (2)
```

Then

```text
||T||<=c.                                                (3)
```

The conclusion and proof remain valid on an arbitrary Hilbert space when
`T` attains its norm.  No assertion is made here for `c<2`, for merely
asymptotically commuting defects, or from any finite prefix of `(1)--(2)`.

The mechanism has three load-bearing coordinates:

```text
all powers:       one first defect is not enough;
uniform bound:    closes the telescoping tail;
commutation:      moves each defect through T^n.         (4)
```

## 2. The proof

Use the convention that the inner product is linear in its first variable.
Put

```text
R_n=c V^*Q^{*n}V,                    E_n=R_n-T^{*n}.    (5)
```

### 2.1 The bounded defect supplies the missing tail

Let `B=sup_n||E_n||`.  Since `Q` is a contraction and `V` is an isometry,

```text
||R_n||<=c,
||T^n||=||T^{*n}||=||R_n-E_n||<=c+B.                   (6)
```

Thus `T` is power bounded.  In particular,

```text
||E_nT^n||<=B(c+B)                                      (7)
```

uniformly in `n`.  This elementary implication is essential: a bound on
`E_n` by itself would not justify discarding a term containing `T^n` unless
one first uses `(5)`.

### 2.2 The singular-vector telescope

Write `kappa=||T||`.  If `kappa<=1`, then `(3)` is immediate.  Suppose
`kappa>1`.  Finite dimensionality gives a unit vector `x` with

```text
T^*T x=kappa^2 x.                                       (8)
```

For `n>=1`, define

```text
m_n=Re <E_nT^n x,x>,
y_n=(R_(n+1)T-kappa R_n)x.                              (9)
```

Using `E_nT=TE_n`, then `(5)` and `(8)`, gives

```text
kappa m_n-m_(n+1)
 =Re <T^n(kappa E_n-E_(n+1)T)x,x>
 =(kappa^2-kappa)||T^{*n}x||^2-Re <T^n y_n,x>.         (10)
```

Completing the square with `a=kappa^2-kappa>0` yields

```text
kappa m_n-m_(n+1)
 =a||T^{*n}x-y_n/(2a)||^2-||y_n||^2/(4a)
 >=-||y_n||^2/(4(kappa^2-kappa)).                       (11)
```

Divide by `kappa` and iterate from `1` through `N`:

```text
m_1>=kappa^(-N)m_(N+1)
     -sum_(n=1)^N kappa^(-n)||y_n||^2/
        (4(kappa^2-kappa)).                              (12)
```

By `(7)`, `m_(N+1)` is uniformly bounded, so the first term tends to zero.

### 2.3 One completion vector forces the sharp bound

Set

```text
w=Q^*VT x-kappa Vx.                                     (13)
```

Equation `(5)` gives

```text
y_n=cV^*Q^{*n}w,                    ||y_n||<=c||w||.    (14)
```

Letting `N` tend to infinity in `(12)` and summing the geometric series
therefore gives

```text
m_1>=-c^2||w||^2/(4kappa(kappa-1)^2).                  (15)
```

On the other hand, `||Tx||=kappa`, contraction of `Q`, and
`R_1=E_1+T^*` give

```text
||w||^2
 <=2kappa^2-2kappa Re <V^*Q^*VT x,x>
 =2kappa^2-(2kappa/c)(m_1+kappa^2).                    (16)
```

Substitution of `(15)` into `(16)` produces

```text
[1-c/(2(kappa-1)^2)]||w||^2
 <=2kappa^2(1-kappa/c).                                 (17)
```

If `kappa>c>=2`, then the right side is negative.  The coefficient on the
left is positive because

```text
(kappa-1)^2>(c-1)^2>=c/2,
2(c-1)^2-c=(2c-1)(c-2)>=0.                              (18)
```

This is impossible.  Hence `kappa<=c`, proving `(3)`.

## 3. Sharpness: a literal unitary `c`-dilation

Let

```text
J= [0 1]
   [0 0],                         T=cJ.                 (19)
```

Then `||J||=1`, so `||T||=c`.  The bound is therefore sharp once `(1)--(2)`
are realized.

For an explicit realization, let `K=ell^2(Z)` with basis `(delta_j)`, let

```text
U delta_j=delta_(j+1),
V e_1=delta_0,                    V e_2=delta_(-1).     (20)
```

The bilateral shift `U` is unitary and `V` is an isometry.  Directly,

```text
V^*U^nV=J^n                     for every n>=1.         (21)
```

For `n=1`, the only returning basis arrow is
`delta_(-1)->delta_0`; for `n>=2`, neither embedded basis vector returns to
the two-dimensional subspace, and both sides of `(21)` are zero.  Taking
adjoints and using `J^2=0` gives

```text
cV^*U^{*n}V-(cJ)^{*n}=0          for every n>=1.        (22)
```

Thus `Q=U` is an explicit unitary contraction, every defect vanishes, and
equality holds in `(3)`.

## 4. Sharp failure boundaries

### 4.1 Commutation cannot be omitted

Fix `M>c`, take `K=H=C^2`, `V=I`, `Q=0`, and `T=MJ`.  Then

```text
E_1=-MJ^*,                         E_n=0 for n>=2.      (23)
```

The defects are uniformly bounded, but

```text
[E_1,T]=M^2(JJ^*-J^*J)=M^2 diag(1,-1) !=0,             (24)
```

while `||T||=M>c`.  This is the minimal-dimensional hostile to deleting the
commutation hypothesis.

### 4.2 Uniformity cannot be omitted

On `H=K=C`, take `V=1`, `Q=0`, and the scalar `T=M>c`.  All defects commute,
but

```text
E_n=-M^n                                                   (25)
```

is unbounded.  Thus commutation without an all-power bound does not control
the norm.

### 4.3 Raw compression is the wrong positive layer

Suppose one tries to use an ordinary compression itself:

```text
T=V^*QV.                                                 (26)
```

Then already at the first iterate,

```text
E_1=cV^*Q^*V-T^*=(c-1)T^*.                              (27)
```

For `c!=1`, the first commutation condition forces `[T^*,T]=0`; hence `T`
is normal.  Since `(26)` is a compression of a contraction, `||T||<=1`, and
the conclusion is trivial.  A genuinely nonnormal application therefore
requires a positive completion map distinct from the target homomorphism.

The unitary compression in `(20)` makes the boundary visible.  With
`T=V^*UV=J` and `c=2`, one has

```text
E_1=J^*,                         [E_1,T]=diag(-1,1).     (28)
```

This is also the minimal Fourier-window hostile: boundary exit and return
appear exactly as nonnormality of the compressed shift.

## 5. Crouzeix provenance and exact citation status

Lorist and Schwenninger's
[*A solution to Crouzeix's conjecture*](https://arxiv.org/abs/2608.03841),
arXiv:2608.03841v1, submitted 2026-08-04, states the case `c=2` of this
lemma and uses it with the numerical-range double-layer potential.  Their
positive unital map has a Stinespring representation

```text
Phi(f^n)=V^*Q^nV,                                        (29)
```

whereas the target is `T=f(A)`.  The adjoint-algebra identity makes

```text
2V^*Q^{*n}V-T^{*n}=alpha(f^n)(A),                       (30)
```

which is uniformly bounded and commutes with `T`.  The case `c=2` of the
present theorem then gives the claimed Crouzeix inequality.

The source is a **very recent v1 preprint**, not a peer-reviewed or stable
published theorem as of this audit.  Accordingly `(29)--(30)` are cited here
as provenance for the mechanism, not imported as established canon.  The
proof of Sections 1--4 is self-contained.  This theorem alone does not prove
Crouzeix's conjecture: that application still needs the double-layer
positivity, completion identity, smooth-domain reduction, and functional
calculus from the cited work.

The scalar/completely bounded distinction is also load-bearing.  The cited
preprint explicitly notes that matrix amplification loses the commutation
used in `(10)`; neither its argument nor THM-3390 proves the complete
constant-`2` statement.

## 6. Typed transfer ledger

The reusable interface is

```text
source:       a contraction dilation of a positive completion Phi;
target:       the norm of a generally nonnormal homomorphic image T;
map:          E_n=c Phi(f^n)^*-T^{*n};
preserved:    every power, uniform boundedness, commutation with T;
destroyed:    matrix-level commutation, coefficient signs, chosen observers;
sidecar:      the full commutator sequence [E_n,T];
test:         cJ sharpness, MJ noncommutation, raw-shift boundary.          (31)
```

Two repo-facing consequences are exact but limited.

1. For Fourier or factorial-moment windows, a raw multiplication-operator
   compression lands in `(27)` and is nontrivial only after adding a
   different positive layer that absorbs boundary flux.  Shared moment or
   kernel syntax does not supply that layer, so no `FC(3)` or `GMC` result
   follows.
2. The two-vertex tournament adjacency is `J`.  Its numerical range is the
   disk of radius `1/2`, and `p(z)=z` gives
   `||p(J)||=1=2 max_(W(J))|p|`; factor `2` is already sharp inside
   tournaments.  The centered matrix `2J-mathbf1 mathbf1^T` is normal, so
   the sharp nonnormality enters through the rank-one ordinary-adjacency
   recovery used by THM-3315.  A Crouzeix envelope does not recover that
   theorem's signed observer numerator, SCC order, or Hamiltonian data.

There is no direct `LRC(14)`, `JC(2)`, `FC(3)`, or tournament classification
consequence.

## 7. Exact companion

The standard-library companion uses rational `2x2` matrices and integer
bilateral-shift coordinates.  It checks the final contradiction signs at
four rational `(c,kappa)` pairs, the unitary dilation and zero defects through
power `64`, the sharp norm, both omission hostiles, and the raw-compression
commutator.  These are controls for the proof, not a finite replacement for
the all-iterate hypotheses.

Reproduce with

```text
python3 04-computation/all_iterate_commuting_completion_norm_bound_thm3390.py
python3 -O 04-computation/all_iterate_commuting_completion_norm_bound_thm3390.py
```

The two transcripts agree byte for byte.  The script has no floating-point
literal, external dependency, or assertion-dependent truth gate.

**QED.**
