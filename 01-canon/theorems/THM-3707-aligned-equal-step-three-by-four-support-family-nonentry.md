---
id: THM-3707
title: "Universal aligned equal-step three-by-four support family nonentry"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every squarefree polynomial
  Sigma over C of degree at least two, no Darboux pair on c^2e=Sigma(b) can
  have reduced supports {a,a+d,a+2d} and
  {a,a+d,a+2d,a+3d}, for any integer a and positive integer d.  A simple-arm
  valuation forces |j-i|d=3, leaving exactly four support placements.  One
  dies by a retained-zero singleton, two compress to forbidden homogeneous
  Darboux pairs, and the sole hard placement is THM-3702.  This closes the
  complete aligned arithmetic-progression family, not arbitrary 3x4 support
  words or JC(2).
source: jc-sparse-direct-search / 2026-08-22
audit: >
  PASS.  Root independently checked the universal simple-arm valuation,
  completeness of all six oriented solutions and four support placements,
  every bucket table, the retained-zero singleton, the sign in the
  low-double homogeneous compression, the literal THM-3702 identification,
  the high scalar singleton, and scaling from a nonzero constant to one.
  Normal and optimized replays byte-match the stored transcript, and hashes
  and documentation checks pass.
depends_on:
  - THM-3572-squarefree-danielewski-affine-modification-and-two-bracket-collapse
  - THM-3702-universal-equal-step-three-by-four-danielewski-nonentry
related:
  - THM-3696-y0-collision-ring-three-branch-conductor-and-graded-modules
  - THM-3699-y0-consecutive-four-weight-three-by-four-nonentry
  - THM-3704-y0-w002-lowest-double-two-orientation-all-scale-nonentry
script: 04-computation/jc2_aligned_equal_step_three_by_four_support_family_thm3707.py
output: 05-knowledge/results/jc2_aligned_equal_step_three_by_four_support_family_thm3707.out
script_sha256: 9d8d8e065ebb1e89e486ae5a030dd6f4d01183a7aad3b7be11927657249fee2b
output_sha256: 76527b7ae913e2478da6f615518acf82873a075393c820e261d9e516d740314a
hash_basis: LF-normalized bytes
---

# THM-3707 -- every aligned equal-step `3 x 4` word is empty

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This turns the difficult
step-three placement of THM-3702 into a complete support-family theorem.  The
first reduction is universal: it uses the simple arms of an arbitrary
squarefree exponent-two Danielewski surface, not the extra three-branch
gluing of the `y=0` collision ring.

## 1. Statement and reduced-support convention

Work over `C`.  Let `Sigma in C[b]` be squarefree of degree at least two and
put

```text
D_Sigma=C[b,c,e]/(c^2e-Sigma(b)),       wt(b,c,e)=(0,1,-2). (1)
```

On `c!=0`, every homogeneous weight-`r` element has a unique expression
`c^r f(b)`.  It is regular in `D_Sigma` precisely when

```text
Sigma^ceil(-r/2) divides f                     if r<0. (2)
```

Scalar weight-zero pieces are deleted before support is counted: they are
bracket-invisible.  Thus every retained weight-zero coefficient below is
nonconstant.

Fix `a in Z` and `d in Z_(>0)`.  Suppose all seven displayed homogeneous
pieces are nonzero and regular, with reduced supports

```text
supp(P)={a,a+d,a+2d},
supp(Q)={a,a+d,a+2d,a+3d}.                          (3)
```

Then

```text
{P,Q} notin C*.                                      (4)
```

The conclusion is unchanged after exchanging the outputs.  Indeed a
nonzero constant bracket can first be scaled to `1`, so it suffices to rule
out `{P,Q}=1` in the orientation `(3)`.

## 2. The universal simple-arm filter

For homogeneous pieces the coefficient rule is

```text
{c^r f,c^s g}=c^(r+s+1)W_(r,s)(f,g),
W_(r,s)(f,g)=s f'g-r f g'.                           (5)
```

Choose a simple zero `beta` of `Sigma`.  A scalar address has complementary
weights.  After reversing its orientation if necessary, write them as
`(-R,R-1)` with `R>=1`.  By `(2)`, its negative coefficient is

```text
f=Sigma^ceil(R/2) F.                                  (6)
```

Up to the harmless sign caused by reversing the two entries, its scalar
coefficient is

```text
J_R(f,g)=R f g'+(R-1)f'g.                             (7)
```

At `b=beta`, equation `(7)` vanishes for `R=1` because `f(beta)=0`.  It
vanishes for every `R>=3` because then both `f(beta)` and `f'(beta)` are
zero.  Therefore only `R=2` can contribute to the value of a scalar row on
a simple arm:

```text
only (-2,1) or (1,-2) can be arm-active.              (8)
```

If `{P,Q}=1`, evaluating its scalar row at `beta` now forces at least one
such address to occur between the supports in `(3)`.

## 3. Exact placement enumeration

Index the three weights of `P` by `i in {0,1,2}` and the four weights of `Q`
by `j in {0,1,2,3}`.  The two orientations in `(8)` give respectively

```text
a+id=-2, a+jd=1       ==> (j-i)d=3,
a+id= 1, a+jd=-2      ==> (i-j)d=3.                   (9)
```

Hence

```text
|j-i|d=3,                                               (10)
```

so `d` is `1` or `3`.  The six oriented solutions of `(9)` collapse to the
following four support placements:

```text
d  a    supp(P)          supp(Q)             arm placements (i,j)
1 -2    {-2,-1,0}        {-2,-1,0,1}         (-2,1):(0,3)
3 -2    {-2,1,4}         {-2,1,4,7}          (-2,1):(0,1); (1,-2):(1,0)
3 -5    {-5,-2,1}        {-5,-2,1,4}         (-2,1):(1,2); (1,-2):(2,1)
3 -8    {-8,-5,-2}       {-8,-5,-2,1}        (-2,1):(2,3).             (11)
```

There are no suppressed reversals in `(11)`: for `d=3`, the first two
placements contain both arm orientations in one support word, while the
last has only the displayed top address.  The companion enumerates `(9)`
before grouping by `(d,a)`, precisely to guard this distinction.

Every addition grid in `(11)` has fibre sizes

```text
1,2,3,3,2,1.                                           (12)
```

We now close its four possible placements.

## 4. The consecutive placement has a retained-zero singleton

For `(d,a)=(1,-2)`, the output-weight-two bucket is the singleton `(0,1)`.
Write these pieces as `p_0(b)` and `c q_1(b)`.  Its required vanishing is

```text
W_(0,1)(p_0,q_1)=p_0'q_1=0.                           (13)
```

The weight-one piece is active, so `q_1!=0`.  Since `C[b]` is a domain,
`p_0'=0`, making `p_0` a scalar.  That contradicts the reduced-support
convention.  This also gives a universal one-row proof of the aligned
consecutive case contained in THM-3699.

## 5. The low-double placement compresses to one homogeneous bracket

For `(d,a)=(3,-2)`, write the low and scalar pieces as

```text
P_-2=c^-2 f,       P_1=c p,
Q_-2=c^-2 g,       Q_1=c q.                           (14)
```

The lowest bucket is the singleton `(-2,-2)`.  Its vanishing says

```text
W_(-2,-2)(f,g)=2(fg'-f'g)=0.                          (15)
```

Thus `(g/f)'=0` in `C(b)`, so `g=lambda f` for some `lambda in C*`.  The
entire scalar bucket contains exactly the two orientations of the arm pair.
Using `(5)` and `(15)`, it compresses without division or a support shear:

```text
W_(-2,1)(f,q)+W_(1,-2)(p,lambda f)
              =W_(-2,1)(f,q-lambda p)=1.             (16)
```

If `q-lambda p=0`, equation `(16)` is already impossible.  Otherwise `(16)`
is a homogeneous Darboux pair of weights `(-2,1)` in `D_Sigma`, forbidden
by the universal homogeneous nonentry theorem THM-3572, Section 5.1.

## 6. The middle placement is exactly THM-3702

For `(d,a)=(3,-5)`, equation `(11)` is literally

```text
supp(P)={-5,-2,1},             supp(Q)={-5,-2,1,4}.    (17)
```

THM-3702 closes this all-degree word on every squarefree exponent-two
Danielewski surface.  This is the only placement in the aligned family whose
scalar fibre has three addresses; its middle-row root peel is not being
replaced by a finite-span computation here.

## 7. The high placement has a scalar singleton

For `(d,a)=(3,-8)`, the only scalar address is the top pair

```text
(P_-2,Q_1).                                            (18)
```

Consequently the bracket of these two homogeneous pieces alone must equal
`1`, again contradicting THM-3572, Section 5.1.

Sections 2--7 exhaust `(9)` and prove `(4)`.

## 8. Scope and reproduction

This theorem closes the complete aligned family in which the three-weight
support is the initial three terms of the four-weight arithmetic progression.
It does not close the other three deletions from a general four-term
progression, a non-arithmetic `3 x 4` word, or `JC(2)`.  Its inheritance is
structural: the universal arm filter turns infinitely many `(a,d)` into four
support placements, and only one of those needs the hard THM-3702 anatomy.

Run

```bash
python3 -B 04-computation/jc2_aligned_equal_step_three_by_four_support_family_thm3707.py
python3 -B -O 04-computation/jc2_aligned_equal_step_three_by_four_support_family_thm3707.py
```

Both modes must agree byte for byte with the frozen transcript.  The exact
companion checks the complete oriented placement enumeration, all bucket
tables, the singleton identities, and the low-double compression.  The
arbitrary-degree simple-arm and homogeneous-nonentry arguments are the proofs
above.  **QED.**
