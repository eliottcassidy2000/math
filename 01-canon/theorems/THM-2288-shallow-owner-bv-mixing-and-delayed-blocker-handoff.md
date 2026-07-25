---
id: THM-2288
title: "Shallow-owner BV mixing and delayed blocker handoff"
status: >
  PROVED + VERIFIED-EXACT, with a CITED word-overlap corollary. In every
  strict first-depth-one scalar profile, one of the two shallow exclusive
  owners has mass at least 5696989/735160140, while the guard-safe
  unit-mask-free residual outside that same blocker has mass at least
  2593/90090. If S is the sum of the nine scalar coefficients, then at
  every time k no earlier than that owner's expiration and satisfying
  13^k >= (2940640560/5696989)S, a set of mass at least
  14772292477/132461154025200 travels from the exclusive owner to a state
  where that owner, the guard absorber, and all five unit masks are
  unavailable. The global cover therefore forces a genuine
  blocker-to-other-blocker handoff on positive Haar mass. A faithful
  base-13 cylinder encoding turns absence of a handoff before n into
  Zakharov's cross-bifix-free hypothesis, but its bound is weaker than the
  native BV estimate and does not select the prescribed expiration time.
  An exact same-root-word family has arbitrarily delayed normalized-owner
  return, proving that pointwise finite root words without terminal phase
  or carry cannot make that selection. No valuation profile is excluded
  and LRC(14) remains open.
source: codex-2026-07-25-word-overlap-lrc-probe
depends_on:
  - THM-2255-valuation-separated-pair-cap-and-exclusive-owner-mass
  - THM-2273-shallow-owner-flow-and-deep-successor-gap-spread
related:
  - THM-2211-carry-regime-root-transducer-and-infinite-autonomous-index
  - THM-2261-expiration-image-surjectivity-and-one-core-carrier-no-go
  - THM-2267-static-owner-coverage-is-flag-and-transition-holonomy-is-a-cut-kernel
  - THM-2268-two-shell-private-owner-trident-and-raw-carry-cocycle-no-go
  - THM-2271-expiration-support-forces-a-weighted-owner-absorber-cut
external:
  - "Dmitrii Zakharov, An isoperimetric inequality for word overlap, arXiv:2602.20143v2."
script: 04-computation/lrc14_shallow_owner_bv_delayed_handoff_thm2288.py
output: 05-knowledge/results/lrc14_shallow_owner_bv_delayed_handoff_thm2288.out
script_sha256: f2f75d778fbea42eeb2ca3fa1084fad0b609f7357413f461bde3dc4abe3fa6c2
output_sha256: 1a3da7db85b9ef452f0670a1e5ff0c22fad0c084eb2d97e25896e59966edf5b0
hash_basis: working-tree bytes (LF)
---

# THM-2288 -- BV mixing forces a delayed blocker-only handoff

Use the scalar five-unit/three-blocker notation

```text
T(x)=13x mod 1,

D_a={x in R/Z:||ax||<1/14},
C_H={x in R/Z:||Hx||>1/7},

A_0=C_H minus union_(i=1)^5 D_(q_i).                 (1)
```

Assume the scalar cover

```text
C_H subset union_(i=1)^5 D_(q_i)
             union union_(r=1)^3 D_(c_r)             (2)
```

almost everywhere and a strict first-depth-one profile

```text
c_1=13u_1,       c_2=13^b u_2,       c_3=13^c u_3,

2<=b<c,                         5<=c<=19,             (3)
```

where `H,q_1,...,q_5,u_1,u_2,u_3` are thirteen-units, `H` is odd, and the
usual distinctness hypotheses hold. Put

```text
S=H+sum_(i=1)^5 q_i+sum_(r=1)^3 c_r.                 (4)
```

Then one shallow label `j in {1,2}` has an exclusive-owner set `E_j`
such that, for every integer `k` satisfying

```text
k>=lambda_j+1,

13^k >= (2940640560/5696989) S,                      (5)
```

the delayed return set

```text
H_(j,k)
 =E_j intersection T^(-k)(A_0 minus D_(c_j))         (6)
```

obeys

```text
measure(H_(j,k))
 >=14772292477/132461154025200
  =0.000111521695440... .                            (7)
```

Every source point in (6) is serviced by `c_j` alone. Almost every target
point is serviced by one of the other two blockers and by none of the guard
absorber, five unit masks, or `c_j`. Thus (7) is a positive-mass,
ancestry-compatible blocker-to-other-blocker handoff at every sufficiently
large time after owner expiration.

The theorem does **not** put this handoff at the prescribed expiration time
`lambda_j+1`. That distinction is the remaining sharp boundary.

## 1. Two positive marginals on one shallow label

Remove the deepest blocker and define THM-2273's two shallow exclusive
pieces

```text
E_1=(A_0 minus D_(c_3))
       intersection D_(c_1) minus D_(c_2),

E_2=(A_0 minus D_(c_3))
       intersection D_(c_2) minus D_(c_1).           (8)
```

They lie in `A_0`, avoid the deepest blocker, and have exactly one blocker
owner. THM-2273 proves

```text
measure(E_1)+measure(E_2)
 >=L_0:=5696989/367580070.                           (9)
```

Choose `j in {1,2}` with

```text
measure(E_j)
 >=e_0:=L_0/2
       =5696989/735160140.                           (10)
```

For this same label define the opposite blocker residual

```text
R_j=A_0 minus D_(c_j).                               (11)
```

This target is not merely the complement of the marked owner. It still
lies in the guard-safe, five-unit-mask-free residual. Consequently the
global cover forces one of the other two blockers at almost every point of
`R_j`.

THM-2273 gives the parity-sharp guard/deep cap

```text
measure(C_H intersection D_q)<=G_d

G_d=
  5/49+5/(49*13^d),            d odd,
  5/49+5/(294*13^d),           d even,              (12)
```

when `nu_13(q)=d`. Over all `d>=1`, its maximum is

```text
G_1=10/91.                                           (13)
```

Since both shallow depths are positive,

```text
measure(R_j)
 >=measure(A_0)-measure(C_H intersection D_(c_j))

 >=961/6930-10/91

 =eta_0:=2593/90090.                                 (14)
```

Equations (10) and (14) are the two positive marginals. Unlike the earlier
one-core image comparison, the second one is already a blocker-only target.

## 2. Exact BV transfer estimate

Let `P` be the normalized Perron operator of `T`:

```text
P f(y)=(1/13)sum_(r=0)^12 f((y+r)/13).               (15)
```

For every circle-BV function,

```text
Var(Pf)<=Var(f)/13.                                  (16)
```

Indeed, apply the triangle inequality along a cyclic partition in the
`y`-circle. For each inverse branch, the lifted partition runs through one
of the thirteen consecutive subintervals of the `x`-circle. Summing the
thirteen branch variations uses at most `Var(f)`, and the normalization in
(15) contributes `1/13`. Iteration gives

```text
Var(P^k f)<=13^(-k)Var(f).                           (17)
```

The operator preserves the integral. A mean-zero circle-BV function `g`
satisfies

```text
||g||_infinity<=Var(g),                              (18)
```

because `g(y)=integral(g(y)-g(z)) dz` at every regular representative
point. Applying (17)--(18) to `f=1_E` gives the uniform mixing estimate

```text
P^k 1_E
 >=measure(E)-13^(-k)Var(1_E)                       (19)
```

almost everywhere.

Every set in (8) is, modulo finitely many endpoints, a finite union of
intervals. Its boundary lies in the combined endpoint bank of the guard,
five unit combs, and three blocker combs. A comb with coefficient `a` has
at most `2a` boundary points, so

```text
Var(1_(E_j))<=2S.                                    (20)
```

Under (5), equations (10), (19), and (20) yield

```text
P^k 1_(E_j)
 >=e_0-2S/13^k
 >=e_0/2.                                            (21)
```

The transfer identity is

```text
measure(E intersection T^(-k)R)
 =integral_R P^k 1_E dmeasure.                       (22)
```

Use (14) and (21):

```text
measure(H_(j,k))
 >=(e_0/2)eta_0

 =14772292477/132461154025200,                       (23)
```

which proves (7).

The bound holds for every `k` satisfying (5), not merely for one selected
time. Requiring `k>=lambda_j+1` places the endpoint after the marked
owner's expiration. It does not say when along the intervening orbit the
first switch occurs.

At the target, discard the null exceptional set of (2). Partition each
remaining point by the least available label among the two other blockers.
One of those two labelled target pieces therefore receives at least half
of (23). The intrinsic cut

```text
{c_j} | {the other two blockers}                     (24)
```

has positive ancestry weight. No guard or unit absorber can pay this cut.

## 3. The faithful Zakharov cylinder map

Zakharov's theorem concerns the full product word space, not a projected
list of orbit labels. There is nevertheless an exact LRC map when one uses
**whole base-13 cylinders**.

Let

```text
Omega={0,1,...,12}.
```

For a word `w=(w_1,...,w_n) in Omega^n`, let `I_w` be its half-open
base-13 cylinder. For a finite union of intervals `U`, first change its
finitely many endpoint memberships arbitrarily, which changes none of the
measures or positive-measure conclusions below, and define its inner
language

```text
L_n(U)={w:I_w subset U},               alpha_n(U)=|L_n(U)|/13^n. (25)
```

Suppose that

```text
U intersection T^(-r)V is empty,
                         0<=r<=n-1.                  (26)
```

Then `L_n(U)` and `L_n(V)` have no directed suffix--prefix overlap. To see
this, suppose a suffix of `w in L_n(U)` of length `h` equals a prefix of
`v in L_n(V)` of length `h`. Concatenate `w` with the final `n-h` digits
of `v`. The resulting phase lies in `I_w subset U`, and after

```text
r=n-h
```

shifts it lies in `I_v subset V`, contradicting (26).

Zakharov's sharp word-overlap theorem therefore gives

```text
alpha_n(U) alpha_n(V)
 <=C_n:=(1/n)(n/(n+1))^(n+1).                       (27)
```

This is the exact source-to-target map:

```text
source object:
  base-13 phase cylinders contained in U and V;

target object:
  two length-n word families over thirteen root digits;

preserved predicate:
  a suffix--prefix overlap constructs a literal T-orbit handoff;

destroyed by pointwise projection:
  terminal phase, endpoint owner, carry, and cylinder thickness;

required sidecar:
  use whole cylinders, or retain the exact endpoint state at every cut.   (28)
```

If `U` has at most `p` boundary points, at most `p` depth-`n` cylinders
meet both `U` and its complement. Hence

```text
alpha_n(U)>=measure(U)-p/13^n.                       (29)
```

For the marked owner's expiration image

```text
B_j=T^(lambda_j+1)(E_j),                             (30)
```

THM-2255 and (10) give

```text
measure(B_j)>=(169/20)e_0
              =b_0:=5696989/87001200.               (31)
```

The image under a finite circle covering has boundary contained in the
image of the source boundary, so both `B_j` and `R_j` have at most `2S`
boundary points. Applying (27)--(29) shows:

```text
if
 (b_0-2S/13^n)_(+)(eta_0-2S/13^n)_(+)>C_n,

then B_j intersection T^(-r)R_j is nonempty
for some 0<=r<n.                                    (32)
```

The overlap cylinders in fact give positive measure. Pulling them back
through (30) produces a positive-mass handoff from `E_j` at total time
`lambda_j+1+r`.

Ignoring the boundary debit only as a numerical diagnostic,

```text
b_0 eta_0=14772292477/7837938108000,
```

and `n=195` is the first integer for which `C_n<b_0 eta_0`.
Equation (32), not that idealized number, is the rigorous finite statement.

The BV estimate in Section 2 is stronger here. It uses the one-dimensional
interval geometry directly, has a logarithmic coefficient-dependent
horizon, gives the explicit mass floor (23), and works at every sufficiently
large time. Zakharov's contribution is the independent language
interpretation and a precise test for any future finite marked-word model.

## 4. Why a finite pointwise root word cannot select expiration

The whole-cylinder condition in Section 3 is essential. There is an exact
hostile family inside every strict profile.

Fix `2<=b<c`, `c>=5`, and take

```text
H=1,
q_i=1+338i,                         1<=i<=5,

c_1=13,             c_2=13^b,      c_3=13^c.        (33)
```

For every integer `R>=c`, put

```text
y_infinity=1/2,
x_infinity=(39+y_infinity)/169=79/338,

y_R=1/2+1/(2*13^R),
x_R=(39+y_R)/169
   =x_infinity+1/(338*13^R).                         (34)
```

All inequalities are strict, and direct arithmetic gives

```text
x_infinity,x_R in C_1,
x_infinity,x_R notin union_i D_(q_i),

x_infinity,x_R in D_(c_1),
x_infinity,x_R notin D_(c_2) union D_(c_3).         (35)
```

Thus both points have the same local `c_1`-exclusive source data. Their
first two forward base-13 digits are also identical:

```text
(floor(13x),floor(13T(x)))=(3,0),

T^2(x_infinity)=y_infinity,
T^2(x_R)=y_R.                                       (36)
```

For the normalized owner core `u_1=1`, let

```text
X_k(y)=1_(D_1)(T^k y).
```

Since thirteen is odd,

```text
X_k(y_infinity)=0                         for every k,

X_k(y_R)=
  0,                                      0<=k<R,
  1,                                      k>=R.     (37)
```

Equivalently, the terminal base-13 tails are

```text
y_infinity:  666666...,
y_R:         66...67000...,
             R-1 copies of 6.                       (38)
```

The same root word `(3,0)` and the same first `h` owner bits therefore
support both permanent nonreturn and return after an arbitrarily chosen
delay `R>h`. This is a local hostile control, not a global scalar cover.

It proves the exact limitation:

```text
finite root digits + finite marked-core bits
do not determine a permanent or fixed-time return;

terminal phase/carry, or an exact whole-cylinder path state,
is the missing coordinate.                          (39)
```

This agrees with THM-2211's stronger all-mask result: exact autonomous
labelled continuation has infinite index. Enriching letters by finite
states does not automatically repair the Zakharov application, because
legal state paths form a constrained transition language rather than the
uniform full product `Omega^n`.

## 5. Scope and next consumer

The new conclusion is:

```text
all 150 strict profiles
 -> one shallow positive exclusive flow
 -> cofinally many positive-mass delayed endpoints
 -> a genuine blocker-only owner cut with exact ancestry.               (40)
```

This strengthens THM-2268's finite private-owner tour and THM-2271's
owner-to-absorber cut in one direction: the target here cannot be paid by
the guard or a unit mask. It is weaker in another: it does not occur at the
fixed expiration map selected by the profile.

A lawful THM-2267 transition graph may use the positive-measure pieces of
`H_(j,k)` as edges across the cut (24), with total weight at least (23).
What is still absent is an inequality charging that delayed switch against
the scalar cover, or a clock theorem transporting it back to the prescribed
expiration/gap-ancestry horizon.

The clean next tests are:

1. retain the first target owner and ask whether one label receives a
   uniform fraction stronger than the automatic half of (23);
2. intersect `H_(j,k)` with THM-2273's named deepest safe-gap ancestry;
3. inside a bounded global coefficient chart, minimize the first legal
   handoff time and compare it with `lambda_j+1`;
4. for any proposed word quotient, verify (28) and the full-product measure
   before invoking (27).

The fifteen repeated-first profiles `(1,1,c)` are outside THM-2273's
strict shallow-flow input and outside this theorem. No strict profile is
excluded, no uniform coefficient-free horizon is asserted, and LRC(14)
remains open.

## 6. Exact companion

The companion checks all `150` strict profiles and `450` exact delayed-tail
instances, every rational constant in (7), (10), (14), and (31), the
depth-one extremum in (12), and the ideal Zakharov threshold `194/195`.
Reproduce with

```bash
python3 04-computation/lrc14_shallow_owner_bv_delayed_handoff_thm2288.py
python3 -O 04-computation/lrc14_shallow_owner_bv_delayed_handoff_thm2288.py
```

Normal and optimized transcripts are byte-identical to the stored output.
QED.
