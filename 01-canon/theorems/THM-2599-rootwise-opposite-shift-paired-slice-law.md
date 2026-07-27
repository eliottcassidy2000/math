---
id: THM-2599
title: "Rootwise opposite-shift target/blocker paired-slice law"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Let
  I_h=(h/13,(h+1)/13), let k and a be unequal positive integers, and put
  P_s(y)=d_L(ky+s/13)u_1(ay-s/13), where L is 1 or 2.  Every physical
  root cell I_h has a shift s for which P_s has positive measure.  For a
  guard target the thirteen paired slices cover pointwise with
  multiplicity at least one and one slice has mass at least 1/169.  For
  an ordinary target, equality of all target/blocker shift components
  would force equality of their essential vector-wall counts 2k and 2a;
  hence one slice is positive, with rational-grid floor 1/(182ak).
  On a positive common wall chamber the whole shift profile is a fixed
  nonempty proper Boolean subset of F_13; primality then forces every one
  of its twelve nontrivial Fourier colours to be nonzero.
  Choosing one common-depth base-13 cylinder inside each of the thirteen
  root chambers embeds a full one-sided 13-shift for a power of T: every
  finite ordered root word has an exact positive same-ancestry paired
  realization by disjoint digit blocks.
  Canonically, k is a thirteen-unit role and a its 13-divisible paired
  blocker, so the unequal-slope hypothesis is automatic.  Combined with
  the proved THM-2604 all-root role selection, the law gives a genuinely later
  same-root target-danger/blocker-safe opposite-shift occurrence.  It does
  not identify the physical root or shift with a THM-2334 relation residue,
  retain the left endpoint mode, complete a full-X endpoint, produce a
  completed THM-2334 target character/current, remove a scalar row, or
  prove LRC(14).
source: lrc-semantic-frontier-2026-07-28
depends_on:
  - THM-2379-anchored-guard-unit-deletion-factor-repair-current
  - THM-2604-unshifted-future-root-accessibility-and-selector-cross-mixing-boundary
related:
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2350-owner-pivot-dual-dipole-normal-form
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2537-cayley-wall-scalarization-and-positive-selector-intertwiner
  - THM-2563-paired-dipole-deep-target-corner-and-partial-bare-boundary
  - THM-2602-commutative-vertex-insertion-and-ordered-transition-curvature-no-go
script: 04-computation/lrc14_rootwise_opposite_shift_paired_slice_thm2599.py
output: 05-knowledge/results/lrc14_rootwise_opposite_shift_paired_slice_thm2599.out
script_sha256: 77b6dd9fb8f21329e3510b1511b10bb3e6b9d7a555226850681a024fcc70d014
output_sha256: 314eb5b9221d5d237e655ab28bcda85bf33bd1f2ca800af0e2acc0bf1059fd22
hash_basis: working-tree bytes (LF)
---

# THM-2599 -- every root admits a lawful opposite-shift transition slice

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The all-root target-access problem and the paired-dipole problem meet at a
smaller object than a completed endpoint current.  A later target-active
first failure asks for a translated **danger** factor.  Its actual blocker
must travel by the opposite translation and remain **safe**.  Surprisingly,
that two-factor transition slice reaches every physical root without any
further role-size threshold.

The mechanism depends on the role type:

```text
guard target:       pointwise translated-tooth surplus;

ordinary target:    equal shift totals + unequal vector-wall degree.       (1)
```

The ordinary argument is a rigidity theorem.  If every paired slice
vanished, the complete target and blocker shift vectors would agree almost
everywhere on the root cell.  Their labelled wall coverings have respective
degrees `2k` and `2a`, which is impossible for unequal slopes.  The
thirteen-adic LRC typing guarantees that inequality, but is not otherwise
used in the local theorem.  Shrinking the surviving slice to one wall
chamber then upgrades positivity to complete primitive Fourier support in
the moving shift coordinate.

## 1. Statement and sign convention

On `R/Z`, put

```text
d_L(x)=1_(||x||<L/14),             u_L=1-d_L,

I_h=(h/13,(h+1)/13),               h in F_13.              (2)
```

Let `k,a` be unequal positive integers and let `L` be one or two.  For
`s in F_13`, define the opposite-shift transition slice

```text
P_s(y)
 =d_L(ky+s/13) u_1(ay-s/13).                            (3)
```

Then, for every `h`,

```text
max_(s in F_13) mu(I_h intersection {P_s=1})>0.           (4)
```

More quantitatively,

```text
L=2  =>  sum_s P_s(y)>=1 a.e.,
          max_s mu(I_h intersection {P_s=1})>=1/169;

L=1  =>  max_s mu(I_h intersection {P_s=1})
          >=1/(182ak).                                    (5)
```

The signs in (3) are THM-2563's first-dipole convention: the `k_a` factor
has `+s/13` and its actual blocker `a` has `-s/13`.  Replacing every `s` by
`-s` gives THM-2350's displayed representative, so no orientation content
depends on that global reversal.  The factor types do matter: (3) is
target-danger/blocker-safe, not THM-2563's moving safe/safe endpoint pair.

## 2. The guard is a pointwise surplus

THM-2379's exact translated-tooth identity, valid away from the finite wall
set, is

```text
sum_(s in F_13)d_L(x+s/13)=2L-d_L(13x).                    (6)
```

The same identity holds with the translation sign reversed.  Put

```text
A_s(y)=d_2(ky+s/13),

B_s(y)=d_1(ay-s/13).                                      (7)
```

For Boolean variables, `A_s(1-B_s)>=A_s-B_s`.  Summing (7), using (6),
and recalling that both remaining indicators take values in `{0,1}` gives

```text
sum_s P_s(y)
 >=[4-d_2(13ky)]-[2-d_1(13ay)]
 =2-d_2(13ky)+d_1(13ay)
 >=1.                                                     (8)
```

Thus the thirteen guard slices cover the entire circle up to null walls.
Integrating (8) over a root cell of length `1/13` and pigeonholing among
the thirteen shifts proves

```text
max_s mu(I_h intersection {P_s=1})>=1/169.                (9)
```

No divisibility, root accessibility, or cover hypothesis entered this
argument; only positivity of the two integer slopes is needed.

## 3. Equal ordinary shift totals force component equality

Now take `L=1` and suppose, toward a contradiction, that every set in (4)
is null.  Write

```text
A_s(y)=d_1(ky+s/13),

B_s(y)=d_1(ay-s/13).                                     (10)
```

Since `P_s=A_s(1-B_s)`, the nullity assumption says

```text
A_s<=B_s a.e. on I_h                         for every s. (11)
```

Both complete shift vectors have exactly the same integrated mass.  By
(6), the substitution `z=13y-h`, and invariance of Lebesgue measure under
multiplication by a positive integer,

```text
sum_s integral_(I_h) A_s(y)dy
 =integral_(I_h)[2-d_1(13ky)]dy
 =2/13-(1/13)integral_0^1 d_1(kz)dz
 =2/13-1/91
 =1/7.                                                    (12)
```

The identical calculation for `B_s` gives

```text
sum_s integral_(I_h) B_s(y)dy=1/7.                       (13)
```

Every term in `B_s-A_s` is nonnegative by (11), while their integral sum
is zero by (12)--(13).  Therefore

```text
A_s=B_s a.e. on I_h                         for every s.  (14)
```

The contradiction is now topological/combinatorial rather than metric.

## 4. The vector-wall degree remembers the slope

The two essential walls of `A_s` occur when

```text
ky+s/13=m+epsilon/14,             epsilon in {-1,+1}.
```

Equivalently,

```text
y=(14n+13epsilon)/(182k),         n=13m-s.               (15)
```

As `s` ranges through `F_13`, `n=13m-s` ranges through every integer
exactly once with its shift label.  For each fixed `epsilon`, (15) is a
lattice of spacing `1/(13k)`.  It has exactly `k` points in `I_h`: neither
endpoint can be a wall, because endpoint equality would imply

```text
14kh=14n+13epsilon
```

or the same equation with `h+1`, impossible modulo fourteen.  The two
`epsilon` lattices do not meet, since a meeting would give

```text
14(n-n')=plus_or_minus 26,                              (16)
```

also impossible.  Each of the resulting `2k` points is a genuine Boolean
jump and belongs to a unique translated component.

For the blocker vector, the corresponding formula is

```text
y=(14n+13epsilon)/(182a),         n=13m+s,               (17)
```

so it has exactly `2a` essential labelled walls in `I_h`, with the same
endpoint and noncollision properties.

Equality almost everywhere of two finite step functions preserves their
essential interior jump sets; applying this componentwise in (14) therefore
forces equality of the total vector-wall counts:

```text
2k=2a.                                                    (18)
```

This contradicts `k!=a` and proves (4) in the ordinary case.

The proof becomes silent at `k=a`; no failure at equality is asserted.
For the canonical application, however,

```text
13 does not divide k,                 13 divides a,       (19)
```

so equality cannot occur.

The root-cell quantifier is sharp: unlike (8), the ordinary result is not
pointwise in `y`.  Take

```text
k=14,                 a=2197,                 y=733/737.   (19a)
```

Then

```text
-ky=ay=56/737 mod 1.
```

The exact target- and blocker-danger shift profiles are both the singleton
`{1}`.  Indeed the active displacement is `9/9581<1/14`, while the two
adjacent displacements are `56/737>1/14` and `746/9581>1/14`.  Hence

```text
P_s(y)=0                         for every s in F_13.       (19b)
```

All displayed inequalities are strict, so the same failure holds on an
open neighbourhood.  This point lies in `I_12`, and another open chamber of
`I_12` has a nonempty paired profile by the theorem.  Thus the wall-degree
argument genuinely finds a different point in the same root; it cannot be
replaced by a universal pointwise deck-covariance claim.

## 5. A uniform rational-grid chamber floor

All walls appearing in (3), for either role type, as well as both endpoints
of `I_h`, lie on the common lattice

```text
(1/(182ak)) Z.                                           (20)
```

For the guard, replace `13epsilon` in the target formula (15) by
`26epsilon`; the same denominator applies.  Thus multiply the target-wall
numerator by `a`, the blocker-wall numerator by `k`, and note that
`h/13=14akh/(182ak)`.  A positive set in (4) is a union of open atoms of
this finite wall arrangement.  At least one such atom has length at least
one lattice step.  This proves the second bound in (5) and also gives the
same chamber floor in the guard case.  The floor is robust, not claimed
sharp.

Strict tooth boundaries are finite and null throughout.  Assigning their
Boolean values differently changes none of (4)--(5).

## 6. One chamber carries all twelve moving-dipole colours

Refine `I_h` by every target and blocker wall occurring in (3).  On each
open chamber `J`, the complete Boolean profile

```text
S_J={s in F_13:P_s(y)=1}                                  (21)
```

is independent of `y in J`.  The proof of (4) supplies a positive chamber
for which `S_J` is nonempty.  It is also proper.  Indeed, (6) shows that at
most `2L<=4` target-danger shifts are active at any point, and `S_J` is a
subset of those shifts:

```text
1<=|S_J|<=4<13.                                           (22)
```

Let `zeta` be a primitive thirteenth root of unity.  For every nonzero
`r in F_13`,

```text
sum_(s in S_J) zeta^(rs) !=0.                             (23)
```

To prove this, put `F_S(X)=sum_(s in S_J)X^s`.  If the left side of (23)
vanished, then the minimal polynomial

```text
Phi_13(X)=1+X+...+X^12
```

would divide the nonzero integer polynomial `F_S`, whose degree is at most
twelve.  Its `0/1` coefficients would then force `F_S=Phi_13`, contrary to
(22).  Thus every one of the twelve primitive moving-shift colours survives
on the **same** positive wall chamber.

There is also a crude uniform magnitude invoice.  The twelve values in
(23) are the Galois conjugates of one nonzero cyclotomic integer, so their
product is a nonzero integer.  Every conjugate has absolute value at most
`|S_J|`.  Therefore

```text
|sum_(s in S_J)zeta^(rs)|>=|S_J|^(-11)>=4^(-11).           (24)
```

The chamber can be chosen with `mu(J)>=1/(182ak)` by Section 5.  Hence the
normalized Fourier transform of

```text
K_J(s)=mu(J)1_(s in S_J)
```

has every nonzero coefficient of magnitude at least

```text
1/[13*4^11*182ak].                                        (25)
```

This is complete spectrum in the moving translation coordinate.  It is not
yet a left-minus-right target coefficient.

## 7. Toothpick self-similarity contains a full root shift

For every `h`, choose one positive chamber `J_h` from Section 6 and one
shift `s_h in S_(J_h)`.  Every nonempty open interval contains a standard
base-thirteen cylinder.  Because there are only thirteen roots, subdivision
gives a common depth `ell` and pairwise disjoint cylinders

```text
C_h=[n_h/13^ell,(n_h+1)/13^ell) subset J_h subset I_h      (25a)
```

for all `h`.  The half-open convention only assigns the null base-expansion
boundaries.  On each `C_h`, the complete profile `S_(J_h)` is constant and

```text
P_(s_h)(y)=1.                                             (25b)
```

Let `(h_0,...,h_m)` be any ordered word in `F_13`, and let `N>=ell`.  Then

```text
E(h_0,...,h_m;N)
 =intersection_(j=0)^m T^(-jN)C_(h_j)                     (25c)
```

fixes `m+1` disjoint blocks of `ell` base-thirteen digits.  Haar product
measure therefore gives the exact specification law

```text
mu(E(h_0,...,h_m;N))=13^(-(m+1)ell)>0.                    (25d)
```

No mixing limit is used.  Each point of (25c) follows the prescribed
physical root word on one orbit and satisfies the actual opposite-shift
target-danger/blocker-safe event (25b) at every stop.

At `N=ell`, let `w_h` be the length-`ell` digit word defining `C_h`.
Subdivide once more if necessary so that no codeword is the all-`12` word.
Concatenation gives

```text
(h_0,h_1,...) -> 0.w_(h_0)w_(h_1)... base 13,             (25e)
```

an injective symbolic copy of the full one-sided thirteen-shift, conjugating
the symbol shift to `T^ell`.  For `N>ell`, fix one filler word in each gap.
This is the toothpick self-similarity hidden inside the local paired law:
one positive tooth in every first digit recursively contains every finite
ordered root itinerary.

In particular, for any `alpha!=0` and `h_0 in F_13`, the eight-stop affine
word

```text
h_j=h_0-j alpha,                    0<=j<=7,               (25f)
```

has exact realization mass `13^(-8ell)`.  Thus any proposed transition
criterion depending only on this ordered physical inverse-root word has a
genuine fibre-product witness.  This does not identify these roots with
THM-2542's chart coordinates, provide its seven-clock transport, make the
root-dependent shifts `s_h` uniform, or attach a relation-address current.

## 8. Canonical later-root and selected-source composition

In the scalar `5+3` branch, every guard/unit role is a thirteen-unit and
each actual blocker is a positive multiple of thirteen.  THM-2604 proves
that on every hypothetical scalar-cover row one may choose a maximal-`nu_7`
graft `q_*` and a distinct, pivot-eligible, target-active role `k_a` whose
danger is accessible from every physical root.  Pair `k_a` with its actual
THM-2309/2350 blocker `a`; then (19) holds automatically.

THM-2537 and THM-2604 give a positive canonical selected-head layer

```text
A^head subset I_t                                             (26)
```

with the lawful role `k_a`.  THM-2599 chooses some `s_t in F_13` for which

```text
B_(t,s_t)
 =I_t intersection {d_(L_a)(k_a y+s_t/13)=1}
      intersection {u_1(a y-s_t/13)=1}                     (27)
```

has positive measure.  Both sets are rational Boolean step sets.  The
elementary inverse-grid mixing argument used at the future-root frontier
therefore gives, for every sufficiently large `N`,

```text
integral A^head(x) 1_(B_(t,s_t))(T^N x)dx>0.               (28)
```

Every point counted in (28) has the same physical first digit `t` at the
old head and the later event, while at the later event

```text
d_(L_a)(k_a T^N x+s_t/13)=1,

u_1(a T^N x-s_t/13)=1.                                    (29)
```

Thus THM-2604 composes with THM-2599 to retain the actual blocker under the
opposite target translation on every scalar-cover row.  The shift may
depend on the root `t`; (28) does not prescribe one global shift for all
heads.

The same composition is stronger on THM-2537's selected **source** packet.
Partition its positive rational Boolean set `S^sel_tau` by the thirteen
physical root cells and choose a positive piece

```text
A^src_t=S^sel_tau intersection I_t.                       (30)
```

Source-side, this piece still retains THM-2537's terminal word, late owner,
old shallow carrier, deepest-comb avoidance, and all twelve old shallow
root colours.  Choose the positive common future chamber `J subset I_t`
from Section 6, with profile `S_J`, and put

```text
rho_N=integral A^src_t(x)1_J(T^N x)dx>0,

K_N(s)=integral A^src_t(x)1_J(T^N x)P_s(T^N x)dx.          (31)
```

For every sufficiently large `N`, inverse-grid mixing gives the first
inequality in (31), and constancy of the profile on `J` gives the exact
factorization

```text
K_N(s)=rho_N 1_(s in S_J).                                (32)
```

Consequently every one of the twelve nontrivial moving-dipole Fourier
coefficients of `K_N` is nonzero.  This puts the selected old source and a
later same-root paired target occurrence on one literal orbit, while
preserving the terminal word and owner as **source-side provenance**.  It
does not move that terminal word to the future endpoint or make the old
source covariant under the target translation.

## 9. The shift label is not the missing relation residue

Equation (29) closes a physical two-factor transition seam, not the
THM-2334 current seam.  Three coordinates must remain distinct:

```text
t:       one physical predecessor/root digit;

s_t:     one coordinate-translation (twist-state) label;

delta_a(u)-delta_a(v):
         the left-minus-right relation-address residue dual to that twist.
                                                               (33)
```

THM-2334 Section 9 proves that the physical predecessor action and the
independent base-coordinate translation action have no known intertwiner.
Moreover, a shift is a Fourier-dual variable, not the residue it probes.
Choosing one positive `s_t` in (27) supplies neither a left endpoint mode
`u` nor nonconstancy in `s`.  Section 6 now proves nonconstancy--and indeed
complete primitive spectrum--for the moving Boolean profile, but its later
complex endpoint refinement can still cancel against an unrecorded left
residue inside a relation-address orbit.

Consequently the valid implication is

```text
all-root target-active access
 + actual opposite blocker shift
 + positive rational cross-mixing
 -> later same-root target-danger/blocker-safe transition
    with all moving-shift colours.                            (34)
```

The invalid promotions are

```text
physical root t                    -/-> left relation residue;
chosen twist state s_t             -/-> left endpoint mode;
full moving-shift spectrum          -/-> completed endpoint current;
target-danger/blocker-safe pair     -/-> THM-2563 safe/safe packet;
later paired occurrence             -/-> future word/owner/full-X packet;
local nonvanishing                  -/-> scalar-row exclusion or LRC(14).
                                                               (35)
```

The exact next object is therefore a charged endpoint-transition measure
that retains **both** endpoint residues over the same physical ancestry,
then tests the difference character.  Adding more untyped root or shift
nonvanishing cannot recover the coordinate discarded in (33).

## 10. Exact companion

The dependency-free referee verifies (6) on every chamber of its complete
`1/182` wall refinement.  It checks the exact `2c` vector-wall law for
every coefficient `1<=c<=80`, every root, and both translation signs.  It
checks all `8,190` nonempty proper Boolean profiles against all twelve
primitive frequencies using exact reduction modulo `Phi_13`.  It then
integrates all thirteen paired slices exactly over `1,976` canonical typed
`(k,a,h,L)` cases, including the guard pointwise cover, both bounds in (5),
and the equal `1/7` ordinary shift totals.  Four valuation-free unequal-
slope banks independently check that the rigidity mechanism is not an
artifact of the LRC divisibility split.  Finally, a physical
`(k,a,L)=(12,13,1)` control finds one positive cylinder in every root at
common depth three, checks their paired profiles, and records the exact
eight-block affine-path mass.

Reproduce with

```bash
python3 04-computation/lrc14_rootwise_opposite_shift_paired_slice_thm2599.py
python3 -O 04-computation/lrc14_rootwise_opposite_shift_paired_slice_thm2599.py
```

Both executions must byte-match

```text
05-knowledge/results/lrc14_rootwise_opposite_shift_paired_slice_thm2599.out
```

after LF normalization.  The general rigidity proof is Sections 3--5 and
the complete-spectrum upgrade is Section 6; the bounded exact bank is an
independent hostile/control surface, not an extrapolation.  An independent
referee rederived the translated-tooth surplus, equal-integral rigidity,
labelled wall bijection and null-boundary treatment, then separately checked
the cyclotomic and disjoint-digit-block corollaries.  **QED.**
