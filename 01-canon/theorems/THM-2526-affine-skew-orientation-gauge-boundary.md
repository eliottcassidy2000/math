---
id: THM-2526
title: "Affine skew no-go and live guard-sheet odd-bank recovery"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  On the thirteen-root
  augmentation module there are two sharp nonzero skew survivors, but each
  spends a different affine gauge.  D_5-D_5^* has rank six and is invariant
  under every multiplicative gauge about its fixed origin, but every nonzero
  translation moves it and its centre average is zero.  The skew part of the
  aligned THM-2521 Radon marginal is a rank-twelve cyclic tournament operator
  H_tau satisfying (I+T_tau)H_tau=T_tau-I; it is translation-invariant but
  changes sign under tau -> -tau and root reflection.  A_tau(D_5+D_5^*) is a
  rank-six skew operator, whereas A_tau(D_5-D_5^*) is symmetric.  The
  Fano-weighted product A_tau H_tau is rank-twelve, skew, and has no
  off-diagonal ties, but lies on exactly the same converse sheet as H_tau.
  No nonzero skew endomorphism or oriented half-set is invariant under the
  full AGL(1,F_13) gauge: its commutant on the augmentation module is scalar,
  and the invariant alternating square is zero.  Thus exterior and
  commutator constructions cannot create chirality after it is quotiented.
  However, the live scalar carrier has not made that quotient: THM-2198's
  physical sheet label s=-Hk selects tau_H=(-H)^(-1), up to harmless
  translation.  On every THM-2522 live depth-one event this turns the even
  translated collision bank into a nonzero, lossless, reflection-odd bank
  H_(tau_H)R_(tau_H)=13A_(tau_H)H_(tau_H)b with all twelve primitive modes.
  Root orientation is therefore already retained on the typed live carrier.
  What remains missing is a Boolean owner/source interpretation and coupling
  to the late owner, not a skew root operator.  No scalar row exclusion or
  LRC(14) proof follows.
source: codex-2026-07-27-affine-skew-orientation-boundary
depends_on:
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
  - THM-2504-endpoint-tournament-no-go-and-root-chart-holonomy
  - THM-2521-k13-drift-k14-potential-module-bridge
  - THM-2522-intrinsic-collision-depth-toothpick-descent-and-late-owner-decoupling
  - THM-2523-chi7-hamilton-energy-split-form-and-fano-boundary
  - THM-2524-translated-chi7-hamilton-polarization-inversion
related:
  - THM-2508-affine-cut-bundle-covariance-and-carry-permutation
  - THM-2525-unit-guard-collision-floor-and-word-owner-cross-cospan-collapse
script: 04-computation/lrc14_affine_skew_orientation_boundary_thm2526.py
output: 05-knowledge/results/lrc14_affine_skew_orientation_boundary_thm2526.out
script_sha256: 327678572ea3e54479d65b54e0995228c02f06b401536f09cc44a721e4493c9d
output_sha256: e48b6861d346d44bc73b0fd8e024223b6e1f5ef49197499c2b27fa004935cd5c
hash_basis: working-tree bytes (LF)
---

# THM-2526 -- the abstract skew quotient is empty, but the live sheet is oriented

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2523 proves that the multiplicative `chi_7` Hamilton operator is
symmetric and real-split, while its three-slope product is the symmetric
`chi_13` Paley operator.  That does not mean the surrounding algebra contains
no skew operators.  It contains two unusually clean ones:

```text
origin carrier:       D_5-D_5^*,
oriented-slope carrier:
    H_tau=sum_(s=1)^6(T_(-2tau s)-T_(2tau s)).              (1)
```

They fail opposite halves of the affine naturality test.  The first remembers
where zero is; the second remembers which of `tau,-tau` is forward.  This
theorem proves that the abstract failure is exact: every attempt to average
either choice into the full affine gauge vanishes, and no tensor, exterior,
or commutator construction can manufacture a gauge-invariant alternating
form.

The live application is stronger than that abstract quotient.  THM-2198
already retains the physical root-sheet coordinate

```text
s=-Hk mod 13.                                               (1a)
```

It therefore lawfully supplies the oriented slope, without choosing an
origin.  Once this inherited sidecar is used, the Fano-weighted product
`A_tau H_tau` is a full-rank skew current with no off-diagonal ties, and the
live even collision profile yields a lossless odd bank.  The algebra and the
live carrier are not missing root chirality.  They are missing a Boolean
owner/source meaning for that signed odd bank and its coupling to the late
owner.

## 1. Gauge and operator conventions

Let

```text
X=F_13,
V=Q^X,
U={p in V:sum_x p_x=0}.                                    (2)
```

Use the standard inner product and pullback operators

```text
(T_a p)_x=p_(x+a),
(M_u p)_x=p_(ux),                  a in F_13, u in F_13^*. (3)
```

Thus

```text
M_u^* T_a M_u=T_(ua).                                      (4)
```

The full predecessor-coordinate gauge is

```text
G=AGL_1(F_13)=F_13 semidirect F_13^*.                      (5)
```

Call an endomorphism **affine-intrinsic** when it commutes with every
pullback in (5).  Call a family `K_tau` **oriented-slope covariant** when

```text
M_u^*K_tau M_u=K_(u tau),
T_a^*K_tau T_a=K_tau.                                      (6)
```

The second notion retains a point of the oriented slope torsor
`F_13^*`; it is not an invariant object on the antipodal quotient
`F_13^*/{+/-1}`.

The exact symmetry ledger proved below is:

| operator | skew? | rank on `U` | translation | reflection | gauge spent |
|---|---:|---:|---|---|---|
| `A_tau`, `C_13` | no, symmetric | `12` | invariant | invariant | none, but no direction |
| `D_5-D_5^*` | yes | `6` | moved; average `0` | invariant | affine origin |
| `H_tau` | yes | `12` | invariant | global converse | oriented slope |
| `A_tau(D_5+D_5^*)` | yes | `6` | moved; average `0` | invariant | affine origin and slope class |
| `A_tau(D_5-D_5^*)` | no, symmetric | `6` | moved | invariant | affine origin |
| `A_tau H_tau` | yes | `12` | invariant | global converse | oriented slope |
| `H_(tau_H)R_(tau_H)` on the THM-2198 carrier | yes, as a bank | lossless on centred collisions | invariant | covariant converse | inherited guard-sheet label |
| full `G`-invariant skew sector | yes | **`0`** | invariant | invariant | impossible without sidecar |

Here reflection means `M_(-1)`.  In the rows involving `tau`, changing a
general multiplicative gauge transports the labelled family as in (6).
The penultimate row is intentionally **covariant**, not invariant: global
time reversal takes it to its negative, while nonvanishing and the paired
positive/negative support survive.  That is lawful because THM-2198 retains
the sheet label rather than quotienting it by converse.

## 2. The full affine skew sector is zero

The action of `G` on `X` is sharply two-transitive.  Consequently its action
on `X x X` has exactly two orbits:

```text
diagonal,                         size 13,
ordered distinct pairs,          size 13*12=156.             (7)
```

An invariant matrix on `V` therefore has one diagonal entry and one common
off-diagonal entry.  Equivalently,

```text
End_G(V)=span_Q{I,J}.                                      (8)
```

If `K in End_G(U)`, extend it by zero on the constant line.  The extension
lies in (8), so `K` is scalar on `U`.  A scalar skew-adjoint operator over
`Q` or `R` is zero.  Hence

```text
{K in End_G(U):K^*=-K}=0.                                  (9)
```

The bilinear version is the same statement.  An invariant alternating form
has one value on all ordered distinct pairs.  But the affine involution

```text
z -> x+y-z                                                  (10)
```

interchanges any prescribed `x!=y`, forcing that value to equal its
negative.  Therefore

```text
(Lambda^2 U^*)^G=0.                                        (11)
```

This also gives the sharp half-set no-go.  A translation-invariant cyclic
orientation chooses

```text
H subset F_13^*,       H disjoint -H,       H union -H=F_13^*. (12)
```

Reflection sends `H` to `-H`, its complementary half.  None of the
`2^6=64` half-systems is reflection-invariant.  Thus there is no equivariant
section of oriented slopes over `F_13^*/{+/-1}`.

Equations (9)--(11) are the general no-go.  Any polynomial tensor operation,
exterior power, contraction, or commutator made from genuinely
`G`-invariant inputs remains `G`-invariant.  If its output is alternating,
it lands in the zero space (11).  A construction from a *covariant family*
may be nonzero, but then its family label is precisely the retained sidecar.

## 3. The even Hamilton--Paley algebra contains no hidden skew part

Let `A_tau` be the THM-2523 operator

```text
A_tau=-sum_(s=1)^6 chi_7(s)
              (T_(2tau s)+T_(-2tau s)).                    (13)
```

Let `C_13` be its symmetric quadratic-character product:

```text
A_1A_2A_4=5C_13.                                          (14)
```

Every `A_tau` is an even real circulant:

```text
A_tau^*=A_tau,
A_(-tau)=A_tau.                                           (15)
```

All `A_tau` commute with one another and with `C_13`.  Reflection is also a
commuting self-adjoint involution.  Hence the real algebra generated by

```text
I,J,{A_tau},C_13,M_(-1)                                   (16)
```

consists entirely of commuting self-adjoint operators.  It has no nonzero
skew element.  In particular, taking differences, products, powers, or
commutators of the Fano-selected Hamilton energies cannot orient a root
edge:

```text
[A_tau,A_sigma]=[A_tau,C_13]=0.                            (17)
```

This is the algebraic reason the Fano selector in THM-2523 does not by
itself choose `tau` against `-tau`.

## 4. A rank-six skew survivor after choosing an origin

Put

```text
D=M_5,                   S=M_(-1).                           (18)
```

Since `5^2=-1 mod 13`,

```text
D^2=S,
D^*=D^(-1)=D^3=SD=DS.                                     (19)
```

Therefore

```text
K_origin=D-D^*=(I-S)D                                    (20)
```

is skew-adjoint.  Multiplication by five has the fixed point zero and the
three four-cycles

```text
(1,5,12,8), (2,10,11,3), (4,7,9,6).                       (21)
```

On each four-cycle, (20) is the skew adjacency of the directed cycle and has
rank two.  It vanishes at zero.  Thus

```text
rank(K_origin|U)=6.                                        (22)
```

Every multiplier `M_u` commutes with `D`, so (20) is invariant under the
entire linear gauge `F_13^*`, including reflection.  This does not make it
affine-intrinsic: it has remembered the fixed point zero.

Conjugate it to the centre `c`:

```text
D_c=T_c^*DT_c,
K_c=D_c-D_c^*.                                             (23)
```

The pullback `D_c` is dilation by five about `c`; its affine map is

```text
x -> c+5(x-c)=5x-4c.                                      (24)
```

For each ordered pair `(x,y)`, exactly one `c` satisfies `y=5x-4c`.
Consequently

```text
sum_c D_c=J=sum_c D_c^*,
sum_c K_c=0.                                               (25)
```

Every nonzero translation moves `K_origin`.  Equation (25) is the exact
failure of origin erasure: the only uniform translation average is zero.

## 5. A rank-twelve skew survivor after orienting a slope

Define the operator already implicit in THM-2521 by

```text
H_tau=sum_(s=1)^6(T_(-2tau s)-T_(2tau s)).                 (26)
```

If `R_tau` is the aligned Radon marginal on the centred potential slice,
then THM-2521 gives

```text
R_tau-(13/2)I=(1/2)H_tau.                                  (27)
```

The first row at `tau=1` is

```text
(0,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1).                       (28)
```

Thus every off-diagonal entry is `+/-1`: `H_tau` is the skew matrix of a
regular cyclic tournament, not merely a partial orientation.

There is an exact telescoping identity in the rational group algebra:

```text
(I+T_tau)H_tau=T_tau-I.                                   (29)
```

Indeed, at `tau=1` the coefficients of `H_1` alternate on the twelve
nonzero powers, so multiplication by `1+X` cancels every interior term and
leaves `X-1`; general `tau` is a relabelling.  Since `-1` is not a
thirteenth root of unity, `I+T_tau` is invertible on `V`.  The kernel of
`T_tau-I` is exactly the constants.  Hence

```text
ker H_tau=Q1,
rank(H_tau|U)=12,
det(H_tau|U)=13                                             (30)
```

The determinant follows from (29): on `U`, the numerator contributes
`Phi_13(1)=13` and the denominator contributes `Phi_13(-1)=1` (the exact
augmentation-basis calculation fixes the displayed positive sign).

Equivalently, its nontrivial Fourier multiplier is the Cayley transform

```text
(zeta^(k tau)-1)/(zeta^(k tau)+1)
 =i tan(pi k tau/13),                                      (31)
```

which never vanishes.  This identifies `H_tau` as a finite cyclic Hilbert
transform, or the Cayley transform of the discrete difference `T_tau-I`:
translation-invariant, skew, and lossless on the augmentation module.  Our
`H_tau` is **twice** THM-2521's skew part
`R_tau-(13/2)I`; the latter has one half of the multiplier in (31).

Its gauge law is equally exact:

```text
T_a^*H_tau T_a=H_tau,
M_u^*H_tau M_u=H_(u tau),
H_(-tau)=-H_tau.                                           (32)
```

In particular reflection takes the tournament to its global converse.
Also

```text
sum_(tau in F_13^*)H_tau=0,                                (33)
```

because the terms pair as `tau,-tau`.  The skew Radon current is therefore
a perfect survivor once an oriented chart is retained and a perfect zero
after that chirality is averaged away.

## 6. What `A_tau D_5` and its commutators actually produce

THM-2523 proves

```text
D^*A_tau D=-A_tau.                                         (34)
```

It follows that `A_tau` anticommutes with both `D` and `D^*`.  Split

```text
S_+=D+D^*,               S_-=D-D^*.                        (35)
```

Here `S_+` is self-adjoint and supported on the reflection-even sector,
while `S_-` is skew-adjoint and supported on the reflection-odd sector.
Both sectors have dimension six inside `U`.  Anticommutation gives

```text
A_tau S_+   skew-adjoint of rank 6,
A_tau S_-   self-adjoint of rank 6.                         (36)
```

In particular `A_tau D` is neither symmetric nor skew.  Its two parts are

```text
A_tau D-(A_tau D)^*=A_tau(D+D^*),
A_tau D+(A_tau D)^*=A_tau(D-D^*).                          (37)
```

This resolves a common commutator trap.  Although

```text
[A_tau,D]=2A_tau D,                                        (38)
```

the right side is not skew because `D` is not self-adjoint.  The genuine
skew commutator is

```text
[A_tau,D+D^*]=2A_tau(D+D^*).                               (39)
```

It still uses the fixed origin, and its translation-orbit average is zero:

```text
sum_c T_c^*[A_tau(D+D^*)]T_c
 =A_tau(2J)=0.                                             (40)
```

The same statements hold with `C_13` in place of `A_tau`, because `5` is a
quadratic nonresidue modulo `13` and hence

```text
D^*C_13D=-C_13.                                            (41)
```

Thus `D_5` exposes an origin-centred skew sector; it does not manufacture an
affine-intrinsic one.

## 7. The strongest Fano-oriented survivor keeps the same converse bit

The oriented Radon operator commutes with every `A_tau`, since both are
circulants.  Therefore

```text
F_tau=A_tau H_tau                                         (42)
```

is skew-adjoint.  Both factors are invertible on `U`, so

```text
rank(F_tau|U)=12.                                          (43)
```

At `tau=1` its first row is

```text
(0,-1,3,-5,7,-7,5,-5,7,-7,5,-3,1).                       (44)
```

Every off-diagonal entry is nonzero.  Moreover its sign is the opposite of
the corresponding entry in (28).  Hence `F_tau` is a lossless weighted
tournament current, but its underlying tournament is exactly the converse
half-set already chosen by `H_tau`.  From (15) and (32),

```text
F_(-tau)=-F_tau.                                           (45)
```

The Fano selector has changed the weights, not supplied a new orientation.

The Paley product behaves the same way.  At `tau=1`,

```text
C_13H_1 first row
 =(0,-1,3,-5,5,-3,3,-3,3,-5,5,-3,1),                     (46)
```

again with no off-diagonal zero.  Equation (14) shows that the three-slope
Fano product and the quadratic-character operator remain on this same
chiral sheet once `H_tau` is inserted.  Since `C_13` itself is even, it
cannot choose that sheet.

One may also combine the two complementary carriers:

```text
[K_origin,H_tau].                                          (47)
```

It is skew and has exact rank twelve for all twelve nonzero `tau`.  But it
spends **both** an origin and an oriented slope.  Reflection negates it, and
averaging its origin over all thirteen centres gives zero by (25).  This is
a nonzero commutator, not an escape from the affine no-go.

## 8. The live scalar carrier already supplies the oriented slope

The abstract no-go in Section 2 applies only after the root sheet has been
quotiented by the full affine gauge.  The live `5+3` carrier has not made
that quotient.  THM-2198 uses the same root disintegration

```text
x_k=(y+k)/13,                    k in F_13,                   (48)
```

and labels it in the reversed guard trivialization by

```text
s=-Hk mod 13.                                                (49)
```

Here `H` is the positive thirteen-unit guard coefficient of the scalar row.
The sign in (49) is explicit: it makes the integer guard window run in its
positive direction.  THM-2198 retains this sheet label, every labelled
mask--sheet incidence, and full integer winding.

Put

```text
u=-H mod 13,
tau_H=u^(-1)=(-H)^(-1) mod 13.                               (50)
```

This is the **chart-slope-one convention**.  To see the factor of two
without ambiguity, a chart endpoint displacement `+/-tau_H r` in the
`k`-coordinate becomes

```text
+/-u tau_H r=+/-r                                           (50a)
```

in the `s`-coordinate, while the one-sided Radon displacement
`-2tau_H r` becomes `-2r`.  If one instead insists that the latter
one-sided displacement itself be `+r`, the equivalent convention is
`tau=(2H)^(-1)`, which is chart slope `-1/2=6` in the `s`-coordinate.
The two choices differ by one fixed multiplicative relabelling.  We use
(50) because `tau` denotes the chart **half-gap** throughout THM-2521--2524.

The coordinate relabelling from a vector indexed by `k` to one indexed by
`s` is

```text
(L_up)_s=p_(u^(-1)s),                  L_u=M_(u^(-1)).        (51)
```

Therefore the baseline oriented operator in the `s`-coordinate pulls back
to

```text
L_u^*H_1L_u=H_(tau_H).                                      (52)
```

This choice needs no root origin.  Replacing the real lift of `y` by another
integer lift translates all `k`, hence all `s`; equation (32) shows that
`H_(tau_H)` is unchanged.  Global time reversal still sends the bank to its
converse, exactly as a covariant oriented current should.  Its nonvanishing
does not depend on which global converse convention is printed.

Now apply this to any one of the `165` live THM-2349 rows.  THM-2522 proves
that every live rational Boolean owner--word event has its first collision
at `L=1`, on the same multiplication-by-thirteen root fibre, while the old
shallow and deep labels remain retained.  Its terminal word factor is common
on every depth-one pair (THM-2522 equation (71a)), so it does not scramble
the inherited THM-2198 sheet coordinate.  Let `b` and `R_tau` be its
centred collision profile and THM-2524 translated Hamilton bank.  Then

```text
b_(-t)=b_t,
b!=0,
R_tau=13A_tau b.                                            (53)
```

First define the **guard-oriented Cayley bank**, and then its Fano-weighted
version:

```text
Xi_H=H_(tau_H)b,

O_H
 =H_(tau_H)R_(tau_H)
 =13A_(tau_H)Xi_H.                                         (54)
```

Let `Jf(t)=f(-t)`.  Since

```text
Jb=b,
JA_tau=A_tau J,
JH_tau=-H_tau J,                                           (55)
```

equation (54) gives

```text
JXi_H=-Xi_H,
JO_H=-O_H.                                                  (56)
```

Thus both banks vanish at zero and have opposite values on each antipodal
coordinate pair.  They are nonzero because `A_(tau_H)` and `H_(tau_H)` are
invertible on `U`.  More strongly, if `zeta=zeta_13`, their nontrivial
Fourier coefficients are

```text
Xi_hat_H(k)
 =[(zeta^(k tau_H)-1)/(zeta^(k tau_H)+1)] B_hat(k),

O_hat_H(k)
 =13 lambda_(tau_H,k) Xi_hat_H(k).                          (57)
```

THM-2524 proves `B_hat(k)>0` for all `k!=0`; THM-2523 proves
`lambda_(tau_H,k)!=0`; and (31) never vanishes.  Hence

```text
Xi_hat_H(k)!=0,
O_hat_H(k)!=0                       for every k in F_13^*.   (58)
```

The transform is lossless.  On `U`, equations (16) and (29) give

```text
b
 =H_(tau_H)^(-1)Xi_H
 =(1/13)H_(tau_H)^(-1)A_(tau_H)^(-1)O_H,

H_tau^(-1)
 =(T_tau-I)^(-1)(I+T_tau),                                 (59)
```

with the displayed factors commuting.  In particular `O_H` is nonzero and
odd, so at least one `t!=0` satisfies

```text
O_H(t)>0>O_H(-t).                                           (60)
```

This is the strongest survivor at the current frontier: every live scalar
row has a lawfully transported, lossless, all-mode odd collision bank
once the already-proved guard sheet is kept.  Root orientation is not the
remaining obstruction.

## 9. The exact remaining Boolean/owner seam

The distinction between (9) and (54) is now precise.

1. **The abstract quotient has no chirality.**  If the THM-2198 sheet label
   is discarded and the full affine gauge is averaged, every skew bank
   vanishes.  This remains the sharp guardrail against calling an unlabelled
   autocorrelation a tournament.
2. **The live typed carrier does have chirality.**  Equation (49) selects
   `tau_H`; Fano and Hamilton operators then give the lossless odd bank
   (54).  No extra arbitrary cyclic half-set is needed.
3. **The bank is signed, not Boolean.**  Each `O_H(t)` is a signed linear
   combination of translated collision intersections.  Positivity in (60)
   does not identify one measurable owner/source event of that measure or
   factor the sign into a Boolean inclusion.
4. **The root direction is not yet a semantic cospan direction.**  Increasing
   guard sheet is a physical cyclic direction, but no theorem identifies it
   with source-to-arrival, old-to-future, or blocker-to-owner order on the
   delayed word cospan.
5. **Autocorrelation still loses predecessor phase.**  Although (59)
   reconstructs `b`, it does not reconstruct the uncontracted predecessor
   vector `p`.  For every skew operator,

   ```text
   <p,Kp>=0.                                                (61)
   ```

   A scalar skew witness needs an ordered cross-pair or a selected vector
   coordinate; it cannot arise from one diagonal self-pairing.
6. **The late owner is not yet coupled to the odd coordinate.**  THM-2522
   may place a positive owner arbitrarily late while retaining the collision
   and old labels, but it does not prove that the owner's Boolean support
   lands on a positive term of (60) on the same ancestry atom.

The next decisive object is therefore narrower than “find an orientation.”
It is a **Boolean refinement of one positive guard-oriented odd coordinate**,
or an ordered source--arrival/old--future cross-cospan whose two legs live on
the retained THM-2198/2522 ancestry sheet and meet the late owner.  The
full-rank no-tie consumer (42)--(44) is already available once that semantic
factorization is supplied.  This theorem does not supply the factorization,
emit a Boolean owner loop, exclude a scalar-cover row, or prove LRC(14).

## 10. Exact dependency-free referee

Run

```bash
python3 04-computation/lrc14_affine_skew_orientation_boundary_thm2526.py
python3 -O 04-computation/lrc14_affine_skew_orientation_boundary_thm2526.py
```

Both executions reproduce

```text
05-knowledge/results/lrc14_affine_skew_orientation_boundary_thm2526.out
```

byte-for-byte.  The referee uses only exact integer/rational arithmetic.  It
checks:

- the two `AGL(1,F_13)` orbits on ordered pairs and all `64` cyclic
  half-systems;
- rank six, three directed four-cycles, all twelve translation failures,
  multiplicative invariance, and zero centre average for `D_5-D_5^*`;
- skewness, complete tournament support, rank twelve, the Cayley identity
  (29), all affine covariance laws, and zero slope average for every
  `H_tau`;
- all `144` Hamilton commutators, the `D_5` anti-isometry, and the two
  rank-six parity sectors in (36)--(37);
- the exact first rows, skewness, rank twelve, and no-tie claims for the
  Fano and Paley oriented products;
- rank twelve and both gauge failures for the mixed commutator (47); and
- vanishing diagonal pairings together with nonzero ordered cross-controls;
- all twelve THM-2198 coordinate transports (50)--(52); and
- all `8,190` nonconstant Boolean predecessor controls, their nonzero odd
  banks, `98,280` rational cyclotomic primitive-mode certificates, and one
  strict antipodal sign pair per control.

The finite referee confirms ranks and identities; the full affine no-go is
the two-transitivity proof in Section 2, while the uniform live application
uses THM-2198's proved sheet law and THM-2522/2524's live collision theorems.
An independent audit rederived the commutant and adjoint arguments, checked
the `tau_H=(-H)^(-1)` chart-half-gap convention against `s=-Hk`, verified
the Fano/Paley rows and live odd-bank inversion, and independently reproduced
the normal/optimized transcript and stored output byte-for-byte.
**QED.**
