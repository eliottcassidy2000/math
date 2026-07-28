---
id: THM-2820
title: "Boolean rigidity, graded translation ranks, and carrier-gauge Hasse boundary"
status: >
  RESERVED / PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT / AWAITING
  INDEPENDENT HOSTILE AUDIT.  Finite pointwise Boolean carrier algebras
  have zero relative derivations and rigid idempotents.  In F_p[C_p] the
  orbit norm is (g-1)^(p-1), fixed by every generator change although the
  cotangent line scales, so norm/Rees data alone select no nonzero first
  jet.  Positively, any rooted group-ring vector of nonzero augmentation
  has normalized Hasse coordinate b=J1/J0 with b(g^aF)=b(F)+a.
  The first normal orbit of a two-axis vector F is the quotient of the
  linear augmentation space by the degree-one annihilator of its leading
  graded piece; the THM-2806 B/P/Q/H/raw-face ranks are 0/1/1/2/2.
  The raw mixed face has base b=(0,0) and all 169 translates recover their
  labels, and its positive integral inverse/push profile is
  2171N+delta_(q1+q2), but it is supported on the 144 joint-absent cells.
  After the THM-2763 carrier gauge, the first-Hasse endpoint sees only a
  nonzero effective transverse vector.  A lawful physical transverse
  translation and common
  non-idempotent raw face or broader two-object sidecar remain missing.
source: root/boolean-norm-cotangent-boundary-2026-07-28
depends_on:
  - THM-2201-cyclic-root-fibre-hasse-jet-transition-carrier
  - THM-2763-carrier-equivariant-endpoint-address-extension-and-gauge-obstruction
  - THM-2806-literal-fixed-sheet-central-allocation-scalar-law-and-endpoint-translation-no-go
  - THM-2813-affine-lift-transvection-and-projective-horn-decoder
related:
  - THM-2771-joint-c7-c13-right-wing-mixed-spectrum-and-commuting-square-no-go
  - THM-2779-bockstein-symplectic-decoder-frame-torsor-and-heisenberg-root-degree-gate
  - THM-2580-hasse-bockstein-carry-tower-and-salem-local-unit-boundary
  - THM-2611-principal-c13-bibundle-lift-torsor-and-holonomy-section-obstruction
  - THM-2658-balanced-lift-helly-circular-arc-gain-nerve-and-wrap-boundary
  - THM-2806-literal-fixed-sheet-central-allocation-scalar-law-and-endpoint-translation-no-go
  - THM-2813-affine-lift-transvection-and-projective-horn-decoder
  - THM-2814-projective-allocation-square-holonomy-and-idempotent-provenance-no-go
  - THM-2792-cyclic-unit-intertwiner-and-positive-naturality-boundary
script: 04-computation/lrc14_boolean_norm_cotangent_boundary_thm2820.py
output: 05-knowledge/results/lrc14_boolean_norm_cotangent_boundary_thm2820.out
script_sha256: e6adbc28ff27712891e48c6a40dd4123c17c0d515d311b69f91ab005f432314c
output_sha256: 2b8886817b1464e44be036dcec067a4df58fffb3dc9487c39ad33d43c7eb55ec
hash_basis: LF-normalized bytes
---

# THM-2820 -- the Boolean norm is top-degree, while a rooted face has an affine jet

**RESERVED / PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT / AWAITING
INDEPENDENT HOSTILE AUDIT.**

THM-2806 leaves two tempting but incompatible pictures.  Pointwise, its
carrier masks are rigid Boolean idempotents and its sole fourfold atom is
flat.  After integral Fourier marginalization, its mixed coordinate is the
unit

```text
D3/v11=144=1 mod13.                                      (1)
```

THM-2813 independently produces the unit normal displacement of one affine
lift.  Equality of the two scalars does not construct a map.

This theorem locates the exact algebra between them.  The full modular orbit
norm survives as the **top** augmentation power, not a first jet.  But the
rooted raw mixed face has nonzero augmentation, so THM-2201's first Hasse
coordinate turns any externally supplied lawful translation into an exact
affine label.  That face is the 144-cell **joint-absent aggregate**, not the
fourfold common atom.  Moreover, THM-2763's address/carrier gauge can absorb
an apparent translation.  The remaining problem is therefore a physical
co-supported, gauge-transverse equivariant family, not another scalar
decoder.

## 1. Pointwise Boolean carriers are formally rigid

Let `R` be a commutative ring, let `X` be finite, and put

```text
A=Map(X,R)=product_(x in X) R                            (2)
```

with pointwise multiplication.  Let `M` be a symmetric `A`-module and

```text
d:A -> M                                                 (3)
```

an `R`-linear derivation.  For each coordinate idempotent `e_x`,

```text
d(e_x)=d(e_x^2)=2e_x d(e_x),

(2e_x-1)d(e_x)=0.                                       (4)
```

But

```text
(2e_x-1)^2=1,                                           (5)
```

so `(4)` forces `d(e_x)=0`.  The coordinate idempotents and constants span
`A` over `R`, hence

```text
Der_R(A,M)=0.                                           (6)
```

The same calculation gives first-order rigidity directly.  In the
commutative square-zero extension `A+epsilon M`, an idempotent lift

```text
e_epsilon=e+epsilon u,              epsilon^2=0          (7)
```

must satisfy

```text
(2e-1)u=0,
```

and therefore `u=0`.  Thus a Boolean selector has no internal infinitesimal
amplitude.

The tempting interpolation

```text
e+epsilon(1-e)                                           (8)
```

exhibits the sharp failure:

```text
(e+epsilon(1-e))^2-(e+epsilon(1-e))
 =-epsilon(1-e).                                        (9)
```

For a delta selector on thirteen atoms, `(9)` is nonzero on exactly twelve.
One can introduce such an amplitude as extra data, but it is not a
deformation through physical idempotent masks.

Discrete atom permutations are not contradicted by `(6)`: they are
automorphisms, not first-order derivations.  THM-2806's moving-sheet orbit is
exactly such a separate construction.

## 2. The modular convolution norm is the top Hasse class

Now let `p` be prime, `k=F_p`, and `C_p=<g>`.  In the convolution algebra

```text
k[C_p]=k[g]/(g^p-1),             eta=g-1,               (10)
```

the freshman's dream gives

```text
k[C_p]=k[eta]/(eta^p).                                  (11)
```

The orbit norm is

```text
N=1+g+...+g^(p-1)
  =(g^p-1)/(g-1)
  =eta^(p-1).                                           (12)
```

So the full orbit is not zero in the modular group ring.  It is the socle,
or top augmentation class.

Let `I=(eta)`.  The cotangent line is

```text
I/I^2=k eta.                                            (13)
```

Changing the cyclic generator by `g->g^u`,
`u in F_p^*`, gives

```text
(1+eta)^u-1 =u eta mod I^2.                             (14)
```

Hence the generator-change group scales `(13)` by `u`.  But it fixes `N`,
both because it permutes the group elements and because

```text
(u eta)^(p-1)=eta^(p-1).                                (15)
```

For `p>2`, the common fixed subspace of `I/I^2` is zero.  Therefore:

> **Norm-to-jet no-go.**  No construction natural under all generator
> changes can take the norm class `N` alone to a nonzero vector of `I/I^2`.

The norm canonically exhibits the filtered line and its top power.  It does
not choose an oriented nonzero first jet.  At `p=2` the automorphism
obstruction is vacuous, while the pointwise idempotent rigidity `(6)--(9)`
still holds; the LRC application is at `p=13`.

The distinction of multiplications is load-bearing.  The coefficient-vector
identification

```text
Map(C_p,k) = k[C_p]                         as k-modules  (16)
```

does not identify pointwise multiplication with convolution.  In
characteristic prime to `p`, a split Fourier transform relates the two
algebras.  In characteristic `p`, the function algebra remains a reduced
product while `(11)` is nonreduced.  No modular Fourier isomorphism is being
claimed.

## 3. The tensor allocation square has only degree zero and top degrees

Apply `(16)` in two rooted axes.  THM-2806's raw Boolean arrays

```text
B=1,
P=delta_0 tensor 1,
Q=1 tensor delta_0,
H=delta_(0,0)                                           (17)
```

become, in `k[C_p x C_p]`,

```text
Phi(B)=N_1N_2,
Phi(P)=1*N_2,
Phi(Q)=N_1*1,
Phi(H)=1.                                               (18)
```

Their exact augmentation bidegrees are

```text
(2p-2,p-1,p-1,0).                                      (19)
```

At `p=13`, this is

```text
(24,12,12,0).                                          (20)
```

The raw Boolean Mobius face is

```text
Omega=B-P-Q+H=(1-delta_0) tensor (1-delta_0).           (21)
```

In the modular convolution coordinates,

```text
Phi(Omega)
 =(eta^12-1)(theta^12-1)
 =1-eta^12-theta^12+eta^12 theta^12.                   (22)
```

Thus its Hasse support is exactly

```text
(0,0),(12,0),(0,12),(12,12),                           (23)
```

with no intrinsic `(1,0)` or `(0,1)` component.  Retaining the full modular
group ring repairs the central-scalar collapse, but it retains top norm
classes, not a first jet.

The pointwise support must be kept separate from the group-ring coordinates.
Equation `(21)` equals one precisely when both rooted coordinates are
nonzero, hence on `12^2=144` joint-absent cells, and it vanishes at the sole
fourfold common atom `(0,0)`.  The common vertex `H=delta_(0,0)` also has a
faithful normalized Hasse coordinate after translation, but it is not the
mixed face.  Thus neither faithful coordinate is yet a nonzero
common-atom mixed contrast.

This is distinct from the integral central filtration

```text
(v00,v10,v01,v11)=(p^2,p,p,1),

(v_p(v00),v_p(v10),v_p(v01),v_p(v11))=(2,1,1,0).       (24)
```

For `p=13`, its Hadamard mixed entry is

```text
D3=(p-1)^2=144=1 mod13.                                (25)
```

The augmentation degrees `(20)`, the `p`-adic valuations `(24)`, and the
scalar residue `(25)` are three different objects.

## 4. Positive survivor: a rooted face converts translation into a jet

THM-2201 supplies the exact positive mechanism.  For

```text
F=sum_(j=0)^(p-1) J_j(F) eta^j                         (26)
```

with `J_0(F)!=0`, define its normalized first Hasse coordinate

```text
b(F)=J_1(F)/J_0(F).                                    (27)
```

Multiplication by `g^a=(1+eta)^a` gives

```text
J_0(g^aF)=J_0(F),
J_1(g^aF)=J_1(F)+aJ_0(F),

b(g^aF)=b(F)+a.                                        (28)
```

Thus `b` is an affine `C_p` torsor coordinate.  Its absolute value depends
on the root origin, but relative differences

```text
b(g^aF)-b(F)=a                                         (29)
```

are origin-gauge invariant.

There is a useful general form.  Put

```text
A_2=F_p[eta,theta]/(eta^p,theta^p)
```

with its total augmentation filtration.  Let `F!=0`, let `f_d` be its first
nonzero total graded piece, and let

```text
T_(a,b)F=(1+eta)^a(1+theta)^b F.
```

Then direct expansion in degree `d+1` gives

```text
gr_(d+1)(T_(a,b)F-F)=(a eta+b theta)f_d.              (29a)
```

Consequently the kernel of the first-normal action is
`Ann_1(f_d)`, its image is `(A_2)_1 f_d`, and its first-normal parameter
space is

```text
(A_2)_1/Ann_1(f_d),                                   (29b)
```

where `(A_2)_1=span(eta,theta)`.  For the five allocation vectors in
`(18),(21)` this gives

```text
             B     P     Q     H     Omega
normal rank  0     1     1     2       2.             (29c)
```

Indeed `B` is killed by both linear directions, `P` and `Q` each retain
only the transverse axis, while `H` and `Omega` have leading piece `1`.
This translation differential is not an internal Boolean tangent:
Section 1 proves that the latter vanishes.  For a general `F`, a direction
killed in `(29a)` may reappear at higher degree.  For these five displayed
vectors, however, the ranks in `(29c)` also give their exact translation
stabilizers.

For two axes use `b=(b_1,b_2)`.  Equation `(22)` has

```text
J_(0,0)=1,                  J_(1,0)=J_(0,1)=0,          (30)
```

so

```text
b(Phi(Omega))=(0,0),

b(g^a h^b Phi(Omega))=(a,b).                           (31)
```

All `169` translations are distinguished exactly as rooted coefficient
vectors.  In fact this is a full-orbit fact, not merely a first-order one.
Writing `N_eta=eta^12` and `N_theta=theta^12`, one has

```text
Phi(Omega)=(N_eta-1)(N_theta-1),
Phi(Omega)^(-1)=(N_eta+1)(N_theta+1),                  (31a)
```

because `N_eta^2=N_theta^2=0`.  Hence `Phi(Omega)` is a unit.  A stabilizer
equation may be multiplied by its inverse, proving that its `C_13^2`
translation orbit is regular.  Every nonzero one-dimensional translation
line is therefore a principal `C_13` torsor.  In the sense of
THM-2611/2792, a physical basepoint and generator intertwiner would determine
the map to THM-2813's off-sheet normal torsor uniquely; without a basepoint
there are thirteen origin choices.

There is also an exact integral shadow of this inverse.  In
`Z[C_p^2]`, let

```text
Theta=(N_eta+1)(N_theta+1)
```

and let `pi_lambda` push coefficients along
`lambda(a,b)=a+b`.  The coefficients of `Theta` are positive weights
`1,2,4`.  Since `N^2=pN` over `Z`, direct convolution gives, for every
`q=(q_1,q_2)`,

```text
pi_lambda(Theta*T_q Omega)
 =p(p^2-2)N+delta_(q_1+q_2).                           (31b)
```

At `p=13` the baseline is `2171N`.  Thus the one exceptional coordinate
recovers `q_1+q_2` exactly over the integers, not merely modulo thirteen.
It is the canonical `ell`-forgotten quotient decoder and kills the marked
gauge direction `(1,-1)` sharply.  By the linear pushforward boundary of
THM-2814, however, `(31b)` remains sourced in the joint-absent aggregate;
positive inverse weights do not turn it into common-atom/root-Cech
physicality.

The normalized coordinate `(27)` requires `J_(0,0)` to be a unit.  Here it
is literally one.  A physically scaled face `w Omega` may be normalized only
after retaining a coefficient line whose reduction `w mod13` is nonzero;
an arbitrary physical scalar must not silently be inverted.

The lawful carrier-address gauge is sharper than merely
forgetting a common origin.  THM-2763 gives

```text
(ell,b_1,b_2)~(ell+sW,b_1+s,b_2-s).                   (32)
```

Choose `z` with `z.W=1`.  The exact gauge-invariant coordinate is

```text
kappa_z=(b_1-z.ell,b_2+z.ell).                        (32a)
```

For a one-parameter family

```text
ell_t=ell_0+tL,        F_t=T_(tq_1,tq_2)F_0,
```

its observable first-normal increment is

```text
kappa_z(t)-kappa_z(0)
 =t(q_1-z.L,q_2+z.L)=t q_eff.                         (32b)
```

Thus the marked motion `(L,q)=(W,(1,-1))` is **pure gauge**:
`q_eff=(0,0)`.  If `ell` is forgotten, only `b_1+b_2` descends, and the
sharp visibility test is `q_1+q_2!=0`.  A faithful rooted Hasse chart does
not by itself prove that a proposed physical motion is transverse to the
carrier gauge.  The vector `z` is a chosen section (canonically `e_0` only
on THM-2763's typed row), so `kappa_z` is not itself section-free.  Also,
`q_eff=0` only says that this first-Hasse endpoint component is flat; unless
`(L,q)` is actually a multiple of `(W,(1,-1))`, the full `(ell,b)` quotient
may still move.

Equations `(30)--(31)` explain the role of `(25)`: `D3/v11=1` is the
nonzero denominator `J_(0,0)` which permits normalization.  It is not the
numerator jet.  The jet appears only after a rooted **translation of the
whole raw face**.

## 5. Conditional bridge to THM-2813

THM-2813's relative affine lifts obey

```text
A_t(y)=y+t13^5(y-7 mod13).
```

At an off-sheet atom with

```text
r=y-7 mod13 !=0,                                       (33)
```

its normal displacement is `tr`.  Suppose a future physical construction
provides:

1. a rooted coefficient aggregate `F_0` with nonzero augmentation;
2. the same endpoint origin, selector, allocation flags, clock, and word;
3. a lawful address path `ell_t=ell_0+trL` and carrier translation
   `F_t=T_(trq_1,trq_2)F_0`; and
4. a nonzero effective transverse vector

   ```text
   q_eff=(q_1-z.L,q_2+z.L)!=0.                         (34)
   ```

Choose any linear functional `lambda` on `F_13^2` with
`lambda(q_eff)=1`.  Then `(32b)` gives the exact decoder

```text
beta_normal(F_t)
 =r^(-1)lambda(kappa_z(F_t,ell_t)-kappa_z(F_0,ell_0))
 =t.                                                    (35)
```

At the residue-eight sheet `r=1`, one such rooted aggregate suffices.  If
the address is discarded, `(34)` must be replaced by the stronger visible
test `q_1+q_2!=0`.

This coefficient-side decoder is not yet a common-atom beta source.
`Omega` is supported on the 144 joint-absent cells, while the common
`H`-vertex is not the mixed face.  A physical allocation consequence
therefore needs either a future non-idempotent common raw face of the kind
isolated in THM-2814 Branch A, or an explicitly broader two-object sidecar
which lawfully couples a rooted aggregate to the common component.
THM-2658 certifies common-component structure, not a nonzero contrast.

Finally, the total-degree augmentation filtration used in `(29a)` is the
modular convolution filtration.  It is not THM-2806's central `13`-adic
Rees valuation filtration `(24)`.  A physical bridge must intertwine those
objects; equality of their displayed scalar residues does not do so.

This is a genuine conditional exit, not a closure.  THM-2806 proves that
literal allocation has zero endpoint translation, so its currently typed
family remains at `b=(0,0)`.  The known marked direction `(W,(1,-1))` is
also killed exactly by `(32b)`.  THM-2813 supplies an oriented normal
generator, but canon contains neither a physical nonzero `q_eff` nor the
required common raw face/two-object coupling.  The generator is external
sidecar data; it is not recovered from the norm by `(12)`.

## 6. Information and failure ledger

| object | retained information | first missing datum |
|---|---|---|
| pointwise Boolean algebra | exact atoms and idempotent masks | every internal tangent vanishes |
| modular orbit norm `N` | top augmentation class `eta^12` | nonzero vector in `I/I^2` |
| central Rees profile | valuations `(2,1,1,0)` and scalar unit `(25)` | rooted first-Hasse numerator |
| rooted raw face `Omega` | regular `C_13^2` translation torsor on 144 joint-absent cells | non-idempotent common raw face or broader coupling |
| normalized Hasse `b` | every supplied rooted translation `(a,b)` | carrier-gauge-effective vector `q_eff` |
| gauge invariant `kappa_z` | every transverse translation modulo `(W,1,-1)` | lawful physical family with `q_eff!=0` |
| THM-2813 normal jet | affine-lift label `t` off the fixed sheet | common physical coefficient/coupling and transverse map |

THM-2814 studies projective four-corner holonomy and provenance.  The
present theorem is nonduplicate: it classifies formal Boolean tangents, the
norm-to-cotangent obstruction, and the positive rooted Hasse translation
coordinate.

It does not identify the two cyclic generators, construct the family in
`(34)`, produce a non-idempotent common raw face, give a root/Cech map,
exclude a relation row, or prove LRC(14).

## 7. Exact companion

Run

```text
python 04-computation/lrc14_boolean_norm_cotangent_boundary_thm2820.py
python -O 04-computation/lrc14_boolean_norm_cotangent_boundary_thm2820.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_boolean_norm_cotangent_boundary_thm2820.out.
```

The dependency-free companion:

1. exhausts all idempotent tangent kernels in `F_p^4` for
   `p=2,3,5,7,11,13`;
2. verifies `N=eta^(p-1)`, every generator automorphism, and every fixed
   cotangent vector at those primes;
3. reconstructs `(20)--(25)` at `p=13`;
4. checks the non-idempotent interpolation defect on all thirteen atoms;
5. computes all `169` translated raw-face Hasse coordinates in `(31)`;
6. verifies `(29a)` on all `845` allocation-state/translation pairs, the
   rank list `(29c)`;
7. checks the unit inverse `(31a)`, regular `169`-element orbit, and all
   `2,197` entries of the integral push profiles `(31b)`;
8. exhausts the `28,561` THM-2763 gauge-invariant triples and `28,561`
   effective-family instances; and
9. verifies the pure-gauge and forgotten-address visibility boundaries.

It uses explicit exception gates, no Python `assert`, no floating point, and
no scratch dependency.

**Awaiting independent hostile audit; not QED.**
