---
id: THM-2820
title: "Boolean rigidity, norm-top cotangent no-go, and rooted Hasse translation coordinate"
status: >
  RESERVED / PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT / AWAITING
  INDEPENDENT HOSTILE AUDIT.  Finite pointwise Boolean carrier algebras
  have zero relative derivations and rigid idempotents.  In F_p[C_p] the
  orbit norm is (g-1)^(p-1), fixed by every generator change although the
  cotangent line scales, so norm/Rees data alone select no nonzero first
  jet.  Positively, any rooted group-ring vector of nonzero augmentation
  has normalized Hasse coordinate b=J1/J0 with b(g^aF)=b(F)+a.
  THM-2806's normalized raw mixed face has base b=(0,0) and all 169
  translates recover their labels.  A lawful common off-sheet translation
  and allocation-to-address identification remain missing.
source: root/boolean-norm-cotangent-boundary-2026-07-28
depends_on:
  - THM-2201-cyclic-root-fibre-hasse-jet-transition-carrier
  - THM-2806-literal-fixed-sheet-central-allocation-scalar-law-and-endpoint-translation-no-go
  - THM-2813-affine-lift-transvection-and-projective-horn-decoder
related:
  - THM-2771-joint-c7-c13-right-wing-mixed-spectrum-and-commuting-square-no-go
  - THM-2779-bockstein-symplectic-decoder-frame-torsor-and-heisenberg-root-degree-gate
  - THM-2806-literal-fixed-sheet-central-allocation-scalar-law-and-endpoint-translation-no-go
  - THM-2813-affine-lift-transvection-and-projective-horn-decoder
  - THM-2814-projective-allocation-square-holonomy-and-idempotent-provenance-no-go
script: 04-computation/lrc14_boolean_norm_cotangent_boundary_thm2820.py
output: 05-knowledge/results/lrc14_boolean_norm_cotangent_boundary_thm2820.out
script_sha256: e38ada505f0ad729b6027aab1282779f3d9ee2e10e00785c292ee3735971783c
output_sha256: 77d05c9c77dd7fd696137b76f6f3a022ebf790f01ce63fc1e4fd9699acbd86f2
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
affine label.  The remaining problem is physical co-support and
equivariance, not another scalar decoder.

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

For two axes use `b=(b_1,b_2)`.  Equation `(22)` has

```text
J_(0,0)=1,                  J_(1,0)=J_(0,1)=0,          (30)
```

so

```text
b(Phi(Omega))=(0,0),

b(g^a h^b Phi(Omega))=(a,b).                           (31)
```

All `169` translations are distinguished exactly.  Under THM-2806's marked
gauge

```text
(a,b)->(a+1,b-1),                                      (32)
```

every absolute coordinate shifts by the same `(1,-1)`, so the relative
differences `(29)` survive.

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

1. one nonzero common raw allocation atom at that same `y`;
2. the same endpoint origin, selector, allocation flags, clock, and word;
3. a lawful equivariant identification of the relevant allocation
   translation line with an injective line

   ```text
   iota:F_13 -> F_13^2;
   ```

4. the covariance

   ```text
   Phi(Omega_t)=g^(iota_1(tr))h^(iota_2(tr))Phi(Omega_0). (34)
   ```

Then `(31)` gives the exact decoder

```text
beta_normal(Omega_t)
 =r^(-1)iota^(-1)(b(Omega_t)-b(Omega_0))
 =t.                                                    (35)
```

At the residue-eight sheet `r=1`, a single atom suffices.

This is a genuine conditional exit, not a closure.  THM-2806 proves that
literal allocation has zero endpoint translation, so its currently typed
family remains at `b=(0,0)`.  THM-2813 supplies an oriented address
generator, and THM-2806 separately supplies a marked carrier-twist
generator, but canon contains no physical map identifying them on a common
off-sheet atom.  The generator is external sidecar data; it is not recovered
from the norm by `(12)`.

## 6. Information and failure ledger

| object | retained information | first missing datum |
|---|---|---|
| pointwise Boolean algebra | exact atoms and idempotent masks | every internal tangent vanishes |
| modular orbit norm `N` | top augmentation class `eta^12` | nonzero vector in `I/I^2` |
| central Rees profile | valuations `(2,1,1,0)` and scalar unit `(25)` | rooted first-Hasse numerator |
| rooted raw face `Omega` | denominator and all top provenance classes | a lawful translated physical family |
| normalized Hasse `b` | every supplied translation `(a,b)` | allocation-to-address line `iota` |
| THM-2813 normal jet | affine-lift label `t` off the fixed sheet | common physical allocation atom |

THM-2814 studies projective four-corner holonomy and provenance.  The
present theorem is nonduplicate: it classifies formal Boolean tangents, the
norm-to-cotangent obstruction, and the positive rooted Hasse translation
coordinate.

It does not identify the two cyclic generators, construct `(34)`, produce a
root/Cech map, exclude a relation row, or prove LRC(14).

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
5. computes all `169` translated raw-face Hasse coordinates in `(31)`; and
6. verifies all `169` relative marked-gauge identities.

It uses explicit exception gates, no Python `assert`, no floating point, and
no scratch dependency.

**Awaiting independent hostile audit; not QED.**
