---
id: THM-3532
title: "Fixed Keller two-sided one-step covariance and all-level conjugacy transport"
status: >
  PROVED + VERIFIED-EXACT.  Independent polynomial source/target changes
  transport the fixed sporadic map's norm, finite divisor and unit, one
  prime-image edge, one-step intrinsic discriminant/different, and Jelonek
  divisor.  Honest polynomial conjugacy transports the complete raw packet,
  image-prime, intrinsic discriminant, reduced-different, and Jelonek towers.
  The five weights must move with the coordinates.  W1/W2 are target
  postcompositions and sharp one-step controls, not self-conjugacies.
source: codex/tame-covariance-keller/2026-08-16
audit: >
  The hostile review found that the first companion used optimization-erased
  Python asserts while comparing normal and optimized transcripts.  The
  repaired companion uses only explicit require/raise gates; normal, -O, and
  stored transcripts agree.  Exact W1/W2 witnesses separate target
  postcomposition from conjugacy, and both nonlinear conjugates plus an
  affine conjugate pass the exact finite control bank.
depends_on:
  - THM-2465-g1-exclusion-package-for-degree-four-twojet-keller
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-2576-composite-jelonek-image-divisor-and-two-component-nonproperness-law
  - THM-3528-fixed-keller-all-level-cleared-norm-polynomiality-and-finite-sheet-defect
  - THM-3529-fixed-keller-complete-packet-finite-sheet-unit
  - THM-3530-fixed-keller-all-level-image-prime-and-component-tower
  - THM-3531-fixed-keller-intrinsic-all-level-discriminant-square-class
  - THM-3533-fixed-keller-newest-prime-reduced-different-and-index-square
  - THM-3535-fixed-keller-full-wreath-and-all-level-linear-primitivity
related:
  - MISTAKE-413
  - MISTAKE-419
  - HYP-9033-discriminant-tower-and-genus-axis-of-the-keller-monoid
script: 04-computation/keller_tame_conjugacy_covariance_audit_20260816.py
output: 05-knowledge/results/keller_tame_conjugacy_covariance_audit_20260816.out
script_sha256: 76878a0190c065afaa85b1ce5c67d83a39122441684e4bb67e5211b51ca0b9c0
output_sha256: a6d5a1db30d33adc2a6e6315e08ed311b49eda8bc80eb43c2dd8ed4141f70b43
hash_basis: LF-normalized bytes
---

# THM-3532 -- two-sided covariance is one-step; conjugacy transports the tower

**PROVED + VERIFIED-EXACT.**

Work over `k=Q`.  Retain the fixed sporadic Keller self-map

```text
F:A^3->A^3,                 det J_F=-2,                 (1)
```

its irreducible target boundary `L`, and the raw tower

```text
P_0=L,
P_(n+1)=L^e_n N_F(P_n),
(e_(n+1),m_(n+1))=(7e_n-2m_n,3e_n-2m_n).              (2)
```

For an automorphism `rho` of `A^3`, define the pushforward of an equation by

```text
rho_#P=P o rho^(-1).                                   (3)
```

Then `V(rho_#P)=rho(V(P))`.  The result below holds for every polynomial
automorphism over `Q`; affine-linear and tame automorphisms are subclasses.

## 1. Independent source and target changes: exact one-step theorem

Let `sigma` and `tau` be polynomial automorphisms and set

```text
G=tau o F o sigma,                 L_G=tau_#L.          (4)
```

Then `G` is generically degree three and Keller, with

```text
det J_G=(det J_tau)(det J_F)(det J_sigma) in Q^*.       (4a)
```

In particular an honest conjugate has the original determinant `-2`.

An old source equation `P` is represented on the new source by

```text
sigma^*P=P o sigma=(sigma^(-1))_#P.                    (5)
```

For every nonzero rational function `P`,

```text
N_G(sigma^*P)=tau_#N_F(P).                             (6)
```

If `P` has the complete packet `A(e,m)`, then

```text
L_G^e N_G(sigma^*P)=tau_#(L^eN_F(P)).                 (7)
```

Thus THM-3528 transports from the source chart `sigma^*A(e,m)` to the target
chart

```text
tau_#A(7e-2m,3e-2m).                                  (8)
```

This is equality of polynomials, not merely equality of divisors or equality
up to a rational scalar.

Let

```text
B=F^*L=L o F                                           (9)
```

be THM-3529's irreducible finite-sheet equation.  The new finite sheet and
its valuation are

```text
G^*L_G=B o sigma,
C_fin,G=sigma^(-1)(C_fin,F),                           (10)
s_(L_G),G(sigma^*P)=s_L,F(P).                          (11)
```

Consequently every transported complete packet is a finite-sheet unit:

```text
s_(L_G),G(sigma^*P)=0.                                 (12)
```

The geometric conclusion survives even when the standard-coordinate
beta/x-degree proof of THM-3529 no longer has the same visible faces.

For every old raw rung, THM-3530's corresponding image edge transports as

```text
closure(G(V(P_n o sigma)))=V(tau_#P_(n+1)),            (13)
deg(V(P_n o sigma)->V(tau_#P_(n+1)))=1.               (14)
```

The source and target isomorphisms preserve absolute irreducibility and the
generic degree.  They also preserve properness, giving

```text
S_G=tau(S_F)=V(L_G).                                   (15)
```

Finally, trace pairings and integral closures commute with the same field
isomorphisms.  The one-step intrinsic square class and normalization
discriminant exponent are therefore

```text
d(G)=tau_#d(F)=[-L_G],                                 (16)
delta_G(V(L_G))=1.                                    (17)
```

Equation `(17)` is THM-3533's reduced different transported to the new
target.  No Jacobian factor appears in `(6)`, `(16)`, or `(17)`.

### Proof of the one-step formulas

Use separate copies `K_t=Q(a,b,c)` and `K_s=Q(x,y,z)` of the rational
function field.  Let

```text
U_sigma(P)=P o sigma,              Psi_tau(Q)=Q o tau^(-1). (17a)
```

The pullbacks defining the two degree-three extensions form the commuting
square

```text
K_t  --F^*-->  K_s
 |              |
Psi_tau       U_sigma
 |              |
 v              v
K_t  --G^*-->  K_s,                                      (17b)

G^* Psi_tau=U_sigma F^*.
```

Both vertical arrows are field isomorphisms and restrict to polynomial-ring
isomorphisms.  Thus `(17b)` is an isomorphism of finite separable field
extensions, not just a bijection on a sample of fibres.  Naturality of the
determinant of multiplication gives

```text
Norm_G(U_sigma P)=Psi_tau(Norm_F(P)),                    (17c)
```

which is `(6)`.

For a new target point `t`,

```text
G(q)=t  iff  F(sigma(q))=tau^(-1)(t).                  (18)
```

The isomorphism `q -> sigma(q)` identifies the complete generic fibres.  The
product of `P o sigma` over the first fibre is therefore the product of `P`
over the second fibre, evaluated at `tau^(-1)(t)`.  This proves `(6)`, and
`L_G=L o tau^(-1)` proves `(7)`.

The fibre argument is a geometric reading of `(17c)`.  The same commuting
square sends the regular finite branch to its inverse image under `sigma`,
proving `(10)`--`(11)`.  Applying `G` to
`sigma^(-1)(V(P_n))` proves `(13)`--`(14)`.  Isomorphisms are proper, so an
escape sequence for `G` is exactly an escape sequence for `F` with its target
limit moved by `tau`, proving `(15)`.  The induced base/source field
isomorphisms intertwine traces; a basis change alters a trace determinant by
a square.  More explicitly, if `B_F` is the integral closure of the old
target polynomial ring in `K_s`, then `U_sigma(B_F)` is the integral closure
for `G` over `Psi_tau(A)`.  The trace Gram determinant, different ideal, and
their height-one valuations are carried by these two ring isomorphisms.
This proves `(16)`--`(17)` and completes the one-step proof.

## 2. Honest conjugacy: all-level covariance

Now take one polynomial automorphism `phi` and specialize

```text
sigma=phi^(-1),             tau=phi,
F^phi=phi o F o phi^(-1),   Psi=phi_#.                 (19)
```

The source and target transports coincide.  Formula `(6)` becomes

```text
N_(F^phi) Psi=Psi N_F.                                 (20)
```

Put

```text
L^phi=Psi L,                 P_n^phi=Psi P_n.          (21)
```

Then `(2)`, `(7)`, and `(20)` give the exact raw recurrence

```text
P_0^phi=L^phi,
P_(n+1)^phi=(L^phi)^e_n N_(F^phi)(P_n^phi).            (22)
```

For every `n>=0`,

```text
P_n^phi is absolutely irreducible;                     (23)
closure(F^phi(V(P_n^phi)))=V(P_(n+1)^phi);             (24)
the restriction in (24) has generic degree one;        (25)
s_(L^phi),F^phi(P_n^phi)=0.                            (26)
```

For every `n>=1`, THM-3531 and THM-3533 become

```text
d_n^phi=[(-1)^n P_(n-1)^phi],                          (27)
delta_n^phi(V(P_(n-1)^phi))=1.                         (28)
```

The local index-square formula also transports.  If `theta` is integral and
primitive at the old newest prime and `theta^phi` is its image, then the
local index lengths agree and

```text
v_(P_(n-1)^phi)(Disc(theta^phi))=1+2i_(theta,n).       (29)
```

In particular, locally maximal primitive elements still attain multiplicity
one, but an arbitrarily chosen literal coordinate need not do so.

THM-3535's all-level observation theorem also transports.  For every nonzero
old constant linear form `ell`,

```text
ell^phi=Psi ell=ell o phi^(-1)                         (29a)
```

is primitive for the conjugated extension at every depth.  If `phi` is
affine-linear, every nonzero standard linear form in the new coordinates is,
up to adding a rational constant, of this shape; hence all standard linear
directions remain primitive.  If `phi` is nonlinear, `(29a)` is generally a
polynomial observation, so THM-3535 gives no automatic statement about every
untransported standard linear form.  The local index of the transported
observation is preserved but remains whatever the old index was.

Since

```text
(F^phi)^r=phi o F^r o phi^(-1),                        (30)
```

THM-3530's Jelonek law transports at every depth:

```text
S_((F^phi)^r)
 =phi(S_(F^r))
 =union_(j=0)^(r-1)V(P_j^phi)
 =V(product_(j=0)^(r-1)P_j^phi).                      (31)
```

It has exactly `r` reduced irreducible components.

## 3. The five-weight chart is a required sidecar

Write

```text
(X,Y,Z)=phi^(-1)(x,y,z).                               (32)
```

On monomials `X^iY^jZ^k`, use the transported weights

```text
lambda=i-k,       beta=i-j-2k,       gamma=i-j-5k,     (33)
```

and retain `k` as the `Z`-exponent.  Equivalently, compose a new polynomial
with `phi` before applying an old weight.  In this chart a packet of grade
`(e,m)` has the five complete faces

```text
X^e(3XZ-2Y)^m,
Y^(3e-2m)Z^e,
Y^(3e-2m)Z^(e-m)(Y^2+27Z)^(2m/3)(Y^2+108Z)^(m/3),
X^(2e-4m/3)Z^(2e-2m/3),
Z^e(27X^2Z+Y^3)^(e-2m/3),                              (34)
```

up to the independent nonzero scalars allowed in THM-3506.  Hence the grade
matrix, grade primitivity, Pell-57 identity, and Cassini identity survive.

They need not survive in the unmoved standard `(x,y,z)` weights.  Nonlinear
tame maps can change every visible extremum; affine maps can change literal
face equations.  Standard total degrees, multidegrees, term counts, and
primitive integer contents are not conjugacy invariants.

For independent `sigma,tau`, the source atlas is `sigma^*A` and the target
atlas is `tau_#A`.  They are the same atlas exactly when

```text
(sigma o tau)_#A=A.                                   (35)
```

This atlas equality alone does not identify the self-iterate tower.

## 4. The universal postcomposition/conjugacy boundary

For `G=tau F sigma`, direct composition gives

```text
G^r=tau F (sigma tau)F ... (sigma tau)F sigma,         (36)
```

with `r-1` inserted copies of `sigma tau`.  In old `F`-coordinates, feeding
the output of `(7)` back as a new source input replaces the raw operator by

```text
P -> (sigma tau)_#(L^eN_F(P)).                         (37)
```

The honest conjugacy condition `sigma=tau^(-1)` removes the insertion and is
therefore the universal, assumption-free all-level case.  For another
two-sided presentation, an all-level identification with the *old* tower
requires an additional proved symmetry or intertwiner for `(36)`--`(37)`.
The words “affine,” “elementary,” and “tame” do not supply one.  This is a
sufficiency boundary for the theorem, not a claim that no exceptional
non-conjugate presentation can have such an intertwiner.

| object | independent `tau F sigma` | honest `phi F phi^-1` |
|---|---|---|
| raw norm packet | one edge, two transported charts | every rung, one chart |
| finite-sheet unit | transported complete inputs | every raw rung |
| image prime | each transported old edge | complete prime tower |
| intrinsic class / newest different | one step | every depth |
| Jelonek set | `S_G=tau(S_F)` | every iterate/component |
| five standard weights | generally lost | recovered only in chart `(32)` |

## 5. W1/W2: exact hostile and positive controls

THM-2465 uses the elementary target automorphisms

```text
T_1(a,b,c)=(a,b,c+a^2),
T_1^(-1)(a,b,c)=(a,b,c-a^2),                           (38)

T_2(a,b,c)=(a+bc,b+c^2,c),
T_2^(-1)(a,b,c)=(a-bc+c^3,b-c^2,c).                   (39)
```

Both have determinant one, and

```text
W_1=T_1 o F=(F_1,F_2,F_3+F_1^2),
W_2=T_2 o F=(F_1+F_2F_3,F_2+F_3^2,F_3).               (40)
```

These are target postcompositions, not conjugacies as self-maps.  With

```text
L_i=L o T_i^(-1),             S_i=S o T_i^(-1).        (41)
```

the one-step identities are

```text
S_(W_i)=V(L_i),
W_i^*L_i=F^*L=B,
N_(W_i)(P)=(N_F(P)) o T_i^(-1),                       (42)
L_i^eN_(W_i)(P)=(L^eN_F(P)) o T_i^(-1).               (43)
```

The source coordinates did not move, so their three coordinate cores are
target pullbacks.  This explains the earlier exact identities

```text
Disc= -4S_i^2L_i=-(det J_(W_i))^2S_i^2L_i.            (44)
```

The standard packet does move.  In the order
`(max lambda,min lambda,min beta,max k,min gamma)`, the seed extrema are

```text
L:   (1,-1,-5,2,-8),
L_1: (6,-1,-5,2,-8),
L_2: (1,-8,-16,8,-40).                                (45)
```

The new maximum face of `L_1` is `27x^6`; the four changed `L_2` faces are
led by `27z^8`.  Composing with `T_i` recovers `L` and `A(1,0)` exactly.

Postcomposition fails already at the second iterate:

```text
W_1^2(0,0,-1)=(1,-3,0),       T_1F^2(0,0,-1)=(0,0,-2),
W_2^2(-1,0,0)=(62,4,0),       T_2F^2(-1,0,0)=(-2,0,0). (46)
```

By contrast, the honest conjugates

```text
C_i=T_i o F o T_i^(-1)                                (47)
```

inherit the full theorem with

```text
P_(n,i)=P_n o T_i^(-1),
S_(C_i^r)=T_i(S_(F^r)),
d_(n,i)=[(-1)^n(P_(n-1) o T_i^(-1))],                 (48)
```

and newest-prime normalization-discriminant multiplicity one.

## 6. Scalar normalization is not invisible in square classes

The exact tower `(21)` uses `L^phi=Psi L`.  If a primitive-integral or other
convention replaces it by

```text
L_tilde=cL^phi,                    c in Q^*,            (49)
```

and rebuilds the raw recurrence, then

```text
P_tilde_n=c^q_nP_n^phi,
q_0=1,                 q_(n+1)=e_n+3q_n.               (50)
```

The first values are

```text
q_n=1,4,19,100,571,3412,20899,...,                    (51)
q_n=n+1 mod 2.                                         (52)
```

Thus `(27)` becomes

```text
d_n^phi=[(-1)^n c^(-q_(n-1))P_tilde_(n-1)].           (53)
```

The constant `[c]` survives exactly for odd `n`.  Deleting it would repeat
MISTAKE-413.  The zero divisor, image prime, generic degree, reduced-different
multiplicity, and Jelonek component count are scalar-insensitive.

## 7. Scope and losses

The theorem preserves the fixed map's fibre degree, norm/trace algebra,
finite-sheet order, prime ancestry, intrinsic square class, newest-prime
different, local index length, and reduced component count.  The transported
packet chart preserves the two grade coordinates.

It does not preserve, without additional sidecars, standard Newton weights,
multidegrees, term counts, named integral normalizations, literal-coordinate
primitivity/separability/local maximality under nonlinear conjugacy, or
self-iteration under independent left/right equivalence.  It proves nothing
for an arbitrary Keller map,
another tame-equivalence class, weighted-lift atoms in a numerical grade,
old-prime positive discriminant multiplicities, `JC(2)`, `DC(2)`, LRC, or a
classification of Jacobian counterexamples.

## Reproduction and validity gate

```text
python -B 04-computation/keller_tame_conjugacy_covariance_audit_20260816.py
python -B -O 04-computation/keller_tame_conjugacy_covariance_audit_20260816.py
```

The repaired companion contains no executable Python `assert`; every
truth-bearing check uses an explicit `require` that raises under both modes.
Normal and optimized transcripts match each other and the stored output.  It
checks `(38)`--`(46)`, the standard/transported packet signatures, all `27`
points of `{-1,0,1}^3` for both nonlinear honest conjugates and one affine
triangular conjugate, and `(50)`--`(52)`.  The finite banks audit the
implementation; Sections 1--4 give the general proof.

**QED.**
