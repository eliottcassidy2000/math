---
id: THM-4299
title: "D-zero square-face elliptic splitting and off-corner Keller extinction"
status: >
  PROVED RELATIVE TO THM-4103/4230/4297 + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED. On the exact-weight-twelve wall D=0 with U*Z*Lambda!=0, the
  genus-seven main factor splits into two j=0 elliptic components. Both are
  Keller-constant. Their unique order-six mutual contact admits an exact
  parameterized-Morse model w^2=z^12-kappa(sigma*z); every positive-genus
  Newton tail has strictly positive good-differential order 9s+5beta and is
  therefore constant. All remaining components are rational, so proper-flat
  degree conservation excludes the stratum. The cubic corner D=Lambda=0,
  the U=0 and Z=0 walls, exact-M=12 seam entry, and JC(2) remain open.
source: root / planar Jacobian signal session, 2026-08-31
depends_on:
  - THM-4103-jc23-theta-boundary-ramification-and-degree-response
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
  - THM-4297-general-lambda-zero-central-and-tail-keller-extinction
related:
  - THM-4290-exact-weight-twelve-deck-equivariant-visible-quotient-exclusion
  - THM-4298-weighted-face-source-normal-unimodular-visibility-transform
primary_script: 04-computation/jc23_d_zero_square_face_elliptic_extinction_thm4299.py
primary_output: 05-knowledge/results/jc23_d_zero_square_face_elliptic_extinction_thm4299.out
primary_script_sha256: f140d0ba7a559d0788a2c310f70c273d22e2204937bac49c88dc44d3ef8b650d
primary_output_sha256: ca836de0ea12fa6af84d528c183f38a6a11693e565d7d4333a3fdccfc4126500
independent_audit_script: 04-computation/jc23_d_zero_square_face_elliptic_extinction_independent_audit_thm4299.py
independent_audit_output: 05-knowledge/results/jc23_d_zero_square_face_elliptic_extinction_independent_audit_thm4299.out
independent_audit_script_sha256: 002058a4d959229e7716fc64a4b25ad25beb8c6aa24ea6946507471612f3d3df
independent_audit_output_sha256: 34795b884a0618ea052855c2a58e6767264a75c9fb08b37642965b3a941fceb2
hash_basis: raw LF bytes
audit: >
  PASS. The primary SymPy certificate checks the square factorization,
  smooth elliptic normalizations, Newton genera, deck characters, all
  contact and arithmetic-genus ledgers, the Morse Hessian, every possible
  m=1,...,11 compact tail, the m=1 genus-five hostile, and the cubic failure
  boundary. A dependency-free sparse-polynomial and integer audit
  independently reconstructs the factorization, polygon, character,
  contact, tail, hostile, and cubic ledgers. Normal, optimized, and fixed-
  hash-seed streams byte-match both frozen outputs. Two independent
  adversarial derivations audited the formal Morse and component-inventory
  arguments.
---

# THM-4299 -- D-zero square-face elliptic splitting and off-corner extinction

**PROVED RELATIVE TO THM-4103/4230/4297 + VERIFIED-EXACT + INDEPENDENTLY
AUDITED. THE EXACT-`M=12` STRATUM `D=0`, `U*Z*Lambda!=0` IS EXCLUDED. THE
CORNER `D=Lambda=0`, THE ENDPOINT WALLS, SEAM ENTRY, AND `JC(2)` REMAIN
OPEN.**

## 1. Statement and inheritance

Retain the exact-weight-twelve reduced `(2,3)` seam and put

```text
D=W^2-4UZ,                         Lambda=U+W+Z.             (1)
```

> **Theorem.** There is no nonautomorphic planar Keller candidate on
>
> ```text
> D=0,                             U*Z*Lambda!=0.             (2)
> ```

Together with THM-4290 and THM-4297, this leaves inside the region `U*Z!=0`
only the codimension-two corner `D=Lambda=0` at exact weight twelve.

The inheritance pass is:

- closest proved mechanism: THM-4297's good-differential extinction followed
  by proper-flat degree conservation;
- canonical hostile: the literal `m=1` genus-five collision tail in Section
  7, which shows that splitting the central curve does not make every tail
  rational;
- corrected near miss: MISTAKE-531 forbids transporting a simple-root
  response packet through a repeated top-edge root, while MISTAKE-521 forbids
  replacing the actual component map by a rational isogeny decomposition;
- least-used decisive sidecar: divide by the now-unit factor `r-1` and retain
  the one-variable critical value `kappa(t)` before reading the wall.

The live concept board was

```text
{square face, deck character, A11 critical value,
 good differential, degree conservation, cubic corner}.                  (3)
```

The decisive connection sends the exact collision germ to its critical-value
series by the parameterized Morse map. It preserves every divisorial
valuation and the good-differential order, destroys the original lower-row
labels, and needs `kappa(t)` as its restoration sidecar. The cheapest hostile
is therefore a nonzero first critical value, not another generic coefficient
sample.

## 2. The discriminant wall is a square face

Work over `C`. Under `(2)`, choose `a,b in C*` so that

```text
U=a^2,                    W=2ab,                    Z=b^2.   (4)
```

Changing one square-root sign realizes either sign of `W`. Put

```text
L=aP^3+bS^2P^2.                                             (5)
```

Then the exact central factor and the remaining wall coordinate are

```text
C=1-UP^6-WS^2P^5-ZS^4P^4=(1-L)(1+L),
Lambda=(a+b)^2.                                             (6)
```

Write `T=SP`. The normalization of the two factors is

```text
E_epsilon: aP^3+bT^2=epsilon,              epsilon in {+1,-1}. (7)
```

The projective cubic

```text
aP^3+bT^2H-epsilon H^3=0                                  (8)
```

is smooth: its three partials are `3aP^2`, `2bTH`, and
`bT^2-3epsilon H^2`, and their common vanishing gives no projective point.
After scaling over `C`, `(8)` is `Y^2=X^3+1`. Thus each `E_epsilon` is a
smooth `j=0` elliptic curve. Independently, its Newton polygon

```text
conv((0,0),(0,3),(2,2))
```

has `(2Area,boundary,interior)=(6,6,1)`.

## 3. Complete contact and genus ledger

At top infinity put

```text
z=1/S,                         r=P/S^2.                    (9)
```

The regular closure factors as

```text
Ctilde=z^12-r^4(ar+b)^2=H_+ H_-,
H_epsilon=z^6-epsilon r^2(ar+b).                          (10)
```

The two elliptics meet at

```text
z=0,                         r=r_0=-b/a.                  (11)
```

Since

```text
d/dr [r^2(ar+b)] at r_0 = b^2/a !=0,                     (12)
```

the quantity `w=r^2(ar+b)` is a local coordinate and their union is

```text
z^12-w^2=0.                                               (13)
```

Hence `(11)` is a two-branch `A_11` contact of intersection length six. The
apparent common point `z=r=0` in the coarse chart does not add an
intersection: in the toric refinement `r=z^3x`, the two strict transforms
meet `z=0` at the disjoint equations

```text
1-bx^2=0,                         1+bx^2=0.               (14)
```

Completeness also follows from the genus ledger. THM-4230's central polygon
has arithmetic genus seven, while

```text
g(E_+)+g(E_-)+I(E_+,E_-)-1=1+1+6-1=7.                    (15)
```

There is no room for another singularity or intersection on the normalized
toric model.

The rational component is `R:r=1`. Since `a+b!=0` under `(2)`, its
intersections with `E_epsilon` are the six distinct roots

```text
z^6=epsilon(a+b).                                        (16)
```

They are transverse and the two six-point sets are disjoint. The complete
special-fibre ledger is therefore

```text
sum genera = 0+1+1=2,
sum intersection lengths = 6+6+6=18,
p_a = 2+18-(3-1)=18,                                    (17)
```

exactly THM-4230's global polygon genus. Thus the only nontransverse point is
the single `A_11` contact `(11)`. Transverse node smoothings and their fixed
resolutions contribute only rational chains. This is also the complete
boundary inventory from THM-4230's four edge schemes: `U`, `Z`, and `Lambda`
remain nonzero, while `D=0` creates precisely the quadratic top-edge double
root `(11)`. No other edge root or component intersection collides.

## 4. Both central elliptics are Keller-constant

There are two independent proofs.

First use the good differential. THM-4103/4297 gives the coefficient-
independent identity

```text
phi^*eta_0=sigma^9 dS/G_P,                               (18)
```

where `eta_0` is nowhere zero on the good elliptic target. At the generic
point of either `E_epsilon`, the special fibre is `R(1-L)(1+L)` and

```text
G_P=unit * R * (partial_P(L-epsilon)) * (L+epsilon).     (19)
```

None of the three displayed factors vanishes identically on `E_epsilon`.
The exact certificate checks

```text
gcd(L-epsilon,L_P)=gcd(L-epsilon,S^2-P)=1.               (20)
```

Thus `G_P` is a unit at the component generic point, `dS` is a relative
generator, and `(18)` has exact vertical order nine. The specialized map
pulls back `eta_0` to zero and is constant in characteristic zero.

For an independent character proof, THM-4290's deck generator is

```text
tau:(S,P)->(xi S,xi^2P),                    xi^12=1.     (21)
```

It sends `L` to `-L` and exchanges the two elliptics. Its square preserves
each and, in `(P,T)` coordinates, acts by

```text
tau^2:(P,T)->(omega P,-T)=[-omega],          omega=xi^4. (22)
```

The squared target action is `[omega^2]`. Equivariant graph specialization
therefore gives, for a component map `m`,

```text
m o [-omega]=[omega^2] o m.                              (23)
```

Pulling back the target invariant differential through a hypothetical
nonconstant `m` would equate the two characters `-omega` and `omega^2`.
Their difference is the unit

```text
(-omega)-omega^2=1                                      (24)
```

because `1+omega+omega^2=0`. Hence the pullback differential is zero and
`m` is constant. This proof does not use an integral Hom-lattice exhaustion.

## 5. Parameterized Morse reduction of the only collision

The exact local source equation used in THM-4297 is

```text
F=(r-1)(Hhat(r,t)-z^12)-t^12/2,               t=sigma z. (25)
```

At `t=0`, its central part is

```text
Hhat(r,0)=r^4(ar+b)^2.                                  (26)
```

The hypothesis `Lambda!=0` is exactly `r_0!=1`, so `r-1` is a unit at the
contact. Divide the proved exceptional multiplicity and write

```text
F=(r-1)J,
J=A(r,t)-z^12,
A(r,t)=Hhat(r,t)-t^12/[2(r-1)].                         (27)
```

At `(r,t)=(r_0,0)`,

```text
A=0,                 A_r=0,                 A_rr=2a^2r_0^4!=0. (28)
```

The parameterized Morse reduction here is elementary and exact. The formal
implicit-function theorem gives a unique critical series `rho(t)` with
`A_r(rho(t),t)=0`. Put `kappa(t)=A(rho(t),t)`. After translating
`r=rho(t)+u`, the absence of constant and linear terms gives

```text
A-kappa=u^2 B(u,t),                         B(0,0)!=0.   (29)
```

Over `C[[u,t]]`, the unit `B` has a square root. The coordinate
`w=u sqrt(B)` therefore turns the complete collision germ into

```text
w^2=z^12-kappa(t),                         t=sigma z.    (30)
```

No lower-row coefficient has been discarded: all of it is retained in the
single critical-value series `kappa`.

## 6. Every collision tail has constant Keller map

Let `gamma=v(w)`. Every Newton face of `(30)` is obtained by requiring the
minimum of

```text
2gamma,                         12beta,
m(s+beta)                                                   (31)
```

to occur at least twice, where the third entry is omitted if `kappa=0`.
If `kappa=0`, or if `m=ord_t(kappa)>=12`, one always has
`m(s+beta)>12beta`; hence only the binomial face `w^2=z^12` remains. In fact

```text
z^12-kappa(sigma z)=z^12*unit
```

with `unit=1+O(sigma,z)`, so a Hensel square root gives two smooth formal
branches. Its fixed `A_11` resolution pieces are rational, and the form order
on them is at least `9s+5beta>0`.

Now suppose

```text
kappa(t)=c t^m+O(t^(m+1)),                 c!=0,
1<=m<12.                                                   (32)
```

For a divisorial valuation centered at the collision, put

```text
s=v(sigma)>0,                         beta=v(z)>0.        (33)
```

There are three face types.

1. If `2gamma=12beta<m(s+beta)`, the face is `w^2=z^12`.
   Its normalization is binomial/rational and its good-form order is
   `9s+11beta-6beta=9s+5beta>0`.
2. If `2gamma=m(s+beta)<12beta`, the face is `w^2=c t^m`.
   Its normalization is a monomial rational curve. Since
   `gamma=m(s+beta)/2<6beta`, its good-form order is strictly larger than
   `9s+5beta`.
3. The only possibly positive-genus face is the three-term equality below.

Terms after the first nonzero `kappa` coefficient have strictly larger
weight. The unique three-term compact face satisfies

```text
12 beta=m(s+beta),
(s,beta)=((12-m)/g,m/g),                 g=gcd(12-m,m).   (34)
```

After removing the square factor at `X=0`, its normalization is

```text
Y^2=X^(12-m)-c                       if m is even,
Y^2=X(X^(12-m)-c)                    if m is odd.         (35)
```

These polynomials are squarefree. The tail genus is

```text
g_m=floor((11-m)/2).                                      (36)
```

Thus `m=1` really can create genus five; the tail inventory is not being
replaced by a rationality assumption.

On the face `(34)`, one has `v(w)=6beta`. Equation `(25)` and `(27)` give,
up to a unit,

```text
phi^*eta_0=-sigma^9 z^10 dz/F_r,
F_r=unit*w                         on J=0.               (37)
```

At the generic point of the normalized tail, `v(d_rel z)>=beta`, so

```text
v(phi^*eta_0)>=9s+11beta-6beta=9s+5beta>0.               (38)
```

The restriction of the good target differential is zero, hence every
positive-genus tail map is constant. The apparent fourth equality
`12beta=m(s+beta)<2gamma` is a point, not another component: on the normalized
curve `(35)`, every nonzero root of `X^(12-m)-c` is a simple branch point.
The point `X=0` is also smooth after the displayed square removal. Thus later
point blowups add only rational components. This Newton-nondegenerate
resolution exhausts every divisorial component above `(11)`.

## 7. A genuine genus-five hostile

Take

```text
(U,W,Z)=(1,2,1),                         Lambda=4,         (39)
alpha_11=0,             beta_11=1,
upsilon_5=xi_10=0.
```

Then `a=b=1`, `r_0=-1`, and the beginning of `A` is

```text
A=r^4(r+1)^2+t r^4+O(t^3).                              (40)
```

The critical series begins `rho=-1+2t+O(t^2)`, and exact substitution gives

```text
kappa(t)=t-4t^2+O(t^3).                                  (41)
```

Thus `m=1`, the compact tail has genus five, and the primitive weights are
`(s,beta)=(11,1)`. Its good-differential order is

```text
9*11+5*1=104>0.                                          (42)
```

This is the positive control: a maximal-genus collision tail exists, but it
cannot carry the Keller map.

## 8. Component extinction and degree contradiction

After finite base change, normalization, and resolution of the exact graph,
every special-fibre component is one of:

1. `E_+` or `E_-`, constant by Section 4;
2. the rational component `R`;
3. a rational chain over one of the twelve transverse contacts;
4. a fixed toric or ordinary point-blowup component, hence rational; or
5. a component over the `A_11` contact, constant by Section 6.

Thus every component map is constant. With the actual pullback of the
relative degree-one origin bundle, THM-4297's proper-flat identity gives

```text
deg(M_generic)=sum_i multiplicity_i * deg(M|X_i)=0.       (43)
```

A nonautomorphic Keller response in the inherited seam is nonconstant and
has positive degree. This contradiction proves the theorem.

## 9. Exact failure boundary: the cubic corner

If `D=Lambda=0` with `UZ!=0`, then `b=-a` and `r_0=1`. Division by `r-1`
in `(27)` is illegal. Put `q=r-1`. The three smooth special branches
`R,E_+,E_-` meet at one point, and every pair has intersection length six:

```text
H_epsilon=z^6-epsilon a(1+q)^2q.                         (44)
```

The complete top contribution to `(25)` is

```text
F_top=a^2q^3(1+q)^4-qz^12-t^12/2.                        (45)
```

Its first face is genuinely cubic,

```text
a^2q^3-qz^12-t^12/2,                                    (46)
```

and the exact correction is

```text
a^2q^4(q+2)(q^2+2q+2).                                  (47)
```

This realizes THM-4297's generated `D=0` task and explains the failure of
its quadratic ladder. The canonical hostile is `(U,W,Z)=(1,-2,1)`. No claim
about extinction at this cubic corner is made here.

If `U=0` or `Z=0`, then `D=0` forces `W=0` and a Newton-hull endpoint is
lost. The main face splits into four or six rational factors respectively,
but the present component inventory cannot be transported across that hull
change. Those endpoint walls remain separate open tasks.

Finally, THM-4298's source-normal transform observes all three top
coefficients if the relevant rows exist. Current proved source-normal canon
does not construct absolute `G` rows through `6,7,8`; visibility is not seam
entry, and this theorem does not repair that missing interface.

## 10. Reproduction and scope

Run

```bash
python3 -B 04-computation/jc23_d_zero_square_face_elliptic_extinction_thm4299.py
python3 -B -O 04-computation/jc23_d_zero_square_face_elliptic_extinction_thm4299.py
python3 -B 04-computation/jc23_d_zero_square_face_elliptic_extinction_independent_audit_thm4299.py
python3 -B -O 04-computation/jc23_d_zero_square_face_elliptic_extinction_independent_audit_thm4299.py
```

The theorem closes exactly `(2)` inside the inherited exact-`M=12` reduced
seam. It does not close `D=Lambda=0`, `U=0`, `Z=0`, prove entry into the
exact-`M=12` seam from arbitrary reduced `(2,3)` data, cross to another cell,
or prove `JC(2)` or `DC(2)`.

**QED.**
