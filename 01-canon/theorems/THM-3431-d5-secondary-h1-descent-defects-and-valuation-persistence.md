---
id: THM-3431
title: "D5 secondary H1 descent defects and valuation persistence"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  The degree-13 LRC
  trivialization induces an explicit
  isomorphism H1(C7;F13) -> H1(C13_deck;F13), sending chart holonomy to the
  deck defect of the cover primitive.  Every selected one-root JC observer
  has a canonical cyclic embedding into H1_(lambda)(K[lambda]) sending it to
  [lambda^(-q_sigma)].  A lawful common forgetting constructed here is the
  lossy valuation-indexed DeathBar record: all additive cross-maps vanish in both
  directions, even when bar lengths agree.  This constructs no semantic
  current, polynomial mate, LRC(14), or JC(2) result.
source: codex2-d5-secondary-defects-2026-08-15
audit: independent graph/group cohomology, cover primitive, action convention, JC cyclic injection and annihilator, affine normalization, both-direction additive no-go, DeathBar scope, q=1 hostile, normal/-O/stored replay, and hash audit CLEAN after removing unproved uniqueness/minimality language
depends_on:
  - THM-2542-seven-chart-cech-holonomy-and-c91-arrival-obstruction
  - THM-3354-inequivalent-h1-carriers-and-typed-obstruction-cospan
  - THM-3422-one-root-nonlinear-integral-hamiltonian-response
  - THM-3427-all-sector-constant-observer-rigidity-and-polynomial-presentation
related:
  - THM-3412-hamiltonian-principal-part-differential-and-prufer-torsion-arms
  - HYP-9031-d5-h1-dictionary-lrc-word-current-vs-jc-flux
script: 04-computation/d5_secondary_h1_descent_defects_thm3431.py
output: 05-knowledge/results/d5_secondary_h1_descent_defects_thm3431.out
script_sha256: 9bcd3f10c1741b7328905b31d8401a205a6916b2a0767c3860a8f5a2b45e0ee7
output_sha256: 41375ee1c22b2ed1c39adf5a6a985896f7fbd976980cb1753d93b69fa39f5b52
semantic_sha256: 79f37de16658ac1e755e0acb5640507dba716756d9b3136f378004a4d29bcc21
hash_basis: LF-normalized bytes
---

# THM-3431 -- D5 secondary H1 descent defects and valuation persistence

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. What the repaired D5 construction is

THM-3354 proves that the original direct coefficient-induced D5 map is
impossible: the LRC class and the JC response classes live on different
sites, with incompatible coefficients and target predicates.  It leaves open
the introduction of a new correspondence space or coefficient-changing
sidecar.

The construction supplies a pair of explicit **secondary trivialization
defects**:

```text
LRC:
  base graph H1 class
    -> deck-group H1 class of its degree-13 cover primitive;

JC:
  distinguished vertical-torsion observer
    -> local-cohomology H1 principal part at its support divisor. (1)
```

Both arrows in `(1)` are explicit and preserve the exact death mechanism.
There is still no arrow between their targets.  Forgetting each target to a
valuation-indexed death interval gives a useful common invariant, but destroys
the class and cannot transport a theorem.

## 2. The LRC deck-defect isomorphism

Orient the seven-cycle `C_7`, let

```text
pi:C_91 -> C_7                                             (2)
```

be its connected cyclic cover of degree thirteen, and write its deck group as
`Delta=C_13=<tau>`, where `tau(j)=j+7` on vertices `j modulo 91`.
Coefficients are the trivial `Delta`-module `F_13`.

For any one-cocycle `g` on `C_7`, the pullback `pi^*g` has total holonomy
thirteen times that of `g`, hence zero in `F_13`.  Choose a zero-cochain
`h` on `C_91` such that

```text
delta h=pi^*g.                                             (3)
```

The primitive is unique up to a constant.  Since `tau` commutes with
`delta`,

We use the forward pullback convention

```text
((tau^m)^*h)(j)=h(tau^m(j))=h(j+7m).
```

With the inverse left-action convention, `chi_g` and `T_L` are replaced by
their negatives; the isomorphism and death law are unchanged.  Under the
forward convention,

```text
chi_g(tau^m)=tau^(m)*h-h                                  (4)
```

is a constant zero-cochain.  It is independent of the choice of `h`.
Equation `(4)` is additive in `m`, so it is a group one-cocycle

```text
chi_g in Z^1(Delta;F_13)=Hom(C_13,F_13).                  (5)
```

It also depends only on the graph-cohomology class of `g`.  Indeed, replacing
`g` by `g+delta f` replaces `h` by `h+pi^*f`, and the descended term
`pi^*f` is deck invariant.

Let

```text
s(g)=sum_(edges of C_7) g.                                (6)
```

Integrating `(3)` along the seven-edge path from `j` to `j+7` gives

```text
chi_g(tau)=s(g).                                          (7)
```

Both `H^1_graph(C_7;F_13)` and
`H^1_group(C_13;F_13)=Hom(C_13,F_13)` are one-dimensional over `F_13`,
and `(6)` identifies the first with `F_13`.  Therefore

```text
T_L:H^1_graph(C_7;F_13) -> H^1_group(C_13;F_13),
T_L([g])=[chi_g],                                        (8)
```

is an isomorphism.

For THM-2542's constant chart cochain `g_k=a!=0`,

```text
s(g)=7a,
h(j)=ja,
chi_g(tau)=h(j+7)-h(j)=7a!=0.                            (9)
```

Thus the old C91 carrier has a genuine secondary deck-group `H^1` class,
while its coefficient object remains `F_13`.  It is not a primitive
`Z/91` class.

More generally, the degree-`k` cyclic pullback is exact precisely when

```text
k s(g)=0
iff 13 divides k                                         (10)
```

for nonzero `[g]`.  Hence degree thirteen is the minimal death operation.

## 3. The JC local-cohomology class

Let `K` have characteristic zero and consider THM-3427's selected one-root
observer for

```text
P=ax+b+c(x-alpha)^e z^d,
lambda=(P-(a alpha+b))/a,
A=K[lambda].                                             (11)
```

For `1<=sigma<=d`, assume

```text
e>1,        d divides sigma(e-1),
q_sigma=sigma(e-1)/d.                                    (12)
```

Put `theta_sigma=[z^(sigma-1)]`.  THM-3422/3427 prove

```text
theta_sigma!=0,
Ann_A(theta_sigma)=(lambda^q_sigma).                      (13)
```

The first local-cohomology module at the vertical divisor is computed by the
one-element Cech complex:

```text
H^1_((lambda))(A)=A[lambda^(-1)]/A.                       (14)
```

There is a canonical injective `A`-linear map on the distinguished cyclic
submodule

```text
iota_sigma:A theta_sigma -> H^1_((lambda))(A),
iota_sigma(f theta_sigma)=[f lambda^(-q_sigma)].          (15)
```

To check well-definedness, `f theta_sigma=f' theta_sigma` exactly when
`f-f'` is divisible by `lambda^q_sigma`, in which case the two Laurent
fractions in `(15)` differ by an element of `A`.  Conversely,
`f lambda^(-q_sigma)` lies in `A` exactly when `lambda^q_sigma|f`, so the
kernel is zero.  In particular,

```text
xi_sigma=iota_sigma(theta_sigma)=[lambda^(-q_sigma)],
Ann_A(xi_sigma)=(lambda^q_sigma).                         (16)
```

Thus `xi_sigma` is an actual local-cohomology `H^1` class retaining exactly
the vertical torsion that generic fibre localization destroys.  The map
`(15)` is canonical because both the generator `theta_sigma` and the
normalized parameter `lambda` are specified.  It does not choose or identify
a splitting of the full sector.

For the wrap character `sigma=d`, the depth is `q_d=e-1`.  THM-3427's
canonical wrap response gives a second reading of the same integer:

```text
H_d(alpha)/[d S'(alpha)]=e-1=q_d                         (17)
```

in the one-root case.  In multiroot cases the response polynomial still
recovers each multiplicity excess, but the displayed constant observers are
nontorsion; `(15)` must not be extrapolated to them.

## 4. Exact additive no-go between the secondary classes

Although `(8)` and `(15)` now place both defects in honestly named
`H^1` objects, the coefficient mismatch remains exact.

The additive group of `H^1_group(C_13;F_13)` has exponent thirteen.  The
additive group of `H^1_((lambda))(A)` is a `K`-vector space and is
torsion-free.  Therefore every additive homomorphism

```text
H^1_group(C_13;F_13) -> H^1_((lambda))(A)                (18)
```

is zero.  Conversely, multiplication by thirteen is surjective on every
characteristic-zero `K`-vector space.  If `f` maps that vector space to an
exponent-thirteen group, then

```text
f(y)=f(13(y/13))=13f(y/13)=0.                            (19)
```

So every additive map in the reverse direction is also zero.  In particular,
`T_L([g])` cannot map to `xi_sigma`, or conversely, through their additive
coefficient structures.  A common ambient direct sum would merely juxtapose
the classes and add no interaction.

## 5. The lawful common forgetting: valuation death bars

Cover degrees form the multiplicative monoid `N_(>0)`.  Pullback of a
nonzero LRC class by degree `k` multiplies its holonomy by `k`; composition
multiplies degrees.  The valuation

```text
v_13(k ell)=v_13(k)+v_13(ell)                            (20)
```

turns this action into an additive death scale.  The class survives at
`v_13=0` and dies at `v_13>=1`.  Along the degree-`13^n` tower its death
barcode is

```text
Bar_L=[0,1).                                             (21)
```

Ordinary numerical degree is not the scale: degree `13` dies while degree
`14` survives.  The annihilator is the principal ideal `13 N_(>0)` in the
multiplicative degree monoid.

On the JC side, multiplication by `lambda^n` composes by addition of
exponents.  Equations `(13)` and `(16)` give

```text
lambda^n theta_sigma!=0 iff 0<=n<q_sigma,
Bar_J=[0,q_sigma),                                      (22)
```

with ring annihilator `(lambda^q_sigma)`.

Define `DeathBar` to retain only the additive index set, the first death
time, and whether death is exact.  Then there are tautological forgetting
maps

```text
(T_L([g]),degree action) -> DeathBar <- (xi_sigma,lambda action). (23)
```

This is the repaired D5 **cospan**.  It preserves composition-to-addition and
minimal death, but destroys the site, coefficient object, class order,
primitive, target predicate, and realization sidecars.  It is therefore an
experiment scheduler, not a cohomology map.

The sharp hostile is `q_sigma=1`: both bars in `(21)--(22)` then have length
one, but `(18)--(19)` still force every additive cross-map to vanish.  Equal
barcodes do not identify classes.

## 6. Exact replay and stopping boundary

The standard-library companion checks all twelve nonzero LRC chart steps,
`1,092` cover-primitive edge cells, `14,196` deck-defect cells, `2,028`
group-cocycle identities, degree death through `200`, and `2,500`
composition/valuation identities.  On the JC side it enumerates `562`
selected observer profiles over

```text
2<=d<=12,       2<=e<=25,       1<=sigma<=d,             (24)
```

checking `7,058` annihilator states and `5,372` cyclic-embedding basis
vectors.  It freezes the degree-14 scale hostile, depths one and two, the
additive-order no-go, and the equal-length barcode hostile.  Normal and
optimized outputs agree exactly with the LF-normalized stored transcript.

The theorem constructs no LRC semantic vertical path, physical current,
coefficient-changing correspondence, polynomial mate, Keller inverse,
LRC(14), or JC(2) conclusion.  It does not identify group, graph, local,
de Rham, or etale cohomology.  **QED.**
