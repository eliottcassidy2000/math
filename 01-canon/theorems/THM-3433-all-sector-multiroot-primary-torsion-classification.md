---
id: THM-3433
title: "All-sector multiroot primary-torsion classification"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Let K have characteristic
  zero, d>=2, a!=0, and P=ax+b+g(x)z^d with nonconstant g.  Over a splitting
  field, root alpha_i of multiplicity e_i>1 contributes a Pruefer arm to
  target sector sigma-1 exactly when d divides sigma(e_i-1) and d divides
  sigma e_j at every other root.  Its exact stage-k thickness is
  k(e_i-1)+sigma(e_i-1)/d.  Every primary torsion element lies in one of
  these arms.  A nonwrap sector has at most one arm; the wrap sector has one
  per repeated root.  Generic sector rank remains deg(rad(g)).  Split
  coordinates are noncanonical, no general torsion-free complement is
  claimed, and this gives no polynomial mate, new Keller case, or open
  JC(2) conclusion.
source: root-jc-all-sector-primary-torsion-2026-08-15
audit: independent selected-root congruence, global-power dying quotient, injective alternative, T-action intertwining, algebraic-closure primary exhaustion and norm/flat descent, nonwrap uniqueness, strict filtration, thickness/gcd character count, wrap regression, exact-artifact, normal/-O, hash, AST/security, and documentation audit CLEAN
depends_on:
  - THM-3418-one-monomial-nonlinear-fiber-keller-classification
  - THM-3422-one-root-nonlinear-integral-hamiltonian-response
  - THM-3419-generic-kummer-response-regular-sector-rank
related:
  - THM-3430-nonlinear-wrap-linearization-and-canonical-vertical-torsion
  - THM-3431-d5-secondary-h1-descent-defects-and-valuation-persistence
  - THM-3427-all-sector-constant-observer-rigidity-and-polynomial-presentation
  - THM-3424-nonlinear-multiroot-unit-observer-rigidity
script: 04-computation/jc_all_sector_multiroot_primary_torsion_thm3433.py
output: 05-knowledge/results/jc_all_sector_multiroot_primary_torsion_thm3433.out
script_sha256: 7bb4db4f6d67436a2739dd79d0579243b8d0f5474a8f739a82b281014c114d8d
output_sha256: 203cd02b837651a608b7834cb745abdb187d0ed38c34df4af2a1b7e0b9088453
hash_basis: LF-normalized bytes
---

# THM-3433 -- all-sector multiroot primary-torsion classification

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Statement and connection contract

Let `K` be a field of characteristic zero, let `d>=2`, and fix

```text
P=ax+b+g(x)z^d in K[x,z],              a in K*,           (1)
```

with `g` nonconstant.  Put `D=Jac(P,-)` and

```text
C=K[x,z]/D(K[x,z])=direct_sum_(sigma=1)^d C_(sigma-1),    (2)
```

where `sigma=d` denotes the wrap target sector.  Let `K'/K` be a finite
normal splitting field and write

```text
g=gamma product_(i=1)^N (x-alpha_i)^e_i,
beta_i=a alpha_i+b,                 e_i>=1,
N=deg(rad(g)).                                             (3)
```

Identify `K'[P]` with `K'[T]`.  For `1<=sigma<=d`, define the selected-root
set

```text
Sigma_sigma={i:
  e_i>1,
  d divides sigma(e_i-1),
  d divides sigma e_j for every j!=i}.                    (4)
```

Then the complete geometric primary torsion is

```text
tors_(K'[T])(C_(sigma-1) tensor_K K')
 ~= direct_sum_(i in Sigma_sigma)
    Pr_(T-beta_i),                                        (5)

Pr_(T-beta_i)
 =K'[T,(T-beta_i)^(-1)]/K'[T].                            (6)
```

Let `F_(sigma,k)` be the image of the depth-`k` THM-3418 stage, whose target
fiber exponent is `sigma-1+kd`.  For a selected root put

```text
t_(i,sigma,k)
 =k(e_i-1)+sigma(e_i-1)/d.                                (7)
```

The filtration on every arm is exact:

```text
F_(sigma,k) intersect tors_(T-beta_i)(C_(sigma-1) tensor K')
 ~=K'[T]/(T-beta_i)^(t_(i,sigma,k)).                      (8)
```

In particular:

- if `sigma<d`, then `|Sigma_sigma|<=1`;
- `Sigma_d={i:e_i>1}` and `(7)` becomes `(k+1)(e_i-1)`;
- a repeated root `alpha_i` occurs in exactly

```text
gcd(d,e_i-1,{e_j:j!=i})                                  (9)
```

  of the `d` target characters;
- every sector still has generic `K(P)`-rank `N` by THM-3419.  Multiplicity
  changes the torsion characters, intercepts, and slopes, not generic rank.

There is one exceptional proof regime, but no exception to `(5)`.  If
`sigma<d` and

```text
d divides sigma e_j for every j,                          (10)
```

all transition maps have a one-dimensional global-power kernel.  Its dying
subdiagram can be quotiented out, and one obtains the stronger full module
formula

```text
C_(sigma-1) tensor_K K'
 ~=direct_sum_(i=1)^N K'[T,(T-beta_i)^(-1)].              (11)
```

Thus this regime is torsion-free.  Outside `(10)`, every nonwrap transition
is injective; formula `(5)` then follows from the unique polynomial
horizontal solution and primary exhaustion below.  No full decomposition of
the torsion-free quotient is asserted in that injective regime.

Descent is sharp.  A selected nonwrap root is unique and hence fixed by the
Galois group, so it is `K`-rational.  Therefore the nonwrap torsion over `K`
is either zero or one `Pr_(T-beta)` with `beta in K`.  No nonsplit orbit can
hide a nonwrap arm.  In the wrap sector, conjugate arms package by irreducible
factors.  If

```text
S=rad(g),                         c=g/S,                   (12)
```

then, noncanonically as a filtered `K[T]`-module,

```text
tors_(K[T])(C_(d-1)) ~= K[x,c^(-1)]/K[x],                (13)
```

where `T` acts by `ax+b`.  The torsion submodule and filtration are
canonical; split primary coordinates and the isomorphism in `(13)` are not.

| item | exact content |
|---|---|
| source | every THM-3418 character-sector directed system |
| target | its full geometric primary torsion and exact finite filtration |
| map | multiplication by the algebraized other-root gauge in the injective regime; quotient by the global-power dying subdiagram in regime `(10)` |
| preserved | `K[T]` action, root support, transition depth, annihilator order, and character label |
| forgotten | the unsplit primary coordinates and, outside `(10)`, the structure of the torsion-free quotient |
| required sidecar | transition-kernel dichotomy plus the valuation-strictness test at every other root |
| sharp hostiles | `(4,2;3,1)` is locally resonant but globally blocked; `(2,1;3,3)` has two local resonances but no arm |

This is a Hamiltonian-cokernel classification.  It is not a coefficient map
from LRC to JC, not a polynomial mate, and not a new Keller stratum.

## 2. Nonwrap transition dichotomy and the dying quotient

Fix `1<=sigma<d` and work over `K'`.  THM-3418 presents the sector as

```text
C_(sigma-1) tensor K'
 ~=colim_(k>=0)(K'[x],L_(sigma+kd)),                      (14)

L_n(q)=(d/(an))gq'-(1/a)g'q.                             (15)
```

Over `K'(x)`, a nonzero transition-kernel element must satisfy

```text
q'/q=(n/d)g'/g,
q=C product_j (x-alpha_j)^(n e_j/d).                     (16)
```

when the displayed root orders are integral.  The leading scalar of `g` is
irrelevant and no `d`th root of it is being adjoined.  For
`n=sigma+kd`, `(16)` is a polynomial exactly when

```text
d divides sigma e_j for every j.                          (17)
```

Thus either every transition has a one-dimensional kernel, namely `(10)`,
or every transition is injective.  There is no stage-dependent third case.

Assume `(10)` and put

```text
delta_j=sigma e_j/d,
H_k=product_j (x-alpha_j)^(k e_j+delta_j).                (18)
```

Then `ker(L_(sigma+kd))=K'H_k`, and the subspaces

```text
B_k=H_k K'[x]                                             (19)
```

form a subdiagram.  If `p in K'[x]`, direct cancellation in `(15)` gives

```text
L_(sigma+kd)(H_k p)
 =gamma d/[a(sigma+kd)] H_(k+1)p'.                       (20)
```

Repeated differentiation eventually kills every polynomial, so
`colim B_k=0`.  Exactness of filtered colimits therefore identifies `(14)`
with the colimit of the quotient diagram `K'[x]/(H_k)`.

CRT now becomes legitimate.  At root `i`, put

```text
G_(i,k)=product_(j!=i)(x-alpha_j)^(k e_j+delta_j).        (21)
```

Modulo `H_k`, multiplication by `G_(i,k)` identifies the `i`th CRT summand
with

```text
K'[y_i]/y_i^(k e_i+delta_i),              y_i=x-alpha_i. (22)
```

Moreover, for every polynomial `p`,

```text
L_(sigma+kd)^g(G_(i,k)p)
 =G_(i,k+1)L_(sigma+kd)^(gamma y_i^e_i)(p).              (23)
```

The dying subdiagram `(19)` is `T`-stable, so the quotient inherits the
Hamiltonian action.  On the `i`th CRT channel that action also intertwines:

```text
(T-beta_i)[G_(i,k)p]_k
 =[a y_i G_(i,k)p]_k+[gG_(i,k)p]_(k+1)
 =[G_(i,k)a y_i p]_k
  +[G_(i,k+1)gamma y_i^e_i p]_(k+1).                    (23a)
```

The quotient in `(22)` is exactly the one-root nonwrap diagram modulo its
own dying global-power subdiagram, so its colimit is the original one-root
sector.  Since `sigma<d`, the two congruences

```text
d|sigma e_i,                   d|sigma(e_i-1)             (24)
```

cannot both hold: their difference would give `d|sigma`.  THM-3422 therefore
identifies every local colimit with the unselected Laurent line
`K'[T,(T-beta_i)^(-1)]`.  Summing the CRT channels proves `(11)`, and hence
`(5)` in the global-power regime.

## 3. The unique horizontal solution in the injective regime

Assume `(10)` fails, so every map in `(14)` is injective.  Extend faithfully
once more to an algebraic closure `Omega/K'`; all stages, transitions, and
Hamiltonian actions commute with this scalar extension.  Let `rho in Omega`,
put

```text
xi=(rho-b)/a,                  n=sigma+kd,                (25)
```

and take a nonzero class `v=[p]_k`.  THM-3418's induced Hamiltonian action is

```text
T[p]_k=[(ax+b)p]_k+[gp]_(k+1).                            (26)
```

Moving the first term to stage `k+1` shows that `(T-rho)v=0` if and only if

```text
M_(n,xi)(p)=0,                                             (27)

M_(n,xi)(p)
 =d g(x-xi)p'+(n+d)gp-n(x-xi)g'p.                        (28)
```

Here injectivity is load bearing: zero in the colimit is zero at stage
`k+1`.  The first-order equation `(28)` has the unique formal logarithmic
solution, up to a scalar,

```text
p=product_j (x-alpha_j)^(n e_j/d)
  (x-xi)^(-(n+d)/d).                                     (29)
```

Here `(29)` records root orders; it is rational only when all combined
orders are integral, and again suppresses the irrelevant leading scalar.
Read its logarithmic residues.  If `xi` is not a root of `g`, the residue
`-(n+d)/d` at `xi` is nonintegral because `sigma<d`, so `(29)` is not
rational.  If `xi=alpha_i`, the root exponents are

```text
at alpha_j, j!=i:   k e_j+sigma e_j/d,
at alpha_i:         k(e_i-1)+sigma(e_i-1)/d-1.            (30)
```

They are nonnegative integers exactly when `i in Sigma_sigma`.  The
condition `e_i>1` is essential: for `e_i=1`, the last exponent is `-1` even
though the divisibility congruence is vacuous.  For selected `i`, a polynomial
generator is

```text
p_(i,k)=G_(i,k)
        y_i^(k(e_i-1)+sigma(e_i-1)/d-1),                  (31)
```

where `(21)` is now formed using the selected-root other-factor exponents.
The exact persistence identity is

```text
L_(sigma+kd)(p_(i,k))
 =-gamma(sigma+(k+1)d)/[a(sigma+kd)] p_(i,k+1),           (32)
```

so these stagewise kernel lines define one nonzero line in the colimit.
Consequently

```text
dim_(Omega) ker(T-rho)=
  1  if rho=beta_i for an i in Sigma_sigma,
  0  otherwise.                                          (33)
```

Finally, two distinct roots cannot both be selected when `sigma<d`.  If
`i` and `j` were selected, root `j` would give both
`d|sigma e_j` and `d|sigma(e_j-1)`, again forcing `d|sigma`.  This proves the
at-most-one-arm assertion before any primary-module classification.

## 4. Embedded arm, filtration strictness, and primary exhaustion

Fix the unique selected root `i`, if it exists, and write

```text
g=gamma y_i^e_i u_i,
h_i=product_(j!=i)(x-alpha_j)^(sigma e_j/d).              (34)
```

The congruences in `(4)` make `h_i` a polynomial.  The subspaces

```text
A_(i,k)=h_i u_i^k K'[y_i] subset K'[x]                   (35)
```

form a subdiagram, because the exact gauge identity is

```text
L_(sigma+kd)^g(h_i u_i^k p)
 =h_i u_i^(k+1)L_(sigma+kd)^(gamma y_i^e_i)(p).          (36)
```

The two terms in `(26)` also intertwine with the one-root
`T-beta_i` action.  Thus exactness of filtered colimits embeds the complete
selected one-root sector from THM-3422.  In particular it supplies one
`Pr_(T-beta_i)` arm with endpoint `(31)`.

The embedding is filtration-strict; this is what turns the one-root lower
bound into the exact thickness `(8)`.  Suppose

```text
L_(sigma+kd)(p) in A_(i,k+1).                             (37)
```

At another root `alpha_j`, set

```text
v=ord_(alpha_j)(p),
r_(j,k)=k e_j+sigma e_j/d.                                (38)
```

If `v<r_(j,k)`, the leading term of `(15)` has order
`e_j+v-1` and nonzero coefficient proportional to

```text
d v-(sigma+kd)e_j.                                        (39)
```

The coefficient in `(39)` could vanish only at `v=r_(j,k)`, while
`e_j+v-1<r_(j,k)+e_j`, the order required by `(37)`.  This is a
contradiction.  Hence `v>=r_(j,k)` for every `j!=i`, so
`p in A_(i,k)`.  The pullback equality just proved makes the induced maps
`K'[x]/A_(i,k) -> K'[x]/A_(i,k+1)` injective.  Therefore

```text
F_(sigma,k) intersect colim A_(i,*)=image(A_(i,k)),       (40)
```

and THM-3422's one-root filtration gives exactly the block in `(8)`, with
length `(7)`.

It remains to exclude extra primary torsion outside the embedded arm.  This
uses `(33)` and a short general lemma.  Let `lambda=T-beta_i`, and let
`Q_i` be the embedded Prüfer arm.  If `lambda^m v=0` with `m` minimal, then

```text
lambda^(m-1)v in ker(lambda)=Q_i[lambda].                 (41)
```

The Prüfer module is `lambda`-divisible.  Choose `q in Q_i` with
`lambda^(m-1)q=lambda^(m-1)v`; then `v-q` is killed by
`lambda^(m-1)`.  Induction gives `v in Q_i`.  If `(33)` has zero kernel,
there can be no primary torsion at that support.  Finally every torsion
element over `Omega[T]` decomposes into its linear primary parts by Bezout.
This proves after scalar extension to `Omega` that the torsion is exactly the
sum of the displayed arms.

There is no hidden descent issue here.  Every element and annihilating
polynomial over `Omega` is defined over some finite extension `L/K'`; a
nonzero norm multiple in `K'[T]` kills it.  Flatness of the corresponding
multiplication kernel shows

```text
tors_(K'[T])(C_(sigma-1) tensor K') tensor_(K') Omega
 =tors_(Omega[T])(C_(sigma-1) tensor Omega).              (41a)
```

Faithful flatness therefore descends the equality of the arm sum and the
whole torsion to `K'`.  This proves `(5)` for every injective nonwrap sector
and exhausts all primary torsion, not merely the visible endpoints.

## 5. Wrap regression and descent

The wrap proof is separate; it does not use provisional THM-3430.  THM-3418
gives

```text
C_(d-1) tensor K'
 ~=colim_(k>=0) K'[x]/g^(k+1),                            (42)
```

with `n=(k+1)d`.  At root `i`, put

```text
u_i=product_(j!=i)(x-alpha_j)^e_j.                        (43)
```

CRT and direct differentiation give the diagram-level identities

```text
K'[x]/g^(k+1)
 =direct_sum_i u_i^(k+1)K'[x]/g^(k+1),                   (44)

L_((k+1)d)^g(u_i^(k+1)p)
 =u_i^(k+2)L_((k+1)d)^(gamma y_i^e_i)(p).                (45)
```

The `T-beta_i` action intertwines as well.  Hence the wrap diagram is the
direct sum of the one-root wrap diagrams from THM-3422.  Every repeated root
contributes one Prüfer arm, every simple root contributes none, and the
stage-`k` primary length is `(k+1)(e_i-1)`.  This proves `(5),(7),(8)` for
`sigma=d` directly.

Torsion commutes with the finite normal extension `K'/K`.  One inclusion is
flatness.  Conversely, if an element after extension is killed by
`p(T) in K'[T]`, it is killed by the nonzero norm multiple
`Norm_(K'/K)(p)(T) in K[T]`; flatness of that multiplication kernel puts it
in the scalar extension of the original torsion.

For nonwrap `sigma`, the selected set is Galois-stable and has size at most
one.  A selected root is therefore fixed, so both `alpha_i` and `beta_i`
lie in `K`; the polynomial gauge `(34)` may be chosen over `K`.  This proves
the stated zero-or-one-arm descent.

For wrap, factor `c=g/rad(g)` over `K`.  The split finite window is

```text
direct_sum_(i:e_i>1)
 K'[T]/(T-beta_i)^((k+1)(e_i-1)).                         (46)
```

The elementary-divisor theorem over `K[T]` identifies its descent,
noncanonically, with `K[x]/c^(k+1)`.  Compatible primary rescalings put the
window inclusions into

```text
K[x]/c^(k+1) -> K[x]/c^(k+2),             [r] |-> [cr].  (47)
```

Taking the colimit gives `(13)`.  This also records exactly what survives
descent: the irreducible primary packets and filtration, not individual
coordinates inside a nonsplit orbit.

## 6. Multiplicity arithmetic, collisions, and boundaries

For fixed repeated root `i`, the conditions `(4)` say that `sigma mod d`
annihilates the integer vector

```text
(e_i-1,{e_j:j!=i}) mod d.                                 (48)
```

The annihilator subgroup has order `(9)`.  This proves the exact character
count.  All selected characters have slope `e_i-1`; their intercept is
`sigma(e_i-1)/d`.  Thus equal generic rank across the `d` sectors coexists
with a character-shifted torsion distribution.  Multiplicity is visible in
the integral boundary filtration even when the generic Kummer packet is the
regular `N`-packet of THM-3419.

This filtration has an exact related-only reading through THM-3431's
`DeathBar` forgetting.  At the first visible stage, selected root `i` has

```text
q_i=sigma(e_i-1)/d,             Bar_(i,sigma,0)=[0,q_i), (49)
```

and every further depth extends the death time by the persistence slope
`e_i-1`:

```text
Bar_(i,sigma,k)=[0,q_i+k(e_i-1)).                        (50)
```

For wrap, `sigma=d`, so the denominator level is `q=k+1` and the length is
`q(e_i-1)` with no indexing ambiguity.  The selected arms therefore define a
multiset of exact valuation bars.  THM-3431's two-way additive no-go remains
in force: forgetting to these bars destroys the coefficient object, site,
class, and target predicate.  No LRC-to-JC or JC-to-LRC map is constructed.

The hypotheses and hostiles are sharp:

- `(d,sigma;e_1,e_2)=(4,2;3,1)` has local resonance at root `1`, but
  `4` does not divide `2e_2`, so it has no arm;
- `(2,1;3,3)` has two local resonances, but each root supplies nontrivial
  monodromy for the other, so it has no arm rather than two colliding arms;
- `(2,1;3,2)` has exactly one selected root and one genuine nonwrap arm;
- the global-power profile `(2,1;2,4)` has transition kernels but is the
  torsion-free Laurent sum `(11)`;
- over `K=Q`, `(x^2+1)^3` is a nonsplit version of the double-local hostile,
  while `x^3(x^2+1)^2` has exactly the rational-root arm;
- because `a!=0`, the values `beta_i=a alpha_i+b` are distinct.  Primary
  collisions would become possible if this hypothesis were removed;
- if every root is simple, every selected set is empty and all sectors are
  torsion-free;
- for one root, `(4),(7),(9)` reduce exactly to THM-3422;
- if `g` is a nonzero constant, every nonwrap transition is a nonzero scalar
  times differentiation and the wrap stages vanish, so `C=0`; if `g=0`,
  then `D=a partial_z` and again `C=0`;
- characteristic zero is load bearing in the differential equations,
  repeated differentiation, norm descent, and all denominators.

MISTAKE-374 remains active: torsion is destroyed by generic localization.
Neither generic vanishing nor a finite annihilator supplies a polynomial
mate.  THM-3418 already settles the sparse Keller classification, so this
theorem gives no new case of `JC(2)` or `DC(2)`.

## 7. Exact replay

The companion checks the congruence classification on every

```text
2<=d<=12,
1<=sigma<d,
one to four roots,
1<=e_i<=6,                                                (51)
```

for `102564` nonwrap sector profiles.  It separately checks the gcd
character count for `54175` repeated-root profiles and exact thickness for
`282055` selected finite windows.  Symbolic polynomial identities cover six
global-power diagrams through four depths and six selected-root gauges
through four depths at the canonical roots `(-3,1,4,8)`.  Exact modular-minor
certificates audit `45` low-grid `(T-rho)` kernels; the explicit horizontal
polynomial supplies the matching null vector in every predicted
one-dimensional case.  The script also checks the strict valuation gate,
finite primary models, wrap regression, nonsplit examples, and support
collisions.  Reproduce with

```bash
python3 04-computation/jc_all_sector_multiroot_primary_torsion_thm3433.py
```

The finite universe is evidence for the coefficient-independent proof, not
an extrapolated cutoff.  **QED.**
