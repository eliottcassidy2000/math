---
id: THM-3271
title: "Universal quartic centered-norm packet and local singleton projector"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every monic separable quartic f over a
  characteristic-zero field, the centered norm of the three roots
  complementary to a tautological root is one canonical cubic element nu_f
  of the quartic etale algebra.  Its explicit polynomial is
  -(20X^3+15a3X^2+(18a2-3a3^2)X+27a1-9a2a3+2a3^3)/27.
  It is translation invariant and scales by u^3 under X->uX+v.  At every
  irreducible 1+3 localization, its singleton value cannot collide with a
  moving value.  In the characteristic polynomial P_nu, that value is the
  unique simple base-field root, so it is extractable without first naming
  the fixed quartic root.  If P_nu=(Z-n0)H, then
  e_fix=H(nu)/H(n0) is the canonical field-level singleton idempotent.  Thus
  the global algebra-valued packet bypasses a coherent sheet choice for the
  local THM-3230 norm, but the local idempotents do not glue in a connected
  quartic field.  Over a DVR with 3 invertible, H(n0)=P_nu'(n0) is the exact
  denominator gate for e_fix to lie in the monogenic order O[nu]; larger-
  order conductor and affine-open incidence remain separate.  No integral
  descent, cross-place scalar section, Keller cofactor, C3/S4 exclusion, or
  JC(2) theorem follows.
source: jc-global-norm-packet-2026-08-02
audit: >
  An independent audit rederived the universal cubic packet, affine cube
  covariance, the degree-three-versus-four-root separation, unique simple
  scalar root, spectral idempotent, exact DVR denominator equivalence, and
  connected-field gluing obstruction.  Fresh normal and optimized replays
  byte-match the archived transcript.
depends_on:
  - THM-3230-marked-c3-trace-centered-norm-and-terminal-prefactor-recovery
related:
  - THM-3037-cusp-braid-s4-lift-dichotomy-and-common-sheet-owner-boundary
  - THM-3064-pointed-cubic-norm-keller-decoder-and-inverse-different-boundary
  - THM-3066-k4-initial-face-product-quotient-blind-to-keller-sheetwise-cofactor
  - THM-3201-c3-local-resolvent-splitting-and-matching-newton-gate
  - THM-3270-gamma3-peripheral-cusp-fixed-sheet-transport-and-holonomy-obstruction
script: 04-computation/jc_universal_quartic_centered_norm_packet_thm3271.py
output: 05-knowledge/results/jc_universal_quartic_centered_norm_packet_thm3271.out
script_sha256: 170e221ce79a101c81f389b83b457bebf7417910a1b6c79d8ad882811aa46820
output_sha256: 141bacfe49a755e5ad2e0b4a1827234cf98a0b08f3a9dad5d8da207f9042cd28
hash_basis: LF-normalized bytes
---

# THM-3271 -- retain all four centered norms before choosing a sheet

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. The universal packet

Let `K` be a characteristic-zero field and let

```text
f(X)=X^4+a_3 X^3+a_2 X^2+a_1 X+a_0 in K[X]             (1)
```

be monic and separable.  Put

```text
A=K[X]/(f),                     xi=X mod f.              (2)
```

For a root `alpha` of `f`, divide

```text
f(X)=(X-alpha)g_alpha(X),
g_alpha=X^3+b_2X^2+b_1X+b_0.                             (3)
```

Synthetic division gives

```text
b_2=a_3+alpha,
b_1=a_2+a_3 alpha+alpha^2,
b_0=a_1+a_2 alpha+a_3 alpha^2+alpha^3.                   (4)
```

The mean of the three complementary roots is

```text
mu_alpha=-(a_3+alpha)/3.                                 (5)
```

Define their centered norm by

```text
n_alpha=-g_alpha(mu_alpha)
        =product_(beta root of g_alpha)(beta-mu_alpha).  (6)
```

Substituting `(4)--(5)` into `(6)` gives one cubic polynomial, independent
of any root choice:

```text
N_f(X)=-1/27 * (
          20X^3+15a_3X^2+(18a_2-3a_3^2)X
          +27a_1-9a_2a_3+2a_3^3).                       (7)
```

Thus

```text
nu_f=N_f(xi) in A                                        (8)
```

is a canonical **four-sheet packet**.  Under the four geometric embeddings
of `A`, its values are exactly the four numbers `n_alpha` obtained by taking
each root in turn as the fixed sheet.  It is defined over `K` before any
sheet is selected.

Its scalar characteristic polynomial is likewise canonical:

```text
P_f(Z)=Norm_(A/K)(Z-nu_f)
      =Res_X(f(X),Z-N_f(X)) in K[Z].                     (9)
```

The pair `(A,nu_f)` retains which value belongs to which tautological sheet.
The polynomial `P_f` alone retains only the multiset of four values.

## 2. Affine-coordinate invariance

Make an affine root-coordinate change

```text
alpha |-> alpha'=u alpha+v,              u in K*, v in K. (10)
```

The complementary mean changes to `mu'=u mu+v`, so every centered
difference changes by the factor `u`.  Therefore

```text
n_(alpha')'=u^3 n_alpha.                                  (11)
```

Equivalently, after the induced isomorphism of quartic algebras,

```text
nu_(f')=u^3 nu_f.                                         (12)
```

Translations disappear exactly, and rescaling contributes a cube.  Hence
the cubeclass of the local centered norm—and, after removing its valuation,
the cubeclass of its leading unit—is invariant under affine changes of the
quartic primitive coordinate.  This is the exact coordinate behavior needed
by THM-3230; no general Mobius-coordinate invariance is asserted.

## 3. Irreducible `1+3` factors are automatically separated

Now let a field extension or completion of `K` be fixed, rename it `K`, and
suppose

```text
f(X)=(X-alpha)g(X),
alpha in K,
g in K[X] irreducible of degree three.                    (13)
```

Then

```text
A=K x L,                         L=K[X]/(g),
nu_f=(n_0,n_L),                  n_0=N_f(alpha).          (14)
```

The singleton value cannot collide with a moving value.  Indeed, if `beta`
is a root of `g` and

```text
N_f(beta)=N_f(alpha),                                     (15)
```

then irreducibility forces `g` to divide the degree-at-most-three polynomial

```text
N_f(X)-N_f(alpha).                                        (16)
```

Its leading coefficient is `-20/27`, so `(16)` is nonzero and has degree
three.  It would therefore equal `c g(X)` for some `c!=0`.  Evaluating at
`alpha` gives `0=c g(alpha)`, impossible because the two factors in `(13)`
are distinct.  Thus

```text
N_f(alpha) != N_f(beta) for every root beta of g.          (17)
```

This degree-three-versus-four-roots argument is the mechanism; cyclic
ordering and explicit radicals play no role.

## 4. Noncircular extraction from the packet polynomial

Define

```text
H(Z)=Norm_(L/K)(Z-n_L).                                   (18)
```

Then `(9)` factors as

```text
P_f(Z)=(Z-n_0)H(Z),                   H(n_0)!=0.           (19)
```

Formula `(18)` is useful for proving the statement, but it is not required
as input.  The singleton value can be extracted directly from the already
global polynomial `P_f`.

Since `[L:K]=3`, the subfield `K(n_L)` has degree either one or three.

- If it has degree three, `H` is the irreducible separable minimal polynomial
  of `n_L`, so `H` has no `K`-rational root.
- If it has degree one, `n_L in K` and
  `H(Z)=(Z-n_L)^3`; equation `(17)` says `n_L!=n_0`.

Consequently

```text
n_0 is the unique simple K-rational root of P_f.           (20)
```

Thus, once the local algebra is known to have irreducible `1+3` type, the
centered singleton norm is found from `P_f` without first naming `alpha`,
factoring out `X-alpha`, or choosing a quartic-sheet label.  After `(20)` is
found, `H=P_f/(Z-n_0)` is also determined from `P_f` alone.

The corresponding spectral projector is

```text
e_0=H(nu_f)/H(n_0)
   =H(nu_f)/P_f'(n_0) in K[nu_f] subset A.                (21)
```

In the decomposition `(14)`, Cayley--Hamilton gives

```text
H(nu_f)=(H(n_0),0),
e_0=(1,0).                                                 (22)
```

Therefore the universal packet not only identifies the centered norm: at
the field level it reconstructs the singleton factor idempotent.  This is a
noncircular construction from `(A,nu_f)` and the local `1+3` type.

## 5. What globalizes and what cannot

Suppose now that `A` is the generic quartic algebra over a global function
field.  The element `nu_f` and polynomial `P_f` in `(8)--(9)` are global.
At every place where `A` has irreducible `1+3` type, equations `(20)--(22)`
canonically produce a local value `n_(0,v)` and a local idempotent `e_(0,v)`.
The local fixed sheet can vary with `v`; no comparison of sheet labels is
needed to define the packet.

The correct gluing statement is

```text
global object:       one rational function nu_f on Spec(A);
local operation:     restrict it to the degree-one prime above a 1+3 place;
base scalar section: not supplied.                        (23)
```

If `A` is a connected quartic field, it has no nontrivial global idempotent.
Hence the different `e_(0,v)` cannot be localizations of one global
`e in A`: such an `e` would split the quartic field.  Equivalently, the
simple local roots `n_(0,v)` are branches of the global spectral polynomial
`P_f`, but need not be values of one rational root of `P_f` over the base.

Thus THM-3270's holonomy obstruction is not contradicted.  It is bypassed by
retaining all four components in the algebra-valued packet, not solved by
constructing one coherent sheet.  Passing from `(A,nu_f)` to the scalar
polynomial `P_f` forgets sheet-value incidence and can reintroduce ambiguity.

## 6. Exact integral denominator gate

Let `O` be a discrete valuation ring with fraction field `K`, assume `3` is
a unit, and suppose `f in O[X]` is monic.  Then `(7)` gives

```text
nu_f in O[xi].                                             (24)
```

At an irreducible `1+3` generic factor, integrality makes `n_0` integral and
synthetic division of the monic `P_f` gives `H in O[Z]`.  Put

```text
Delta_nu=H(n_0)=P_f'(n_0).                                (25)
```

Then

```text
e_0 belongs to the monogenic order O[nu_f]
  iff Delta_nu is a unit of O.                             (26)
```

For the forward algebra, `(21)` proves sufficiency.  For necessity, let `M`
be the minimal polynomial of `n_L`.  It has degree one or three, and
`H=M^(3/deg M)`.  A polynomial in `O[nu_f]` can equal `(1,0)` exactly when
the ideals

```text
(Z-n_0), (M(Z)) in O[Z]                                   (27)
```

are comaximal.  Their resultant is `M(n_0)`, which is a unit exactly when
`H(n_0)` is a unit.  This proves `(26)`.

The scope of `(26)` is the monogenic packet order.  If `Delta_nu` is not a
unit, a larger quartic order or its conductor may still contain `e_0`.
Conversely, even an integral singleton idempotent in the finite normalization
does not put its boundary trace into the affine Zariski-main open; THM-3037's
owner hostile remains.  Therefore

```text
field projector
 -> monogenic integrality (Delta_nu unit)
 -> full-order/conductor membership
 -> affine-open incidence
 -> cofactor/inverse-different compatibility                         (28)
```

are distinct gates.

## 7. Application to the THM-3230 centered norm

At a tame pure-`C3` completion from THM-3230, the quartic algebra has the
form `(13)--(14)`, with the unique degree-one factor being the inertia-fixed
sheet.  Its centered norm is precisely `n_0`, by `(3)--(7)`.  Therefore the
leading coefficient `N_h` used there is obtained by projecting the one
global packet `nu_f` with `(21)`.  The recovered local cubeclasses `[K]` and
`Lambda` no longer require a coherent **combinatorial** mark across all
places, provided the full quartic algebra and its tautological element are
retained.

This does not make those local classes one base-field function.  Their
valuations occur at different degree-one primes of `Spec(A)`, and signed
divisor comparison still needs a geometric transport law.  Nor does
`nu_f` provide the sheetwise cofactor or inverse different required by
THM-3064.

## 8. Scalar-packet collision hostile

The irreducible `1+3` hypothesis is load-bearing for `(20)`.  Take

```text
f(X)=(X^2-9)(X^2-1)=X^4-10X^2+9.                         (29)
```

Formula `(7)` becomes

```text
N_f(X)=-20X(X-3)(X+3)/27.                                (30)
```

At the four roots `(3,-3,1,-1)`, the packet values are

```text
(0,0,160/27,-160/27),                                    (31)
```

and

```text
P_f(Z)=Z^2(27Z-160)(27Z+160)/729.                         (32)
```

Thus the scalar characteristic polynomial can identify two distinct sheets.
The algebra-valued element `(A,nu_f)` still retains them.  Without a known
irreducible `1+3` local type, neither a repeated value nor an arbitrary
simple scalar root is a fixed-sheet certificate.

## 9. Exact companion and scope

Run

```bash
python3 04-computation/jc_universal_quartic_centered_norm_packet_thm3271.py
python3 -O 04-computation/jc_universal_quartic_centered_norm_packet_thm3271.py
```

Both modes byte-match the stored transcript.  The companion checks the
symbolic synthetic identity `(7)`, all four rootwise identities, all four
affine cube-covariance identities, the cyclic control

```text
f=(X-2)(X^3-3X+1),              disc(X^3-3X+1)=81,
n_0=-1,
Delta_nu=424/729,
e_0=X^3/3-X+1/3,                                        (33)
```

the idempotent congruences on both factors, `360` irreducible bounded
`1+3` separation controls, and the collision hostile `(29)--(32)`.  Every
truth-bearing check uses an explicit exception and remains live under
optimized execution.

The theorem proves no integral projector when `(25)` is nonunit, no larger-
order conductor membership, no affine section, no cross-place scalar gluing,
no forbidden cubeclass, no Keller cofactor, and no exclusion of pure `C3`,
`A4`, `S4`, a degree-four Keller map, `JC(2)`, or `DC(2)`.

**QED.**
