---
id: THM-2590
title: "Boolean Bockstein and theta-selector incidence spectrum"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE UNDER AUDIT.  On the fixed primitive
  coefficient carrier of THM-2585, the thirteen septimal target-slice factors
  span the full six-dimensional algebra over F_13.  The section map has rank
  six and kernel dimension seven; 1,581 of its 1,716 six-column minors are
  nonzero.  Its restriction to the 8,192 Boolean subsets has 8,184 images,
  with the empty set as the unique zero fibre and one reduced collision
  circuit, up to adjoining a common subset of three free indices.  Thus every
  nonempty Boolean sum has nonzero first Bockstein in all six labelled owner
  colour.  Formally inserting every one of the 1,312 admissible THM-2586
  theta-zero selector patterns into the septimal coefficient rows never gives
  zero in any labelled owner colour; 1,272 patterns are units, and the actual
  priority selector is a unit
  in all twelve displacement cells.  This last statement is only a finite
  coefficient-incidence splice across distinct sigma={a} and sigma={b}
  carriers.  It is not a physical carrier composition, semantic endpoint,
  row exclusion, or LRC(14) conclusion.
source: wild-holotopy-mining-2026-07-28
depends_on:
  - THM-2585-saturated-normalized-target-projector-and-bockstein-noncommutation
  - THM-2586-depth-five-arrival-to-future-root-diagonal
related:
  - THM-2532-cyclic-tournament-cayley-algebra-and-chi7-even-quotient
script: 04-computation/lrc14_boolean_bockstein_selector_incidence_thm2590.py
output: 05-knowledge/results/lrc14_boolean_bockstein_selector_incidence_thm2590.out
script_sha256: c5ac1a25ed4fafa09cbb838b368fd4be990020e8d09fbbcbaf6cd7ccce9681da
output_sha256: 34674303ad72cca08fe608486dffdee2aca679faa0a69d6f3f77c1c2f69d60f5
hash_basis: LF-normalized bytes
---

# THM-2590 -- Boolean Bockstein and theta-selector incidence spectrum

**RESERVED / PROVISIONAL PROOF CANDIDATE UNDER AUDIT.**

THM-2585 proves that every one of its thirteen literal target-shift sections
has nonzero first Bockstein in every owner colour.  That does not by itself
control cancellation after several sections are selected or after a root
choice is allowed to vary with the seven owner-clock cells.  The exact finite
answer is stronger than pointwise nonvanishing:

```text
nonempty Boolean set of target sections -> all six owner Bocksteins nonzero;

any admissible THM-2586 theta-zero selector pattern
  -> nonzero formal owner polynomial;

the actual THM-2586 priority selector
  -> a unit, hence all six owner colours survive.                 (1)
```

The first line lives wholly on THM-2585's one common primitive carrier.  The
last two lines are deliberately only coefficient-incidence statements.  They
do not identify THM-2585's `sigma={a}` carrier with THM-2586's
`sigma={b}`, depth-five physical packet.

## 1. The exact section map

Work in the septimal algebra

```text
R_7=F_13[z]/(Phi_7(z)),

Phi_7(z)=1+z+...+z^6
        =(z^2+3z+1)(z^2+5z+1)(z^2+6z+1).                 (2)
```

Use the thirteen rows printed in THM-2585:

```text
Y_q(z)=sum_(ell=0)^6 Y_(ell,q)z^ell,       0<=q<13.       (3)
```

In the power basis `1,z,...,z^5`, reduction modulo `Phi_7` sends a seven-row
`(y_0,...,y_6)` to

```text
bar y=(y_0-y_6,...,y_5-y_6).                            (4)
```

Define the linear section map

```text
L:F_13^13 -> R_7,

L(c_0,...,c_12)=sum_(q=0)^12 c_q bar Y_q.                (5)
```

The exact row reduction of its `6 x 13` matrix is

```text
[1 0 0 0 0 0 |  2  0 11  3  2  1  9]
[0 1 0 0 0 0 |  9  5  2  0  6  7  4]
[0 0 1 0 0 0 |  0 12 12  2  0  7  2]
[0 0 0 1 0 0 | 11  2 12 11  6 11  6]
[0 0 0 0 1 0 |  2  1  6 12  5  8  2]
[0 0 0 0 0 1 | 10  5  9 12  9  0  0].                  (6)
```

In particular, the determinant of columns `q=0,...,5` is `11 mod 13`.
Therefore

```text
rank(L)=6,                    dim ker(L)=7.               (7)
```

For completeness, a kernel basis, with coordinates ordered `q=0,...,12`, is

```text
(11,4,0,2,11,3, 1,0,0,0,0,0,0),
(0,8,1,11,12,8, 0,1,0,0,0,0,0),
(2,11,1,1,7,4, 0,0,1,0,0,0,0),
(10,0,11,2,1,1, 0,0,0,1,0,0,0),
(11,7,0,7,8,4, 0,0,0,0,1,0,0),
(12,6,6,2,5,0, 0,0,0,0,0,1,0),
(4,9,11,7,11,0, 0,0,0,0,0,0,1).                        (8)
```

The companion independently checks every vector in (8).  It also exhausts
all `binom(13,6)=1716` six-column minors:

```text
nonzero minors =1581,                 zero minors =135.   (9)
```

Thus the rank certificate is highly redundant; it is not an accident of one
chosen basis.

## 2. The Boolean cube has no nonempty zero fibre

Restrict (5) to the vertices of the Boolean cube:

```text
B:{0,1}^13 -> R_7,

B(1_S)=sum_(q in S)bar Y_q.                              (10)
```

Exhaustion of all `8192` subsets gives

```text
number of images =8184,
maximum fibre size=2,
number of double fibres=8,
B^(-1)(0)={empty set}.                                   (11)
```

The eight collisions all reduce to the single signed circuit

```text
Y_3+Y_4+Y_6+Y_7+Y_9+Y_10
  =Y_1+Y_5+Y_8+Y_12                       in R_7,          (12)
```

with an arbitrary common subset of `{0,2,11}` adjoined to both sides.  There
are no other Boolean collisions.

THM-2585 gives, on its fixed common primitive carrier,

```text
beta(D^(kappa,q))=Omega Y_q(zeta_7^kappa).                (13)
```

The target sections are literal disjoint `s=-q` slices before Fourier
reduction, so their coefficient sum is lawful on that carrier.  Equations
(11) and (13) imply

```text
beta(sum_(q in S)D^(kappa,q))
 =Omega (sum_(q in S)Y_q)(zeta_7^kappa),                  (14)
```

For every `kappa in F_7^*`, substitution

```text
sigma_kappa:R_7 -> R_7,               z -> z^kappa       (14a)
```

is an automorphism, with inverse `sigma_(kappa^(-1) mod 7)`.  Therefore a
nonzero class remains nonzero under every labelled owner substitution.
There is one further zero-divisor issue to settle: the factor `Omega` in
(14).  Put `u=zeta_13-1`.  In characteristic thirteen,

```text
Phi_13(1+u)=u^12,

Omega=sum_(m=0)^12 m(1+u)^m=-u^11,                       (14b)

M/13M = R_7 tensor_F13 F_13[u]/(u^12).
```

Thus multiplication `P -> Omega P=-P tensor u^11` embeds `R_7` into the
socle summand `R_7 tensor F_13 u^11`.  In particular,
`Omega P=0` if and only if `P=0`, even when `P` is a zero divisor of `R_7`.
Equations (11)--(14) imply that, for every nonempty `S`, the right side of
(14) is nonzero for **all six** `kappa in F_7^*`.  This is the first
assertion in (1).

Nonzero is load-bearing: it is not always a unit.  From (2), a class is a
unit exactly when none of the three displayed quadratics divides it.  The
complete factor-profile census is

| divided factors, recorded by middle coefficients | subsets |
|:---|---:|
| none | 8011 |
| `(3)` | 59 |
| `(5)` | 55 |
| `(6)` | 63 |
| `(3,5)` | 1 |
| `(3,6)` | 1 |
| `(5,6)` | 1 |
| `(3,5,6)` | 1 |

The last row is the empty subset.  Hence among the `8191` nonempty subsets,
`8011` are units and `180` are nonzero zero divisors.  A quadratic factor
kills one CRT field component.  It does **not** kill a labelled owner colour:
the six maps (14a) permute the three CRT components and preserve nonzeroness
of the whole ring element.  A nonunit can instead be annihilated by a nonzero
downstream coefficient supported in complementary CRT components.  The
collision (12) and these `180` zero divisors are the sharp hostile controls
against injectivity and against replacing "nonzero" by "unit".

## 3. Every admissible theta-selector splice is nonzero

THM-2586 proves that its two theta-zero rails have root values `0` and `6`.
For `s in F_13^*`, let

```text
h=(h_0,...,h_6) in {0,6}^7                              (15)
```

choose any positive rail edge.  Its exact two zero sets say precisely

```text
s=6:       h_4=h_5=h_6=0,
s=7:       h_4=h_5=h_6=6,
s notin {6,7}: no further restriction.                   (16)
```

Thus there are

```text
10*2^7+2*2^4=1312                                       (17)
```

labelled admissible pairs `(s,h)`.  Now make only the naked label
substitution `q=h_ell`, identifying two copies of the set `F_13` while not
identifying their carriers, and define the formal incidence splice by

```text
Z_h(z)=sum_(ell=0)^6 Y_(ell,h_ell)z^ell in R_7.           (18)
```

This is an exact operation on the printed coefficient table, not yet an
operation on a common physical packet.  Exhausting (16) gives

```text
labelled patterns                1312,
distinct classes                   32,
zero classes                         0,
unit labelled patterns            1272,
nonunit labelled patterns           40,
distinct units                       31,
distinct nonunits                     1.                 (19)
```

The only nonunit class is

```text
(9,5,8,4,0,0)                                             (20)
```

in the six-term power basis.  It is divisible by `z^2+3z+1` and by neither
other factor, so it is nonzero but not invertible.  All six owner
substitutions (14a) remain nonzero; their power-basis support sizes are
`4,6,6,4,4,6` for `kappa=1,...,6`.  It occurs exactly when

```text
s notin {6,7},
h_0=h_6=0,                h_2=h_3=h_4=6,                 (21)
```

with the inert choices `h_1,h_5` arbitrary.  Their inertness is visible
directly in THM-2585's rows:

```text
Y_(1,0)=Y_(1,6)=9,        Y_(5,0)=Y_(5,6)=4.              (22)
```

Consequently every admissible splice has a nonzero Bockstein class in all six
labelled owner colours.  All but the forty patterns in (21) are additionally
units, so they cannot be annihilated by any nonzero downstream coefficient.

## 4. The actual priority selector is everywhere a unit

THM-2586's priority rule chooses its root-`0` edge whenever positive and
uses root `6` only in the three cells

```text
(s,ell)=(7,4),(7,5),(7,6).                               (23)
```

For eleven displacements, (18) is therefore exactly `Y_0`, whose
multiplication norm is `7 mod 13` by THM-2585.  At `s=7`, the selected row is

```text
h=(0,0,0,0,6,6,6),

bar Z_h=(6,2,2,6,10,10),

Norm_(R_7/F_13)(Z_h)=10 mod 13.                           (24)
```

Both norms are nonzero.  Hence the formal priority-selector splice is a unit
for all twelve nonzero displacements and all six owner colours survive.

This isolates the exact remaining issue: if a lawful common-carrier map
realized the coefficient substitution (18), Bockstein cancellation would not
be the obstruction.  The theorem does not construct that map.

## 5. Scope and failure boundary

The following distinctions are part of the statement.

1. Sections in Section 2 all come from THM-2585's single globally cleared
   `sigma={a}` coefficient tensor.  No division is performed on independent
   representatives.
2. Sections 3--4 read THM-2586's proved selector incidence but do not pull its
   `sigma={b}`, depth-five Boolean fibres back to the THM-2585 tensor.  The
   common speed row does not make the words, clocks, sheets, or primitive
   normalizations identical.
3. A nonzero class in `R_7` survives every labelled owner automorphism (14a).
   A unit additionally survives multiplication by every nonzero downstream
   coefficient; (20) is the exact selector zero-divisor hostile.
4. The rank-six surjection in (7) is a signed linear correction statement.
   It does not preserve positivity, Boolean realizability, semantic owner
   labels, old-head/future-root identity, or higher 13-adic carry data.
5. Nothing here constructs a Cech transition correction, a THM-2334 relation
   current, a positive terminal endpoint, an all-165 theorem, a row exclusion,
   or LRC(14).

## 6. Exact verification

The dependency-free companion uses only Python's standard library.  It
checks the irreducible factorization (2), row reduction (6), determinant and
all minors,
all seven kernel vectors, all `8192` Boolean subsets, every collision fibre,
all `1312` selector patterns, all `57,018` labelled owner substitutions, the
unique nonunit selector class, and the priority norms.  Its immutable
transcript contains `83,159` explicit exact checks.

Reproduce with

```bash
python 04-computation/lrc14_boolean_bockstein_selector_incidence_thm2590.py
python -O 04-computation/lrc14_boolean_bockstein_selector_incidence_thm2590.py
```

Both modes must byte-match

```text
05-knowledge/results/lrc14_boolean_bockstein_selector_incidence_thm2590.out
```

on the declared LF-normalized hash basis.

QED for the provisional coefficient-incidence theorem, pending independent
hostile audit and status promotion.
