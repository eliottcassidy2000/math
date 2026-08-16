---
id: THM-3508
title: "Level-two sporadic Keller three-coordinate primitive discriminant square class"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENT CLEAN-ROOM AUDIT.  For the fixed
  sporadic Keller map at level two, each source coordinate x,y,z is a
  primitive element of the generic degree-nine function-field extension.
  Their monic minimal-polynomial discriminants, and every saturated
  degree-nine coordinate eliminant obtained from them by target-field
  scaling, therefore have the common full square class [H].  The proof uses
  exact nonzero 9x9 power-basis determinants at the lawful etale target
  (1,1,1), a direct split-fibre replay over F_41, the trace-form basis-change
  identity, and THM-2582's exact x-coordinate class [H].  This is fixed-map,
  level-two scope; it gives no higher-level, exact multiplicity, arbitrary
  Keller, JC(2), DC(2), or LRC consequence.
source: codex/primitive-coordinate-audit/2026-08-16
audit: >
  Independent regular-representation implementation using python-flint;
  imports no repository computation module, candidate script, precomputed
  coordinate eliminant, or precomputed power-basis matrix.  A second direct
  finite-field route evaluates the original map equations on all of F_p^3.
depends_on:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-2576-composite-jelonek-image-divisor-and-two-component-nonproperness-law
  - THM-2582-odd-block-discriminant-tower-and-composite-jelonek-square-class
related:
  - THM-2546-integral-coordinate-dichotomy-and-parity-lens-scope
  - THM-3494-weighted-lift-primitive-coordinate-discriminant-atlas
  - MISTAKE-413
script: 04-computation/keller_level_two_three_coordinate_primitive_independent_audit.py
output: 05-knowledge/results/keller_level_two_three_coordinate_primitive_independent_audit.out
script_sha256: ecf1821e9d918904d7b878931aa9ff26137e8ef4eac97234c6aec48cc5ecae5f
output_sha256: 2d659f01f404745a73f0b03fb3556bba0bf99cf897a11ea2607bd6fc4b705311
hash_basis: LF-normalized bytes
---

# THM-3508 -- all three level-two coordinate views carry `[H]`

**PROVED + VERIFIED-EXACT + INDEPENDENT CLEAN-ROOM AUDIT.**

## 1. Statement

Let `F` be the fixed sporadic Keller map of THM-2473.  Write

```text
(A,B,C)=F(x,y,z),              (a,b,c)=F(A,B,C),
K=Q(a,b,c),                    M=Q(A,B,C),
L=Q(x,y,z).                                               (1)
```

The generic degrees are

```text
[M:K]=3,                       [L:M]=3,
[L:K]=9.                                                  (2)
```

Then each of the three actual source coordinates is primitive:

```text
K(x)=K(y)=K(z)=L.                                        (3)
```

Let `m_x,m_y,m_z in K[T]` be their monic minimal polynomials.  If `H` is
THM-2576's second irreducible image polynomial, then

```text
[Disc(m_x)]=[Disc(m_y)]=[Disc(m_z)]=[H]
                                      in K*/K*2.          (4)
```

Equivalently, let `E_x,E_y,E_z` be any **saturated, actual degree-nine**
coordinate eliminants, so that each is a nonzero element of `K` times the
corresponding monic minimal polynomial.  Then their discriminants also obey
(4).  The saturation and degree clauses exclude an unsaturated resultant
with chart factors or a repeated polynomial attached to a nonprimitive view.

## 2. Inheritance and the missing gate

THM-2582 already proves the exact full field square class

```text
[Disc(m_x)]=[H].                                         (5)
```

It does not prove that `y` or `z` generates `L/K`.  THM-3494 supplies the
general trace-form mechanism once primitivity is known, but its fixed-map
composition hostile explicitly leaves this level-two primitivity gate open.
The closest corrected near miss is MISTAKE-413: an odd divisor carrier does
not determine a full field square class if a constant unit was discarded.
No unit is discarded here; (5) is THM-2582's exact class, and the remaining
comparisons are literal squares in `K`.

The canonical hostile is an intermediate coordinate such as `A`.  It
generates at most the degree-three field `M`, not the degree-nine field `L`.
Section 5 computes this failure exactly.

## 3. A lawful exact specialization reconstructed from the map

Use the suggested final target

```text
(a,b,c)=(1,1,1).                                         (6)
```

The outer first-coordinate core is

```text
g(t)=25t^3+t-2=(5t-2)(5t^2+2t+1),
Disc(g)=-67600 != 0.                                     (7)
```

Thus

```text
R=Q[t]/(g)                                               (8)
```

is a rank-three finite etale algebra.  It is **not** a cubic field.  This
disconnectedness is harmless: the specialization argument needs a nonzero
matrix determinant in a rank-nine algebra, not a connected number field.

For a target `(a,b,c)` and a first source coordinate `X`, put `s=XY` and

```text
D=-3X^2ac+X^2b^2c-X^2b+2Xbc-2X+c,

N=-3X^3abc+4X^3a-6X^2ac+X^2b^2c-X^2b
  +2Xbc-2X+c.                                           (9)
```

The inverse section is

```text
s=-N/D,                  Y=s/X,
Z=[X(2-3s)-c]/X^3.                                     (10)
```

The independent companion constructs (8) from scratch, applies (9)-(10) to
`X=t`, and verifies all of

```text
g(t)=0,
D s+N=0,
3aX^2-(1+s)(bX-s)=0,
aX^3+c(1+s)^3-X(1+s)(2+s)=0,
F(t,Y,Z)=(1,1,1).                                      (11)
```

The relevant denominators are units:

```text
Norm_R(D)=676/625,                  Norm_R(t)=2/25.      (12)
```

Call the reconstructed middle point `q=(t,Y,Z)`.  The inner first-coordinate
core is

```text
L(q)X^3+(4-3q_2q_3)X-2q_3.                            (13)
```

The script computes directly, by the determinant of multiplication in `R`,

```text
Norm_R(L(q))=190265288239/320
             =H(1,1,1)/(64*25),
H(1,1,1)=951326441195.                                 (14)
```

In particular `L(q)` is a unit.  Normalize (13) to a monic cubic and define
the rank-nine algebra

```text
S=R[X]/(X^3+[(4-3q_2q_3)/L(q)]X-2q_3/L(q)).             (15)
```

Applying (9)-(10) again in `S` reconstructs the inner `y,z`.  A second exact
substitution verifies

```text
F(x,y,z)=q,                    F(q)=(1,1,1),
F^2(x,y,z)=(1,1,1).                                    (16)
```

The new denominators are again units.  Their exact norms are

```text
Norm_S(D)=2667250410249863301596360358125
          /6854604479748118669088,

Norm_S(x)=2344000/190265288239.                         (17)
```

Equations (7), (12), (14), and (17) make the inverse tower well-defined at
(6).  Section 4's nonzero degree-nine discriminants add separability and
complete the lawful-specialization audit.

## 4. Three exact rank-nine power bases

Use the fixed tower basis

```text
e=(1,t,t^2,x,tx,t^2x,x^2,tx^2,t^2x^2)                 (18)
```

of `S/Q`.  For `xi in {x,y,z}`, let `P_xi` be the `9x9` matrix whose columns
are the coordinates of

```text
1,xi,xi^2,...,xi^8                                     (19)
```

in (18).  FLINT computes the matrices, ranks, determinants, regular
representation characteristic polynomials, minimal polynomials, and
discriminants over exact `Q`.  The result is

| coordinate | `rank(P_xi)` | `det(P_xi)` | minimal degree | minimal-polynomial coefficient hash |
|---|---:|---|---:|---|
| `x` | `9` | nonzero | `9` | `6cff42de...f591d89` |
| `y` | `9` | nonzero | `9` | `9eb35cd8...df91a0` |
| `z` | `9` | nonzero | `9` | `5e08a76c...95291d` |

The full exact rational determinants and discriminants are printed in the
matching stored output.  No floating-point rank decision is made.  For every
coordinate the matrix minimal polynomial equals its degree-nine
characteristic polynomial and has nonzero discriminant.

This already proves that the three specialized power bases span all of `S`.

## 5. Hostile: a degree-three intermediate coordinate is nonprimitive

Embed the outer coordinate `t=A` into the rank-nine algebra (15).  Its exact
minimal polynomial is

```text
t^3+(1/25)t-2/25,                                      (20)
```

while multiplication by `t` on `S` has characteristic polynomial equal to
the cube of (20).  Its nine-column power matrix satisfies

```text
rank(1,t,...,t^8)=3,                  determinant=0.     (21)
```

Thus `t` is primitive only for the intermediate rank-three algebra and is
nonprimitive for the full rank-nine algebra.  A naive degree-nine eliminant
is `(20)^3` and has zero discriminant.  This is the first failed implication
if one tries to compare coordinate discriminants without proving
primitivity.

## 6. Why one specialization proves generic primitivity

On the generic open set where the two cubic leading coefficients and the
inverse denominators are units, the tower basis corresponding to (18) is a
`K`-basis of `L`.  For each source coordinate `xi`, form the same power matrix
`P_xi` over `K`.  Its determinant is a rational function of `(a,b,c)`.

The specialization map at (6) is defined on every entry and sends this
generic matrix to the exact matrix of Section 4.  Since its specialized
determinant is nonzero, the generic determinant cannot be the zero rational
function.  Hence

```text
1,xi,...,xi^8
```

are linearly independent over `K`.  Because `[L:K]=9`, they are a basis and
`K(xi)=L`.  Applying this to `x,y,z` proves (3).

This is a nonvanishing argument, not an assertion that every specialization
is connected or primitive.  The split factorization in (7) is therefore a
boundary audit, not a defect.

## 7. Basis change gives one discriminant class

For a primitive coordinate `xi`, let `v_xi` be the determinant taking the
tower basis (18) to its power basis (19).  The trace Gram matrices are related
by congruence, so

```text
Disc(m_xi)=v_xi^2 Disc(e).                             (22)
```

Consequently

```text
Disc(m_y)/Disc(m_x)=(v_y/v_x)^2,
Disc(m_z)/Disc(m_x)=(v_z/v_x)^2.                       (23)
```

The exact companion verifies (23) as rational identities at (6), not merely
as square/nonsquare tests.  It also verifies separately that each of

```text
Disc(m_x)/H(1,1,1),
Disc(m_y)/H(1,1,1),
Disc(m_z)/H(1,1,1)                                    (24)
```

is an exact positive rational square and prints all three square roots.

Generically, (22) and THM-2582's exact identity (5) give (4).  This uses the
full class `[H]`, including its constant unit; it does not infer a field class
from divisor parity.

Finally, if an actual coordinate eliminant is

```text
E_xi=c_xi m_xi,                     c_xi in K*,        (25)
```

then its degree is nine and

```text
Disc(E_xi)=c_xi^(2*9-2) Disc(m_xi)
           =c_xi^16 Disc(m_xi)
           =(c_xi^8)^2 Disc(m_xi).                     (26)
```

Thus monic normalization, denominator clearing, and any other scalar change
of a saturated degree-nine eliminant preserve (4).  Equation (26) is the
required eliminant-scaling parity audit.

## 8. Orthogonal direct finite-field replay

As a second route, the companion ignores (9)-(15) and evaluates the original
three displayed coordinate formulas for `F` directly on every point of
`F_p^3` for the prime bank

```text
p=11,13,17,19,23,29,31,37,41,43.                       (27)
```

At `p=41`, target `(13,0,11)` has a split fibre of exactly nine enumerated
source points.  The counts of distinct source-coordinate values are

```text
(#x,#y,#z)=(9,9,9),                                    (28)
```

and the corresponding Vandermonde determinants are

```text
(1,14,12) in F_41*.                                    (29)
```

The three middle points are

```text
(36,6,10),              (8,12,6),              (38,23,13),
```

The old Jelonek polynomial has values `20` at the final target and
`37,36,8` at these three middle points.  All four level-one targets are
therefore off `V(L)` in this good odd characteristic, so the three displayed
middle points and their three children each exhaust the corresponding cubic
fibres.

Their first coordinates occur on the nine source points with
multiplicities `(3,3,3)`.  Thus this route independently reproduces both the
three primitive views and the degree-three block hostile directly from the
map equations.  It is a redundant control; the characteristic-zero proof is
Sections 3-7.

## 9. Scope, losses, and boundaries

- The theorem concerns the one fixed map `F^2` and the generic
  characteristic-zero function field.  It does not assert that every source
  coordinate of every finite polynomial map is primitive.
- Primitivity does not make `F^2` proper or invertible.  THM-2576 still gives
  `S_(F^2)=V(LH)`, while THM-2582 and this theorem say that the odd
  discriminant class sees `H`, not the full reduced nonproperness set.
- Equality of discriminant square classes preserves only the permutation
  sign character.  It forgets labelled sheets, full wreath monodromy,
  projection-collision divisors, basis volumes, and exact multiplicities.
- No exact positive multiplicity of `H`, no boundary saturation theorem for
  an arbitrary raw resultant, and no level-three-or-higher three-coordinate
  primitivity law follows.
- There is no consequence for `JC(2)`, `DC(2)`, `GMC(2)`, or LRC(14).

The connection is therefore precisely typed:

```text
primitive coordinate power basis
    --trace-form congruence--> common discriminant square class;

preserved: permutation sign and [H];
destroyed: labelled roots, basis volume and exact divisor multiplicity;
sidecar: unsquared basis determinant plus saturated eliminant;
hostile: the rank-three intermediate coordinate (20)-(21).             (30)
```

## 10. Reproduction

Run

```bash
python3 -B 04-computation/keller_level_two_three_coordinate_primitive_independent_audit.py
python3 -B -O 04-computation/keller_level_two_three_coordinate_primitive_independent_audit.py
```

The ordinary and optimized transcripts agree byte-for-byte after newline
normalization.  Every gate uses `require`, not `assert`, so optimized mode
executes the same checks.  The script imports only python-flint and the Python
standard library.  It reconstructs the two inverse stages, verifies the map
equations, computes all rational matrices and discriminants, runs the hostile,
and performs the direct finite-field enumeration.  The LF-normalized hashes
are pinned in frontmatter.

**QED.**
