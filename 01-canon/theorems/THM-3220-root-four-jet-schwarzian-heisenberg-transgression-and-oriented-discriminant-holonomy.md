---
id: THM-3220
title: "Root four-jet Schwarzian--Heisenberg transgression and oriented-discriminant holonomy"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Tangent-to-identity four-jets form the exact Heisenberg central extension
  of their abelian three-jet quotient.  At a selected common simple root the
  first coordinate lifts THM-3215's affine shear, the Schwarzian is the
  second coordinate, and a covariant fourth jet fills the half-symplectic
  triangle.  Quadratic factorial germs lie on a twisted cubic: their central
  commutator is an oriented cubic discriminant and its two-root deck line has
  the exact quadratic Frobenius character.  Odd-prime resets return as a
  primitive central carry.  No root selector or physical carrier is supplied.
source: root/multiscale-newton-flag/low-child-flag-extension/2026-08-02
audit: >
  The assertion-independent exact companion pins promoted THM-3215, proves
  the symbolic four-jet law and quadratic commutator, checks the exact
  THM-2779 coordinate dictionary, the characteristic-two D8 boundary, 216
  strict transition triangles, 110 positive-Witt brackets, 4,023 quadratic
  Vandermonde triples, 4,018 direct nonzero order-p coefficient jets, 68
  Frobenius-character cases, six full divided carries, coordinate weights,
  and the nonconstant-unit and higher-jet hostiles.  Two independent hostile
  proof audits accepted every algebraic sign and scope boundary.  A second
  implementation pins the primary script and output, rebuilds the laws from
  raw truncated series, checks nonlinear coordinate cancellation, 300
  quadratic-extension Frobenius eigenvectors, 492 direct order-p controls,
  and 24 raw full carries.  Both normal/-O/stored replays pass byte-for-byte.
depends_on:
  - THM-3215-arbitrary-degree-root-jet-hamiltonian-affine-dihedral-holonomy-and-p-fold-carry
related:
  - THM-2591-theta-zero-selector-cech-coboundary-and-c91-holonomy-no-go
  - THM-2779-bockstein-symplectic-decoder-frame-torsor-and-heisenberg-root-degree-gate
  - THM-3175-first-witt-hensel-wedge-and-infinitesimal-pluecker-covariance
  - THM-3188-quadratic-character-pre-reset-holonomy-and-exterior-flag-rigidity
  - THM-3206-heterogeneous-factorial-exterior-reflection-groupoid-and-fixed-plane-holonomy
script: 04-computation/root_four_jet_schwarzian_heisenberg_thm3220.py
output: 05-knowledge/results/root_four_jet_schwarzian_heisenberg_thm3220.out
script_sha256: f399b6e86b639e9b27859dcbc9d2fcc6b44f4fb24aae51706e205a39e38c9e98
output_sha256: 1255b5ceebdf2d6b26ee92d9cbc3901e686728ba7ee79fdf6dbe7be52934e743
independent_script: 04-computation/root_four_jet_schwarzian_heisenberg_independent_audit_thm3220.py
independent_output: 05-knowledge/results/root_four_jet_schwarzian_heisenberg_independent_audit_thm3220.out
independent_script_sha256: b3b4044295ef875b00eaa8b06569db2f074b17dfe9e1654c320db519e130c3f0
independent_output_sha256: d2365ee5cd397f530c5757d6bcd83ff323284cd947d4eaa069cfe71d824739fd
hash_basis: LF-normalized bytes
---

# THM-3220 -- root four-jet Schwarzian--Heisenberg transgression and oriented-discriminant holonomy

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3215 proves that the selected-simple-root two-jet connection is flat.
The first noncommutative layer occurs two orders later.  It is not an
analogy: the four-jet group is an exact Heisenberg central extension, and its
central coordinate is an oriented discriminant on the factorial quadratic
curve.

The theorem carefully separates three statements:

```text
full selected-root transition groupoid:     still strictly flat;
three-jet projection:                       has a central area invoice;
fourth jet:                                 explicitly fills that invoice. (1)
```

A genuine global holonomy requires a branch/deck seam or another obstruction
to choosing the fourth coordinate globally.

## 1. The tangent four-jet Heisenberg group

Let `R` be a commutative ring in which `2` is a unit.  Work modulo `u^5` with
tangent-to-identity germs

```text
phi(u)=u+alpha u^2+beta u^3+gamma u^4.                    (2)
```

Define the logarithmic coordinates

```text
A=alpha,
B=beta-alpha^2,
C=gamma-(5/2)alpha beta+(3/2)alpha^3.                     (3)
```

The convention for multiplication is ordinary functional composition, with
the outer germ written on the left.  If `phi_o` is outside `phi_i`, direct
substitution gives

```text
X(phi_o composed with phi_i)
 =(A_i+A_o,
   B_i+B_o,
   C_i+C_o+(A_i B_o-B_i A_o)/2).                          (4)
```

Thus `(A,B)` is an abelian quotient and `C` is central.  The inverse has
logarithmic coordinate `(-A,-B,-C)`.  For the functional commutator

```text
[phi,psi]=phi psi phi^(-1) psi^(-1),
```

equation `(4)` gives

```text
X([phi,psi])
 =(0,0,B_phi A_psi-A_phi B_psi).                          (5)
```

This proves that the tangent four-jet group is the class-two Heisenberg
extension of the `(A,B)` plane.

### Exact dictionary with the finite Heisenberg convention

THM-2779 uses

```text
(x,y,z)(x',y',z')=(x+x',y+y',z+z'-yx').                  (6)
```

The exact isomorphism, with the same outer-left composition order, is

```text
(x,y,z)=(B,A,C-AB/2).                                     (7)
```

Under `(7)`, equation `(5)` becomes THM-2779's determinant commutator.  This
is an algebraic group identification only.  It does not map a polynomial
root jet to THM-2779's LRC endpoint carrier.

### The sharp characteristic-two boundary

Formula `(3)` is unavailable in characteristic two, but the raw coefficient
group `(2)` still has exactly eight elements over `F_2`.  Direct composition
gives

```text
order census:       1^1, 2^5, 4^2;
center:             {(0,0,0),(0,0,1)};
(alpha,beta,gamma)^2=(0,0,alpha(alpha+beta)).              (8)
```

Hence the group is `D_8`.  With quotient coordinates
`(x,y)=(alpha,alpha+beta)`, its square form is `q(x,y)=xy`, exactly the
characteristic-two refinement in THM-2779.  The odd-characteristic
Heisenberg logarithm must not be transported through this boundary.

## 2. Selected-root normalization and the Schwarzian flag

Let `f` be a polynomial or formal germ with a selected simple root `a`.
Normalize

```text
phi_f(u)=f(a+u)/f'(a)
        =u+alpha_f u^2+beta_f u^3+gamma_f u^4+O(u^5).     (9)
```

This normalization is unchanged by multiplying `f` by a constant unit.  If
`6` is invertible, write

```text
kappa_f=f''/f',
S(f)=f'''/f'-(3/2)(f''/f')^2                              (10)
```

at `a`.  Substituting the Taylor coefficients in `(3)` gives

```text
A_f=kappa_f/2,
B_f=S(f)/6,
C_f=(S(f)'-kappa_f S(f))/24.                              (11)
```

Thus `B` is the Schwarzian coordinate and `C` is its first covariant
derivative.  The coefficient definition `(3)` remains valid in residue
characteristic three; equation `(11)` is only the derivative presentation.

The THM-3215 hostile `f(u)=u`, `g(u)=u+u^3` has equal two-jets but

```text
A_g-A_f=0,                  B_g-B_f=1.                    (12)
```

The Schwarzian layer is therefore the first exact repair of that local
higher-jet collision.

## 3. Strict transitions and the central triangle transgression

For selected-root germs `f_i`, define

```text
H_(j<-i)=phi_j composed with phi_i^(-1),
X(H_(j<-i))=(A_ji,B_ji,C_ji).                             (13)
```

Equations `(3)--(4)` give

```text
A_ji=A_j-A_i,
B_ji=B_j-B_i,
C_ji=C_j-C_i-(A_i B_j-B_i A_j)/2.                        (14)
```

The first coordinate is the literal lift of THM-3215's affine transition:

```text
R(kappa_j)R(kappa_i)=T(kappa_i-kappa_j)=T(-2A_ji).        (15)
```

For a composable triangle, the full germs telescope strictly, while the
central coordinates obey

```text
C_ki=C_ji+C_kj+(A_ji B_kj-B_ji A_kj)/2.                  (16)
```

The last term is the half-symplectic area of the projected two edges.  It is
the curvature invoice of the abelian three-jet projection, and `C` is its
explicit transgression.

### The area term is load bearing

Take vertex logarithms

```text
X_0=(0,0,0),           X_1=(1,0,0),           X_2=(0,2,0).
```

They are realized by

```text
phi_0=u,
phi_1=u+u^2+u^3+u^4,
phi_2=u+2u^3.                                                (17)
```

Then

```text
C_10=0,                 C_21=-1,              C_20=0.     (18)
```

Naive central addition misses by one.  The area term in `(16)` is exactly
`+1`, so the strict transition identity closes.  This is both a positive
control and the sharp refutation of an uncorrected additive central law.

### Coordinate covariance

Let `y` be another source coordinate and put `rho=dy/dx(a)`.  The normalized
root transitions satisfy

```text
H_(j<-i)^(y)=D_rho H_(j<-i)^(x) D_rho^(-1),
D_rho(u)=rho u.                                            (19)
```

The nonlinear part of the source change cancels between the two normalized
germs.  Hence

```text
(A_ji,B_ji,C_ji)^(y)
 =(rho^(-1)A_ji,rho^(-2)B_ji,rho^(-3)C_ji)^(x).           (20)
```

The area in `(16)` has the same cubic weight as `C`.  Equivalently,
`A dx`, `B dx^2`, and `C dx^3` are intrinsic transition tensors.

## 4. Factorial quadratics form a twisted cubic

Let

```text
q_i(t)=v_i t^2-t+s_i
```

share the selected affine root `x`, and put

```text
lambda_i=q_i'(x)=2v_i x-1,            r_i=v_i/lambda_i.   (21)
```

The normalized root germ is exactly

```text
phi_i(u)=u+r_i u^2.                                       (22)
```

Therefore

```text
X_i=(r_i,-r_i^2,(3/2)r_i^3).                             (23)
```

In the THM-2779 coordinates `(7)`, this is

```text
(x_i,y_i,z_i)=(-r_i^2,r_i,2r_i^3).                       (24)
```

After adjoining the constant coordinate, `(24)` is a linear permutation and
rescaling of `[1:r:r^2:r^3]`, the twisted cubic.

For three root charts, the central area in `(16)` is

```text
A_ji B_kj-B_ji A_kj
 =-(r_j-r_i)(r_k-r_j)(r_k-r_i).                          (25)
```

Its square is exactly

```text
Disc((T-r_i)(T-r_j)(T-r_k)).                              (26)
```

Thus the fourth-order central fill is an oriented discriminant square root.
It is a unit exactly when the three normalized curvatures are pairwise
separated by units.

Two quadratic germs already have a nontrivial functional commutator:

```text
[u+r u^2,u+s u^2]
 =u+r s(s-r)u^4                         (mod u^5).         (27)
```

The square of its central coefficient is
`Disc(T(T-r)(T-s))`.  At the common factorial root `x=0`, the two choices
`v=-1,-2` give `r=1,2`, hence the primitive central jet

```text
u+2u^4.                                                      (28)
```

## 5. The positive-Witt bracket and its resonance wall

The Heisenberg center is the first case of a degree-independent law.  Let

```text
Phi(u)=u+c u^m+O(u^(m+1)),
Psi(u)=u+d u^n+O(u^(n+1)),                m!=n.            (29)
```

Direct substitution, followed by the formal inverses, gives

```text
Phi Psi Phi^(-1) Psi^(-1)
 =u+(m-n)cd u^(m+n-1)+O(u^(m+n)).                         (30)
```

This is the positive-Witt bracket.  The sign is fixed by the functional
commutator order in `(30)`; derivations themselves satisfy

```text
[u^m d/du,u^n d/du]=(n-m)u^(m+n-1)d/du.                  (31)
```

The case `(m,n)=(2,3)` is `(5)`.  In characteristic dividing `m-n`, the
leading bracket dies.  This resonance is a sharp obstruction to extrapolating
the unit four-jet commutator blindly through the whole osculating tower.

## 6. The discriminant deck is a genuine character seam

For one genuine quadratic `q(t)=vt^2-t+s`, assume `v(1-4sv)` is a unit.  On
the discriminant double cover put

```text
Delta=1-4sv,                 delta^2=Delta.                (32)
```

The two simple roots have

```text
q'(x_+)=delta,       q'(x_-)=-delta,
r_+=v/delta,         r_-=-v/delta.                        (33)
```

Their ordered central commutator is

```text
Gamma=r_+ r_-(r_--r_+)=2v^3/delta^3,
Gamma^2=4v^6/Delta^3.                                    (34)
```

The deck involution swaps the ordered branches and sends `Gamma` to
`-Gamma`.  Over `F_p`, with `p` odd and `Delta!=0`, Euler's criterion gives

```text
delta^p=chi(Delta)delta,
Gamma^p=chi(Delta)Gamma.                                  (35)
```

This is a genuine global branch-seam holonomy: the central four-jet is the
quadratic-character eigenline on the discriminant cover.  Equation `(35)`
matches the character in THM-3188 at the representation level.  There is no
proved comparison map from this root-jet line to THM-3188's Gauss--Manin
exterior state, so no exterior conclusion is asserted.

## 7. Odd-prime reset and first central carry

Let `O` be a mixed-characteristic DVR with odd residue characteristic `p`.
Suppose the logarithmic coordinates `(A,B,C)` of an integral four-jet lie in
`O`.  Since the alternating form of a vector with itself is zero, `(4)` gives

```text
X(phi^n)=(nA,nB,nC).                                      (36)
```

Returning to ordinary coefficients gives

```text
alpha_n=nA,
beta_n=nB+n^2A^2,
gamma_n=nC+(5/2)n^2AB+n^3A^3.                            (37)
```

If `(A,B,C)` is nonzero modulo `p`, equations `(36)--(37)` show that the
reduced jet has exact order `p`.  Moreover every coefficient of `phi^p-u`
is divisible by `p`, and

```text
(alpha_p,beta_p,gamma_p)/p
 ==(A,B,C)                                      (mod p).  (38)
```

Thus the mod-`p` reset returns as the full primitive logarithmic carry.  For
the factorial control `(28)`, it is the central carry `(0,0,2)` for every odd
prime.

## 8. Sharp scope boundaries

### 8.1 The full selected-branch Cech class is still zero

Equation `(16)` is an exact transgression.  On one globally selected simple
root branch, the full transition germs compose strictly and every triangle
closes.  Calling the area alone a nonzero Cech class would discard its
displayed fourth-coordinate filling.  Genuine global holonomy in `(35)`
appears only because the discriminant deck obstructs a single ordered branch.

### 8.2 Four jets do not recover an arbitrary formal germ

The distinct germs

```text
f(u)=u,                         g(u)=u+u^5                (39)
```

have identical four-jets.  Arbitrary formal recovery needs the full
osculating tower.  Finite-degree reconstruction is treated separately by the
next osculating theorem; it is not hidden in the present claim.

### 8.3 A nonconstant common unit changes the higher jet

For `f=u`, `g=u+u^2`, the original `B` transition is `-1`.  Multiplying both
by the common unit `1+t u` changes it to

```text
B_(g<-f)=-1-t.                                             (40)
```

Therefore the four-jet invariant belongs to the actual normalized polynomial
germs, not merely to their common zero divisor modulo nonconstant units.

### 8.4 Multiple roots and infinity remain separate charts

At a multiple root, division by `f'(a)` is impossible.  At infinity,
polynomials of different degrees require a reciprocal coordinate and common
homogeneous trivialization.  Neither boundary is covered by `(9)`.

## 9. Cross-frontier connection contract

THM-2591 proves that an LRC vertex root selector changes transition data only
by a Cech coboundary and that the missing object must be a genuinely
transition-dependent mixed-square two-cell.  The present theorem supplies an
exact algebraic model of that mechanism type:

```text
source:       selected-simple-root polynomial four-jets;
map:          projected edge pair |-> half-symplectic central area;
target:       Heisenberg center / oriented discriminant line;
preserved:    transition order, coordinate weight, branch character;
destroyed:    physical root chart, owner/current provenance, positivity;
needed:       a carrier map identifying this center with a lawful mixed cell. (41)
```

The exact group dictionary `(7)` does not construct that map.  Likewise the
character equality `(35)` does not identify the THM-3188 exterior carrier.
No LRC row, factorial staircase row, Jacobian stratum, or Gaussian-moment
conclusion follows from these algebraic identifications alone.

## 10. Exact evidence

Run

```text
python 04-computation/root_four_jet_schwarzian_heisenberg_thm3220.py
python -O 04-computation/root_four_jet_schwarzian_heisenberg_thm3220.py
python 04-computation/root_four_jet_schwarzian_heisenberg_independent_audit_thm3220.py
python -O 04-computation/root_four_jet_schwarzian_heisenberg_independent_audit_thm3220.py
```

and compare LF-normalized bytes with the declared output.  The companion is
assertion-independent and contains no floating literals or randomness.  It
pins the promoted THM-3215 theorem; proves `(4)` and `(27)` symbolically;
checks the exact THM-2779 coordinate dictionary and `D_8` boundary; replays
`216` strict transition triangles, `110` Witt brackets, `4,023` Vandermonde
triples, `4,018` direct coefficient-jet order computations, `68` Frobenius
character cases, six full carries, all displayed coordinate weights, and the
two sharp hostiles `(39)--(40)`.

Two independent hostile audits rederived the outer-left group law, the exact
THM-2779 dictionary and characteristic-two square form, the Schwarzian and
covariant fourth coordinate, strict transition transgression, source-change
weights, every quadratic orientation sign, the Witt sign convention, deck
character, and divided prime carry.  Both found the carrier/no-selector scope
honest.  A fresh immutable normal replay byte-matches the stored transcript;
the candidate's recorded normal and optimized replays and LF hashes remain
unchanged.

The independent companion imports no primary implementation.  Its generic
raw-series engine separately derives composition and inversion; checks `216`
strict triangles and `72` genuinely nonlinear source-coordinate changes;
reconstructs the characteristic-two group, twisted cubic, oriented
discriminant, and positive-Witt sign; computes the deck action and Frobenius
eigencharacter in `300` explicit quadratic algebras; and replays `492`
nonzero order-`p` jets and `24` full first carries.  It pins the immutable
hashes of the primary script and output, and its normal, optimized, and stored
outputs agree byte-for-byte.

**QED.**
