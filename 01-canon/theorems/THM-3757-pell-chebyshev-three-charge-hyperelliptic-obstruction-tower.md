---
id: THM-3757
title: "Pell-Chebyshev three-charge hyperelliptic obstruction tower"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Over an algebraically
  closed characteristic-zero field, an explicit Chebyshev tower of
  three-charge polynomials Q_n=X+chi_n(XT)+T psi_n(XT), n>=1, has no critical
  point.  Its torus critical resultant is exactly -psi_n, and every root of
  psi_n is absorbed with multiplicity.  Nevertheless no Q_n has even a
  rational constant-Jacobian mate.  At n=1 the generic time form has two
  nonzero logarithmic residues; for n>=2 it is a nonzero holomorphic
  differential on a positive-genus hyperelliptic generic fibre.  This is a
  counterexample-component obstruction tower, not a counterexample to JC(2).
source: root + jc_quartic_c3_construct / 2026-08-23
audit: >
  PASS.  An independent hostile audit rederived the Chebyshev quotient and
  transport identities, Pell sign, exact resultant, absorbed-root and axis
  smoothness, generic-function-field equality and squarefreeness, Jacobian
  sign, positive-genus holomorphic obstruction, and conic infinity residues.
  Normal, optimized, and frozen output agree; script/output/semantic hashes
  and CHECKS=167 match.
depends_on: []
related:
  - THM-3551-one-ray-planar-jacobian-mate-no-go
  - THM-3598-danielewski-rational-exact-polar-graph-family-and-classification
  - THM-3741-radial-two-charge-keller-component-classification
script: 04-computation/jc2_pell_chebyshev_three_charge_hyperelliptic_thm3757.py
output: 05-knowledge/results/jc2_pell_chebyshev_three_charge_hyperelliptic_thm3757.out
script_sha256: 8679eee8131322243342e5748a36069bd5cee0216d6d12f7eb4d4cc71d8ea36e
output_sha256: a3add7c55f261bb1f0d90e3b24b44f50ff261a0ff5086feff7af7ef4c59a0c10
semantic_sha256: 551870d0b77f4ba55104af7d7b6a3112748a4d2646c3634239f5b89cca2e6234
hash_basis: raw LF bytes
---

# THM-3757 -- Pell cancellation buys smoothness but exposes hyperelliptic debt

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  THM-3741 proves
that two opposite radial charges cannot produce a planar Keller pair.  Adding
a charge-zero profile changes the critical equation from a derivative into a
polynomial Pell expression.  The construction below solves that smoothness
problem at every depth.  It then fails at a strictly later gate: the forced
time form on the generic fibre is not rationally exact.

Work over an algebraically closed field `k` of characteristic zero.  Denote
the Chebyshev polynomials of the first and second kinds by
`C_n(tau),U_n(tau)`, normalized by

```text
C_n(tau)^2-(tau^2-1)U_(n-1)(tau)^2=1.                 (1)
```

For every integer `n>=1`, put `tau=2z-1` and define

```text
G_n(z)=C_n(tau)+(tau-1)U_(n-1)(tau),
H_n(z)=C_n(tau)+(tau+1)U_(n-1)(tau),                  (2)

S_n(z)=U_n(tau)/[2(n+1)]-U_(n-1)(tau)/(2n),           (3)
chi_n(z)=integral_0^z G_n(s) ds,
psi_n(z)=(z-1)S_n(z)^2,
F_n(z)=z psi_n(z).                                     (4)
```

Finally set

```text
Q_n(X,T)=X+chi_n(XT)+T psi_n(XT).                      (5)
```

Then:

1. `Q_n` has no critical point in `k^2` for every `n>=1`.
2. There is no `P in k(X,T)` with `J(P,Q_n) in k*`; in particular there is
   no polynomial Jacobian mate.
3. The generic `Q_n`-fibre is hyperelliptic.  At `n>=2` its genus is at least
   `floor(n/2)>=1`, and the mate equation would make its nonzero holomorphic
   differential exact.  At `n=1` the same differential instead has two
   nonzero residues at infinity.

Thus `(5)` is an infinite tower of smooth, genuinely three-charge
noncoordinates, not a Jacobian-conjecture counterexample.

## 1. Chebyshev transport produces a supported Pell remainder

The standard recurrences give the quotient identities

```text
(tau+1)G_n=C_(n+1)+C_n,
(tau-1)H_n=C_(n+1)-C_n.                               (6)
```

Substituting `(2)` into `(1)`, or expanding `(6)`, yields

```text
zG_n(z)^2-(z-1)H_n(z)^2=1.                            (7)
```

The first-order transport operator

```text
L=(tau^2-1)d/dtau+tau                                  (8)
```

satisfies `L(U_j)=(j+1)C_(j+1)`.  Formula `(3)` therefore gives

```text
L(S_n)=[C_(n+1)-C_n]/2=(tau-1)H_n/2.                  (9)
```

Returning to `z`, equation `(9)` is exactly

```text
(2z-1)S_n+2z(z-1)S_n'=(z-1)H_n.                      (10)
```

Consequently

```text
F_n'=(z-1)S_nH_n.                                     (11)
```

Combining `(7)` and `(11)` produces the decisive identity

```text
F_n'^2-F_nG_n^2
 =(z-1)S_n^2[(z-1)H_n^2-zG_n^2]
 =-(z-1)S_n^2
 =-psi_n.                                             (12)
```

This is not a zero remainder.  It is a nonconstant remainder whose entire
zero divisor is deliberately placed on the profile boundary.

Two endpoint values will control the affine axes.  Since
`U_j(1)=j+1` and `U_j(-1)=(-1)^j(j+1)`, formula `(3)` gives

```text
S_n(1)=0,                    S_n(0)=(-1)^n,
psi_n(0)=-1.                                           (13)
```

Here `S_n(0)` means evaluation at `z=0`, hence `tau=-1`.

## 2. The critical resultant is exactly the absorbed profile

The following calculation explains why `(12)` is the right Pell equation.
For an arbitrary charge-zero derivative `G=chi'`, profile `psi`, and

```text
Q=X+chi(XT)+Tpsi(XT),             z=XT,
F=zpsi,                                               (14)
```

the two critical equations on the torus become, after eliminating
`T=z/X`,

```text
X^2+zGX+z^2psi'=0,
GX+F'=0.                                               (15)
```

Their exact resultant in `X` is

```text
Res_X(15)=F'^2-FG^2.                                   (16)
```

For the tower, `(12)` makes this resultant `-psi_n`.  If
`psi_n(z)!=0`, equations `(15)` therefore have no common root.  At a root
of `psi_n`, one also has `psi_n'=0`: every root of `S_n` occurs squared,
and the extra root `z=1` is also a root of `S_n` by `(13)`.  Hence
`F_n'=psi_n+zpsi_n'=0` there.  The original derivatives reduce to

```text
(Q_n)_T=XG_n,                    (Q_n)_X=1+TG_n.        (17)
```

On the torus, if `G_n!=0` the first entry is nonzero, while if `G_n=0` the
second is one.  Thus the absorbed resultant roots are all safe rather than
hidden critical points.

On `T=0`, equation `(5)` gives `(Q_n)_X=1`.  On `X=0`, equations
`(4)` and `(13)` give `(Q_n)_T=F_n'(0)=psi_n(0)=-1`.
Together with the torus analysis, this proves global smoothness.

The preservation/loss ledger of this eliminant is therefore precise:

```text
source       the two torus gradient equations;
target       R(z)=F'(z)^2-F(z)G(z)^2;
preserved    every critical point with XT!=0 and F(XT)!=0;
lost         which profile vanishes at a root of R;
sidecar      axis values and the multiplicities of psi;
decisive test R=-psi together with psi|psi' on its radical.          (18)
```

## 3. The generic fibre is hyperelliptic

Let `Lambda` be transcendental over `k` and work on the generic fibre
`Q_n=Lambda`.  Put

```text
Y=X-Tpsi_n(z)=X-F_n(z)/X.                             (19)
```

Since

```text
Lambda-chi_n(z)=X+F_n(z)/X,                            (20)
```

one obtains

```text
Y^2=Delta_n(z),
Delta_n(z)=[Lambda-chi_n(z)]^2-4F_n(z).                (21)
```

Conversely `X=[Lambda-chi_n+Y]/2` and `T=z/X`, so `(21)` is the generic
function field, not merely a quotient of it.

The polynomial `Delta_n` is squarefree in `k(Lambda)[z]`.  Indeed, if
`alpha` were a common root of `Delta_n` and its derivative, writing
`W=Lambda-chi_n(alpha)` would give

```text
W^2=4F_n(alpha),                 G_n(alpha)W+2F_n'(alpha)=0. (22)
```

Squaring the second equation and using the first forces
`F_n'^2-F_nG_n^2=0`, hence `psi_n(alpha)=0` by `(12)`.  Then
`F_n(alpha)=0`, so `(22)` gives `W=0`.  But `alpha` is algebraic over `k`
as a root of the fixed nonzero polynomial `psi_n`; the equality
`Lambda=chi_n(alpha)` contradicts the transcendence of `Lambda`.

Now `G_n` has degree `n`, so `chi_n` has degree `n+1`.  The coefficient
`-2Lambda chi_n` in `(21)` cannot cancel against a term independent of
`Lambda`; hence

```text
d_n=deg_z Delta_n >=n+1.                                (23)
```

The smooth projective model of `(21)` has genus
`floor((d_n-1)/2)`, which is at least `floor(n/2)`.  In fact the exact
controls through depth eight find `d_n=n+1`, but the lower bound `(23)` is
all that the proof needs and no all-depth equality is asserted here.

## 4. The mate equation is a forbidden exact differential

The other reason for choosing `Y` is the exact identity

```text
J(Q_n,z)=X(Q_n)_X-T(Q_n)_T=Y.                          (24)
```

In generic-fibre coordinates `(Lambda,z)`, a rational equation
`J(P,Q_n)=c in k*` therefore becomes

```text
partial P/partial z=-c/Y,
dP=-c dz/Y.                                            (25)
```

For `n>=2`, equation `(23)` gives `d_n>=3`.  On the smooth projective
hyperelliptic curve `(21)`, the differential `dz/Y` is nonzero and
holomorphic: at a finite branch point use `Y` as a local parameter, and the
standard infinity calculation for a squarefree polynomial of degree at
least three has no pole.  A derivative of a rational function cannot be a
nonzero holomorphic differential.  Any pole of the function would create a
pole of its derivative; without poles the function is constant.  Thus
`(25)` is impossible.

At `n=1`, the definitions specialize to

```text
G_1=4z-3,             S_1=z-1,
chi_1=2z^2-3z,        psi_1=(z-1)^3,                   (26)
```

and

```text
Delta_1=Lambda^2+(6Lambda+4)z-(4Lambda+3)z^2.          (27)
```

This is a squarefree conic over `k(Lambda)`.  Over its algebraic closure it
has two points at infinity, and `dz/Y` has opposite nonzero residues
`+/-1/sqrt[-(4Lambda+3)]` there.  A rational derivative has zero residue at
every point, so `(25)` is again impossible.

For comparison, the same boundary has the exceptional split fibre

```text
Q_1+1=[1+T(XT-1)][X+(XT-1)^2].                         (28)
```

Its two disjoint components make the logarithmic debt visible before
compactification.  Equations `(25)--(28)` prove that no member of the tower
has even a rational mate, and hence none is a coordinate.  **QED.**

## 5. Exact verification and search consequence

Reproduce the exact audit surface with

```bash
python3 -B 04-computation/jc2_pell_chebyshev_three_charge_hyperelliptic_thm3757.py
python3 -B -O 04-computation/jc2_pell_chebyshev_three_charge_hyperelliptic_thm3757.py
```

The assertion-free companion verifies `(1)--(24)` through eight Pell depths,
direct Groebner smoothness through depth four, squarefreeness and the exact
generic degrees `2,...,9`, the split conic boundary, and 23 bounded mate
systems in the first two depths.  These are hostile controls; the uniform
Chebyshev, resultant, and differential proofs establish every `n>=1`.

This tower sharpens the counterexample search in two ways.  First, a
nonconstant critical eliminant supported on the repeated profile locus can
indeed make all three charges smooth; constant-eliminant searches miss this.
Second, smoothness must be followed immediately by the generic-fibre
cohomology test.  Here the exact Pell cancellation that removes affine
critical points simultaneously exposes a nonzero `H^1` time-form class.  A
viable construction must cancel that class as well, before attempting to
clear its remaining boundary poles polynomially.
