---
id: THM-3229
title: "Hasse--Pluecker simple-root contact gcd flag and degree termination"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For two polynomials, the Hasse--Pluecker numerators
  Omega_m=f^[m]g'-g^[m]f' cut the derivative-normalized contact order at
  every common simple root, with the exact THM-3221 orientation
  c_m=-Omega_m/(f'g').  After squarefree gcd and simple-root saturation,
  successive gcds with Omega_m form a nested root-free divisor flag; the
  quotient G_(m-1)/G_m is exactly the first-contact-m stratum.  For
  nonproportional degrees bounded by D the flag terminates at G_D=1, and
  deg Omega_m<=deg(f)+deg(g)-m-1.  Hasse derivatives are load bearing in
  small characteristic.  The flag does not order/select roots and is not
  the resonant factorial-moment PRS.
source: root/multiscale-newton-flag/low-child-flag-extension/2026-08-02
audit: >
  The assertion-independent exact companion pins promoted THM-3221; builds
  one four-root pair with prescribed first contacts 2,3,4,5; verifies 52
  Hasse--Pluecker orientation identities, ten prescribed contact gates,
  thirteen exact nested gcd stages through degree fourteen, thirteen degree
  bounds and constant-scaling laws, the proportional boundary, the F3
  x versus x+x^3 Hasse-only third contact, and the multiple-root saturation
  boundary.  Normal/-O/stored replay agrees with the LF-normalized hashes
  below.  An independent immutable audit rederived the Hasse--Pluecker
  orientation and degree bound, the squarefree simple-root saturation, every
  nested gcd stratum, and degree-bound termination.  It also checked the
  proportional, multiple-root, and small-characteristic boundaries and
  accepted both exact replays and the pinned hashes.
depends_on:
  - THM-3221-selected-root-osculating-separation-and-minimal-jet-prime-carry
related:
  - THM-3215-arbitrary-degree-root-jet-hamiltonian-affine-dihedral-holonomy-and-p-fold-carry
  - THM-3217-universal-resonant-degree-prs-wall-atlas-and-fixed-offset-exception-set
  - THM-3167-inverse-different-three-gate-target-shear-descent-and-full-marked-jet-no-go
  - THM-2022-gmc2-frobenius-lowest-balanced-face
script: 04-computation/hasse_pluecker_contact_gcd_flag_thm3229.py
output: 05-knowledge/results/hasse_pluecker_contact_gcd_flag_thm3229.out
script_sha256: f8db4569812fcd9abc964088ade2ed787b59fc602edf949141ee07812754cb52
output_sha256: 1177bc28bf575dfdf835abc88a994733bc2644463f8f8960f717a39d6ed24b43
hash_basis: LF-normalized bytes
---

# THM-3229 -- Hasse--Pluecker simple-root contact gcd flag and degree termination

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3221 proves that a selected common simple root has a canonical first
separating normalized coefficient.  The missing elimination question is
whether one can see the **set of roots at each contact depth** without first
choosing and ordering those roots.  Hasse derivatives give an exact answer:
successive exterior derivative minors cut a finite nested gcd flag.

The flag is root-selection-free as a divisor.  It is not a root selector.

## 1. Hasse--Pluecker numerators

Let `K` be a perfect field, let `f,g in K[x]` be nonzero polynomials, and
write Hasse derivatives by

```text
f(x+u)=sum_(m>=0) f^[m](x)u^m,
g(x+u)=sum_(m>=0) g^[m](x)u^m.                           (1)
```

Thus `f^[1]=f'`; when `m!` is invertible,

```text
f^[m]=f^(m)/m!.                                          (2)
```

For `m>=2` define the Hasse--Pluecker polynomial

```text
Omega_m(f,g)=f^[m]g'-g^[m]f'.                            (3)
```

It is alternating in `f,g` and scales under constant target rescaling by

```text
Omega_m(lambda f,mu g)=lambda mu Omega_m(f,g).            (4)
```

If `d_f=deg(f)` and `d_g=deg(g)`, then

```text
deg Omega_m <= d_f+d_g-m-1,                              (5)
```

with the zero polynomial allowed when the right side is negative.  This is
immediate from `deg f^[m]<=d_f-m` and `deg g'<=d_g-1`.

## 2. Exact selected-root orientation

Let `a` be a common simple root:

```text
f(a)=g(a)=0,       f'(a)g'(a)!=0.                         (6)
```

The normalized germs are

```text
phi_f(u)=f(a+u)/f'(a),       phi_g(u)=g(a+u)/g'(a).       (7)
```

Their order-`m` coefficient difference is

```text
c_m(a)
 =[u^m]phi_g-[u^m]phi_f
 =g^[m](a)/g'(a)-f^[m](a)/f'(a)
 =-Omega_m(f,g)(a)/(f'(a)g'(a)).                         (8)
```

The minus sign matches THM-3221's transition convention
`H_(g<-f)=phi_g composed phi_f^(-1)`.  It also matches THM-3215 after noting
that its ordinary second Pluecker numerator is `2 Omega_2` whenever `2` is
invertible.

If

```text
Omega_2(a)=...=Omega_(m-1)(a)=0,
Omega_m(a)!=0,                                           (9)
```

then THM-3221 gives

```text
H_(g<-f)(u)=u+c_m(a)u^m mod u^(m+1).                     (10)
```

Thus `(3)` is a polynomial numerator for every selected-root first-live jet.

## 3. The simple-common-root divisor

Let

```text
C=sqfree(gcd(f,g)).                                       (11)
```

Remove every factor on which `f'g'` vanishes.  Equivalently, define the
monic squarefree polynomial `G_0` whose geometric roots are

```text
X_0={a in Kbar:f(a)=g(a)=0, f'(a)g'(a)!=0}.              (12)
```

In a chosen affine chart one may compute it by

```text
G_0=C/gcd(C,f'g').                                       (13)
```

Put `G_1=G_0`, and for `m>=2` define the nested monic gcd flag

```text
G_m=gcd(G_(m-1),Omega_m).                                (14)
```

Since `G_0` is squarefree, every `G_m` is squarefree and divides its
predecessor.

## 4. Exact contact strata

For every `m>=2`, the roots of `G_m` are exactly the simple common roots at
which the normalized germs agree through order `m`:

```text
roots(G_m)
 ={a in X_0:c_2(a)=...=c_m(a)=0}.                        (15)
```

This follows inductively from `(8)` because `f'(a)g'(a)` is invertible on
`X_0`.  Consequently

```text
roots(G_(m-1)/G_m)
 ={a in X_0:first normalized contact is exactly m}.      (16)
```

The quotient in `(16)` is squarefree.  The collection of these quotients is
a disjoint root-free divisor partition: it says how many roots have each
first contact and gives their defining polynomial, but it does not choose an
ordering or a section of those roots.

## 5. Degree-bound termination

Let

```text
D=max(deg f,deg g).                                       (17)
```

If `f,g` are not proportional, then

```text
G_D=1.                                                    (18)
```

Indeed, a root `a` of `G_D` would make every normalized coefficient through
degree `D` equal.  Since both normalized germs are polynomials of degree at
most `D`, coefficientwise equality gives

```text
f(a+u)/f'(a)=g(a+u)/g'(a)                                (19)
```

as polynomials.  Hence `g'(a)f=f'(a)g`, contradicting nonproportionality.

Conversely, if `g=lambda f`, then every `Omega_m` vanishes and

```text
G_m=G_0 for every m.                                      (20)
```

Thus proportionality is the exact nontermination boundary.  No uniform jet
depth independent of the supplied degrees is asserted.

## 6. Coordinate and elimination meaning

The polynomial presentation `(11)--(14)` uses the affine coordinate `x`, but
the root strata do not.  Under an invertible local source-coordinate change,
normalized coefficient differences transform triangularly; once the lower
differences vanish, the first live one scales by the weighted factor from
THM-3221.  Therefore the vanishing conditions `(15)` are intrinsic on the
simple-root locus.

Because all objects in `(14)` are univariate polynomials, the flag is
computable by squarefree factorization followed by any exact gcd or
subresultant/PRS implementation.  This is a conceptual elimination target
for higher root jets.  It is **not** THM-3217's resonant factorial-moment PRS:
that theorem applies pseudo-division to reciprocal moment jets on a degree
line, whereas `(14)` applies gcd elimination to derivative-normalized
coefficient germs.  A bridge must identify the two input pairs before their
walls can be compared.

## 7. Small characteristic makes Hasse derivatives essential

Over `F_3`, take

```text
f(x)=x,       g(x)=x+x^3.                                (21)
```

At the common simple root zero, the normalized germs first differ at order
three.  Hasse derivatives give

```text
Omega_2(0)=0,       Omega_3(0)=-1=2 in F_3,
c_3(0)=1=-Omega_3(0).                                    (22)
```

But the ordinary third derivative of `x^3` is `3!=0` in `F_3`.  Replacing
Hasse derivatives with ordinary derivatives would erase the exact contact
which carries the characteristic-three reset.

## 8. Sharp boundaries

1. **Multiple roots.**  If `f'(a)g'(a)=0`, normalization `(7)` is illegal.
   The saturation in `(12),(13)` is load bearing; multiple-root strata need
   their own divided-power or homogeneous treatment.
2. **Infinity.**  A root at infinity requires a homogeneous chart.  The
   affine polynomial `G_m` does not silently include it.
3. **Root ordering.**  The divisor flag gives finite root sets, not a
   canonical root label, monodromy section, or owner.
4. **Nonconstant local units.**  Constant target scalings preserve the flag;
   multiplying one polynomial by a nonconstant unit generally changes its
   normalized contact tower.
5. **Physical survival.**  No Newton/Wick channel, LRC carrier, or Jacobian
   owner is proved to retain a selected quotient `G_(m-1)/G_m`.

## 9. Frontier consequences

### Gaussian moments

THM-2022 already proves NC2/GMC(2).  If a future Gaussian-moment construction
produces a pair of radial coefficient polynomials and a common-simple-root
condition, `(14)` gives a finite, degree-adaptive way to partition all such
roots by first live jet.  It does not prove that any one stratum survives the
whole moment sum.

### Jacobian and PRS work

The flag is compatible with bounded-degree algebraic elimination, but it
does not contradict THM-3167's fixed-finite-jet global-owner no-go: the depth
is degree dependent, multiple roots are removed, and no polynomial owner is
selected.  Its most concrete possible use is as a target specification for
a fraction-free PRS atlas: a proposed scalar wall should be checked against
the elimination divisor of one `Omega_m`, not merely matched by degree or
numerical roots.

## 10. Connection contract

```text
source:      two finite-degree polynomials in one affine chart;
map:         squarefree simple gcd, then successive Hasse--Pluecker gcds;
target:      nested contact divisors G_m and exact first-contact strata;
preserved:   simple-root set, normalized contact depth, constant scaling;
destroyed:   root ordering, branch monodromy, nonconstant-unit provenance;
sidecar:     root selector/owner plus a physical whole-layer survivor map.
```

## 11. Exact companion

The assertion-independent companion

```text
04-computation/hasse_pluecker_contact_gcd_flag_thm3229.py
```

pins promoted THM-3221 and checks exactly:

```text
one four-root bank with prescribed contacts 2,3,4,5;
52 Hasse--Pluecker orientation identities;
10 lower/equality contact gates;
13 nested gcd stages and termination at D=14;
13 polynomial degree bounds and 13 constant-scaling laws;
the proportional, multiple-root, and F3 Hasse-only boundaries. (23)
```

Ordinary and optimized runs byte-match

```text
05-knowledge/results/hasse_pluecker_contact_gcd_flag_thm3229.out
```

and the LF-normalized hashes are pinned in the frontmatter.

QED.
