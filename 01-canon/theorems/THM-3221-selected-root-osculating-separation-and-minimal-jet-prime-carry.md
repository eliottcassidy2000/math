---
id: THM-3221
title: "Selected-root osculating separation and minimal-jet prime carry"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  At one selected common simple root, the first unequal coefficient of two
  derivative-normalized polynomial germs is exactly the first nonidentity
  coefficient of their transition germ.  It is a weighted cotangent tensor
  under arbitrary nonlinear source-coordinate change.  Its minimal jet
  quotient is additive, has exact order p in residue characteristic p, and
  returns as the undivided coefficient after division by p.  A root,
  derivative anchor, and the normalized tower through the supplied degree
  reconstruct a finite-degree polynomial exactly.  No global root selector,
  fixed degree-independent jet, or whole-layer moment survivor is supplied.
source: root/multiscale-newton-flag/low-child-flag-extension/2026-08-02
audit: >
  The assertion-independent exact companion pins promoted THM-3215; checks
  27 common-prefix transition inversions and 27 genuinely nonlinear source
  charts; checks 44 minimal-quotient compositions and coordinate weights,
  561 nonzero exact-order-p classes, 231 divided carries, 30 finite-degree
  reconstructions, ten fixed-depth hostiles, the same-divisor nonconstant-unit
  hostile, and the F3 higher-Nottingham-tail hostile.  An independent hostile
  audit rederived every mechanism and scope boundary.  Normal/-O/stored
  replay passes with LF-normalized hashes below.
depends_on:
  - THM-3215-arbitrary-degree-root-jet-hamiltonian-affine-dihedral-holonomy-and-p-fold-carry
related:
  - THM-3220-root-four-jet-schwarzian-heisenberg-transgression-and-oriented-discriminant-holonomy
  - THM-3167-inverse-different-three-gate-target-shear-descent-and-full-marked-jet-no-go
  - THM-2022-gmc2-frobenius-lowest-balanced-face
  - THM-2591-theta-zero-selector-cech-coboundary-and-c91-holonomy-no-go
script: 04-computation/selected_root_osculating_separation_prime_carry_thm3221.py
output: 05-knowledge/results/selected_root_osculating_separation_prime_carry_thm3221.out
script_sha256: 415c09c6e6e6bdf3f2e5a9b07987418698b9915a9b9b2ab8d8131246580d8312
output_sha256: b9c046b58fd5bc26fa20092388d6c561ff3537e56a347d8a151c0fcaa352b31b
hash_basis: LF-normalized bytes
---

# THM-3221 -- selected-root osculating separation and minimal-jet prime carry

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3215 extracts the universal two-jet reflection at a selected simple root,
and THM-3220 identifies the first simultaneous nonabelian four-jet layer.  The
degree-independent statement between them is not that one fixed truncation
sees every polynomial.  It is that the **first live osculating layer** is
canonical, additive, and prime-resettable, while the depth required for full
separation grows with the supplied degree.

This distinction is the point of the theorem.

## 1. Selected-root normalization

Let `A` be a commutative ring.  Let `f,g` be polynomials or formal germs with
a selected common root `a`, and suppose `f'(a),g'(a)` are units.  In a local
source coordinate `u`, define

```text
phi_f(u)=f(a+u)/f'(a),       phi_g(u)=g(a+u)/g'(a).       (1)
```

Both series are tangent to the identity.  Their selected-root transition is

```text
H_(g<-f)=phi_g composed phi_f^(-1).                       (2)
```

Suppose `m>=2` is their first unequal normalized coefficient.  Thus, modulo
`u^(m+1)`,

```text
phi_f=P+a_m u^m,       phi_g=P+b_m u^m,
P=u+O(u^2),            deg(P)<m,                         (3)
```

and `c_m=b_m-a_m` is nonzero.  Then

```text
H_(g<-f)(u)=u+c_m u^m mod u^(m+1).                       (4)
```

Indeed, `phi_g=phi_f+c_m u^m mod u^(m+1)`.  Substitution of
`phi_f^(-1)=u+O(u^2)` gives

```text
phi_g composed phi_f^(-1)
 =u+c_m(phi_f^(-1))^m
 =u+c_m u^m                    mod u^(m+1).               (5)
```

No lower coefficient appears in `(4)`, and no division beyond the two unit
first derivatives is used.

## 2. Nonlinear coordinate covariance

Let `v` be another source coordinate.  Write

```text
v=psi(u),       psi'(0)=rho in A^*,       u=chi(v).       (6)
```

The derivative normalization in `(1)` gives the exact identities

```text
phi_f^v(v)=rho phi_f^u(chi(v)),
phi_g^v(v)=rho phi_g^u(chi(v)).                            (7)
```

Although `chi` may be nonlinear, it cancels from the transition:

```text
H^v=L_rho composed H^u composed L_rho^(-1),
L_rho(u)=rho u.                                           (8)
```

Consequently the first separating coefficient obeys

```text
c_m^v=rho^(1-m)c_m^u.                                    (9)
```

Thus

```text
c_m (du)^(m-1)                                            (10)
```

is independent of the chosen nonlinear source chart.  It is the intrinsic
first osculating defect on the selected simple-root branch.  Multiplying an
entire germ by a constant unit does not change `(1)`; a nonconstant unit
generally does.

## 3. The minimal layer is additive

Let `F^m` denote tangent-to-identity germs congruent to `u` modulo `u^m`.
The coefficient of `u^m` identifies the associated-graded quotient with the
additive group of `A`:

```text
F^m/F^(m+1)  ->  (A,+),
u+c u^m      |-> c.                                      (11)
```

Explicitly,

```text
(u+c u^m) composed (u+d u^m)
 =u+(c+d)u^m                      mod u^(m+1),

(u+c u^m)^n
 =u+n c u^m                      mod u^(m+1).             (12)
```

This is the abelian first-live-layer shadow of the nonabelian four-jet group
in THM-3220.  Retaining several live layers simultaneously creates the
Heisenberg cross term; isolating the first nonzero layer removes it exactly.

## 4. Prime reset and divided carry

Let `R` be a mixed-characteristic DVR with residue characteristic `p`, and
let `H in F^m(R)` have first coefficient `c_m` a unit.  Reducing `(12)` modulo
the maximal ideal shows that the class of `H` in

```text
(F^m/F^(m+1))(k)                                         (13)
```

has exact order `p`: none of the iterates `1,...,p-1` kills its nonzero
additive coefficient, while the `p`-th iterate does.  Before reduction,

```text
[u^m](H^p-u)=p c_m,

(1/p)[u^m](H^p-u)=c_m.                                   (14)
```

so the first divided return is again a unit modulo `p`.

Equation `(14)` is a statement only in the minimal quotient.  It requires no
nilpotence-class-`<p` hypothesis, but it makes no claim that the full
untruncated germ has order `p`.

## 5. Degree-adaptive finite separation

Now let `A` be a characteristic-zero field, or any ring where the displayed
Taylor coefficients are defined, and suppose `f` is a polynomial of degree
at most `D`.  Write

```text
phi_f(u)=u+a_2u^2+...+a_Du^D.                             (15)
```

Then the selected root, the derivative anchor `d=f'(a)`, and the normalized
tower through degree `D` reconstruct `f` exactly:

```text
f(x)=d((x-a)+a_2(x-a)^2+...+a_D(x-a)^D).                 (16)
```

Hence two degree-`<=D` polynomials with the same selected root and the same
normalized tower through `D` are scalar multiples.  If their derivative
anchors also agree, they are equal.

This is a genuine arbitrary-finite-degree mechanism: the required depth is
the supplied degree, not a universal constant.  It therefore gives local
separation for arbitrary radial coefficient polynomials once a simple root
branch and derivative anchor have been chosen.

## 6. Sharp boundaries

### 6.1 Higher Nottingham tails survive the reset

Over `F_3`, let

```text
H(u)=u+u^2.
```

Then

```text
H^3(u)=u                  mod u^3,
H^3(u)=u+u^5              mod u^6.                       (17)
```

Thus the minimal visible layer resets exactly while a higher resonant tail
survives.  No full-germ exponent claim may be read into `(14)`.

### 6.2 No uniform fixed depth

For every `D`, the germs

```text
u,       u+u^(D+1)                                        (18)
```

agree through depth `D` and separate at the next coefficient.  Therefore no
fixed finite jet identifies arbitrary unbounded-degree polynomials or formal
germs.  The **full** infinite formal tower does identify a normalized formal
germ coefficientwise; what fails is any one finite uniform truncation.

### 6.3 The derivative anchor is load bearing

The normalized tower is unchanged by constant scalar multiplication.  It
recovers the scalar-proportionality class, not an absolute polynomial, until
`f'(a)` is supplied.  Conversely, the same zero divisor is not enough: `u`
and `u(1+u)` have the same local zero divisor and the same derivative at zero,
but their normalized second jets differ.

### 6.4 Locality is not ownership

Multiple roots do not admit normalization `(1)`, and a root at infinity needs
a homogeneous chart.  More importantly, this theorem does not choose a root
consistently across Newton/Wick channels.  It preserves a branch only after
that branch has been supplied.

## 7. Interaction with the live frontiers

### Gaussian moments

THM-2022 already proves NC2 and GMC(2) for arbitrary finite exact support by
the lowest balanced face and Frobenius amplification.  The present theorem is
not a competing proof.  It answers a narrower structural question left by
the early termwise-factorial approaches: arbitrary radial coefficient
polynomials can be separated **locally** by a degree-adaptive osculating
tower.  To turn that local fact into a moment proof one would still need both
a canonical common-root selection and a whole-layer noncancellation theorem.

### Jacobian finite-jet no-go

There is no conflict with THM-3167.  That theorem rules out a uniform fixed
finite-jet decision of global polynomial ownership on its supplied hostile
class.  Equation `(16)` instead assumes a selected simple root, retains a
derivative anchor, and uses a truncation depth growing with the degree.  It
can serve bounded-degree strata without deciding the global owner gate.

### Holonomy

THM-3220 shows that the simultaneous four-jet connection is a flat
Heisenberg extension with a central area transgression.  THM-3221 supplies
the complementary filtration view: at the first order where two branches
separate, the holonomy lands in one weighted additive line and has an exact
prime reset.  A branch seam can obstruct global gluing; the local associated
graded class itself is flat.

## 8. Connection contract

```text
source:      finite-degree germs at one selected common simple root;
map:         first nonidentity normalized transition coefficient c_m;
target:      (T_a^*)^(m-1), with additive prime reset and divided carry;
preserved:   selected branch, derivative normalization, first live order;
destroyed:   global root choice, channel/owner provenance, higher tail;
sidecar:     root selection plus whole-layer noncancellation or owner routing.
```

## 9. Exact companion

The companion

```text
04-computation/selected_root_osculating_separation_prime_carry_thm3221.py
```

is assertion-independent and uses only exact rational/integer arithmetic.  It
pins promoted THM-3215 and verifies:

```text
27 common-prefix transition inversions;
27 nonlinear source-coordinate cancellations;
44 minimal-layer compositions and 44 coordinate weights;
561 nonzero exact-order-p classes;
231 divided carries;
30 finite-degree reconstructions;
10 fixed-depth hostiles;
the nonconstant-unit and F3 higher-tail boundaries.        (19)
```

Ordinary and optimized runs byte-match

```text
05-knowledge/results/selected_root_osculating_separation_prime_carry_thm3221.out
```

and the LF-normalized hashes are pinned in the frontmatter.

An independent hostile audit separately rederived the first unequal
coefficient formula, cancellation of an arbitrary nonlinear source chart,
the cotangent weight, associated-graded additivity, exact quotient order and
divided carry, degree-adaptive reconstruction, and every stated hostile.  It
confirmed that `(14)` is deliberately quotient-only and that neither a full
Nottingham exponent nor a global root selector is smuggled into the result.
Fresh normal and optimized runs each byte-match the stored transcript and the
declared LF-normalized hashes.

QED.
