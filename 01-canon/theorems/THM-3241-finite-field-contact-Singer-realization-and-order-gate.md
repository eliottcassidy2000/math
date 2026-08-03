---
id: THM-3241
title: "Finite-field contact Singer realization and multiplicative-order gate"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Every prescribed nonzero element of an irreducible finite-field contact
  algebra can be realized as the uniform first-live normalized coefficient
  of an explicit polynomial pair.  Multiplication by that coefficient is a
  Singer carrier exactly when the usual prime-divisor order tests pass; each
  test has a root-free resultant form.  Over F13, S=x^2-2 and
  g=S+S^2(1+10x) realize c_2=1+2x, whose multiplication matrix is exactly the
  order-168 Singer matrix of THM-3234.  The element alpha^49 has the same
  primitive norm and still generates F169, but has order 24, sharply proving
  that norm and contact-spectrum generation do not imply Singer transitivity.
source: root/multiscale-newton-flag/low-child-flag-extension/2026-08-02
audit: >
  The assertion-independent exact companion verifies eight prescribed-unit
  contact realizations in irreducible degree-two and degree-three algebras;
  enumerates the full 169-class deformation atlas, all 169 Singer
  intertwining identities, and, as a scoped formal pullback control, all
  371,293 affine-Heisenberg state/action pairs;
  checks the explicit p=13 helper, contact
  numerator/denominator, order 168,
  norm 6 of order 12, all three prime-divisor resultant gates, and the exact
  THM-3234 matrix; and proves that alpha^49 has the same norm, generates F169,
  has order 24, and gives seven punctured orbits.  Normal and optimized runs
  byte-match the stored transcript and LF-normalized hashes below.  An
  independent hostile audit rederived the prescribed-contact formula and
  deformation isomorphism, the resultant order gate, all explicit F13
  arithmetic and the order-24 same-norm hostile, and confirmed that the
  Heisenberg pullback is only a formal deformation reparametrization.
depends_on:
  - THM-3232-root-free-contact-stratum-norm-and-discriminant-power
  - THM-3234-singer-owner-compactification-and-pointed-heisenberg-carrier-gate
related:
  - THM-3228-four-jet-heisenberg-minimal-faithful-permutation-carrier-gate
  - THM-3227-selected-root-residue-contact-trie-primitive-carry-and-delayed-resplitting
  - THM-2022-gmc2-frobenius-lowest-balanced-face
script: 04-computation/finite_field_contact_singer_gate_thm3241.py
output: 05-knowledge/results/finite_field_contact_singer_gate_thm3241.out
script_sha256: a4de526002b926caaf1f184fddb3e73e2c79bc06ebc6c914ff2e4925d02688c7
output_sha256: e12ee1ec196ab7ae64bc8ae3bde21814636f419089b883deefbedc7550a8b103
hash_basis: LF-normalized bytes
---

# THM-3241 -- finite-field contact Singer realization and multiplicative-order gate

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3232 turns first-live root jets into units of finite etale contact
algebras.  THM-3234 independently identifies a Singer multiplier on the
`169`-point affine plane as the exact abstract carrier compatible with its
`168`-owner cycle.  The missing algebraic question is whether a contact jet
can itself be such a multiplier.

It can.  In fact every nonzero field element can occur.  Multiplicative
primitivity is an additional and exactly testable gate.

## 1. Prescribed contact-unit realization

Let

```text
K=F_q,
S in K[x] monic irreducible of degree r,
A=K[x]/(S)=F_(q^r),                                      (1)
```

and fix `m>=2` and `a in A^*`.  Since `S` is separable, `S'` is a unit in
`A`.  Let `H in K[x]`, of degree less than `r`, be the unique representative
of

```text
H == a/(S')^(m-1)                              in A.      (2)
```

Define

```text
f=S,                    g=S+S^m H.                       (3)
```

Then the complete simple common-root divisor is `S`, every one of its roots
has first normalized contact exactly `m`, and the THM-3232 contact class is

```text
c_m=a in A^*.                                             (4)
```

### Proof

The common gcd in `(3)` is exactly `S`, because

```text
g=S(1+S^(m-1)H).                                         (5)
```

At a root `z` of `S`, both derivatives equal `S'(z)`.  The difference
`g-f=S^mH` begins at order `m`, and its order-`m` Hasse coefficient is

```text
(S'(z))^m H(z).                                          (6)
```

After derivative normalization this becomes

```text
(S'(z))^(m-1)H(z)=a(z),                                  (7)
```

while all lower differences vanish.  Since `a` is a unit, `(7)` is nonzero
at every conjugate root.  This proves `(4)` and exact first contact.

Thus there is no algebraic realizability restriction on a nonzero contact
coefficient inside an irreducible finite-field stratum.

## 2. The complete truncated deformation carrier

Fix `f=S` and retain only the order-`m` transverse deformation modulo
`S^(m+1)`.  Its classes are

```text
D_m={g_H=S+S^mH mod S^(m+1):H in A}.                    (8a)
```

The calculation `(6)--(7)` says that the contact map is

```text
D_m -> A,
g_H |-> c_m=(S')^(m-1)H.                                (8b)
```

Because `S'` is a unit, `(8b)` is an `F_q`-linear isomorphism.  Thus the
deformation slice has exactly `q^r` states:

```text
H=0       <-> contact is delayed beyond m,
H!=0      <-> first contact is exactly m.                (8c)
```

Multiplication by any `a in A^*` intertwines across `(8b)`: multiplying `H`
by `a` multiplies `c_m` by `a`.  Hence a primitive `a` makes the `q^r-1`
exact-contact deformation classes one Singer orbit and fixes the unique
delayed-contact class.  This is a genuine algebraic deformation carrier,
not merely the underlying set of the field viewed after the fact.

## 3. The Singer order gate

Multiplication by `a` is an `F_q`-linear permutation of `A`.  It fixes zero,
and every orbit in `A^*` has length `ord(a)`.  Put

```text
N=q^r-1.                                                  (8)
```

Then multiplication by `a` is a Singer carrier--one fixed zero and one
regular orbit on all nonzero vectors--if and only if

```text
ord(a)=N.                                                 (9)
```

Equivalently, for every prime `ell` dividing `N`,

```text
a^(N/ell)!=1.                                            (10)
```

For an arbitrary THM-3232 presentation

```text
a=n/d,             n=-Omega_m,          d=f'g'           (11)
```

in the irreducible algebra `A`, `(10)` has the root-free resultant form

```text
Res_x(S,n^(N/ell)-d^(N/ell)) !=0
                 for every prime ell|N.                  (12)
```

Indeed, `d` is a unit and irreducibility of `S` makes the resultant in `(12)`
zero exactly when the displayed class vanishes in `A`.  Powers may be reduced
modulo `S` before taking the resultant.

The ordinary norm is only the projection

```text
Norm_(A/F_q)(a)=a^[N/(q-1)] in F_q^*.                    (13)
```

If `a` is Singer, `(13)` has order `q-1`; the converse need not hold.

## 4. Exact p=13 realization of the THM-3234 carrier

Work over `F_13` with

```text
S=x^2-2,                 A=F_169,
alpha=1+2x.                                                 (14)
```

The element `2` is a nonsquare modulo `13`, so `S` is irreducible.  At
contact order two,

```text
H=alpha/S'=alpha/(2x)=1+10x.                             (15)
```

Therefore

```text
f=S,
g=S+S^2(1+10x)                                           (16)
```

have the uniform first-contact class

```text
c_2=alpha=1+2x.                                          (17)
```

In the basis `(1,x)`, multiplication by `alpha` is

```text
[1 4]
[2 1],                                                    (18)
```

because `x^2=2`.  This is exactly the order-`168` matrix in THM-3234, not
merely a group of the same cardinality.  The three prime-divisor tests for

```text
168=2^3*3*7                                               (19)
```

use exponents `84,56,24`; the reduced gate classes and their resultants are

```text
ell=2:  11,          resultant 4,
ell=3:   8,          resultant 12,
ell=7: 9+11x,        resultant 8.                        (20)
```

All are nonzero, so `ord(alpha)=168`.  Its norm is

```text
Norm(alpha)=6,                 ord_(F_13^*)(6)=12,        (21)
```

matching the determinant of `(18)` and THM-3234's scalar gauge.

Consequently the contact algebra itself supplies the exact abstract
decomposition

```text
A={0} disjoint union A^*,
|A^*|=168,              alpha acts regularly on A^*.     (22)
```

Under `(8b)`, this is exactly one delayed-contact deformation class plus all
`168` exact-order-two deformation classes.  The companion enumerates all
`169` classes, all `169` Singer intertwining identities, and all
`13^5=371,293` state/group-element Heisenberg pullbacks.

## 5. Sharp hierarchy hostile

Let

```text
beta=alpha^49=7x.                                        (23)
```

Then

```text
Norm(beta)=6=Norm(alpha),
beta notin F_13,
F_13[beta]=F_169,                                        (24)
```

but

```text
ord(beta)=24.                                            (25)
```

Thus `beta` is nonzero, its norm is even a primitive base-field unit, and it
passes the full algebra-generation/contact-spectrum gate, yet multiplication
by `beta` has

```text
168/24=7                                                 (26)
```

punctured orbits instead of one.  It too is realized by `(3)`: at `m=2`, its
helper is the constant `H=10`.

Section 1 with prescribed `a=1` in any degree `r>1` proves that a nonzero
contact unit need not generate the contact field.  Equations `(23)--(26)`
prove the second strict separation.  Hence the hierarchy

```text
nonzero contact unit
  < generator of the contact field
  < multiplicatively primitive Singer element             (27)
```

is strict.  Neither the THM-3232 norm nor the open root-reconstruction gate
can replace the order tests `(10)--(12)`.

## 6. Relation to the Heisenberg carrier

Once a basis of `A` is chosen, its additive group is `F_q^r`.  For `q=13`
and `r=2`, `(18)` supplies the Singer part of THM-3234's
Singer--Heisenberg construction.  Combining it with the standard affine
Heisenberg action then generates `AGL_2(F_13)` exactly as proved there.

The new content here is the polynomial-contact deformation realization of
the Singer carrier.  There is one useful but strictly formal related
observation.  After choosing a basis and the derivative trivialization
`(8b)`, any `F_q`-affine map `L:A->A` can be transported to the coefficient
reparametrization

```text
H |-> (S')^(1-m) L((S')^(m-1)H).                         (28)
```

The companion checks all `371,293` state/action pairs for THM-3234's affine
Heisenberg formulas.  This does **not** supply THM-3240's exact-address
clutch: no LRC owner, endpoint, or root packet is identified with an element
of `A`, and `(28)` is not an intrinsic transformation of one fixed physical
pair.  Hence `(22)` remains an algebraic deformation-moduli carrier, not a
lawful physical owner map.  The delayed-contact zero is fixed by the Singer
action but moved by formal Heisenberg translations, matching THM-3234's
fixed-head obstruction only at the abstract carrier level.

## 7. Consequences and boundaries

1. **Irreducibility is load bearing.**  For a reducible squarefree `S`, the
   contact algebra is a product of fields.  Multiplication can be tested
   componentwise, but there is no single punctured Singer orbit of size
   `q^r-1`.
2. **The prescribed value must be a unit.**  If `a=0`, the order-`m`
   coefficient vanishes and first contact occurs later.
3. **Norm is insufficient.**  The exact `alpha/beta` pair has identical
   primitive norm but different orbit decompositions.
4. **Field generation is insufficient.**  `beta` generates the full degree-
   two field and has a squarefree contact spectrum, but is not Singer.
5. **No characteristic-zero lift is automatic.**  The finite order and
   carrier live in the residue algebra; a Hensel or monodromy lift needs a
   separate theorem.
6. **No physical LRC map is supplied.**  The realization `(16)` is a new
   polynomial pair over `F_13`, not an extracted canonical LRC packet.

For Gaussian-moment work, the realization theorem is also a no-go against
an algebra-only shortcut: arbitrary coefficient data can realize subgroup
or Singer behavior at will, so nonzero Frobenius norm does not force maximal
orbit mixing.  THM-2022 already proves GMC(2); no new arbitrary-radial
channel statement is asserted.

## 8. Connection contract

```text
source:      an irreducible finite-field contact algebra and prescribed unit;
map:         g=S+S^m a/(S')^(m-1), then multiplication by c_m;
target:      an explicit linear carrier on the q^r algebra elements;
preserved:   contact depth, prescribed unit, norm, exact multiplicative order;
destroyed:   root ordering and all physical LRC/JC provenance;
sidecar:     lawful identification of packet states plus compatible affine
             Heisenberg operations from the same source geometry.
```

## 9. Exact companion

The assertion-independent companion

```text
04-computation/finite_field_contact_singer_gate_thm3241.py
```

uses exact finite-field polynomial arithmetic.  It checks eight prescribed
contact realizations over `F_5,F_7,F_13` in degrees two and three, enumerates
the full `169`-class deformation atlas and `371,293` Heisenberg pullbacks,
reconstructs `(15)--(18)`, verifies all order/resultant/norm claims, and
compares the Singer element with the exact `beta` hostile.  Normal and
optimized runs
byte-match

```text
05-knowledge/results/finite_field_contact_singer_gate_thm3241.out
```

and the LF-normalized hashes are pinned in the frontmatter.

QED.
