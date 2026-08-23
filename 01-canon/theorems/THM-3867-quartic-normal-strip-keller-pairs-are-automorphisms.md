---
id: THM-3867
title: "Quartic normal-strip Keller pairs are automorphisms"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Over every
  characteristic-zero field, a polynomial Keller pair in k[s,z] whose two
  transverse z-degrees are at most four is a polynomial automorphism.  After
  a constant target change, the (4,1) and (4,2) rows reduce by target shears
  to the cubic theorem.  The only genuinely new (4,3) row admits a complete
  depressed Kummer packet.  Its three local denominator channels all fail
  because canceling the C-arm pole leaves a nonzero A-arm residual.  Constant
  Kummer scale is contradictory as well.  Arbitrary polynomial Keller pairs,
  rational or infinite normal series, and the planar Jacobian conjecture
  remain open.
source: root / quartic_normal_scout / planar JC quartic strip lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE PROOF/REPLAY AUDIT PASSED.  The independent companion
  reconstructs all eight buckets by coefficient convolution, both SL2 target
  charts, every q/v/beta zero stratum, both degree-drop shears, the h^4/h^3
  UFD packet, the depressed integrations, and the load-bearing outer-B factor
  in the constant row.  Its 4,872 gates cover 1,152 simultaneous-polynomiality
  dominance profiles, W-pole/W-zero/W-identically-zero and constant-h strata,
  nonclosed-field and irreducible-DVR edges, both nonzero terminal residues,
  and two rational hostiles.  Normal and optimized executions of both
  companions byte-match their frozen LF outputs.
depends_on:
  - THM-3861-cubic-normal-strip-keller-pairs-are-automorphisms
related:
  - THM-3856-quadratic-normal-strip-keller-pairs-are-automorphisms
  - THM-3843-russell-arm-birational-immersion-and-forced-self-identification
  - THM-3846-formal-arm-darboux-lift-and-algebraization-gate
  - THM-3860-russell-higher-normal-rational-lifts-and-vertical-pole-barrier
script: 04-computation/jc2_quartic_normal_strip_keller_thm3867.py
output: 05-knowledge/results/jc2_quartic_normal_strip_keller_thm3867.out
script_sha256: 37b6fa92c3063f281eef2d2cb3075d1dfd3400c700365d21dffedf71bf2d06d1
output_sha256: 966a1da1b5d20847b3d2aac2be087814d29c12efa9cca4cf305966b0f7402b91
semantic_sha256: 3b14427454f046f289065ff37da208ff731cf7b8e699af0a550568ceceea4380
independent_script: 04-computation/jc2_quartic_normal_strip_keller_independent_audit_thm3867.py
independent_output: 05-knowledge/results/jc2_quartic_normal_strip_keller_independent_audit_thm3867.out
independent_script_sha256: 6b9b84fdea04ad66ba06f977d80aa2a29c5ead05cc0b6e1094a786f86046300a
independent_output_sha256: 7c50c4f9f3f28f50e04c613f5aecbbd17925e347038386f811e073c8fc663c9a
independent_semantic_sha256: 737ccc63a339f0defde8e0d4995b2582529821f19fb194bccf691a57a69059f0
hash_basis: raw LF bytes
---

# THM-3867 -- quartic transverse depth is still triangular

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Let `k` be a
field of characteristic zero.  Suppose

```text
A,C in k[s,z],                  deg_z A,deg_z C <=4,
J_(z,s)(A,C)=A_z C_s-A_s C_z=lambda in k*.                    (1)
```

Then `(A,C)` is a polynomial automorphism.

More precisely, after a constant target `SL_2(k)` change, every genuinely
quartic pair lies in exactly one of the following rows.

1. In row `(4,1)`, there are `rho,beta in k*`, `d_0,...,d_3 in k`, and
   `b in k[s]` such that

   ```text
   C=b+beta z,
   A=rho C^4+d_3 C^3+d_2 C^2+d_1 C-(lambda/beta)s+d_0.       (2)
   ```

2. In row `(4,2)`, there are `rho,sigma in k` and a polynomial `L` with
   `deg_z L<=1` such that

   ```text
   A=rho C^2+sigma C+L,                                      (3)
   ```

   and `(L,C)` is a quadratic normal-strip Keller automorphism.

3. The apparent row `(4,3)` is empty in `k[s,z]`.

The proof of the last assertion is the new content.  Its obstruction is not
a single uncancelled principal part.  There are genuine cancellation
channels, but the two arm coordinates demand incompatible leading
cancellations.

## 1. The eight exact coefficient buckets

Write

```text
A=a+alpha z+u z^2+p z^3+r z^4,
C=b+beta z+v z^2+q z^3+t z^4,                               (4)
```

with all ten coefficients in `k[s]`.  For `m=0,...,7`, the coefficient of
`z^m` in the Jacobian is

```text
E_m=sum_(i+j=m+1) (i a_i b_j'-j a_i' b_j).                  (5)
```

In expanded form, the eight buckets are exactly

```text
E_0 = alpha b'-a'beta,                                      (6)

E_1 = alpha beta'-alpha'beta+2u b'-2a'v,                    (7)

E_2 = alpha v'-2alpha'v+2u beta'-u'beta+3p b'-3a'q,         (8)

E_3 = alpha q'-3alpha'q+2u v'-2u'v+3p beta'-p'beta
      +4r b'-4a't,                                          (9)

E_4 = alpha t'-4alpha't+2u q'-3u'q+3p v'-2p'v
      +4r beta'-r'beta,                                    (10)

E_5 = 2u t'-4u't+3p q'-3p'q+4r v'-2r'v,                   (11)

E_6 = 3p t'-4p't+4r q'-3r'q,                              (12)

E_7 = 4(rt'-r't).                                          (13)
```

Thus `E_0=lambda` and `E_1=...=E_7=0`.  These are identities in the
ordinary polynomial ring; no completion or saturation is being used.

If `(r,t)=(0,0)`, THM-3861 applies.  Otherwise `(13)` says that `(r,t)`
has constant direction over `k`: the constant field of `k(s)` is `k`.
A constant target `SL_2(k)` matrix therefore normalizes

```text
t=0,                         r!=0.                           (14)
```

## 2. The `(4,1)` and `(4,2)` rows are removable

First suppose `q=0`.

If `v=0` and `beta=0`, equation `(9)` reduces to `4r b'=0`; then `(6)`
cannot equal `lambda`.  Hence the `(4,1)` row has `beta!=0`.  Equation
`(10)` becomes

```text
4r beta'-r'beta=0,
```

so `(r/beta^4)'=0` and

```text
r=rho beta^4,                  rho in k*.                    (15)
```

The target shear

```text
A_1=A-rho C^4                                                     (16)
```

preserves the Jacobian and lowers `deg_z A_1` to at most three.  Applying
THM-3861 to a surviving cubic row, and THM-3856 on its lower-degree edges,
gives exactly `(2)` and its polynomial inverse.

Now suppose `v!=0`.  Equation `(11)` is

```text
4r v'-2r'v=0,
```

whence

```text
r=rho v^2,                      rho in k*.                    (17)
```

The target shear

```text
A_1=A-rho C^2                                                     (18)
```

again preserves the Jacobian and lowers `deg_z A_1` to at most three.
THM-3861 excludes a genuine `(3,2)` survivor, so in fact
`deg_z A_1<=2`.  The top quadratic Wronskian then makes the `z^2`
coefficient of `A_1` a constant multiple `sigma` of `v`.  Subtracting
`sigma C` leaves the coordinate `L` in `(3)`, and THM-3856 proves that
`(L,C)` is an automorphism.  This exhausts `q=0`, including all zero
components.

It remains to assume

```text
q!=0.                                                             (19)
```

## 3. Kummer scaling and depression of the `(4,3)` row

Equation `(12)` now reads

```text
4r q'-3r'q=0.                                              (20)
```

At each irreducible prime, `(20)` gives
`3 ord(r)=4 ord(q)`.  Unique factorization therefore supplies

```text
r=R h^4,                      q=Q h^3,
R,Q in k*,                    0!=h in k[s].                 (21)
```

Put `y=h(s)z`.  In `k(s)[y]`, write the two profiles as

```text
mathcal A=R y^4+P y^3+U_0 y^2+L_0 y+a_0,
mathcal C=Q y^3+V y^2+B_0 y+b_0.                           (22)
```

The moving scale does not add a hidden connection term.  Indeed, for any
two profiles `F,G`, the two Euler terms cancel and

```text
J_(z,s)(F(hz,s),G(hz,s))=h J_(y,s)(F,G).                  (23)
```

Consequently

```text
J_(y,s)(mathcal A,mathcal C)=delta:=lambda/h.              (24)
```

The highest remaining coefficient in `(24)` is

```text
4R V'-3Q P'=0.                                             (25)
```

Set

```text
g=V/(3Q),                         x=y+g.                    (26)
```

Translation by a function of `s` also preserves the bracket.  It depresses
the cubic profile and turns `(25)` into the assertion that the new cubic
coefficient of the quartic is constant.  Thus

```text
A_*=R x^4+D x^3+U x^2+Lx+a,
C_*=Q x^3+Bx+b,                    D in k.                  (27)
```

All coefficients in `(27)` lie in `k(s)`.  The five remaining buckets of
`J_(x,s)(A_*,C_*)=delta` are

```text
x^4: 4R B'-3Q U'=0,                                      (28)
x^3: 3D B'+4R b'-3Q L'=0,                               (29)
x^2: 3D b'-3Q a'-B U'+2U B'=0,                          (30)
x^1: -B L'+L B'+2U b'=0,                                (31)
x^0: L b'-B a'=delta.                                   (32)
```

Introduce the constants

```text
K=4R/(3Q),       H_0=D/Q,       I=2R/(9Q^2).             (33)
```

Successive integration of `(28)--(30)` gives, exhaustively,

```text
U=KB+E,
L=Kb+H_0B+F,
a=H_0b+I B^2+(2E/(3Q))B+G,       E,F,G in k.              (34)
```

The factor `B` multiplying `a'` in `(32)` is load-bearing.  Substitution of
`(34)` gives

```text
(KB+2E)b'+(Kb+F)B'=0,                                    (35)

delta=(Kb+F)b'-[B(KB+2E)/(3Q)]B'.                         (36)
```

Omitting the outer `B` in the second term of `(36)` changes the autonomous
equation and invalidates the valuation audit.  The exact companion checks
this transcription directly against the expanded Jacobian.

Put

```text
W=KB+2E.                                                   (37)
```

When `W` is not identically zero, `(35)` integrates to

```text
Wb+FB=N,                         N in k.                    (38)
```

If

```text
S=KN+2EF,                                                    (39)
```

then

```text
Kb+F=S/W.                                                   (40)
```

Equations `(36)--(40)` factor the constant bucket completely:

```text
lambda/h=-(3Q/(16R^2))
          [W^2-2EW+4R S^2/W^3]W'.                          (41)
```

Equivalently, with

```text
Phi(W)=W^3/3-EW^2-2R S^2/W^2,                              (42)
```

equation `(41)` says

```text
lambda/h=-(3Q/(16R^2)) Phi(W)'.                            (43)
```

There is a separate stratum that must not be divided away.  If
`W` is identically zero, then `B` and `U` are constant and `(36)` becomes

```text
lambda/h=(Kb+F)b'.                                         (44)
```

Equations `(21),(27),(34),(38)--(44)` are the complete rational
`(4,3)` packet.

## 4. The two original arm coefficients

The translation in `(26)` means that the original profiles in the variable
`y` are obtained by evaluating `(27)` at `x=y+g`.  In particular, their
two arm coefficients are

```text
b_0=C_*(g)=Qg^3+Bg+b,
a_0=A_*(g)=Rg^4+Dg^3+Ug^2+Lg+a.                            (45)
```

The original hypothesis says `a_0,b_0 in k[s]`.  Hence both expressions in
`(45)` are regular at every finite prime.  This simultaneous regularity,
rather than regularity of a single transformed coefficient, closes every
denominator channel below.

## 5. Nonconstant Kummer scale: three channels, three contradictions

Suppose `h` is nonconstant.  Fix an irreducible `pi|h`, and put

```text
H_pi=ord_pi(h)>0.                                           (46)
```

For every nonzero rational function `f` with `ord_pi(f)!=0`, separability in
characteristic zero gives

```text
ord_pi(f')=ord_pi(f)-1.                                    (47)
```

The derivative also preserves the DVR when `ord_pi(f)>=0`.

### 5.1 `W` has a pole

Assume first that `W` is not identically zero and

```text
ord_pi(W)=-m<0.                                             (48)
```

In `(41)`, the term `W^2W'` is uniquely most singular, so

```text
H_pi=3m+1.                                                  (49)
```

Moreover `B` has pole order `m`, while `(40)` makes `b` regular.  If `g`
has no pole, the term `I B^2` is the unique leading pole of `a_0`.  If `g`
has pole order `d>0`, the only candidates for the largest pole order in
`a_0` are

```text
2m          from I B^2,
m+2d        from K B g^2,
4d          from R g^4.                                    (50)
```

Every omitted term has strictly smaller order than the relevant maximum.
If `d<m/2`, the first entry in `(50)` is unique; if `d>m/2`, the last is
unique.  Regularity therefore forces

```text
d=m/2.                                                      (51)
```

At this order, regularity of `b_0` in `(45)` forces the leading cancellation

```text
B+Qg^2=0.                                                   (52)
```

But the leading order of `a_0` is then

```text
I B^2+K B g^2+R g^4
 = [I-K/Q+R/Q^2]B^2
 = -R B^2/(9Q^2),                                          (53)
```

which is nonzero.  Thus the best possible cancellation of the `C` arm
leaves a forced pole in the `A` arm.

### 5.2 `W` has a zero

If

```text
ord_pi(W)=n>0,                                              (54)
```

then `(41)` can have the pole `lambda/h` only when `S!=0`.  Its
`S^2/W^3` term gives

```text
H_pi=2n+1.                                                  (55)
```

Now `B` is regular, while `(40)` makes both `b` and `L` have pole order
`n`.  Regularity of `b_0` forces `g` to have pole order `n/3` and

```text
Qg^3+b=0                                                    (56)
```

at leading order.  Equation `(40)` then says

```text
S/W=-KQ g^3                                                 (57)
```

at the same order.  The unique order `4n/3` part of `a_0` is

```text
Rg^4+(S/W)g=(R-KQ)g^4=-R g^4/3,                            (58)
```

again nonzero.  If `ord_pi(W)=0`, the right side of `(41)` is DVR-regular,
so it cannot equal `lambda/h`.  Thus `(48)` and `(54)` exhaust nonzero `W`.

### 5.3 `W` is identically zero

Finally suppose `W=0` identically.  Write

```text
Y=Kb+F.                                                     (59)
```

Since `Y'=Kb'`, equation `(44)` is `lambda/h=YY'/K`.  It forces

```text
ord_pi(Y)=ord_pi(b)=-n,
H_pi=2n+1.                                                  (60)
```

The same arm argument makes `g` have pole order `n/3` and forces
`Qg^3+b=0`.  Hence `Y=-KQg^3` at leading order, and `(58)` repeats
verbatim:

```text
Rg^4+Yg=-R g^4/3!=0.                                       (61)
```

This closes the stratum that division by `W` would have lost.  Every
possible behavior at a prime factor of nonconstant `h` is contradictory.
Therefore

```text
h in k*.                                                    (62)
```

## 6. Constant Kummer scale is also impossible

When `h` is a unit, the scaling and translation in `(22),(26)` are
polynomial, so `B,b` and all coefficients in `(27)` lie in `k[s]`.

If `W=0` identically, `(44)` is a product of polynomials equal to a nonzero
constant.  Thus `Kb+F` is a unit, making `b` constant, while `b'` would also
have to be a unit.  This is impossible.

Suppose `W!=0`.  If `S!=0`, equations `(38),(40)` give

```text
W(Kb+F)=S in k*.                                            (63)
```

Both factors are units.  Hence `B,b` are constant and `(36)` gives
`delta=0`, a contradiction.  If `S=0`, equation `(40)` makes `b` constant,
and the load-bearing factor `B` in `(36)` gives

```text
delta=-[B W/(3Q)]B'.                                        (64)
```

A nonzero constant product would make `B` a unit, hence constant, and then
`B'=0`.  This is the final contradiction.

Thus the `(4,3)` row is empty.  Together with Sections 1--2 and THM-3861,
this exhausts every quartic and zero-component branch in `(1)`.

## 7. Sharp rational hostiles

The simultaneous-arm step is load-bearing.  The exact companion contains
two positive controls in which the `C`-arm pole cancels but the residual
predicted above survives in the `A` arm.

For the pole channel, work over `k=Q(d)` with `d^2=-3/4` and set

```text
C=s^21 z^3+3d s^13 z^2-(3/2)s^5 z,

A=s^28 z^4+4d s^20 z^3-(7/2)s^12 z^2-d s^4 z
  -1/(16s^4).                                               (65)
```

Direct differentiation gives

```text
J_(z,s)(A,C)=3/8.                                           (66)
```

This is the packet with `R=Q=1`, `h=s^7`, `W=s^-2`, `E=S=0`, and
`g=d/s`.  The relation `d^2=-3/4` makes `b_0=0`; the only nonpolynomial
coefficient is

```text
a_0=-1/(16s^4)=-B^2/9.                                     (67)
```

For the identically-zero `W` channel, set

```text
C=s^21 z^3+3s^13 z^2+3s^5 z,

A=s^28 z^4+4s^20 z^3+6s^12 z^2+(8/3)s^4 z
  -1/(3s^4).                                                (68)
```

Then

```text
J_(z,s)(A,C)=-4.                                            (69)
```

Here `h=s^7`, `g=1/s`, and the `C` arm again vanishes, while
`a_0=-g^4/3` is exactly the residual in `(61)`.  Neither `(65)` nor `(68)`
is a polynomial Keller map, so neither has a Jacobian-conjecture
consequence.  They are sharp hostiles against deleting either arm or the
`W=0` stratum.

## 8. Russell and planar-Jacobian boundary

THM-3846 identifies the completed Russell arm with `k[s][[z]]`, and
THM-3843 forces every hypothetical global Russell Darboux pair to be
noninjective on that arm.  The restriction of a polynomial automorphism of
`A2_(s,z)` to `z=0` is injective.  The present theorem therefore excludes
polynomial canonical normal expansions through transverse degree four.

This does not truncate a general global element of the Russell surface,
which can have an infinite canonical `z`-series.  THM-3860 also shows that
rational higher-normal flexibility exists beyond the polynomial strip.
Transverse degree five is only the next bounded polynomial cell; arbitrary
Keller maps and the planar Jacobian conjecture remain open.  **QED.**

## 9. Reproduction

```bash
python3 -B 04-computation/jc2_quartic_normal_strip_keller_thm3867.py
python3 -B -O 04-computation/jc2_quartic_normal_strip_keller_thm3867.py
python3 -B 04-computation/jc2_quartic_normal_strip_keller_independent_audit_thm3867.py
python3 -B -O 04-computation/jc2_quartic_normal_strip_keller_independent_audit_thm3867.py
```

Each primary execution must byte-match
`05-knowledge/results/jc2_quartic_normal_strip_keller_thm3867.out`; each
independent execution must byte-match
`05-knowledge/results/jc2_quartic_normal_strip_keller_independent_audit_thm3867.out`.
