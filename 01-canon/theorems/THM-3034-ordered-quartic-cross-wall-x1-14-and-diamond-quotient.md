---
id: THM-3034
title: "Ordered quartic cross-wall X1(14) identification and diamond quotient"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE AUDIT.
  The positive ordered quartic cross-wall is explicitly Q-isomorphic to
  X_1(14); its even flank cycle is translation by the visible rational
  three-torsion point, and its normalized quotient is the natural degree-three
  diamond cover to X_0(14).  An odd flank move exchanges the two sign sheets
  and is not a diamond C2 or Fricke involution.  No Keller or LRC consequence
  is asserted.
source: codex-ordered-quartic-x1-14-diamond-quotient-2026-08-01
depends_on:
  - THM-2998-quartic-star-complement-cross-wall-real-locus-and-c3-x0-14-quotient
related:
  - THM-2950-three-conjugate-pair-v-four-torsor-and-quartic-resolvent-frame
  - THM-2968-quartic-edge-and-oriented-cycle-s4-complements
  - THM-2996-prime-modular-affine-defect-trichotomy-and-spherical-quartic-uniqueness
script: 04-computation/quartic_cross_wall_x1_14_diamond_thm3034.py
output: 05-knowledge/results/quartic_cross_wall_x1_14_diamond_thm3034.out
script_sha256: 39019d1a0efb8e2f8bc343d1db2e48a1842bfb2125fd1e0ba126ce0370a5ca2d
output_sha256: f939f2201e6b3c1db08911e32b26f72570fd80314185cc803d689fbcc3122d10
hash_basis: LF-normalized bytes
---

# THM-3034 -- ordered quartic cross-wall as `X_1(14)`

**PROVED CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE AUDIT.**

## 1. Inheritance and statement

[THM-2998](THM-2998-quartic-star-complement-cross-wall-real-locus-and-c3-x0-14-quotient.md)
proves that the positive ordered cross-wall

```text
C_+ : ABC=(A-B)(A-C)(B-C)                              (1)
```

is a smooth genus-one curve.  The even flank cycle

```text
rho[A:B:C]=[B:C:A]                                     (2)
```

acts freely, and its quotient is the standard `X_0(14)` model

```text
E_0 : V^2+UV+V=U^3+4U-6.                               (3)
```

The theorem left the ordered source curve unnamed.  It is the other natural
level-`14` modular curve:

```text
E_1 : v^2+uv+v=u^3-u.                                  (4)
```

More precisely:

1. `C_+` is explicitly `Q`-isomorphic to `(4)`.  The negative sign sheet
   `C_-` is another `Q`-isomorphic copy, obtained by an odd flank swap.
2. Under the stated normalization, `rho` is translation by
   `T=(0,0)`, where `T` has exact order three.
3. Quotienting by `<T>` gives the explicit degree-three map `(14)` below to
   `(3)`.  As an unpointed cover, this is the natural diamond/forgetful cover

   ```text
   X_1(14) -> X_0(14),          (E,P) -> (E,<P>).        (5)
   ```

4. An odd flank permutation exchanges `C_+` and `C_-` and conjugates `rho`
   to `rho^{-1}`.  It is not a second deck transformation of one sheet and is
   neither a proved diamond `C_2` nor a Fricke/Atkin--Lehner involution.

The external modular naming is pinned to the LMFDB records
[14.a5 / `X_1(14)`](https://www.lmfdb.org/EllipticCurve/Q/14/a/5) and
[14.a6 / `X_0(14)`](https://www.lmfdb.org/EllipticCurve/Q/14/a/6).
Those records give exactly the models `(4)` and `(3)`.  Every change of
variables and quotient formula used here is independently checked over `Q`.

## 2. Direct projection and birational equivalence

On the chart `B!=0`, write

```text
tau=C/B,                         xi=A/B.
```

Then `(1)` is

```text
xi tau=(xi-1)(xi-tau)(1-tau),                         (6)
```

and

```text
eta=2(tau-1)xi-tau^2+tau+1                            (7)
```

satisfies the exact identity

```text
eta^2-(tau^4-6tau^3+7tau^2-2tau+1)
 =4(tau-1)[xi tau-(xi-1)(xi-tau)(1-tau)].             (8)
```

Thus projection gives the THM-2998 quartic

```text
eta^2=tau^4-6tau^3+7tau^2-2tau+1.                     (9)
```

On its dense `tau!=0` chart, put

```text
chi=(eta+1-tau)/tau^2,
u=(chi+1)/2,
v=(chi^2-1)tau/4.                                     (10)
```

Direct reduction by `(9)` gives `(4)`.  Conversely, where `u(u-1)!=0`,

```text
tau=v/[u(u-1)],
eta=(2u-1)tau^2+tau-1,
xi=(eta+tau^2-tau-1)/[2(tau-1)].                      (11)
```

Substitution recovers `(6)` and `(9)`, and both round trips are the identity
on the common dense chart.  Since both projective curves are smooth, the
birational equivalence extends uniquely to an isomorphism.

The exact short-model change is

```text
X=36u+3,                    Y=108(2v+u+1),
Y^2=X^3-675X+13662.                                      (12)
```

The generalized invariants are

```text
(c4,c6,Delta,j)=(25,-253,-28,-15625/28),               (13)
```

matching the independently named `X_1(14)` record.

## 3. The cyclic flank phase is rational three-torsion

Choose `O_+=[0:1:0]` as origin on `C_+`.  On the chart above, pulling `(2)`
through `(7)--(10)` gives

```text
u'=(v/u)^2+v/u-u,
v'=-(v/u+1)u'-1.
```

These are exactly the addition formulas for `T=(0,0)` on `(4)`.  The tangent
at `T` is `v=-u`; its intersection polynomial with `(4)` is `-u^3`.  Hence

```text
2T=(0,-1)=-T,                       3T=O.                (14a)
```

Equivalently, the rational coordinate orbit is

```text
[0:1:0] -> [1:0:0] -> [0:0:1] -> [0:1:0]
O        -> T         -> -T        -> O.                 (14b)
```

The `3`-division polynomial is

```text
psi_3(u)=u(3u^3+u^2-3u+3).
```

The cubic factor is irreducible over `Q`, so `{O,T,-T}` is the unique
rational order-three subgroup.  This identifies the action subgroup, not
merely the source curve's `j`-invariant.

## 4. The normalized quotient is the diamond map

The quotient functions are

```text
U=(u^3-u+1)/u^2,
V=[v(u^3+u-2)+u^2-u-1]/u^3.                            (14)
```

Reducing by `(4)` gives exactly `(3)`.  In the short coordinates

```text
X =36u+3,                  Y =108(2v+u+1),
X'=36U+3,                  Y'=108(2V+U+1),
```

this is the normalized Velu quotient

```text
X'=(X^3-6X^2-1287X+50544)/(X-3)^2,
Y'=Y(X^3-9X^2+1323X-97227)/(X-3)^3,                   (15)
```

from `(12)` to

```text
Y'^2=X'^3+5805X'-285714.                               (16)
```

The target invariants are

```text
(c4,c6,Delta,j)=(-215,5291,-21952,9938375/21952),      (17)
```

the standard `X_0(14)` values already identified in THM-2998.

The coarse diamond deck group is

```text
(Z/14Z)^*/{+-1}
 = {{1,13},{3,11},{5,9}}
 = C_3.                                                 (18)
```

The natural forgetful map therefore has degree three and a rational free
`C_3` deck action.  By the uniqueness of the rational order-three subgroup in
Section 3, its kernel is `<T>`.  Hence `(14)` and the THM-2998 cyclic quotient
are the modular cover `(5)`, with the generator determined only up to
inversion (`<3>` versus `<5>`) until a modular cusp labelling is chosen.

## 5. Pointed and unpointed quotients are different data

THM-2998's symmetric quotient coordinates are not the pointed Velu map
`(14)`.  At `O_+`, its symmetric coordinates are `t=s=0`, and its equations
(33)--(35) give

```text
O_+ -> (9,-33) in E_0(Q),                              (19)
```

not the point at infinity.  The point `(9,-33)` satisfies `(3)`.  Thus the
displayed THM-2998 quotient is the same unpointed `C_3` cover followed by a
target translation and possibly target negation.  A compatible choice of
source and target cusps is the sidecar needed to identify literal pointed
coordinate formulas.

The first failed implication is therefore

```text
same degree-three quotient cover
  does not imply that two displayed quotient formulas preserve the same
  origins or name the same generator.                               (20)
```

## 6. Where the binary and ternary faces really meet

An odd flank transposition `sigma` changes the Vandermonde sign, so

```text
sigma:C_+ -> C_-,                   sigma rho sigma^-1=rho^-1.       (21)
```

After transporting origins and a generator orientation between the two
copies, it may be represented by elliptic negation

```text
(u,v) -> (u,-v-u-1),                                       (22)
```

which sends `T` to `-T`.  Only the two-sheet groupoid relation `(21)` is
canonical.  It is false to call `(22)` an extra diamond involution on one
sheet without that synchronization.

This gives a precise version of the recurring `2`/`3` picture.  On the six
flank orderings, `sigma` and `rho` generate

```text
<sigma,rho | sigma^2=rho^3=(sigma rho)^2=1> = S_3.       (23)
```

Thus the quartic ordering data see the finite quotient

```text
PSL_2(Z)=C_2*C_3 -> S_3=PSL_2(F_2),                     (24)
```

not a faithful modular-group action.  Together with the quartic affine
extension `V_4 normal S_4` and `S_4/V_4=S_3`, this explains why the binary
`V_4` frame and ternary cyclic quotient repeatedly co-occur.  The additional
relation `(sigma rho)^2=1` is exactly the information collapse: the
`Gamma(2)`/Farey-origin kernel, affine root origin, and physical owner are no
longer visible.

## 7. Connection contract and stopping boundary

```text
source:     one ordered sign sheet of the marked quartic cross-wall;
map:        the explicit birational map (7)--(11), followed by (14);
target:     X_1(14) and its natural degree-three quotient X_0(14);
preserved:  sign sheet, cyclic flank phase, rational three-torsion kernel;
destroyed:  odd sheet sign after quotient, pointed target origin, named
            diamond generator, marked quartic root, affine scale;
sidecars:   compatible cusp/origin labels for a pointed modular map, plus an
            independent graph owner/current for any physical application;
hostiles:   (19), the sheet exchange (21), and the extra S_3 relation (23).
                                                                    (25)
```

This theorem supplies no universal modular interpretation of a quartic, no
Keller/Jelonek owner, no LRC carrier or current, and no ledger decrement.

## 8. Exact companion

Run

```text
python 04-computation/quartic_cross_wall_x1_14_diamond_thm3034.py
python -O 04-computation/quartic_cross_wall_x1_14_diamond_thm3034.py
```

Both modes byte-match

```text
05-knowledge/results/quartic_cross_wall_x1_14_diamond_thm3034.out.
```

The companion checks the plane-cubic projection, both birational directions,
the exact flank-translation formulas, the division polynomial, the normalized
quotient and both short models, the modular invariants, the diamond classes,
the target-origin hostile, and the odd/even permutation signs.  It contains no
truth-bearing Python `assert`.
