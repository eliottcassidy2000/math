# Ordered quartic cross-wall as `X_1(14)` and its diamond quotient

**Status:** PROOF-COMPLETE SCRATCH CANDIDATE / FINITE-EXACT COMPANION / NOT
CANON.  This is an input to reserved THM-3034, not a proved dependency.

## 1. Inheritance and exact claim

THM-2998 proves that the positive ordered cross-wall

```text
C_+ : ABC=(A-B)(A-C)(B-C)
```

is a smooth genus-one curve, that

```text
rho[A:B:C]=[B:C:A]
```

acts freely with order three, and that the quotient is the standard
`X_0(14)` model

```text
E_0 : V^2+UV+V=U^3+4U-6.                         (1)
```

Its projection from `[1:0:0]` gives

```text
eta^2=tau^4-6tau^3+7tau^2-2tau+1.                (2)
```

The missing identification is exact:

```text
E_1 : v^2+uv+v=u^3-u                             (3)
```

is birational over `Q` to `(2)`, the action `rho` is translation by the
rational point `T=(0,0)` of order three, and the normalized quotient by
`<T>` is the explicit degree-three isogeny `(12)` below.  The LMFDB records
`E_1` as its model of `X_1(14)` and `E_0` as its model of `X_0(14)`:

- <https://www.lmfdb.org/EllipticCurve/Q/14/a/5>
- <https://www.lmfdb.org/EllipticCurve/Q/14/a/6>

Consequently the THM-2998 quotient is, as a cover over `Q`, the natural
diamond/forgetful quotient

```text
X_1(14) -> X_0(14),
(E,P)   -> (E,<P>),
```

whose deck group is

```text
(Z/14Z)^x/{+-1} = {{1,13},{3,11},{5,9}} = C_3.   (4)
```

This statement identifies the deck **subgroup**.  Without a chosen modular
cusp labelling it does not distinguish `rho` as the generator `<3>` rather
than its inverse `<5>`.

## 2. Direct projection and birational map

On the chart `B!=0`, write

```text
tau=C/B,       xi=A/B.
```

Then the plane cubic becomes

```text
xi*tau=(xi-1)(xi-tau)(1-tau),                    (5)
```

and

```text
eta=2(tau-1)xi-tau^2+tau+1                       (6)
```

satisfies

```text
eta^2-(tau^4-6tau^3+7tau^2-2tau+1)
 =4(tau-1)[xi*tau-(xi-1)(xi-tau)(1-tau)].        (7)
```

For a point of `(2)`, put

```text
x=(eta+1-tau)/tau^2,
y=(x^2-1)tau+x+3.                                 (8)
```

Exact completion of the quadratic in `tau` gives

```text
y^2=2x^3+7x^2+4x+3.                              (9)
```

Now set

```text
u=(x+1)/2,
v=(y-x-3)/4=(x^2-1)tau/4.                        (10)
```

Substitution gives `(3)`.  Conversely, on `u(u-1)!=0`,

```text
tau=v/[u(u-1)],
eta=(2u-1)tau^2+tau-1,
xi=(eta+tau^2-tau-1)/[2(tau-1)],                 (11)
```

and `[A:B:C]=[xi:1:tau]`.  Thus `(6)--(11)` are an explicit
birational equivalence, not merely a match of `j`-invariants.

The intermediate cubic `(9)` transforms to the short model by

```text
X=18x+21=36u+3,
Y=54y=108(2v+u+1),
Y^2=X^3-675X+13662.
```

This independently recovers the exact short model printed for LMFDB
`14.a5`.

## 3. The cyclic flank action is translation by `(0,0)`

Choose

```text
O_+=[0:1:0]
```

as origin on `C_+`.  It is the quartic point `(tau,eta)=(0,1)` and maps to
the point at infinity of `(3)`.  Under `rho`,

```text
tau' = xi/tau,             xi'=1/tau,
eta'=2(tau'-1)xi'-tau'^2+tau'+1.                  (12a)
```

Substituting `(12a)` into `(8)--(10)` gives exactly the addition formulas
for `T=(0,0)`:

```text
u'=(v/u)^2+v/u-u,
v'=-(v/u+1)u'-1.                                  (12b)
```

The tangent at `T` is `v=-u`; its intersection polynomial with `(3)` is
`-u^3`.  Hence `2T=-T=(0,-1)` and `3T=O`.  Translation by nonzero torsion is
fixed-point-free, recovering the geometric freeness in THM-2998.

The coordinate orbit is consistent with this normalization:

```text
[0:1:0] -> [1:0:0] -> [0:0:1] -> [0:1:0]
O        -> T         -> -T        -> O.
```

In particular, the free cyclic action is not merely an abstract order-three
automorphism of a curve with the same `j`; it is the rational translation
subgroup `<(0,0)>` on the standard `X_1(14)` model.

## 4. Exact normalized quotient

The normalized quotient by `<T>` is

```text
U=(u^3-u+1)/u^2,
V=[v(u^3+u-2)+u^2-u-1]/u^3.                       (13)
```

Direct substitution gives `(1)`.  In the short coordinates

```text
X =36u+3,                 Y =108(2v+u+1),
X'=36U+3,                 Y'=108(2V+U+1),
```

the same map is the normalized Velu map

```text
X'=(X^3-6X^2-1287X+50544)/(X-3)^2,
Y'=Y(X^3-9X^2+1323X-97227)/(X-3)^3,               (14)
```

from

```text
Y^2 =X^3-675X+13662
```

to

```text
Y'^2=X'^3+5805X'-285714.
```

The exact invariants are

```text
E_1: (c4,c6,Delta,j)=(25,-253,-28,-15625/28),
E_0: (c4,c6,Delta,j)=(-215,5291,-21952,9938375/21952).
```

The LMFDB identifies these equations with `X_1(14)` and `X_0(14)`, records
`E_1(Q)_tors=Z/6Z`, and its isogeny-class matrix records degree three between
`14.a5` and `14.a6`.  Independently, `(4)` shows that the natural forgetful
map has deck group `C_3`.  Because `j(E_1)` is nonintegral, `E_1` is non-CM;
any free rational order-three curve automorphism is translation by rational
three-torsion.  The visible subgroup `{O,T,-T}` is the unique such subgroup
inside `E_1(Q)_tors`.  Therefore the geometric quotient by `rho` and the
diamond quotient are the same quotient cover (with the generator determined
only up to inversion).

## 5. The target-origin normalization caveat

THM-2998's symmetric quotient coordinates are not themselves the pointed
Velu map `(13)`.  At the chosen source origin, `t=s=0`; equations (33)--(35)
of THM-2998 give

```text
O_+ -> (9,-33) in E_0(Q),                          (15)
```

not the point at infinity.  Thus the displayed THM-2998 map is the same
unpointed quotient cover as `(13)`, followed by a target translation and
possibly target negation.  Translating `(15)` to infinity produces the
pointed isogeny.  Likewise, identifying the cover with the **named** modular
forgetful morphism requires a compatible choice of source and target cusps;
the quotient subgroup itself requires no such choice.

This is the first failed implication to keep visible:

```text
same degree-three quotient cover
  does not mean the displayed affine formulas preserve the chosen origins.
```

## 6. The odd `C_2` is only a sheet exchange

For an odd flank permutation `sigma`, the Vandermonde changes sign, so

```text
sigma : C_+ -> C_- ,
sigma rho sigma^{-1}=rho^{-1}.                    (16)
```

It is not an automorphism of one ordered sign sheet.  Consequently it is
neither an extra diamond involution (the diamond deck group is already
`C_3`) nor a proved Fricke/Atkin--Lehner involution.

There are two equivalent bookkeeping conventions.

1. Transport the origin across `sigma`; then `sigma` may be represented as
   the identity between two copies of `E_1`, while it reverses which
   translation is called `rho`.
2. Synchronize both copies so that the labelled `rho` is translation by the
   same `T`; then `sigma` is represented by elliptic negation

   ```text
   (u,v) -> (u,-v-u-1),                             (17)
   ```

   which sends `T` to `-T`.

Only the groupoid relation `(16)` is canonical.  Calling `(17)` a modular
`C_2` without the synchronization sidecar would reintroduce exactly the
orientation loss that the ordered cross-wall was meant to retain.

## 7. Scope

This identifies the ordered quartic wall and its quotient as modular curves;
it does **not** construct a universal modular interpretation of a quartic,
a Keller/Jelonek owner, an LRC carrier, or a physical current.  The preserved
data are the ordered sign sheet and its cyclic flank phase.  The destroyed
data are the chosen affine origin after quotient, odd sheet sign, marked
colliding root, and every external physical owner coordinate.

The exact companion is

```text
python .scratch/thm3034_ordered_quartic_x1_14_diamond_referee.py
python -O .scratch/thm3034_ordered_quartic_x1_14_diamond_referee.py
```

and both modes must byte-match

```text
.scratch/thm3034_ordered_quartic_x1_14_diamond_referee.out.
```

The frozen LF-normalized hashes are

```text
script  bce505591d0fad9c61b0ff5edcec7974d3be94e326d8be78824d4c4089999771
output  b20bb786bff577aaacac08a4ea186a5916eed595ae0d7bb1f53b528d23ced52d
```
