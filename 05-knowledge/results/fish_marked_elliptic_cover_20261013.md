# A marked elliptic obstruction inside the fish, and a necessary genus-nine cover

2026-10-13. **PROVED, scoped, with the recorded exact 2-descent inputs:**
the necessary cover, the exclusion of `Q(zeta11,sqrt5)`, and the exclusion
of a rational central parameter in `Q(zeta55)`. **VERIFIED:** PARI's
unconditional rank bounds, checked torsion and finite-field counts.
**FINITE-EXACT:** graph and algebra controls. Self-audited; no independent
peer audit or formalization. The fish's placement in the full field,
six-colour forcing and Hadwiger–Nelson remain **OPEN**.

**UPDATE 2026-10-16 (opus-2026-10-16-S38):** the open test `C3(F_4)` of section 5 is decided: **`C3(F_4)` is empty**
([THM-4662](../../01-canon/theorems/THM-4662-the-fish-curve-has-no-point-over-the-real-quartic-field-so-the-fish-and-ei21-are-not-in-q-zeta55-or-any-odd-carrier.md); a two-descent over `Q(sqrt5)`:
`rank E(F_4) = 2` and the class `(33, 33)` of `C3` is not a descent class). With (5), the fish has **no** placement in
`Q(zeta55)`, nor in any conjugation-stable odd-degree extension of `Q(zeta5, sqrt(-11))` (THM-4661 4(b)).

## 1. Inheritance and the exact question

The [incoming fish analysis](ei19_fish_placement_q_zeta55_20261012.md)
uses the House of Graphs 51404 labels `0,...,22`. Its direction relations
leave a one-parameter real family; its finite parameter scan finds no
placement in `L=Q(zeta55)`. Absence of a support-three **linear** relation
does not remove consequences of the unit-norm equations. We apply those
equations to the already identified inscribed parallelogram.

The [six-colour extension theorem, THM-4614](../../01-canon/theorems/THM-4614-three-colourable-contexts-cannot-change-the-six-colour-port-relation-in-q-zeta55.md)
requires a four-chromatic exterior before the eight original ports can be
forced. The fish is one candidate for that exterior, not itself a
six-colour forcing context. The companion independently checks its abstract
chromatic number is four.

Closest mechanism: the [rhombus square-class and elliptic gates](rhombus_overlap_elliptic_20261011.md).
Hostile: an elliptic point that fails its auxiliary square-root lift.
Corrected near miss: treating a rank or one quotient curve as the entire
placement. Least-used sidecar: the field of the **whole tuple** of square
roots. The golden field enters through that field test, not a decimal match.

Here a placement is injective and sends every graph edge to a unit segment;
additional unit distances are allowed. The exclusion therefore also applies
to induced faithful placements. Translate the hub `v22` to zero, keeping
all vertex labels and the unit scale.

## 2. A forced central rhombus

Every noncollapsed equilateral four-cycle in the plane is a parallelogram.
The five actual cycles

```
(0,10,8,20), (2,8,10,19), (1,11,9,21),
(3,9,11,18), (18,19,21,20)
```

give

```
v0+v2 = v19+v20 = v18+v21 = v1+v3.
```

The four distinct points `v0,v1,v2,v3` are on the unit circle about zero.
Two distinct unordered pairs on this circle can have the same sum only
when that sum is zero: for a nonzero midpoint the chord is uniquely
determined. Consequently write

```
v0=a, v2=-a, v1=b, v3=-b,
v18=c, v21=-c, v19=d, v20=-d.
```

The central edges give `|c-d|=|c+d|=1`, hence

\[
 \operatorname{Re}(c\bar d)=0,\qquad |c|^2+|d|^2=1.
\]

Put `D=|d|^2`, so `0<D<1`. Neither endpoint is allowed by injectivity.
Let `F=L^+`, and `beta=sqrt(-11)`. Every conjugation-anti-invariant element
of `L` is `beta` times an element of `F`.

The actual marked positions supply elements `z0,z1,z2,z3` of `F` satisfying

\[
\begin{split}
11z_0^2&=D(1-D),\\
11z_1^2&=D(4-D),\\
11z_2^2&=(1-D)(D+3),\\
11z_3^2&=(1+2D)(3-2D).                         \tag{1}
\end{split}
\]

For the first equation take
`z0=(c*bar(d)-bar(c)*d)/(2 beta)`.
The edges `(0,20)` and `(1,21)` give
`Re(a*bar(d))=-D/2` and `Re(b*bar(c))=-(1-D)/2`.
Thus take the corresponding skew products divided by `beta` for `z1,z2`.
Finally `v10` is at unit distance from both `a,d`; their squared separation
is `1+2D`. For `u=v10-a`, `v=v10-d`, take
`z3=(u*bar(v)-bar(u)*v)/beta` and use
`(u*bar(v)-bar(u)*v)^2=(1+2D)((1+2D)-4)`.

All embeddings commute with conjugation. Thus `D` lies between zero and
one at **every** real embedding of `F`. Equations (1), their simultaneous
square roots, and this total-positivity condition are necessary. They are
not claimed sufficient for all 23 vertices, or even for distinctness of
all vertices in the smaller subframework.

## 3. The genus calculation and elliptic quotient

Let `C3` be the smooth projective normalization of the first three equations
of (1), and `C4` that of all four. Over an algebraic closure their square
classes have branch sets

```
{0,1}, {0,4}, {1,-3}, {-1/2,3/2}.
```

They are independent: every nonempty product has a root of odd valuation.
The first three give a degree-eight cover of the `D`-line, with four branch
values and ramification index two at each. Riemann–Hurwitz gives
`2g(C3)-2=8(-2)+4*4=0`, so `g(C3)=1`.
The fourth gives a degree-sixteen cover with six branch values, whence

\[
 2g(C4)-2=16(-2)+6\cdot8=16,\qquad g(C4)=9.       \tag{2}
\]

Equivalently `C4 -> C3` is a double cover ramified at sixteen points above
`D=-1/2,3/2`. The fifteen character quotients have seven genus-one pieces,
one genus-two piece and seven genus-zero pieces, whose genera sum to nine.
This concerns the **necessary subframework cover**, not the still unknown
genus of the complete fish realization curve.

Taking `y=11z1z2` gives the elliptic quotient

\[
 y^2=D(1-D)(4-D)(D+3).
\]

For `D != 0`, the exact change

\[
 X=12/D-3,\qquad Y=12y/D^2
\]

gives

\[
             E:\quad Y^2=X(X-9)(X+7).              \tag{3}
\]

Discarding the `z_i` would forget the crucial lifting obstruction.

## 4. The golden-field descent excludes a whole proper CM field

Put `K0=Q(sqrt5)`. Exact PARI 2-descent returns

```
ellrank(E)             = [0,0,0,[]]
ellrank(elltwist(E,5)) = [0,0,0,[]].
```

These are algebraic lower/upper rank bounds, not analytic-rank estimates
or applications of BSD. See the primary
[PARI ellrank documentation](https://pari.math.u-bordeaux.fr/dochtml/html-stable/Elliptic_curves.html#ellrank).
The trace/anti-trace decomposition therefore gives `rank E(K0)=0`.

The rational points generated by `(-3,12)` and `(0,0)` give a subgroup
`Z/4 x Z/2` of order eight. The good primes 11 and 19 split in `K0`, and
direct counting gives `#E(F11)=8`, `#E(F19)=16`. Prime-to-residue-characteristic
torsion injects under reduction; using both primes bounds **all** torsion
over `K0` by eight. Hence `E(K0)=E(Q)` consists of those eight points.
Their possible `D` coordinates on the quartic are

```
0, 1, 4, -3, 1/2, infinity.
```

In fact **`C3(K0)` is empty**, including the compactification:

- A real finite point needs `0<=D<=1` from its first equation.
- At `D=0` or `D=1`, a remaining equation requires `sqrt(3/11)` in `K0`.
- At `D=1/2`, the first equation requires `sqrt(1/44)` in `K0`.
- At infinity, `z_i/D` would have square `-1/11`, impossible in a real field.

Neither `3/11` nor `1/44` is a rational square or five times one.
This is the precise golden-field obstruction.

The degree-eight map `C3 -> P1` supplies a rational divisor of degree eight.
If `C3` acquired a point over any finite **odd-degree** extension of `K0`,
its trace would supply an odd-degree divisor. An integer combination gives
a divisor of degree one. Riemann–Roch on this genus-one curve then gives
a `K0`-rational point, contradicting the preceding paragraph. Therefore

\[
 C3(K')=\varnothing\quad\text{for every finite odd-degree }K'/K0. \tag{4}
\]

In particular the CM field

\[
 M=Q(\zeta_{11},\sqrt5),\qquad [M:Q]=20,
\]

has real subfield `M^+=Q(zeta11+zeta11^-1,sqrt5)` of degree five over `K0`.
It contains `sqrt(-11)`, so the construction of (1) works there unchanged.
Equation (4) proves that **the fish has no injective unit-edge placement in
`M`**. This excludes an infinite field, not a finite generator scan.

There is a second consequence in the full `L`. Its real Galois group is
cyclic of order twenty, so its only quadratic subfield is `K0`.
If `D` were rational, all three square roots in the first three equations
would already belong to `K0`, giving a forbidden point of `C3(K0)`.
Thus **every fish placement in `L`, if one exists, has nonrational `D`**.
Neither result excludes the full quartic real direction available in `L`.

## 5. The first existence test descends to a quartic field

There is a sharper reduction before attempting a degree-twenty real-field
calculation. Define

\[
 F_4=Q(\tau),\qquad
 \tau^2=\frac{55+11\sqrt5}{2},\qquad
 \tau^4-55\tau^2+605=0.
\]

This is the real quartic subfield of `Q(zeta5,sqrt(-11))`; one can take
`tau=sqrt(-11)(zeta5-zeta5^-1)`. Since `[F:F4]=5`, the same degree-eight
divisor argument proves

\[
                  C3(F)\ne\varnothing
                  \quad\Longleftrightarrow\quad
                  C3(F_4)\ne\varnothing.              \tag{5}
\]

This is an **existence** reduction. Riemann–Roch may produce a different
point, so it does not preserve a supplied `D`, ports, or a lift to `C4`.
It must not be used to replace a marked shape by another shape. The
genus-nine cover has no corresponding degree-one-divisor implication:
there a degree-one divisor need not have a nonzero section. This isolates
where the extra quintic field could still matter.

For an explicit arithmetic model of `C3`, put `r=z1/z0` and
`s=11(r^2-1)z0`. Then

```
D=(r^2-4)/(r^2-1),       s^2=33(r^2-4).
```

The latter conic has the rational point `(r,s)=(7/2,33/2)`.
The line of slope `t` through it gives the second intersection

\[
 r=\frac{7t^2-66t+231}{2(t^2-33)},\qquad
 D=\frac{11(t-11)^2(t-3)^2}
 {(3t^2-22t+55)(5t^2-66t+297)}.
\]

Writing `H=(3t^2-22t+55)(5t^2-66t+297)`, the remaining equation becomes

\[
 Y^2=22(7t^4-154t^3+1342t^2-5082t+7623),\qquad
 z2=\frac{4(t^2-33)}{11H}Y.                            \tag{6}
\]

The factors recover the actual square roots; the smooth projective model
takes care of the parametrization's exceptional points.

**FINITE-EXACT probe:** for every primitive triple
`-30<=a,b<=30`, `1<=d<=12`, `b!=0`, test `t=(a+b sqrt5)/d`.
There are 36,518 such parameters. None makes the right side of (6) a
square in `F4`. This uses exact arithmetic: an element `q` of `Q(sqrt5)`
is square in `F4` iff either `q` or `q/((55+11sqrt5)/2)` is square in
`Q(sqrt5)`, and both tests reduce to rational square tests and a norm.
This is a bounded search, not emptiness of `C3(F4)`; in particular it does
not examine arbitrary quartic-field values of `t`.

A separate first probe in the missing quartic direction uses PARI's exact
number-field square test on
`t=(a+b*tau+c*tau^2+e*tau^3)/d`, with all four coefficients between
`-2` and `2`, `1<=d<=3`, primitive five-tuples, and `(b,e)!=(0,0)`.
None of these 1,728 genuinely quartic parameters lifts. Positive and
negative square controls are recorded. This small second box also leaves
`C3(F4)` undecided; it identifies no universal height cutoff.

## 6. What the attached papers contribute, and the next tests

The supplied low-Selmer-corank manuscripts distinguish rank from the full
integral arithmetic data. Here that distinction becomes concrete: even a
complete elliptic point list does not give a graph placement until the
marked square roots lift. The rank calculations above use established
2-descent, not the manuscripts' headline claims. The intake and cross-topic
map are in [the companion synthesis](marked_switches_intake_20261013.md).

The next exact tasks are the explicit genus-one test `C3(F4)` in (6), and
then the marked ramified double cover `C4(F)`,
with the totally positive parameter and the remaining fish edges retained.
A point on (3) alone is insufficient. A route that stays in `M` is closed;
a candidate must use the extra quadratic real direction of `L/M`.
No original pentagon/rhombus port has been moved or recoloured: the work
screens possible exteriors before attempting to attach them.

The concurrent [finite abelian Cayley realization criterion](unit_cayley_realization_20261013a.md)
uses the same unit-four-cycle identity but retains different global data.
There the closed generator walks force an internal product of cyclic
subgroups; here the labelled central chords force simultaneous square
classes. The fish is not being replaced by a Cayley graph. In particular,
large chromatic numbers of finite-field norm-circle graphs do not supply
missing fish placements: that packet excludes their complete injective
plane realizations except at `q=3`. This rules out a shortcut while leaving
the non-Cayley, marked-cover test above intact.

## 7. Reproduction

[Script](../../04-computation/experiments/marked_switches_20261013.py),
[output](marked_switches_20261013.json): run normally and with `-O`.
The companion checks the actual graph edges and midpoint relations, abstract
3-colour rejection and a 4-colour witness, an injective real thirteen-vertex
subframework control, square-class independence, both genus calculations,
the birational identities and 36,518 exact quartic tests, exact descent bounds, a separately implemented
rational group law, and independent finite-field point counts.
The real control lives in a larger field and is not an `L`-placement.

Connection contract: fish positions -> labelled norm/skew-product tuple ->
genus-nine necessary cover -> elliptic quotient. Each arrow preserves the
displayed equations; the last loses three square-root choices and the
unrepresented edges. Keeping those coordinates is what produces the field
exclusion rather than a false placement claim.
