# LRC14 reflected levels: the full two-D cone closes

**Proof candidate + independently hardened exact referee, 2026-08-01.**  This
note is intentionally not a canon promotion.  It gives a complete candidate
proof, with a frozen exact referee, that the reflected `k=1` sufficient family
of [THM-2941, the critical-seven-slot scalar wall and balanced
boundary](../01-canon/theorems/THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary.md)
closes
throughout

```text
D >= 6,                    m >= 2D.
```

Combined with the audited reflected assembly in THM-2941 and its
arbitrary-level robust-edge-eight bank, the
remaining *certificate-failure* locus is confined to the same `561` bodies
with `D>=6` and `1<=m<2D`.  This is neither a physical-survivor census nor a
proof of LRC(14).

## 1. Inheritance and the missing coordinate

Let `E={e_1,...,e_6}`, `L=14*lcm(E)`, and attach levels `q_i` with

```text
m=min q_i,                 D=max q_i-min q_i.
```

The closest proved mechanism is THM-2941's audited cross-determinant
transport.  Its
primitive phase floor sees the reduced pair `(P,Q)`, while its loss sees the
label placement `(a,b)` and actual common scale `g`.  The robust-edge-eight
bank closes `2,442/3,003` bodies for arbitrary positive levels, leaving `561`
bodies having at most seven robust edges.

The canonical near miss is to call the zero-channel graph a matching.  It is
not.  At `D=12,m=24`, the distinct levels

```text
32,24,36,27
```

have zero edges `32--24--36--27`.  Thus the correct object must retain both
the reduced-ratio graph and the saturation flag `m=2D`.  Repeated levels are
a separate coordinate: the universal same-level graph is `K6` on every one
of the `561` bodies, so an uncertified word on those bodies has six distinct
levels.

## 2. Sharp phase hierarchy in the C2 cone

For an unequal pair, write its levels as `gP,gQ`, with `P<Q` coprime.  The
cone gives

```text
P >= 2(Q-P),               Q/P <= 3/2.                 (1)
```

The exact primitive fibre from THM-2941 has

```text
F(P,Q) = 1/49 + c(P mod 14,Q mod 14)/(PQ),
c >= -12/49.                                                   (2)
```

The complete `PQ<40` bank subject to `(1)` is

```text
(2,3):0, (3,4):0, (4,5):1/70,
(5,6):2/105, (5,7):1/49.                                  (3)
```

For `PQ>=40`, `(2)` is at least `1/70`; no admissible coprime channel has
`PQ=40`.  Hence the only zeros are `(2,3),(3,4)`, every nonzero channel is
at least `1/70`, and equality is unique at `(4,5)`.  Repeating the argument
at `PQ=110`, with the complete thirteen-channel finite bank below that
cutoff, gives

```text
F(P,Q) >= 1/55 off (2,3),(3,4),(4,5).                   (4)
```

This hierarchy, rather than a single coarse floor, is what makes `C=2`
possible.

## 3. The exact zero graph

Put an edge between distinct levels when their reduced ratio is `4/3` or
`3/2`.

- The `3:4` edges form a matching.  If one level had neighbours on both
  sides, their quotient would be `16/9>3/2`, contradicting `(1)`.
- A `2:3` edge realizes the full quotient `3/2`.  It therefore exists only
  when `m=2D`, and its endpoints must be the unique minimum `2D` and maximum
  `3D`.
- The endpoint edge can meet one `3:4` edge at either end.  Consequently the
  only nonmatching component is a subpath of the hostile `P4` above; all
  other `3:4` edges are isolated matching edges.

In particular the zero graph is triangle-free, but it need not have maximum
degree one.

## 4. Two quality graphs and a complete projective CSP

For labels `a<b` in a body, let `delta=b-a`.  At the worst cone corner
`D_0=6,m_0=12`, define a quality edge at floor `f` by the strict invoice

```text
f - 2D_0(2 delta+b)/(m_0 L-b)
  - sum_(e in E) e/[7(12L-e)] > 0.                       (5)
```

Write `E70` and `E55` for the graphs obtained from `f=1/70` and `f=1/55`.
They satisfy `E70 subset E55`.  Formula `(5)` is uniform over the full cone.
Indeed cross-determinant transport gives

```text
2|eta| <= 2[m delta+Db]/(mL-b).
```

For fixed `D` this decreases in `m`, so its maximum is at `m=2D`; then
`D/(2DL-b)` decreases in `D`.  Singleton debt also decreases coordinatewise.

Any failure must therefore label every `E70` edge by a zero ratio and every
edge of `E55\E70` by one of

```text
2:3, 3:4, 4:5, and their reverses.                       (6)
```

This is a finite exact graph-homomorphism problem.  In each connected
component of `E55`, normalize one vertex to one and grow along a spanning
tree.  Every new value is forced to be its parent times one of the ratios in
`(6)`; all already exposed edges are then checked.  This enumerates every
component realization, not a sample.  Different components have free
relative scales.  They fit distinctly in an interval of quotient `3/2`
exactly when at most one component already spans the full quotient; otherwise
generic rational component scales avoid the finitely many collisions.  Put
the unique full-span component first, if it exists.  A later normalized
component with values `x` may collide with an old value `y` only at the single
scale `y/x`.  When it has `r` vertices and `u` values have already been
placed, there are therefore at most `ru<=r(6-r)<=9` forbidden scales.  Its
admissible scale interval is nondegenerate because it does not span `3/2`.
The referee tries `102` distinct rational scales in that interval, checks the
forbidden set explicitly, and then checks that every scaled value lies in
`[1,3/2]`.  Thus generic placement is a proved finite step, not a heuristic.
Clearing denominators and dilating preserves all ratios, distinctness, and
the global quotient bound, so integrality loses no witness.

The exact body ledger is

```text
561 residual bodies;
539 have an E70 triangle and close because the zero graph is triangle-free;
551 have no realization of the full two-tier CSP;
10 retain a projective obstruction.                              (7)
```

The ten are

```text
(1,2,3,4,6,8)      (1,2,3,4,6,12)
(1,2,3,4,8,12)     (1,2,3,6,8,12)
(1,2,4,5,8,10)     (1,2,4,6,8,12)
(1,2,4,6,9,12)     (1,3,4,6,8,12)
(1,3,4,6,9,12)     (2,3,4,6,8,12).                         (8)
```

Thus `(7)` is not merely a triangle heuristic: it is channel-exhaustive and
keeps disconnected-component scales.

## 5. Eight located selectors

Eight bodies in `(8)` have a fixed selected `E55` edge; on seven it is
already an `E70` edge.  Channels outside `(6)` close by `(4)` and `(5)`.
For each of the six oriented exceptional channels, the referee stores one
body-safe cell and computes its exact primitive skeleton.  If the actual pair
is `g(P,Q)`, then

```text
eta_g = g(Qa-Pb)/(PgL-a),
g/(PgL-a) - (g+1)/(P(g+1)L-a)
  = a/[(PgL-a)(P(g+1)L-a)] > 0.                           (9)
```

Because every level is at least twelve, it is enough to check
`g_0=6,4,3` for the `2:3,3:4,4:5` channels.  The primitive skeleton is
scale-independent, `|eta_g|` decreases by `(9)`, and the six singleton debts
are bounded by setting every level to twelve.  All `48` located margins are
strictly positive.  The weakest is

```text
170554943656361282 / 8980213849475477245
```

on body `(2,3,4,6,8,12)`, slots `(0,1)`, reverse channel `(5,4)`, cell `290`.

The transport homotopy has its own audited gate.  Put `delta=b-a` and

```text
B(D,m)=[m delta+D b]/[mL-b].
```

The two exact finite differences are

```text
B(D,m)-B(D,m+1)
  =b(DL+delta)/[(mL-b)((m+1)L-b)] > 0,

B(D,2D)-B(D+1,2(D+1))
  =b(2delta+b)/[(2DL-b)(2(D+1)L-b)] > 0.               (9a)
```

Hence the full infinite cone is worst at `m=2D,D=6`.  Over all `3,003`
bodies and all label pairs, the resulting exact maximum is

```text
|eta| <= 17/167 < 1
```

on `E=(1,2,3,4,6,12)`, labels `(1,12)`.  If `s in [0,1]` is the homotopy
parameter, then its moving slope obeys the explicit guard

```text
q+s eta >= q-|eta| >= 12-17/167 = 1987/167 > 1.        (9b)
```

Finally the exact transport floor is
`c^{-1}(skeleton-2|eta|)`, where `c=1-a/(pL)` lies in `(0,1)`.  Every located
row first checks

```text
skeleton-2|eta| > singleton debt > 0,
```

then checks `c^{-1}>1`, and only then drops the favourable factor.  The
prefactor inequality is therefore never applied to a zero or negative
bracket.

## 6. The two analytic-tail bodies

The remaining bodies have no edge even in the conservative `E55` graph:

```text
H =(1,2,3,4,6,12),       L=168,  selected labels (1,2);
H2=(1,3,4,6,8,12),       L=336,  selected labels (3,4).  (10)
```

For `H`, put `s=min(P,Q)` and `r=Q/P`.  On
`r in [2/3,3/2]`, convexity gives the exact endpoint maximum

```text
2|r-2| <= 8/3.                                           (11)
```

Also `P,Q>=s`, so `PQ>=s^2`; and the positive common scale gives
`Pg>=P>=s`, so `168-1/(Pg)>=168-1/s`.  The phase-correction bound and exact
cross-determinant therefore give the increasing lower envelope

```text
1/49 - 12/(49s^2) - (8/3)/(168-1/s) - debt_H(12).       (12)
```

At `s=10`, `(12)` is the positive rational

```text
805396046466515008 / 9456367360728129366525.
```

There are exactly `28` oriented primitive channels with `s<10`.  This finite
enumeration is complete: for the threshold `t`,
`max(P,Q)<=3s/2<2t` whenever `s<t`, so the bank is precisely the coprime
oriented pairs in the finite square `1<=P,Q<2t` satisfying the cone ratio.
The referee freezes that exact list and independently gets the same list from
the oversized square `1<=P,Q<=4t`.  Direct exact invoices close `23`; the five
exceptions are

```text
(2,3),(3,2),(3,4),(4,3),(5,4).
```

Located cells `12,129,90,136,142` close them respectively.

For `H2`, the numerator endpoint check is

```text
2|3r-4| <= 4,                                            (13)
```

and `336-3/(Pg)>=336-3/s`.  The analogous increasing envelope is

```text
1/49 - 12/(49s^2) - 4/(336-3/s) - debt_H2(12).          (14)
```

At `s=6`, it equals

```text
43121847752669072 / 90608914671529412235 > 0.
```

The finite head has ten oriented channels, frozen and oversized-audited by the
same rule.  Six close generically; located cells `24,284,311,290` close
`(2,3),(3,2),(3,4),(4,3)`.  All inequalities are strict.  Finally, if `A` is
the appropriate endpoint numerator bound, the two successive loss drops are

```text
12/49 [1/s^2-1/(s+1)^2] > 0,

A s/(Ls-a)-A(s+1)/(L(s+1)-a)
  =Aa/[(Ls-a)(L(s+1)-a)] > 0.                           (15)
```

Thus `(12)` and `(14)` increase for every subsequent integer `s`, rather than
only at the checked endpoints.  The referee checks the first exact increments
as positive rationals and freezes them in its semantic image.

This closes the last two rows of `(8)` and proves the candidate theorem.

## 7. Structural connections and honest stopping boundary

The incoming THM-3045 K4-edge clutch
(`01-canon/theorems/THM-3045-k4-edge-isotypic-binary-ternary-integral-clutch.md`)
gives a useful but noncanonical picture.  Identify the six levels with the six
edges of `K4`.  Declare each
disjoint `3:4` edge to be an opposite-edge pair and complete to the three-pair
opposition involution.  The unique `2:3` edge cannot be one of those pairs, so
it is an octahedral adjacency.  The hostile is therefore exactly

```text
opposition -- adjacency -- opposition.
```

This explains the P4 anatomy and the recurring binary/matching quotient.
It is not a proof dependency: the identification is chosen after seeing the
zero graph, and a body-quality triangle is not canonically a K4 star or face.
No overlap certificate is preserved by pretending otherwise.

An incoming finite audit of the `m=1`, levels-in-`1..10` box gives the
complementary signal.  On `H`, one exact word
`(10,4,6,3,2,1)` defeats the best high-edge-only invoice, while a common-cell
two-star Hunter tree has positive aggregate margin.  This lies outside
`m>=2D`.  It suggests the correct continuation: the C2 cone admits a
one-pair selector because its ratio alphabet collapses to `(6)`; below the
cone, retain common-cell tree credit rather than minimizing every high edge
independently.

## 8. Exact referee and audit boundary

Run

```bash
python3 04-computation/lrc14_j7_reflected_c2_full_cone_closure_thm2941.py
python3 -O 04-computation/lrc14_j7_reflected_c2_full_cone_closure_thm2941.py
```

The referee pins the phase, universal repeated-level, low-phase, and
robust-edge-eight sources, the edge-eight transcript, and the edge-eight
semantic digest `af12592e...`.  It freezes the complete 561-body semantic
image, both quality graphs, every projective component realization, all
`48+5+4` located controls, both finite oriented heads, both analytic-tail
seams, the homotopy floor, and the component-scaling bound.  Fresh ordinary,
optimized, `PYTHONHASHSEED=1`, and optimized `PYTHONHASHSEED=8675309` replays
are byte-identical to one another and to the stored transcript.

```text
source:   f7eb06c50b11fd810320146a2cded2590790fdbcf6db58cff2a358470fcfd0b6
output:   7449cb0ab6b969fa4653d98582423df262329e087cc73d8ea33499cb23ec8422
semantic: e052683e114b067648c61e5451001a4871d9940eb465eb4bf2522537d1cbe665
body:     1f382f27bde53c7afac3ae9cfcd5b79b14ab2b63242b8c3cd9a3ee37976c7a61
located:  1ef0c2811df1fe4426ddd97ae142b97325af0bbce67fd7cca8f69389d5755f32
```

The candidate has passed this independent logical hardening but still needs a
canon promotion decision and integration into THM-2941's assembled dependency
graph.  After that promotion, the reflected certificate-failure frontier is
the `561`-body wedge `D>=6,1<=m<2D`, not LRC(14) itself.
