# Independent all-parity norm-five closure: different network and physical extrema

**Status: PROVED ELEMENTARY + FINITE-EXACT; independent proof/head audit PASS.**
The formulas and complete finite head below give an independent proof of
the signed coefficient-magnitude `(1,2,2)` family. Incoming
[THM-4441 — signed-122 sharp ray closure](../../01-canon/theorems/THM-4441-lrc14-signed-122-sharp-ray-closure.md)
was **RESERVED / UNPROVED EMPTY STUB** at origin commit `d0f9383875` and
was promoted **PROVED** at `058a8ded98`, with the same sharp constants and
a larger native replay. Its concurrent owner retains that namespace. This
note is independent co-discovery/referee evidence and promotes no second
ID. The historical reservation was never a proved dependency.

For every primitive sorted distinct positive ternary-unit triple `w`
admitting a signed `(1,2,2)` relation, the sharp conclusions are

```text
min_i E_i(w)<=46/665<6/77,
mu(F_w)<=51/770<46/665,
```

with unique primitive equality triples `(2,19,20)` and `(1,11,20)`,
respectively. The network and physical extrema are **different**. This is
a local circuit closure, not an arbitrary-body entry or synchronization
theorem.

## 1. Inheritance and the missing residue coordinate

The incoming
[THM-4437 — all-parity reduction to three low circuits](../../01-canon/theorems/THM-4437-lrc14-all-parity-network-reduction-to-three-low-circuits.md)
is now **PROVED** and leaves this family explicit. The
[incoming parity probe](lrc14_parity_empty_core_sep06.md) has finite norm-five
peak `46/665` at `(2,19,20)`, but does not supply the all-height statement.
The additive family is separately closed; no additive theorem is claimed
again here.

The closest proved mechanism is
[THM-4386 — zero-defect incidence](../../01-canon/theorems/THM-4386-lrc14-canonical-component-relation-and-zero-defect-incidence.md),
Lemma 2.1. Its statement is odd-typed; its actual modular argument is
reproved without parity below. The analytic sidecar is the complete-ray
quadrature and concave-roof error in the audited
[boundary-discrepancy note](lrc14_sum_ray_boundary_discrepancy_overnight_hexagon_sep05.md),
Section 5, recovering THM-4373's period-three primitive.

The canonical hostile to naive locking is precise. For

```text
w=(2,19,20), d=(1,2,-2), C=(-8,4,-3),
w dot C=0, d cross C=w,
14|C_i|<3(sum(w)-w_i) for every i.
```

This integer geometric carrier has affine defect one; the uncolored short
relation does **not** lock, since `5*(3/14)>1`. It is excluded only by the
actual owner condition `C_3=0 mod3`. The corrected claim keeps that missing
coordinate rather than treating norm five as an uncolored light relation.

The live board is: owner-live relation; complete affine defect; parameter
product; zero-boundary projection; central cube section; physical versus
network extremum. The source-to-target map preserves every live multiplier,
strict boundary, and physical interval. Dropping owner residues admits the
hostile above; dropping physical overlap confuses the two sharp extrema.

## 2. Complete ray confinement for all parities

Let `d` be the signed primitive `(1,2,2)` relation, `d dot w=0`. Every
coordinate of both `d` and an actual live carrier `C` is a ternary unit.
The triples `(w_i d_i)` and `(w_i C_i)` consist of nonzero residues summing
to zero modulo three, so each is constant. Thus `d,C` are proportional
modulo three and `d cross C=0 mod3`.

Both lie in the plane perpendicular to primitive `w`, so
`d cross C=delta w` with `delta` an integer. This integrality follows by
applying an integer Bezout functional taking value one on `w`. Since `w`
has ternary-unit coordinates, `3|delta`.

The parity-free cube/Helly identity supplies `C=w cross e` with every
`|e_i|<3/14`. Hence `delta=d dot e` and

```text
|delta|<3||d||_1/14=15/14<3.
```

Therefore `delta=0`. Primitivity of `d` makes `C=k d` with integer `k`;
the owner gate is precisely `3∤k`. The full dictionary is

```text
Omega={k d:0<|k|<B, 3∤k},
B=min_i 3(sum(w)-w_i)/(14|d_i|).                       (1)
```

Conversely every multiplier in (1) satisfies every strict roof and owner
gate, so this is equality of full sets, including empty supports. No ray
or affine defect has been silently discarded. This reproduces the modular
gain in THM-4386; it also works through norm fourteen, but only norm five
is consumed here.

## 3. Two exhaustive parameter families

Positive speeds force one sign of the relation to be opposite the other
two. If the singled coefficient has magnitude one, the row is

```text
I: (a,b,c)=(a,b,2(a+b)), a<b, gcd(a,b)=1, h=ab,
d=(2,2,-1).                                          (2)
```

If the singled coefficient has magnitude two, write the relation as
`A+2p=2q`, so `q=p+d0`, `A=2d0`, with `p,d0>0`. The row is

```text
II: sort(2d0,p,p+d0), gcd(d0,p)=1, h=d0*p.             (3)
```

Keep only distinct ternary-unit speeds in either family. The relation in
the unsorted order `(A,p,q)` is `(1,2,-2)`. These two cases exhaust every
signed `(1,2,2)` relation up to permutation and simultaneous sign. The gcd
conditions are exact: in II, `gcd(2d0,p,p+d0)=gcd(d0,p)`.

## 4. A zero-boundary projection and the network product cutoff

Write `r=3/14` and `q0=3/(7 max(w))`. For any coordinate put

```text
f_i(t)=min(q0,[r(w_j+w_k)-|d_i|t]/(w_jw_k)).
```

The complete network is `E_i=2 sum_(0<k<B,3∤k) f_i(k)`. Choose the
coordinate attaining the tight roof in the following way.

In family I, choose the speed `c=2(a+b)` (coefficient magnitude one).
Then

```text
B=3c/28, slope=1/(ab)=1/h,
R_c=(4/3) integral_0^B f_c
   =3/49-3ab/[98(a+b)^2] <=3/49.                      (4)
```

In family II, choose speed `q=p+d0`, even when it is not the largest
speed. Let `c=max(A,q)`, `A=2d0`. The other two roof bounds are larger,
and

```text
B=3(A+p)/28, slope=2/(Ap)=1/h,
R_q=(3/49)[(A+p)/c-Ap/c^2] <=3/49.                    (5)
```

For the final inequality, if `c=A` the bracket is exactly one. If `c=q`,
put `t=p/q`; then `1/2<=t<1` and the bracket is `2-3t+2t^2<=1`.
Both chosen roofs vanish at `B`, so there is no boundary jump. Their
plateau ends at `B-q0*h`, a point in `[0,B]`.

Let `P` denote the inherited period-three primitive with `-1/3<=P<=0`.
The exact selected-coordinate formula is

```text
E_selected=R_selected+(2/h)[P(B)-P(B-q0*h)].
```

Thus both families satisfy

```text
min_i E_i <= E_selected <=3/49+2/(3h).                (6)
```

For `h>=84`, the last expression is at most
`61/882=46/665-1/83790`, strictly below the claimed sharp constant.
It remains to check only the complete parameter product head `h<84`.

## 5. A constant physical coarea term and a smaller product cutoff

The continuous physical roof is `lambda(t)=min_i f_i(t)`. Its bulk term
is independent of the speed ratio; it depends on the coefficient normal.
To see this without guessing density, choose a real vector `v` with
`w cross v=d`. Then `(y,k) -> yw+kv` parametrizes `d^perp`, with area
Jacobian `||d||`. The literal overlap length for real carrier `k d` is
the length of the allowed `y` interval. Consequently

```text
integral_(k in R) lambda(|k|) dk
 = area([-r,r]^3 intersect d^perp)/||d||.              (7)
```

Signs and coordinate permutations of `d` leave this central cube section
unchanged. For coefficient magnitudes `(1,2,2)`, eliminate a coefficient-two
coordinate. The other two coordinates satisfy `|x/2+y|<=r`, so

```text
area(section)/||d||
 =(1/2) integral_(-r)^r (2r-|x|/2) dx=7r^2/4.
```

The density of retained integer multipliers is `2/3`, giving

```text
(4/3) integral_0^B lambda(t) dt=7r^2/6=3/56.           (8)
```

This is a continuous identity, not yet a discrete bound. The missing
arithmetic remainder is supplied by the concave-roof theorem. Its maximum
slope is exactly `1/h`: in I, the largest `|d_i|w_i` is `c`; in II, it is
`2q`. Hence

```text
|mu(F_w)-3/56|<=2/(3h).                               (9)
```

For `h>=53`, the upper side is at most
`589/8904=51/770-41/489720`, strictly below the physical maximum.
The same `h<84` head contains every remaining physical obligation.

## 6. Complete exact heads and distinct equality consumers

Enumerate positive parameter pairs of product below 84, apply the exact gcd,
distinctness and ternary-unit conditions, and in I keep `a<b`. This produces
41 family-I rows and 82 family-II rows, 123 distinct primitive triples.
The independent implementations compare **every full raw dictionary**, all
three network projections, and actual physical mass with literal six-sheet
intervals. All cases pass the claimed bounds, with unique equalities

```text
w=(2,19,20): E=(48/665,11/140,46/665), H=173/2660;
w=(1,11,20): E=(51/770,3/35,3/35),     H=51/770.
```

The strict product tails (6) and (9) prove the all-height conclusions and
exclude every larger-product equality. The first row maximizes the best
network projection, but not actual overlap; the second maximizes actual
overlap. In particular no all-three-coordinate `46/665` claim is made.

Common ternary-unit dilation preserves physical Haar mass, so the physical
result extends after primitive reduction, with equality precisely at
dilates of `(1,11,20)`. A body-specific Haar consumer can use `51/770`, but
that number is not a floor for an arbitrary actual body. The earlier body
counterexamples and exact endpoint compiler remain relevant; this circuit
closure alone does not prove entry or LRC(14).

## Reproduction and concurrent ownership

```bash
python3 -B 04-computation/lrc14_norm_five_overnight_hexagon_sep05.py
python3 -B -O 04-computation/lrc14_norm_five_overnight_hexagon_sep05.py
```

The [source](../../04-computation/lrc14_norm_five_overnight_hexagon_sep05.py)
uses independent raw and literal engines from the earlier one-ray source.
It separately integrates the full physical lower envelope by exact
trapezoids, checking (8) against the coarea proof. The
[output](lrc14_norm_five_overnight_hexagon_sep05.out) records the complete
head, uncolored affine-defect hostile, and positive large-product controls.
Normal and optimized runs match: 1,680 primary and 3,555 inherited
raw/literal gates. Independent `observer_collision` written proof and
complete-head audit: **PASS**. The referee checked modular locking, both
families, zero-boundary roofs, the coarea Jacobian, strict product tails,
distinct equalities and the uncolored affine-defect hostile, then separately
replayed all 123 rows and every gate.

The concurrent owner has now proved THM-4441 for this exact signed-122
target. This independent product-head and coarea certificate is a referee
and alternative proof of that result; no new ID is claimed here. Root also
read Sections 1--6 and independently audited the modular lock, exhaustive
families, zero-boundary projection, coarea, distinct extrema and strict
product tails: **PASS**.

Frozen raw-LF SHA256:

```text
source 781f94c37f4f437a0023a03c8ef4acda1071fe6a41e34243482b4381f27cc4b8
output bd62f6ceb568a566945dc4314c8d25a3575e969e0d2ea9de89bd36874348557f
```
