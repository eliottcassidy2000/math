---
id: THM-1183
title: The general-d continuum ceiling is exactly a six-order torus X-ray; three disjoint directed paths reduce 2/21 to one explicit 1/21 inequality
status: PROVED REDUCTION, OPEN FINAL INEQUALITY — the functional identity, universal mirror symmetry, six-order decomposition, three-path domination, proportional equality, and an exact d<=40 sweep are proved.  The displayed all-integer three-path inequality is not yet proved, so THM-1174's uniform ceiling remains open.
source: codex-2026-07-18 (exact analytic audit of THM-1172/1173/1174)
depends_on: [THM-1172, THM-1173, THM-1174]
script: /tmp/lrc14_continuum_ceiling_exact_codex_s1183.py (draft; repo intentionally untouched)
---

# THM-1183 — exact carrier and the remaining three-path inequality

Let `0<a<b<c` be the drift rates and put

```text
P_d(u) = {0, {a u}, {b u}, {c u}} subset R/Z,
M_d(u) = largest circular spacing of P_d(u).
```

Endpoints have measure zero throughout, so closed/open conventions do not
affect any measure statement.

## 1. Exact functional form

The continuum tooth for drift `q`, in the normalized `[0,1]` gap, is

```text
T_q(u) = [(7/6){-q u}-1/6, (7/6){-q u}] intersect [0,1].       (1)
```

Apply `x -> (6/7)x` and restore the base tooth
`[-1/7,0]` on the unit circle.  The four teeth are the oriented arcs

```text
A_v(u) = [{-v u}-1/7, {-v u}],       v in {0,a,b,c}.            (2)
```

If `Delta_1,...,Delta_4` are the circular spacings of their right endpoints,
the surviving components have circle lengths `(Delta_i-1/7)_+`.  Since
`max Delta_i >= 1/4`, the longest normalized survivor is therefore exactly

```text
F_d(u) = (7/6)(M_d(u)-1/7).                                    (3)
```

Consequently

```text
F_d(u) <= 1/6  iff  M_d(u) <= 2/7.                              (4)
```

If all four spacings are at most `2/7`, each is also at least
`1-3(2/7)=1/7`.  Thus the bad set is exactly

```text
B_d = {u: every circular spacing of P_d(u) lies in [1/7,2/7]}. (5)
```

This also proves the mirror symmetry for **every** drift triple, not only
`(1,2,3)`:

```text
F_d(1-u)=F_d(u),                                                (6)
```

because `P_d(1-u)=-P_d(u)` is a reflection of the circle.

## 2. Six cyclic orders, exactly

Write `I=[1/7,2/7]` and `f=1_I`.  For a cyclic order
`pi=(0,v1,v2,v3)`, define its signed edge differences

```text
w_pi=(v1, v2-v1, v3-v2, -v3),       sum_i w_pi,i=0.             (7)
```

Then its bad chamber is exactly

```text
E(w)={u: {w_i u} in I for i=1,2,3,4}.                           (8)
```

Indeed the four displayed fractional parts lie in `I`, and their sum is an
integer in `[4/7,8/7]`, hence is `1`; they are therefore precisely the four
positive cyclic spacings in that order.  The six order chambers are disjoint
apart from endpoints and cover `B_d`.

Take these three representatives modulo reversal:

```text
wA=( a, b-a, c-b, -c)       order 0,a,b,c,
wB=( a, c-a, b-c, -b)       order 0,a,c,b,
wC=( b, a-b, c-a, -c)       order 0,b,a,c.                      (9)
```

Reflection pairs each with its reverse, so

```text
mu(B_d)=2[mu E(wA)+mu E(wB)+mu E(wC)].                          (10)
```

Equation (10) is an exact rational evaluator, not a grid approximation.

## 3. A three-path majorant

Delete the first edge of each cycle in (9).  Define

```text
P1={u: { (b-a)u},{ (c-b)u},{-c u} all lie in I},
P2={u: { (c-a)u},{ (b-c)u},{-b u} all lie in I},
P3={u: { (a-b)u},{ (c-a)u},{-c u} all lie in I}.                (11)
```

These are the directed Hamiltonian paths

```text
a -> b -> c -> 0,
a -> c -> b -> 0,
b -> a -> c -> 0.                                               (12)
```

Each cycle event is contained in its path event.  Moreover `P1,P2,P3` are
pairwise disjoint off endpoints: three increments in `I` have sum in
`[3/7,6/7]<1`, so each event forces the displayed clockwise order, and the
three orders differ.  Therefore

```text
mu(B_d)/2 <= mu(P1)+mu(P2)+mu(P3).                              (13)
```

The entire analytic ceiling is now reduced to the following precise claim:

```text
PATH(a,b,c):  mu(P1)+mu(P2)+mu(P3) <= 1/21
              for every 0<a<b<c in Z.                          (14)
```

`PATH(a,b,c)` implies `mu(B_d)<=2/21` immediately by (13).  This is the
remaining unproved inequality; none of the exact computations below is a
substitute for its proof.

For checking, if `q>0`,

```text
C(q)  = union_{n=0}^{q-1} [(7n+1)/(7q),(7n+2)/(7q)],
C(-q) = union_{n=0}^{q-1} [(7n+5)/(7q),(7n+6)/(7q)].            (15)
```

Thus every term in (8) or (11) is an intersection of three or four explicit
rational interval lists.  Formula (15) is a finite exact decision procedure
for any proposed counterexample.

## 4. The proportional family, exactly

For `(a,b,c)=m(1,2,3)`, put `x={m u}`.  Chamber A has

```text
x in I, {-3x} in I  iff  x in [5/21,2/7],
```

and hence measure `1/21`; chambers B and C contain opposing `x,-x`
conditions and have measure zero.  Therefore

```text
mu(B_{m,2m,3m})=2/21                                           (16)
```

for every `m>=1`.  This proves attainment.  Uniqueness of attainment for all
integer triples is **not** claimed here.

## 5. Exact finite evidence and corrections to THM-1174

The companion `Fraction` evaluator checks all `9880` triples
`1<=a<b<c<=40`:

```text
actual mu(B)>2/21:                    0
three-path majorant >1/21:            0
actual equality rows:                 exactly (m,2m,3m), 1<=m<=13
largest nonproportional actual value: 4/105
```

This corrects two points in THM-1174 without proving its final ceiling:

1. mirror symmetry is universal by (6);
2. the reported values `0.0960,0.0967,...` on proportional rows are grid
   endpoint artefacts.  Their exact value is always `2/21` by (16).

## 6. Torus/Kakeya and tournament/carrier audit

In `T^3`, `u -> (a u,b u,c u)` is a rational geodesic.  For order
`0,a,b,c`, the bad chamber is the `1/7`-scaled Kuhn orthoscheme

```text
conv{(1,3,5),(2,3,5),(2,4,5),(2,4,6)}/7,                      (17)
```

and the other five are its coordinate permutations.  Hence `mu(B_d)` is the
one-dimensional toral X-ray/Kakeya density of six explicit tetrahedra.  A
direction-free volume argument cannot see the extremal resonance: their
ambient Haar volume is only `1/343`, while the AP geodesic has X-ray density
`2/21`.

The faithful pairwise carrier has runner vertices and observable

```text
theta_ij(u)={(v_j-v_i)u};   i -> j iff theta_ij(u) in I.        (18)
```

The switch/gauge `u -> 1-u` reverses every edge.  On a bad configuration the
relation is exactly a directed `C4`: score histogram `(1,1,1,1)`, no directed
triangles, one SCC of size four, four directed Hamiltonian paths, and zero
edge changes inside a bad component.  The three path carriers in (12) are
the spanning-tree bases obtained by deleting one edge from three reversal
representatives of that cycle.  Increasing velocity order is the fixed tie
Hamiltonian path at wall coincidences.

Assumption challenge:

- runner-only winding tournaments preserve cyclic orientation but destroy
  the metric window `[1/7,2/7]`, so they cannot decide (14);
- the four spacing vertices preserve the predicate exactly and expose the
  simplex, but forget runner ownership;
- seven fixed sectors/residues destroy within-sector slack and are only a
  discretization;
- wall events preserve the exact interval evaluator but obscure the torus
  symmetry;
- Fourier modes encode the rational geodesic resonance, while the three
  directed path obligations are the smallest currently faithful proof
  carrier.

The sharply narrowed obstruction is therefore (14), equivalently the X-ray
bound for the union of the three path relaxations.  Proving (14), or finding
one exact triple violating it via (15), is the next analytic step.
