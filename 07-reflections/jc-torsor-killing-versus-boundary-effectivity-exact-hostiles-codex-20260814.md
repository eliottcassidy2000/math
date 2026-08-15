# JC effectivity: killing the torsor does not kill the boundary poles

**Date:** 2026-08-14  
**Status:** PROVED elementary filtration and valuation-lattice statements for
the THM-3383 terminal monomial family and THM-2655's explicit non-Keller
`V4` quotient; VERIFIED-EXACT companion; unnumbered.  This supplies a
candidate payload for the reserved THM-3397 stub but does not promote that
stub, exclude `A4/S4`, classify G1, or prove JC(2).

## 1. Inheritance and connection contract

THM-3383 is the closest proved mechanism: a finite residue torsor decodes
rationally while one target coordinate remains outside `C[x,y]`.  THM-2655
is the canonical hostile: the even-sign-change `V4` quotient realizes the
required Kummer/class-group carrier sharply, but its induced quartic map has
nonconstant Jacobian.  The corrected near miss is “a killed `H^1` or class
group element has an effective polynomial section.”  The least-used sidecar
is the full vector of boundary valuations.

The exact map is:

| field | connection |
|---|---|
| source | a cyclic or `V4` quasi-etale torsor over a normal invariant ring |
| target | the polynomial total-space cover together with selected rational target views |
| map | pull back Kummer characters/divisor classes and rational functions |
| preserved | function-field identities, finite group action, discriminant/class carrier |
| destroyed | sign and magnitude of valuations along individual boundary components |
| required sidecar | the oriented valuation vector and the effective cone |
| cheapest decisive test | `d/a=yz/x`, `d/b=xz/y`, `d/c=xy/z`: every class dies upstairs and every pair product is polynomial, but no individual ratio is |

This is the JC-side version of compatibility/closure/effectivity.  It is not
an LRC-to-JC `H^1` map.

## 2. Exact terminal denominator filtration

Use THM-3383's notation

```text
n=g-ae in {+1,-1},       v=ut,       L=1+cv=xy.         (1)
```

Let `s` denote the decoded target which is polynomial and `m` the missing
target:

```text
n=+1:  s=u, m=t, d=ae=g-1;
n=-1:  s=t, m=u, d=g=ae-1.                              (2)
```

Both orientations have the same normal form in `C(x,v)`:

```text
s=x^e v^(d+1),              m=x^(-e)v^(-d),
sm=v.                                                        (3)
```

The target polynomial ring has the exact Laurent decomposition

```text
C[x,y]
= direct_sum_(k>=0) x^k C[v]
  direct_sum_(r>=1) x^(-r)L^r C[v].                     (4)
```

Therefore, for every `q>=1` and `f(v) in C[v]`,

```text
m^q f(v) in C[x,y]
iff v^(dq)L^(eq) divides f(v).                           (5)
```

There are two independent requirements.  The factor `v^(dq)` clears the
pole already present in the Laurent coefficient, and `L^(eq)` is the target
boundary divisibility required by `(4)`.  They are coprime because
`L=1+cv` and `c` is nonzero.

The filtration ideals are consequently

```text
I_q=(v^(dq)L^(eq))=(v^dL^e)^q.                          (6)
```

They are exactly multiplicative, not merely submultiplicative.  Equality is
attained by

```text
m^q v^(dq)L^(eq)=x^(-eq)L^(eq)=y^(eq),                  (7)
```

while removing one available `v` factor or one `L` factor leaves a literal
pole.  Thus the entire denominator filtration is generated in degree one by
the boundary vector `(d,e)`.

Let

```text
S=C[x,y]^(mu_e)=C[v,x^e,y^e].                           (8)
```

The polynomial map `Spec C[x,y] -> Spec S` is the cyclic quotient cover;
over the regular locus it kills the `mu_e` torsor.  Nevertheless the missing
target `m=x^(-e)v^(-d)` is invariant, lies in `Frac(S)`, and remains outside
`C[x,y]` by `(5)`.  This is a complete cyclic hostile to

```text
torsor trivial upstairs  =>  decoded target regular upstairs. (9)
```

## 3. The three-view `V4` hostile

Take THM-2655's even-sign-change quotient

```text
A=C[a,b,c,d]/(d^2-abc),
a=x^2, b=y^2, c=z^2, d=xyz,
B=C[x,y,z].                                             (10)
```

The finite map `Spec B -> Spec A` has group `V4`; it is quasi-etale and
becomes a connected torsor on the regular locus.  Its three nonzero Kummer
characters are represented by the squareclasses of `a,b,c`, with

```text
[a]+[b]+[c]=0                                           (11)
```

because `abc=d^2`.  Their class-group images are the three nonzero elements
`P_a,P_b,P_c` of `Cl(A)=(Z/2)^2`.  Pullback kills all of them:

```text
a=x^2, b=y^2, c=z^2,
P_a -> div(x), P_b -> div(y), P_c -> div(z).             (12)
```

Now retain the three symmetric invariant rational views

```text
rho_a=d/a=yz/x,
rho_b=d/b=xz/y,
rho_c=d/c=xy/z.                                         (13)
```

Relative to the coordinate boundary divisors `H_x,H_y,H_z` upstairs, their
valuation vectors are the rows of

```text
R = [ -1  1  1 ]
    [  1 -1  1 ] = J-2I.                               (14)
    [  1  1 -1 ]
```

Each row has one distinct pole, so no `rho_i` lies in `B`.  Nevertheless

```text
rho_a rho_b=c,
rho_a rho_c=b,
rho_b rho_c=a,
rho_a rho_b rho_c=d,                                    (15)
```

so every pair and the triple are polynomial.

## 4. The valuation matrix recovers the `V4` lattice exactly

The matrix in `(14)` has

```text
det R=4,                    SNF(R)=(1,2,2).              (16)
```

Hence `Z^3/RZ^3=(Z/2)^2`.  This is not a numerical echo: it is exactly the
class relation lattice from THM-2655.  If `S=(1,1,1)` and `e_i` are the
coordinate basis vectors, then

```text
row_i(R)=S-2e_i.                                        (17)
```

The class relation lattice is generated by

```text
2e_1,2e_2,2e_3,S.                                      (18)
```

Every row in `(17)` belongs to `(18)`.  Both lattices have index four, so
they coincide.  Thus the same three principal ratio divisors which exhibit
the boundary poles generate the full `V4` class relation lattice.

The missing predicate is the effective cone `N^3 subset Z^3`.  Class-group
closure sees vectors modulo the lattice; polynomiality asks whether the
actual oriented vector lies in the cone.  Every row of `R` lies outside it,
while each pair sum and the triple sum lie inside:

```text
r_a+r_b=2e_3,
r_a+r_c=2e_2,
r_b+r_c=2e_1,
r_a+r_b+r_c=S.                                         (19)
```

This is the exact distinction between a closed abstract carrier and an
effective polynomial realization.

## 5. Why this is not a tournament

If the vertices are the three views `rho_a,rho_b,rho_c` and the binary
observable is “the product is polynomial,” every pair is a tie-positive
edge.  That complete pair shadow loses the negative singleton valuation
which is the whole obstruction.  Orienting the edges supplies no missing
information.

The natural finite object is instead a Boolean effectivity complex:

```text
singletons: non-effective;
pairs:      effective;
triple:     effective.                                 (20)
```

Its sidecar is the signed valuation matrix `(14)`.  This is another precise
instance where a tournament with missing/both-way edges is only a quotient
of the useful carrier.

## 6. Relation to the three `disc=-4(square)^2 L` cubics

THM-2546 proves for the fixed sporadic grade-three map that all three
coordinate cubics have

```text
disc=-4(square)^2 L.                                   (21)
```

This synchronizes their branch divisor and is genuine fixed-map structure.
It does not classify a degree-four Keller map.  A quartic and its cubic
resolvent share a discriminant, but the transfer to a polynomial Keller
resolvent still has to preserve the boundary-effective data.  Equations
`(13)--(19)` show that the abstract `V4` kernel, all three Kummer classes,
and the full class lattice can be present and can trivialize on a polynomial
cover while three natural resolvent-side views remain nonpolynomial.

The determinant `4` in `(16)` is the exact order/index of this `V4` lattice.
It is not the factor `-4` in `(21)`, not a Jacobian determinant, and not a
degree-four counterexample.  The proved composite discriminant law in
THM-2582 further warns that even the grade-two square class changes under
composition.  Any exact classification of the Keller monoid still needs at
least

```text
degree grade + monodromy blocks + Jelonek divisor tower
+ oriented boundary-valuation/effectivity data.          (22)
```

## 7. An explicit JC-side cohomology cospan

Within this `V4` hostile one can write the classes explicitly:

```text
W=Hom(V4,C2) = <[a],[b]>,       [c]=[a]+[b],
W -> H^1_et(A_reg,mu_2) -> Cl(A)[2],
pullback to B: [a],[b],[c] -> 0.                         (23)
```

The valuation vectors of the associated rational views land instead in
`Z^3`, and polynomiality is membership in `N^3`.  Thus the lawful diagram is
a typed cospan from the Kummer carrier and from rational functions into a
common boundary valuation lattice.  There is no direct identification of
this `W` with the LRC word-current class, Berggren ancestry, or Hamiltonian
flux.

## 8. Reproduction and scope

Reproduce with

```text
python 04-computation/jc_torsor_effectivity_boundary_probe_20260814.py
python -O 04-computation/jc_torsor_effectivity_boundary_probe_20260814.py
```

The ordinary and optimized payloads agree.  The semantic digest is
`90e478b47b925203d81978048dcc25d76b97258bbc413982016db7a13d594a18`.
The companion checks `135` normalized terminal cells and `675` filtration
gates, the exact `V4` invariant ratios, all subset-product effectivity laws,
the determinant/Smith data, and equality of the valuation and class relation
lattice indices.

The result is sharp for the two named families.  It proves no general
boundary-effectivity theorem for arbitrary terminal modules, no physical
quartic resolvent embedding, and no `A4/S4`, G1, JC(2), or LRC exclusion.
