---
id: THM-3935
title: "Linear-conic resolvent has class group Z/3 and a unique cubic character"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  The normal quadratic resolvent of the THM-3932 fold-three sextic has scalar
  units and class group Z/3.  Its natural Cardano divisor is the generator,
  not one of several local three-torsion directions.  Consequently the
  natural cyclic cover is the only nontrivial smooth-locus (equivalently,
  codimension-one unramified) C3 character up to inversion, and it recovers
  the already monogenic THM-3932 cubic.
source: jc_degree6_one_place / post-THM-3932 global resolvent audit, 2026-08-23
depends_on:
  - THM-3932-infinity-component-linear-conic-torus-sextic-fold-classification
related:
  - THM-3912-even-degree-split-boundary-cusp-three-torsion-sieve
  - THM-3922-affine-plane-completion-free-boundary-class-group-obstruction
  - THM-3931-degree-two-pole-cubic-principal-ramification-no-atlas
script: 04-computation/jc2_linear_conic_resolvent_class_unique_character_thm3935.py
output: 05-knowledge/results/jc2_linear_conic_resolvent_class_unique_character_thm3935.out
script_sha256: 61c59e636d9122bf6118ec59df962169b5bee80e171d975ac73ed29cd6361096
output_sha256: 4a68ab996bd24301194c351572a7b8cd846ca9f51318fb0e1c086d0ebee36c0a
semantic_sha256: df63f1e9b60d0385c3f736b11c51493520a4115efb27206586208398bb481798
hash_basis: raw LF bytes
---

# THM-3935 -- the local four-line choice collapses to one global character

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**
Work over an algebraically closed field `k` of characteristic zero.  Put

```text
q=t^3-X^2,
H=q^2-4X^3,
B=k[t,X,W]/(W^2-H),                 S=Spec(B).             (1)
```

Then `B` is a normal domain with unique singular point `(t,X,W)=(0,0,0)`,
and

```text
B*=k*,                              Cl(B)=Z/3.             (2)
```

If

```text
D+=(X,q+W),                         D-=(X,q-W),            (3)
```

then either `[D+]` or `[D-]` generates `Cl(B)`, and

```text
div(X)=D++D-,
div(q+W)=3D+,
div(q-W)=3D-.                                             (4)
```

In particular, the order-three class exhibited in THM-3932 is the entire
global three-primary supply.  If `Sreg` is the smooth locus, then

```text
Pic(S)=0,                 H^1_et(S,mu_3)=0,
Pic(Sreg)=Z/3,            H^1_et(Sreg,mu_3)=Z/3.           (5)
```

The two nonzero smooth-locus characters are inverse generators of the same
cyclic cover.  That cover is precisely the natural Cardano cover

```text
z^3=(q+W)/2,        v^3=(q-W)/2,        zv=X,             (6)
u=z+v,
u^3-3Xu-(t^3-X^2)=0.                                     (7)
```

Consequently there is no second normal generically `S3` cubic field, and no
second normal nonmonogenic cubic order, with this exact squarefree branch and
quadratic resolvent: after quadratic base change its cyclic character must be
`(6)` or its inverse, so its Galois closure and cubic subfield are the ones in
`(7)`.  The normal integral closure in that cubic field is the already
monogenic THM-3932 order.

This conclusion is deliberately about normal cubic completions with this
exact branch/resolvent.  It does not classify nonnormal suborders, altered
branch multiplicities, or cubic constructions having merely the same local
analytic singularity.

## 1. Normality and the exact vertical-prime ledger

The equation in `(1)` is irreducible.  Its partial derivatives are

```text
2W,
4X(3X-X^2+t^3),
-6t^2(t^3-X^2).                                          (8)
```

Together with the surface equation they vanish only at the origin.  Indeed,
the `X=0` row forces `t=0`; in the `t=0`, `X!=0` row equation `(8)` forces
`X=3`, which does not lie on `H=0`.  In the remaining `t!=0`, `q=0` row,
the `X` derivative becomes `12X^2`, so `X=0` and then `q=t^3=0`, again a
contradiction.  Thus this irreducible hypersurface is regular in codimension
one and Cohen--Macaulay, hence normal by `R1+S2`.

Let `K=k(t)` and let `C=S_K`.  Every closed fibre of `S -> A1_t` is
irreducible.  With `T=lambda^3`, its quartic is

```text
H_T=X^4-4X^3-2TX^2+T^2.                                 (9)
```

If `(9)` were a square, it would equal `(X^2+bX+c)^2`.  The `X^3`, `X`,
`X^2`, and constant coefficients successively force

```text
b=-2,            c=0,            T=-2,            0=T^2,
```

a contradiction.  Hence the only vertical height-one prime over
`t=lambda` is

```text
V_lambda=(t-lambda),              div(t-lambda)=V_lambda. (10)
```

The Weil localization sequence, with its unit terms retained, is

```text
B* -> Gamma(C,O_C)* -> direct_sum_lambda Z[V_lambda]
   -> Cl(B) -> Pic(C) -> 0.                              (11)
```

All generators in the middle group already have principal preimages by
`(10)`.  Therefore

```text
Cl(B) ~= Pic(C).                                         (12)
```

This is the point at which no hidden vertical class or deleted fibre may be
discarded: `(9)` proves that the displayed ledger contains every vertical
prime, and `(10)` kills each of them individually.

## 2. The generic quartic and its two infinity sections

The discriminant of the quartic `H` in `X` is

```text
-256 t^12(16t^3+27),                                    (13)
```

so `C/K` is a smooth genus-one curve with the two rational points at infinity
distinguished by `W/X^2 -> +1` and `-1`.  Call them `O+` and `O-`, choosing
`O+` as the group origin.

The exact birational transformation starts with

```text
A=W+X^2-2X-t^3-2,                 V=XA,
x0=A/2,                           y0=(V-A-t^3-2)/2,       (14)
x0=t^2 x-1,                       y0=t^3 y.
```

It gives the global minimal Weierstrass model

```text
E/K: y^2=x^3+t x^2+1/4.                                 (15)
```

Let

```text
Q=(0,-1/2).
```

Direct addition on `(15)` gives

```text
2Q=(-t,1/2),
3Q=(t^-2,(t^3+2)/(2t^3)).                               (16)
```

The two finite points over `X=0` map as

```text
(X,W)=(0,+t^3) |-> Q,
(X,W)=(0,-t^3) |-> 2Q.                                  (17)
```

For the negative infinity branch, put `epsilon=1/X`.  Expanding the negative
square root in `(14)` gives

```text
A=(2t^3+4)epsilon+O(epsilon^2),
V=2t^3+4+O(epsilon),                                    (18)
```

which maps that branch to `3Q`; the positive branch maps to the zero section.
Therefore

```text
O--O+=3Q.                                                (19)
```

## 3. The Mordell--Weil group is exactly ZQ

For `(15)`,

```text
c4=16t^2,
c6=-8(8t^3+27),
Delta=-(16t^3+27).                                      (20)
```

There is good reduction at `t=0` and one `I1` fibre at each of the three
roots of `16t^3+27`.  At infinity put `s=1/t`,
`x_inf=s^2x`, and `y_inf=s^3y`.  The integral model is

```text
y_inf^2=x_inf^3+s x_inf^2+s^6/4,                        (21)
```

with valuations

```text
v_s(c4)=2,                 v_s(c6)=3,
v_s(Delta)=9.                                             (22)
```

It is minimal and has type `I3*`, hence root lattice `D7`.  The Euler invoice
is `3+9=12`, so the relatively minimal surface is a rational elliptic surface
with `rho=10`.  Shioda--Tate now yields

```text
rank E(K)=10-2-rank(D7)=1.                               (23)
```

There is no torsion.  Torsion injects into the `I3*` component group, which
has order four.  A rational two-torsion point would give a root in `k[t]` of

```text
x^3+t x^2+1/4.                                           (24)
```

Monicity makes that root polynomial; degree comparison leaves `x=-t+b`,
whose `t^2` coefficient forces `b=0` while its constant coefficient remains
`1/4`.  Thus `(24)` has no rational root.

Finally `Q` is primitive.  In the `s`-chart it is

```text
(x_inf,y_inf)=(0,-s^3/2).
```

Three successive point blowups in `(21)` have exceptional equations

```text
y_1^2=0,                y_2^2=0,
y_3^2=1/4,                                               (25)
```

and `Q` meets the terminal spin component `y_3=-1/2`.  The corresponding
diagonal entry of the inverse `D7` Cartan matrix is `7/4`, so Shioda's height
formula applies with `(Q.O+)=0`: the sections are disjoint at finite fibres,
and at infinity `Q` meets a nonidentity spin component.  Hence

```text
<Q,Q>=2-7/4=1/4.                                         (26)
```

The trivial lattice `U+D7` has absolute discriminant four, while the Neron--
Severi lattice of a rational elliptic surface is unimodular.  With torsion
zero, the Mordell--Weil determinant is therefore `1/4`.  Equation `(26)`
attains it, proving

```text
E(K)=ZQ.                                                  (27)
```

## 4. Class group, units, and the Cardano generator

The affine quartic is exactly

```text
C=E minus {O+,O-}.
```

Using `Pic(E)=Z[O+] + E(K)` and `(19),(27)`,

```text
Pic(C)=E(K)/<O--O+>
      =ZQ/<3Q>
      =Z/3.                                               (28)
```

Together with `(12)`, this proves the class-group assertion in `(2)`.  The
points in `(17)` show that `(3)` maps to `2Q` and `Q`, respectively, and the
principal identities in `(4)` follow from

```text
(q+W)(q-W)=4X^3.                                         (29)
```

Thus `[D+]=-[D-]` is the unique nonzero cyclic subgroup and has exact order
three.

A unit on `C` can have divisor only on `O+` and `O-`.  A nonconstant one would
make a nonzero multiple of `O--O+=3Q` principal, contradicting the infinite
order of `Q`.  Hence `C*=K*`.  Returning to the unit end of `(11)`, an element
of `K*` having valuation zero at every `V_lambda` is a scalar.  Therefore
`B*=k*`, completing `(2)`.

## 5. Cartier failure at the origin and the sole smooth-locus character

The distinction between `S` and its smooth locus is essential.  Let

```text
I=(X,q+W),                                                (30)
```

the divisorial ideal of `D+`.  The quotient `B/I=k[t]` is Cohen--Macaulay.
Since `B` is a two-dimensional Cohen--Macaulay hypersurface, the depth lemma
makes `I` maximal Cohen--Macaulay and hence reflexive.  At the singular origin
`O`, however, its two displayed generators have independent linear initials

```text
X, W in m_O/m_O^2.
```

Equivalently, their classes in `I_O/m_O I_O` are independent, so `I_O` needs
two generators and is not principal or invertible.  Thus the generator of
`Cl(B)=Z/3` is not Cartier.  Its inverse cannot be Cartier either, since that
would make the generator Cartier.  Therefore

```text
Pic(S)=0.                                                 (31)
```

Together with `B*=k*` and cube-divisibility of `k*`, the Kummer sequence on
`S` gives

```text
H^1_et(S,mu_3)=0.                                        (32)
```

Now let `U=Sreg`.  Normal Hartogs gives `Gamma(U,O_U)=B`, while restriction
across the codimension-two singular set gives `Pic(U)=Cl(B)`.  The Kummer
sequence on `U` therefore gives

```text
H^1_et(U,mu_3)=Cl(B)[3]=Z/3.                              (33)
```

The deck involution `W |-> -W` swaps `D+` and `D-`, so it acts by `-1` on
this group.  The generator is represented by `(q+W)/2`; its inverse is
represented by `(q-W)/2`.  Their product is `X^3`, and adjoining their cube
roots gives exactly `(6)-(7)`.

Now take a normal generically `S3` cubic completion having exactly this
squarefree discriminant and quadratic resolvent.  Over `U`, quadratic base
change kills the tame transpositions along the branch and leaves a connected
codimension-one-unramified `C3` cover.  Since `U` is regular, purity makes its
normalization finite etale.  By `(33)` its character is the Cardano character
or its inverse, and those two choices define the same cyclic extension.

For completeness, let `sigma` generate its `C3` group and let `g` lift the
quadratic deck involution.  The `S3` descent relation is
`g sigma g^-1=sigma^-1`.  Since `g^2` lies in `C3` and is fixed by conjugation
with `g`, inversion forces `g^2=1`.  The three lifts `sigma^b g` are conjugate:
conjugation by `sigma^c` replaces `g` by `sigma^(2c)g`, and two is invertible
modulo three.  Their fixed cubic fields are therefore `K`-isomorphic.  Thus
the Galois closure and cubic field are forced.  THM-3932 proves that the normal
integral closure in this field is the monogenic order `(7)`.  A different
nonprincipal resolvent ideal cannot create a second normal nonmonogenic cubic
order with this exact branch.

Thus the Cardano cover is quasi-etale over `S`: it is etale in codimension one
and a genuine torsor over `Sreg`, but it does not extend to an etale torsor
across `O`.  Equation `(32)` is the exact obstruction to confusing the local
Weil class with a Cartier class on the singular affine surface.

## 6. Local abundance versus global scarcity

At the origin the weighted initial equation is

```text
W^2+4X^3-t^6=0,                                          (34)
```

the simple elliptic `E~8` singularity.  Its weighted exceptional curve is the
`j=0` elliptic curve `W^2=t^6-4X^3` in `P(2,1,3)`, and the deck involution acts
as `[-1]`.  Locally its four lines in `E[3]` suggest four possible cyclic
directions.  Equations `(2)` and `(33)` show what the local picture forgets:
only the natural Cardano line globalizes to an affine Weil class.  The other
three local directions die in the global vertical/infinity ledger.

## 7. Exact companion and hostile boundaries

The assertion-free exact companion

```bash
python3 04-computation/jc2_linear_conic_resolvent_class_unique_character_thm3935.py
python3 -O 04-computation/jc2_linear_conic_resolvent_class_unique_character_thm3935.py
```

checks 53 identities: irreducibility and singular support; the uniform fibre
square obstruction; the quartic-to-Weierstrass change; the finite and infinity
fibre invariants; the two infinity addresses; the multiples of `Q`; the three
blowups and `D7` inverse-Cartan correction; absence of two-torsion; the two
independent local generators of `D+`; and the Cardano product/recovery
identities.  Both modes byte-match the frozen output.

The proof would fail, rather than silently persist, in any of the following
nearby settings:

1. a reducible closed fibre could contribute additional vertical primes to
   `(11)`;
2. deleting only one infinity section would retain an infinite degree class;
3. a torsion section or a second Mordell--Weil generator would enlarge
   `(28)`;
4. a nonnormal cubic suborder need not equal the unique integral closure;
5. changing branch multiplicities can make the quadratic base-changed cover
   ramified on `U`, outside `(33)`.

These are scope boundaries, not unresolved gaps in the stated exact model.
