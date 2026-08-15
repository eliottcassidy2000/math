---
id: THM-3397
title: "Torsor killing versus effective boundary valuations"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  A divisorial pole of a rational function on a
  normal affine variety survives every finite dominant normal pullback: each
  valuation is multiplied by a positive ramification index.  For the cyclic terminal initial ring
  R_e=k[v,X,Y]/(XY-(1+cv)^e), the missing rational coordinate
  w=X^(-1)v^(-m) has the exact coefficient-denominator filtration
  {h in k[v]: h w^q is polynomial}=
  (v^(mq)(1+cv)^(eq)).  Its v=0 pole when m>0 and its oriented X-boundary
  pole survive the polynomial cyclic cover even though that cover kills the
  generator of Cl(R_e)=Z/e and is etale-locally split on the regular locus.
  In the sharp V4 quotient R_0=k[a,b,c,d]/(d^2-abc), the three nonzero
  two-torsion classes are C3-compatible and all die on the polynomial
  V4 cover, but every nonzero constant combination of 1/a,1/b,1/c has a
  separated boundary pole; even their invariant sum remains nonpolynomial.
  Thus killed finite monodromy, class compatibility, and a UFD total space
  do not imply polynomial effectivity.  No A4/S4 or JC(2) branch is closed.
source: root/jc-effectivity-next-2026-08-14
audit: independent valuation, class-group, conductor, and finite-pullback line audit; exact normal/optimized replay; one-step boundary hostiles
depends_on:
  - THM-3383-terminal-monomial-cone-polynomiality-fork
  - THM-2655-quartic-keller-resolvent-v4-quasietale-torsor-and-kummer-class-group-gate
related:
  - THM-3072-a4-flag-three-c2-tomography-and-edge-cycle-cospan
  - THM-3386-linear-z-canonical-divergence-minimal-polynomial-collision-law
script: 04-computation/jc_torsor_effective_boundary_thm3397.py
output: 05-knowledge/results/jc_torsor_effective_boundary_thm3397.out
script_sha256: dbb402ba0fc1bfd9ae02f7997ffe5ee13b860cc255430da5e5f1a26de85e3179
output_sha256: 6befc8ec076622ff99eb4fe2848f9226b28c5798a38d49e9a08a7e9e750ab26d
semantic_sha256: a63ab86be7767bc5e160c13d4b0d5e228ef347dc27ee502849c75f2627d0f09f
hash_basis: LF-normalized bytes
---

# THM-3397 -- killing the cover does not make a rational section effective

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Inheritance and result

[THM-3383](THM-3383-terminal-monomial-cone-polynomiality-fork.md)
separates rational decoding from polynomial decoding on a complete terminal
monomial family.  Its two surviving coordinates are the valuations at
`v=0` and `1+cv=0`.  The canonical hostile is one decoded target with a
reciprocal in both coordinates.  The corrected near miss is localization:
removing those boundaries makes the coordinate regular but proves nothing
about its extension across them.

[THM-2655](THM-2655-quartic-keller-resolvent-v4-quasietale-torsor-and-kummer-class-group-gate.md)
supplies the closest topological mechanism and its sharp hostile.  A quartic
`V4` kernel returns as an etale torsor on the full resolvent normalization;
the even-sign quotient realizes all three compatible two-torsion classes on
a common polynomial cover.  The least-used sidecar is not another class: it
is the signed divisor of the rational function one wants to promote.

This theorem makes the distinction exact in two complementary models.
The terminal model computes the entire denominator filtration, while the
`V4` model shows that even simultaneous three-view compatibility and a
polynomial total space do not remove separated poles.

## 2. Universal finite-pullback pole persistence

Let `A subset B` be a finite inclusion of normal domains, with fraction
fields `K subset L`.  Let `P` be a height-one prime of `A`, and let `Q` be
a height-one prime of `B` above it.  The local rings at these primes are
discrete valuation rings, and the restricted valuation satisfies

```text
ord_Q(f)=e(Q/P) ord_P(f)       for every f in K^*,       (0)
```

where the ramification index `e(Q/P)` is a positive integer.  Therefore

```text
ord_P(f)<0  implies  ord_Q(f)<0.                         (0a)
```

Every prime downstairs has a prime above it because the morphism is finite
and dominant.  Thus a rational function with a divisorial pole cannot become
regular after any finite dominant normal pullback.  In ring language,

```text
B intersect K=A.                                        (0b)
```

Indeed an element of the intersection is integral over `A` and lies in `K`,
so normality of `A` puts it back in `A`.

This lemma concerns the pullback of the **same** rational function.  A new
function assembled upstairs can cancel principal parts, but that cancellation
is additional algebraic data; finite monodromy killing alone does not provide
it.  The next two models compute exactly what such a cancellation would have
to pay.

## 3. The cyclic terminal quotient

Let `k` be an algebraically closed field of characteristic zero, let
`e>=2`, and let `c in k^*`.  In `S=k[x,y]` put

```text
v=(xy-1)/c,             L=1+cv=xy,
X=x^e,                  Y=y^e.                          (1)
```

Let `mu_e` act by

```text
zeta.(x,y)=(zeta x,zeta^(-1)y).                         (2)
```

An invariant monomial has exponent difference divisible by `e`.  Removing
the smaller of its two exponents writes it as a monomial in `L,X,Y`.
Consequently

```text
R_e=S^(mu_e)
   =k[v,X,Y]/(XY-L^e).                                  (3)
```

The ring is normal, either as a finite-group invariant ring or directly as
the `A_(e-1)` hypersurface.  Its unique singular point is

```text
o=(L,X,Y)=(0,0,0).                                      (4)
```

On `U=Spec(R_e)\{o}`, the quotient map

```text
pi: A^2\{(0,0)} -> U                                    (5)
```

is a finite etale `mu_e` torsor: the action `(2)` is free away from the
origin.

## 4. Divisor class killing and what it forgets

There are two height-one primes over `L=0`:

```text
D_X=(X,L),                 D_Y=(Y,L).                   (6)
```

After inverting `L`, equation `(3)` eliminates `Y` and gives the Laurent
UFD `k[L^+-1,X^+-1]`.  Nagata's sequence therefore generates the class
group by `(6)`.  The principal-divisor relations are

```text
div(L)=D_X+D_Y,             div(X)=eD_X,                (7)
```

The units of `k[L^+-1,X^+-1]` are exactly `k^*L^aX^b`, so their
boundary divisors are generated by the two rows in (7).  Thus Nagata's
kernel contains no further relation, and

```text
Cl(R_e)=Z/e,                [D_X] a generator.          (8)
```

The pullbacks of the two boundary primes are

```text
pi^*D_X=(x=0),              pi^*D_Y=(y=0).              (9)
```

They are principal on the polynomial total space.  Thus the cover kills the
whole class-group generator.  It also tautologically kills its own torsor
monodromy after base change to the total space.  Neither fact controls the
signs of a chosen rational function's boundary orders.

The ordinary ring conductor does not replace that missing information.  In
fact

```text
(R_e:S)={r in R_e : rS subset R_e}=0.                   (9a)
```

If a nonzero `r` belonged to this conductor, then `rx` would lie in `R_e`,
so `x=(rx)/r` would lie in `Frac(R_e)`.  This contradicts the nontrivial
`mu_e` action, or equivalently the degree-`e` field extension
`Frac(S)/Frac(R_e)`.  Conductors are effective for birational orders inside
one fraction field; a nontrivial quotient cover is a different use case.

## 5. Exact two-boundary denominator filtration

Let

```text
H=(v) in Spec(R_e),          w_m=X^(-1)v^(-m),           (10)
```

where `m>=0`.  The divisor orders at the three displayed primes are

```text
                 H       D_X       D_Y
ord(w_m)        -m        -e         0.                 (11)
```

Indeed `X` is a unit at `H`; at `D_X`, `Y` is a unit and
`X=L^e/Y`, while at `D_Y`, both `X` and `v=-1/c` are units.

For every integer `q>=1`, define the coefficient-clearing ideal

```text
I_q(m)={h(v) in k[v] : h(v)w_m^q lies in S}.             (12)
```

Then the exact answer is

```text
I_q(m)=(v^(mq)L^(eq))=I_1(m)^q.                         (13)
```

For necessity, regularity along `H` forces `v^(mq)` to divide `h`.
Regularity along `D_X` forces a zero of order at least `eq` at
`v=-1/c`, hence `L^(eq)` divides `h`.  The factors `v` and `L` are coprime.
Conversely, if `h=v^(mq)L^(eq)k(v)`, relation `(3)` gives

```text
h w_m^q
 =k(v) X^(-q)L^(eq)
 =k(v)Y^q in S.                                         (14)
```

This also proves minimality independently at each boundary: deleting one
factor from either positive required exponent leaves a pole of order one.

Pulling `(11)` to the polynomial cover does not improve it:

```text
div(pi^*w_m) has orders (-m,-e,0)
on (xy=1), (x=0), (y=0).                                (15)
```

In particular `pi^*w_m` is not polynomial for any `m>=0`.  The torsor is
etale-locally split at the generic points of all three curves, but a local
section cannot change `(11)`: `w_m` is a rational function from the base,
so pulling it through a section returns the same function.  Local section
existence and boundary effectivity are logically independent gates.

### 5.1 Exact specialization to THM-3383

On THM-3383's rational-decoding locus, put `n=g-ae=+/-1`.  Its unique
missing decoded target is precisely

```text
n= 1:       t=X^(-1)v^(-ae)=w_(ae),
n=-1:       u=X^(-1)v^(-g) =w_g.                        (16)
```

Thus `(13)` upgrades the qualitative nonmembership there to all powers:

```text
{h in k[v]: h t^q in S}=(v^(aeq)L^(eq))       if n= 1,
{h in k[v]: h u^q in S}=(v^(gq)L^(eq))        if n=-1.  (17)
```

The filtration is the effective initial-module sidecar discarded by the
field decoder.  The class group sees only linear equivalence; a rational
function's divisor is already principal, so it cannot record whether the
principal divisor has a negative part.

## 6. A simultaneous three-view `V4` hostile

Now use THM-2655's sharp quotient.  Let

```text
T=k[x,y,z],
V4={even sign changes of (x,y,z)},
R_0=T^V4=k[a,b,c,d]/(d^2-abc),                          (18)

a=x^2,                   b=y^2,
c=z^2,                   d=xyz.
```

Let `C3` cycle `x,y,z`; together with `V4` it gives the standard
`A4=V4 semidirect C3` action.  The three boundary primes

```text
P_a=(a,d),               P_b=(b,d),
P_c=(c,d)                                                     (19)
```

satisfy

```text
Cl(R_0)=(Z/2)^2,
[P_a]+[P_b]+[P_c]=0,                                    (20)
```

and `C3` cycles its three nonzero classes.  This is the complete standard
two-dimensional `F_2[C3]` compatibility packet, not merely three unrelated
classes.  On the regular locus the polynomial quotient is a connected etale
`V4` torsor.  Its pullback sends the three classes to the principal
coordinate divisors `(x=0),(y=0),(z=0)`, so all three classes die
simultaneously on the UFD `T`.  For the same unequal-fraction-field reason
as `(9a)`, the cover conductor `(R_0:T)` is zero; it contains no hidden
boundary repair.

Nevertheless consider the three conjugate rational functions

```text
f_a=1/a,                  f_b=1/b,          f_c=1/c.     (21)
```

At the generic point of `P_a`, one has

```text
ord_(P_a)(a)=2,          ord_(P_a)(b)=ord_(P_a)(c)=0,   (22)
```

and cyclically.  Therefore, for constants `lambda_a,lambda_b,lambda_c`,

```text
lambda_a/a+lambda_b/b+lambda_c/c is regular in R_0
iff lambda_a=lambda_b=lambda_c=0.                       (23)
```

If `lambda_a!=0`, the first summand has the unique order `-2` term at
`P_a`; the other two have order zero and cannot cancel it.  The same argument
at `P_b,P_c` proves `(23)`.

In particular the fully `C3`- and `S3`-invariant orbit sum is

```text
1/a+1/b+1/c=(ab+ac+bc)/d^2,                             (24)
```

and has a pole of order two on every prime in `(19)`.  Pulling it to the
common polynomial cover gives

```text
x^(-2)+y^(-2)+z^(-2),                                  (25)
```

which has the same three separated poles.  Thus all of the following hold
at once:

```text
the three class views satisfy their A4 compatibility relation;
the common finite monodromy is killed on one polynomial UFD cover;
the invariant combination is still not polynomial.                     (26)
```

This is a divisor-side three-view hostile.  It is compatible with
[THM-3072](THM-3072-a4-flag-three-c2-tomography-and-edge-cycle-cospan.md),
but does not identify its flag quotient tables with these divisors.  The
common object here is the explicit `V4` polynomial cover in `(18)`.

## 7. Effectivity gate and transfer boundary

For a normal affine domain, a rational function is regular exactly when all
height-one orders are nonnegative.  Passing to `Cl` quotients the divisor
group by principal divisors and therefore discards precisely this sign
information.  Trivializing a torsion line bundle makes a divisor principal;
it does not make a specified rational trivialization regular or nowhere
vanishing.  The zero cover conductors `(9a)` and `(R_0:T)=0` likewise show
that the order-conductor mechanism applies only after a birational
same-field model has been supplied.  The missing coordinate here is

```text
effective divisor / principal-part data,
not finite monodromy or the divisor class alone.                         (27)
```

The recent Lu--Yang Banach preprint uses torsion killing, degree, polynomial
powers, and unique factorization inside its own structured norm problem.
Equations `(13)` and `(23)` show exactly why only the first of those moves
does not transfer to a rational inverse chart: a separate theorem must
force the boundary orders to be nonnegative.  This does not question the
preprint's claimed result and imports none of it into canon.

The comparison with THM-3386 is also typed.  Generic localization kills its
canonical divergence torsion but leaves a minimal-polynomial annihilator;
here a finite cover kills class/torsor data but leaves the Rees filtration
`I_q=I_1^q`.  Both are integral sidecars, but no map between the two modules
is asserted.

Nothing here realizes a physical quartic resolvent, proves that every
terminal module is monomial, imposes the Keller equation on `(18)`, or
excludes any `A4`, `S4`, `JC(2)`, or `DC(2)` branch.

## 8. Exact controls

The standard-library companion uses exact integer and rational arithmetic.
It checks `896` terminal parameter cells, both one-step boundary hostiles,
all divisor Smith forms for `2<=e<=9`, the `V4` invariant-monomial parity
criterion, the `(1,2,2)` divisor Smith form, all `124` nonzero coefficient
triples in `{-2,-1,0,1,2}^3`, and all `300` forced pole places.  These are
positive and hostile controls for the literal formulas; normality, etaleness,
and the all-parameter valuation arguments are proved above rather than
replaced by enumeration.

Reproduce with

```text
python3 04-computation/jc_torsor_effective_boundary_thm3397.py
python3 -O 04-computation/jc_torsor_effective_boundary_thm3397.py
```

The script uses no floating point, external package, shell, network, or
assertion-dependent truth gate.  Artifact and semantic hashes are pinned in
the frontmatter.

**QED.**
