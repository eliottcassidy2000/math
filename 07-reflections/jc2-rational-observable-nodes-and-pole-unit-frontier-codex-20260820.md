# Planar-JC overnight frontier: smooth selected observables, nodal full-field observables, and the first surviving cells

**Status:** synthesis of **PROVED**, **VERIFIED-EXACT**, **REFUTED**, and
explicitly **OPEN** cells.  This note proposes counterexample architectures
but does not assert that any exists.  Current truth lives in the cited canon
files, not in this reflection.

## 1. Inheritance pass

The anchor was the rational-completion / nonlinear-descent frontier after
THM-3553, THM-3558, THM-3559, and THM-3560.  Concurrent work added four
load-bearing gates during the session: the smooth Danielewski completion
THM-3561, balanced pole-unit nonentry THM-3562, the fixed nongraph node-cycle
obstruction THM-3563, and target-graph resonance/classification
THM-3564/3565.

- **Closest positive mechanism:** THM-3561 gives an everywhere-etale
  polynomial map from `A2` to the smooth non-`A2` surface
  `c^2e=b(b+4)`, of generic degree three, with an exact triple collision.
  A polynomial Darboux pair on that surface would be a planar
  counterexample.  THM-3569 has now extended the homogeneous and two-by-two
  gates through the complete two-by-three cell.  Any solution needs at
  least six nonconstant pieces; `3x3` and `2x4` are the first live shapes.
- **Canonical hostile:** THM-3554 retains an etale collision only on
  `G_m x A1`; filling the divisor restores Kummer ramification.  THM-3567
  gives the sharper all-degree separated hostile: the *full target-field*
  observable ring is nodal, and degree at least three contracts companion
  lines.
- **Corrected near miss:** a smooth ring can arise when one intersects the
  polynomial target ring with source polynomials, as in THM-3561.  It does
  not follow that the larger field intersection is smooth.  Conflating
  `C[a,c] intersect C[x,q]` with
  `C(a,c) intersect C[x,q]` was the dangerous missing type.
- **Least-used sidecars:** on rational completions, record both the selected
  target-ring intersection and the full target-field intersection.  On
  Faber responses, retain the common pure-`q` pole value and the marked
  divisor at infinity.  On target graphs, retain the infinity valuation of
  a possible core-cubic root before expanding coefficients.

The live concept board is:

| object | representation | retained predicate | destroyed information | cheapest decisive test |
|---|---|---|---|---|
| smooth Danielewski completion | Poisson ring `C[b,c,e]/(c^2e-b(b+4))` | polynomiality, etaleness, triple collision | affine-plane target | complete `3x3` and `2x4` support cells |
| separated rational seed | `P=f(x), Q=y/f'(x)` | unit rational Jacobian, arbitrary collision degree | smoothness and quasi-finiteness | full-field intersection and companion fibres |
| nonlinear target graph | `c+phi(a,b)=0` | coordinate target and collision value | reducible source pullback | first resonance `deg_a(phi)=4` |
| balanced Faber response | pole divisor and pure-`q` unit | normalized nonsplit chart | unbalanced infinity owner | projective pole-unit identity |
| degree passport | `H^2,H^3` at `(72,108)` | global necessary degree shape | construction mechanism | reject candidates before large coefficient search |

## 2. Exact progress

### 2.1 THM-3561: a smooth selected-observable completion

For

```text
D=1+x^2q,
a=q/D^2,
c=xD(D+2),
```

THM-3561 proves `Jac(a,c)=-3` and an exact rational triple collision.  The
functions in the **polynomial target ring** whose pullbacks extend across
`D=0` form

```text
C[a,c] intersect C[x,q]
 =C[b,c,e]/(c^2e-b(b+4)),
b=ac^2=(D-1)(D+2)^2,
e=a(b+4)=q(D+3).                                      (1)
```

This surface is smooth.  The polynomial completion `Phi:A2->Y` is
everywhere etale of generic degree three, keeps the collision, and misses
exactly `(-4,0,0)`.  It is not a planar counterexample because
`chi(Y)=2` and `Pic(Y)=Z`, so `Y` is not `A2`.

The residual problem is exact rather than philosophical: find

```text
P,Q in C[Y]                    with {P,Q}=kappa!=0.     (2)
```

The symplectic form is globally exact, so de Rham theory does not kill
`(2)`.  THM-3569 proves that every solution needs at least six nonconstant
weight pieces and leaves exactly the first shapes
`(3,3),(2,>=4),(>=4,2)`.  Thus the next bounded cells are `3x3` and `2x4`,
not an unbounded blind coefficient expansion.

### 2.2 THM-3567: the separated full-field completion is exactly nodal

For `f'` squarefree with distinct critical values, put

```text
P=f(x),                     Q=y/f'(x),
Delta(T)=product_(f'(alpha)=0)(T-f(alpha)).             (3)
```

Then `Jac(P,Q)=1`, and each regular value gives a `deg(f)`-fold collision at
`Q=0`.  If `K=k(P,Q)`, THM-3567 proves the full-field equality

```text
K intersect k[x,y]
 =k[P,A,B]/(A^2-Delta(P)B),
A=Delta(P)Q,                 B=Delta(P)Q^2.             (4)
```

The mechanism is the sharp local semigroup

```text
Delta(P)^ceil(n/2) divides the coefficient of Q^n.      (5)
```

Each critical value produces one ordinary node.  For `deg(f)>=3`, every
critical fibre also has a regular companion root and its whole source line
contracts to that node, so the completion is not quasi-finite.  A quadratic
hostile shows the exact boundary: for
`(x,y)->(x^2,2xy,y^2)`, points over `Delta=0` away from the node can still be
etale.  Thus `Delta!=0` is a sufficient etale open, not always the maximal
one.

The rational cubic control

```text
f=x(x-3)(x-8),                 f'=(x-6)(3x-4)           (6)
```

has collision `(0,0),(3,0),(8,0)->(0,0)`, critical values
`-36,400/27`, full-field surface

```text
A^2=(P+36)(27P-400)B,                                  (7)
```

and collapsed companion lines `x=-1,25/3`.  Normal, optimized, and stored
exact replays agree.

There is no conflict between `(1)` and `(4)`: `(1)` is a polynomial-ring
intersection for a nonseparated seed, while `(4)` is the full-field
intersection for a separated seed.

### 2.3 THM-3562: every balanced resonant pole passport is dead

The independent audit reconfirms the promoted theorem.  The pure-`q` face
forces one common nonzero value

```text
E(beta_j)=e0!=0                         at every pole.   (8)
```

The balanced first integral then gives

```text
p_j W'(beta_j)=q!=0                     for every j.    (9)
```

For at least two poles, `(9)` contradicts the Lagrange identity
`sum_j 1/W'(beta_j)=0`.  With one pole, balance makes its order even and the
response split.  This closes all balanced passports in the normalized
nonsplit polynomial exact-square-prefix chart.  It does not close
unbalanced infinity divisors, nonpolynomial prefixes, or other charts.

## 3. Failure anatomy

The separated failure is not caused by a poor generator choice.  The first
false implication is

```text
rational Keller collision
  -> the full field of polynomial pullbacks has a smooth completion.   (10)
```

The local order pair `(pole order 1, critical-value contact 2)` forces
`(5)`, and `(5)` is exactly the nodal suspension `(4)`.  Increasing the
degree adds nodes.  Source or target polynomial automorphisms cannot change
the intersection singularity type, and a rational symplectic conjugation
regular along every critical divisor preserves the same parity semigroup.

THM-3563 closes a different tempting shortcut.  On its fixed nongraph
normalization plane, all target functions become even on a Laurent double
curve, while polynomial exactness produces an odd inverse-cube difference.
The three-node cycle forces incompatible values.  Therefore a cosmetic
node cycle is not a surviving architecture; a new normalization must break
the even double-curve packet or carry a genuinely new sidecar.

For target graphs, THM-3564 says reducibility is possible only when
`deg_a(phi)=1 mod 3`.  THM-3565 classifies the complete degree-one resonance
and THM-3568 proves neither component in that family is `A2`; genuinely
quadratic graphs are irreducible.  THM-3570 then replaces ad hoc factor
searches by a universal rational Pell-conic compiler.  The first polynomial
degree not already excluded is still the quartic resonance, but its search
variable is now the compiler parameter `q`, not an arbitrary cubic root.

## 4. Counterexample architectures that still survive

These are **OPEN proposals**, ordered by the cost of a decisive hostile.

### A. Passport-aware `3x3` and `2x4` Darboux search on `Y`

Use THM-3561's weight decomposition and begin exactly where THM-3569 stops.
Give the outputs either three weights each or two and four weights, subtract
constants, and use

```text
{c^r f(b),c^s g(b)}
 =c^(r+s+1)(s f'(b)g(b)-r f(b)g'(b)).                 (11)
```

Rather than bounding arbitrary coefficients, enumerate the integer support
patterns whose pairwise sums can cancel every nonzero output weight and
leave weight zero.  Impose the required powers of `S=b(b+4)` before solving
the Wronskian equations.  Any solution gives a literal planar Keller
counterexample after pullback through `Phi`; a contradiction closes an
entire Newton cell.

The cheapest bank is support weights in `[-6,6]` with minimal allowed
`b`-degrees.  Every row must be replayed symbolically and then measured in
ordinary source degree.  This route survives all current Danielewski gates
because `3x3` and `2x4` are their proved equality boundary.

### B. Pell-compiled quartic target graph

THM-3570 proves that a nonzero target graph is reducible exactly when some
`q in C(a,b)^*` satisfies

```text
phi q^3-4q^2-4bq-16a=0,
phi=4(q^2+bq+4a)/q^3.                                 (12)
```

Search rational numerator/denominator pairs for `q` for which `(12)` is a
polynomial of `a`-degree four and is genuinely mixed in `a,b`.  Polynomiality
is an exact divisibility problem and should be imposed before expanding the
source pullback.  A positive compiler row must then pass the Jelonek Euler
gate, have a component carrying two collision points, and make that component
an actual coordinate `A2`.  The order is compiler divisibility, Euler
section, component topology, coordinate recognition, then planar Jacobian.
This is the first all-degree graph search that uses every current factor gate.

### C. Two-variable critical-divisor-changing rational seed

The separated semigroup `(5)` persists under conjugations regular at its
critical divisors.  A surviving seed must change those valuations.  Start
with a denominator such as

```text
D=1+x^2q+x q^2                                           (14)
```

or a symplectic birational transform with a prescribed zero/pole on the
critical divisor, and solve the rational constant-Jacobian equation before
asking for collisions.  For every candidate compute two rings separately:

```text
k[a,c] intersect k[x,q],
k(a,c) intersect k[x,q].                               (15)
```

The first may be smooth while the second is singular.  The candidate
survives only if a smooth selected subring supports a constant Poisson pair,
or if the full-field ring is regular and regular companion components go to
infinity rather than contract.  Formula `(14)` is a design seed, not a
Keller claim; its cheapest test is the exact Jacobian PDE and divisor
semigroup.

### D. Projective pole-unit defect for unbalanced responses

Homogenize the response divisor and keep infinity as a marked pole/zero.
Re-derive `(8)` as a projective value law.  If finite poles still share one
value with zero infinity correction, THM-3562's Lagrange contradiction
extends immediately.  If a correction survives, its exact residue is the
only escape coordinate and should be solved before coefficients are
enumerated.

The first test is the existing asymmetric `R=8` hostile, followed by the
first response bank whose degree map can reach height `108`.  The state is
`(finite pole passport, infinity order, common-value defect)`, not a raw
coefficient tuple.

### E. Broken-involution nongraph normalization

THM-3563 kills every target projection on one sharp constant-different
plane because its Laurent double curve identifies `h` with `-h`.  Move to a
nearby nonlinear normalization on which the target packet distinguishes
those two branches, for example by retaining one odd conductor observable
that is not polynomial on the closed slice but becomes polynomial after an
affine modification.  The required sidecar is the sign/holonomy lost in the
old target ring.

The cheapest hostile is not a longer node cycle.  It is a two-branch test:
compute the normalization, conductor, units, class group, and whether every
candidate target observable remains even.  If evenness persists, the route
falls back into THM-3563 immediately.

## 5. The `(72,108)` compiler

A collision is not a degree passport.  Up to transpose, the cited current
sub-`125` classification leaves only

```text
(deg F,deg G)=(72,108),            gcd=36, exponents=(2,3).  (16)
```

The leading and first subleading rows obey

```text
F_72=c H^2,                   G_108=d H^3,             (17)
deg(H)=36,
rad(H) | F_71,                H rad(H) | G_107.        (18)
```

For squarefree `H`, write `F_71=H C_35`; then `G_107` is divisible by
`H^2`.  Any rational-completion or target-graph proposal must compile its
conductor data into this one degree-36 object and the shared subleading
quotient.  Padding a low-degree collision does not do so.

A concrete cross-lane hostile is therefore:

1. take a two-factor squarefree model `H=L1 L2` as the local conductor
   shadow;
2. impose `(17)--(18)` exactly;
3. encode the local choices as a `3x3` or `2x4` Danielewski support, or as a
   Pell-compiled quartic graph factor;
4. solve every homogeneous constant-Jacobian row; and
5. reject the architecture if its collision disappears under the required
   coordinate recognition.

Only after this two-factor hostile survives should the object be scaled to
a labelled 36-factor conductor packet.  The labels must retain branch sign,
source component, and infinity owner; an unlabelled dual graph loses exactly
the information used by THM-3562/3563.

## 6. Exact status and next tests

```text
PROVED / VERIFIED-EXACT:
  THM-3561 smooth selected polynomial-target-ring completion;
  THM-3562 balanced resonant pole-unit nonentry;
  THM-3563 obstruction on its fixed nongraph normalization plane;
  THM-3564/3565 target-graph resonance and complete degree-one factors;
  THM-3568-reducible-target-graph-component-euler-no-go closes both such components;
  THM-3569 closes the complete two-by-three Darboux cell;
  THM-3570 gives the universal Pell-conic graph factor compiler;
  THM-3567 separated full-field nodal completion and line contractions.

REFUTED:
  smoothness of the separated full-field observable surface;
  global etaleness of that full-field completion;
  raw node-cycle projection on the THM-3563 slice;
  nonresonant target-graph degrees and coordinate factors in degree one;
  Danielewski support shapes through two-by-three;
  treating balanced nonentry as only a finite-pole equality check.

OPEN:
  Danielewski `3x3` and `2x4` Darboux cells;
  polynomial quartic rows of the Pell compiler and component recognition;
  two-variable critical-divisor-changing rational seeds;
  unbalanced projective pole-unit defects;
  broken-involution nongraph normalizations;
  any construction entering (72,108), and JC(2).
```

The cheapest next exact computations are the Danielewski `3x3`/`2x4`
support-pattern classifier and the polynomiality locus of the Pell parameter
`q` in quartic graph rows.  These attack the first equality boundaries of
two independent proved no-go mechanisms and can return either a literal
counterexample candidate or a whole-cell obstruction.
