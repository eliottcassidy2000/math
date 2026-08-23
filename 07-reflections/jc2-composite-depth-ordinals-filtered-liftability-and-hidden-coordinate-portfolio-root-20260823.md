# Composite depth, ordinal schedulers, and filtered liftability in the planar Jacobian frontier

**Session:** root / 2026-08-23  
**Portfolio:** Anchor = sextic normal strips; Niche = alternative Russell
seed quotients; Wildcard = operation-compatible natural addresses and
cross-problem filtration hostiles.

**Truth boundary.**  `JC(2)`, `LRC(14)`, and Mahler's `3/2` problem remain
**OPEN**.  THM-3412, THM-3745, THM-3756, THM-3801, THM-3811, THM-3848,
THM-3868, THM-3871, THM-3872, and THM-3874 are used only at their current
proved scope.  The new computational companions in this session are
**FINITE-EXACT SCOUTS**, not theorem dependencies; the cubic and sextic main
packets additionally have independent finite-exact audits.  During the
session THM-3877, THM-3879, THM-3880, THM-3881, and THM-3882 were promoted to
**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED**.  They enter below
only at their typed scopes: THM-3881 still leaves its general arbitrary-degree
norm equation open, and THM-3882 closes a rational-dual architecture rather
than `JC(2)`.

## 1. Outcome first

The session produced four useful advances without pretending to close the
planar Jacobian conjecture.

1. The degree-six top row has an exact composite arithmetic sieve.  After
   top-direction normalization, `(6,1)`, `(6,2)`, and `(6,3)` are target
   shears into THM-3871.  Only `(6,4)` and `(6,5)` are genuinely sextic.
2. On the principal balanced pole face, `(6,5)` has fourteen apparent
   branches after arms and two conserved equations; its next nonexact
   Jacobian bucket kills all fourteen.  The `(6,4)` face instead collapses
   to a nonreduced length-six cusp whose leading constant bracket vanishes.
   It needs a strict transform, so there is no sextic theorem yet.
3. The proposed “first genuine cubic seed residual” was a section artifact.
   The full formal source-addition action has a polynomial inverse at every
   grade.  Constant density survives only when the source determinant
   cocycle is retained.  The raw affine cubic cell is nevertheless empty,
   and its unique two-row near-solution is the cubic truncation of the old
   nonsquare Catalan/Hensel resummation.
4. The Pythagorean and Jacobian task universes share the *same* triangular
   natural-number address, but their arithmetic filters disagree at the
   first Pythagorean hole.  This gives an exact positive use of ordinals as a
   scheduler and an exact negative against treating the scheduler as an
   invariant.

The root ordinal/GCD companion is
[`jc2_normal_strip_ordinal_gcd_sidecar_probe_root_20260823.py`](../04-computation/jc2_normal_strip_ordinal_gcd_sidecar_probe_root_20260823.py).
It performs more than twenty-five million exact valuation-ray checks in both
ordinary and optimized mode.

## 2. Inheritance pass and live board

The closest proved mechanism was THM-3871's quintic normal-strip closure:
normalize the constant top target direction, factor the top Wronskian in
`k[s]`, depress, and keep every lower conserved bucket.  Its canonical
hostile is the arms-only point

```text
(X,Y,Z)=(1,3/10,-23/10),
```

which cancels both arms while the two conserved residuals remain `476/25`
and `321/50`.  The corrected near miss was “recovery failure implies Keller
escape”: the alternative quadratic seed has genuine nonintegrality but a
nonconstant boundary Jacobian.  The least-used relevant sidecars were the
full source-addition ideal, rather than one normal-form section, and the
filtered transition equations below the top bucket.

| Live object | Representation | Predicate | Destroyed coordinate | Cheapest hostile |
|---|---|---|---|---|
| normal strip | twelve Jacobian buckets | constant nonzero bracket | lower extension equations after top scalarization | THM-3871 arms-only point |
| sextic top cell | `(n,j)` plus Kummer ray | target-shear reducibility | `gcd(n,j)` under ordinal rank | consecutive ranks 13 and 14 |
| alternative seed | raw coefficients modulo source additions | constant density and effectivity | determinant cocycle | `(x,u)->(x+x^3,u)` |
| hidden-control recovery | monogenic eliminant or rank-two order | integrality/nonmonogenicity | index/different and deleted divisors | THM-3868 versus THM-3801 |
| triangular address | `T(n-2)+j` | enumeration | native filter and operation | shared rank 8 |
| filtered response | `gr(ker)->ker(gr)` | liftability | lower recurrence/transition maps | THM-3412 block `[4]` versus `[3,1]` |

Every successful pull changed at least two board entries: the gcd sieve
shrunk the sextic target; the full quotient erased the supposed cubic
invariant; and the filtered-kernel hostile explained why conserved buckets
must remain attached to every leading arm solution.

## 3. A triangular edge address for normal-strip tasks

For every potential new top cell

```text
n>=2,                         1<=j<n,
```

define the one-based task address

```text
rho(n,j)=T(n-2)+j.                                           (1)
```

The cells in shell `n` occupy the consecutive interval

```text
[T(n-2)+1,T(n-1)].                                         (2)
```

Thus `(1)` bijects the full triangular task cone with the positive natural
numbers.  Through depth `N`, the cells are literally the edges

```text
(n,j) <-> {j,n} in K_N,                                    (3)
```

so their count is `T(N-1)=binom(N,2)`.  Adding depth `N` is adding one vertex
and its `N-1` incident edges.  This is also why the user's signed identity

```text
T(z+2)-T(z-2)=4z+2                                         (4)
```

measures a four-shell edge increment.  The `r`th odd square can be stored as
the ordinal `r` because

```text
(2r-1)^2=8T(r-1)+1.                                        (5)
```

Equations `(1)--(5)` are exact address laws.  They do not yet say which
tasks are equivalent or which mathematical predicate an edge carries.

There is a second literal count connection.  THM-3745 proves that the
normalization defect of an ordinary `N`-fold monomial-plane branch is

```text
delta=T(N-1)=binom(N,2).                                    (6)
```

Under `(3)`, both sides enumerate unordered pairs.  The map preserves the
complete-graph incidence carrier and nothing else: it does not send a
normal-strip bracket equation to a conductor equation.  It is a lawful
shared scheduler, not a JC reduction.

At `N=14`, the same address set has `91` elements, just as the pairwise
runner graph in `LRC(14)` has `91` edges.  The LRC edge carries a speed ratio,
phase, owner, and covering predicate; the JC edge carries transverse
exponents and coefficient buckets.  No predicate-preserving map between
them is known.

## 4. The top Kummer theorem and the missing gcd sidecar

Let

```text
A=w(s)z^n+lower,                  C=c(s)z^j+lower,
0<j<n.                                                     (7)
```

The highest coefficient of the Jacobian is

```text
n w c'-j w'c.                                             (8)
```

If `(8)=0`, put `g=gcd(n,j)`.  At every irreducible prime `p` of `k[s]`,
logarithmic residues give

```text
n ord_p(c)=j ord_p(w).                                    (9)
```

Since `n/g` and `j/g` are coprime, there are `R,Q in k*` and `h in k[s]`
such that

```text
w=R h^(n/g),                    c=Q h^(j/g).              (10)
```

Conversely `(10)` satisfies `(8)`.  The equality boundary includes constant
`h`; zero top coefficients belong to the preceding degree cell.

A constant target shear removes the `A`-top exactly when `j|n`.  Indeed this
is equivalent to `j/g=1`, in which case

```text
w/c^(n/j)=R/Q^(n/j) in k*.                                (11)
```

Thus the edge `(n,j)` is colored by its primitive Kummer ray
`(n/g,j/g)` and by the divisibility predicate `j|n`.  The exact number of
genuinely new rows in shell `n` is

```text
n-tau(n),                                                 (12)
```

where `tau(n)` is the number of positive divisors of `n`.  Prime depths have
many new rows; composite depths can have unexpectedly few.

For `n=6`, the ordinal/gcd table is

| rank | cell | `gcd(6,j)` | primitive ray | operation |
|---:|---:|---:|---:|---|
| 11 | `(6,1)` | 1 | `(6,1)` | subtract a scalar `C^6` |
| 12 | `(6,2)` | 2 | `(3,1)` | subtract a scalar `C^3` |
| 13 | `(6,3)` | 3 | `(2,1)` | subtract a scalar `C^2` |
| 14 | `(6,4)` | 2 | `(3,2)` | genuine sextic row |
| 15 | `(6,5)` | 1 | `(6,5)` | genuine sextic row |

The first three land in depth at most five and are closed by THM-3871.
Ranks 13 and 14 are consecutive, but one is a target shear and the other is
new.  Ordinal successor therefore does not preserve the mathematical task
type; `gcd(n,j)` is load-bearing.

## 5. The exact Pythagorean cross-hostile at rank eight

THM-3756 uses the ambient Pythagorean address

```text
a(r,s)=T(r-2)+s,                      r>s>=1.             (13)
```

This is literally `(1)` under the identity `(r,s)=(n,j)`.  The map preserves
shell order, complete triangular packing, and round-trip decoding.

It does **not** preserve the native arithmetic filter.  Pythagorean
primitivity asks for

```text
gcd(2r-1,2s-1)=1,                                         (14)
```

whereas the Jacobian top packet uses ordinary `gcd(n,j)` and the shear test
`j|n`.  The first failure is exact:

```text
rank 8 <-> (5,2).
```

On the Pythagorean side, `(14)` is `gcd(9,3)=3`, so this is the first hole;
its ambient decoded triple is `(27,36,45)=9(3,4,5)`.  On the JC side,
`2` does not divide `5`, so `(5,2)` is a genuine quintic row.  THM-3871
closes it by the unrelated local arm residue `3R/8`.

This is the desired exact interpretation of “treat the selected odd square
as `n`.”  The ordinal can schedule both objects.  To act on either object,
one must restore its native filter and operation word.  Shared syntax alone
transfers no theorem.

## 6. Sextic balanced faces: one closure and one strict-transform cusp

The sextic companion and detailed report are:

- [`jc2_sextic_normal_strip_composite_frontier_scout_20260823.py`](../04-computation/jc2_sextic_normal_strip_composite_frontier_scout_20260823.py);
- [`jc2-sextic-normal-strip-composite-frontier-scout-20260823.md`](../05-knowledge/results/jc2-sextic-normal-strip-composite-frontier-scout-20260823.md); and
- [`jc2-sextic-normal-strip-composite-frontier-independent-audit-20260823.md`](../05-knowledge/results/jc2-sextic-normal-strip-composite-frontier-independent-audit-20260823.md),
  a 244-gate independent replay by radial homotopy, Buchberger reduction, and
  quotient multiplication matrices.  It passes with certificate
  strengthening and no formula correction.

They verify all twelve universal buckets and complete depressed rational
packets in the two surviving rows.

### `(6,4)`

The top packet is `w=Rh^3,c_4=Qh^2`.  After a quadratic scale extension and
depression, six high rows integrate, two more become conserved polynomials,
and the constant row remains.  On the principal balanced pole face, the arm
and conserved equations reduce at their unique support to

```text
6v^2-u^3,
u^2(2u+3v),
u^4.                                                       (15)
```

The local algebra has basis

```text
1,u,u^2,u^3,v,uv                                           (16)
```

and length six.  Its reduced support is the common-power cusp

```text
A=T^3,                         C=T^2,                       (17)
```

whose bracket is zero.  The normalized leading constant bracket belongs to
the ideal `(15)`, so it cannot pay a nonzero constant on this face.  Because
the support is nonreduced, this is not a full contradiction: dividing the
exceptional multiplicity may expose a first surviving jet.  The quadratic
scale extension also requires an involution/parity descent audit.

The first strict-transform dividend is already exact.  On the primitive
simple-pole chart, successive arm and conserved coefficients force

```text
D=F=I=0,                         u=O(t^2), v=O(t^4).        (17a)
```

Thus three odd sextic constants die before the next cusp-weighted chart.
The surviving variables are the `u_2,v_4` jets, the regular target-arm term,
the even constants `E,G,A0`, and both conserved values.  Higher pole order
can insert intermediate uniformizer jets, so `(17a)` still is not an
all-valuation closure.

### `(6,5)`

Here the top packet is `w=Rh^6,c_5=Qh^5`, so no quadratic extension is
needed.  Arms plus two conserved rows give a reduced apparent scheme of
length fourteen.  The next `x^1` row is a genuinely nonexact one-form and
satisfies

```text
(arms,conserved rows,source arm,x^1 row)=(1).             (18)
```

It kills all fourteen apparent branches.  The constant leading bracket is
already nonzero at every one.  This closes only the **principal balanced
face**.  Nonbalanced valuation fans and zero, unit, identically-zero, and
constant-scale channels remain open in both rows.

## 7. The cubic residual disappears under the correct quotient

The cubic companion and full loss ledger are:

- [`jc2_russell_cubic_seed_quotient_constant_density_scout_20260823.py`](../04-computation/jc2_russell_cubic_seed_quotient_constant_density_scout_20260823.py);
- [`jc2-russell-cubic-seed-quotient-density-cocycle-scout-20260823.md`](jc2-russell-cubic-seed-quotient-density-cocycle-scout-20260823.md); and
- [`jc2-russell-cubic-seed-constant-density-independent-audit-20260823.md`](jc2-russell-cubic-seed-constant-density-independent-audit-20260823.md),
  which independently passes 66 exact gates without repair.

For the fixed normalized seed

```text
A_0=u^2-x,
C_0=u^3-u-(3/2)ux,                                         (19)
```

the associated-graded source-addition action at every normal grade is

```text
M(u)=[[-1,       2u],
      [-(3/2)u, 3u^2-1]],                    det M=1.       (20)
```

Its inverse is polynomial.  Hence every cubic coefficient residual relative
to a chosen quadratic section is killed by a cubic source addition.  There
is no cubic residual invariant on the full cosmetic quotient.

The quotient is too coarse for the Keller predicate: a source addition
changes the source determinant.  For example `(x,u)->(x+x^3,u)` preserves
the arm and first normal row but multiplies the seed density by `1+3x^2`.
The lawful quotient object is therefore

```text
(formal source orbit, determinant cocycle, complete divisor data), (21)
```

or a proved determinant-one effective subgroup.

In raw affine coefficient space, the first two nonconstant Jacobian rows
force

```text
p=3/4, q=9u/8, r=-9/8, s=-27u/16.                         (22)
```

This is the fixed seed precomposed by

```text
chi_3=x-(3/4)x^2+(9/8)x^3.                                (23)
```

The next density bucket is `135/16`, so no constant-J seed occurs.  The
coefficient is structural: `(23)` is the cubic truncation of

```text
(1+(3/2)chi)chi'=1,
chi=(2/3)(sqrt(1+3x)-1).                                  (24)
```

The first omitted term is `-(135/64)x^4` and cancels the failed bucket.
The exact sheet is algebraic and nonsquare in the rational function field,
reconnecting the calculation to THM-3846 rather than escaping it.

Restoring the entire polynomial quadratic fibre and affine relative cubic
residuals still ends at the same contradiction.  On the known nonintegral
quadratic pole chart, every nonzero affine cubic correction creates a pole
and cannot change the already nonconstant first Jacobian bucket.

## 8. The direct filtered-kernel connection

The independent cross-frontier replay and its typed connection ledger are in
[`jc2-sextic-filtered-kernel-principal-parts-cross-frontier-scout-20260823.md`](jc2-sextic-filtered-kernel-principal-parts-cross-frontier-scout-20260823.md).

THM-3412 supplies more than an analogy.  Its filtered Hamiltonian complex has
a natural symbol map

```text
gr_q^F(ker Dbar) -> ker(gr_q Dbar).                         (25)
```

For `f=x,g=x^3,q=2`, the finite kernel is one `P`-Jordan block `[4]`.
The graded kernel has blocks `[3,1]`; both have total dimension four, but the
liftable-symbol space has dimension two.  The deficit is two.  A dimension
count or top-symbol solve therefore cannot detect liftability.

For a depressed JC packet, the corresponding first-order operator is

```text
L_(A,C)(delta A,delta C)
 ={delta A,C}+{A,delta C}.                                  (26)
```

Filter `(26)` by normal depth and, independently, by every boundary-prime
pole order.  Passing to the top observer preserves the Wronskian/Kummer row
and leading arm cancellations.  It destroys the lower coefficient equations
needed to lift that symbol to a full Keller deformation.  THM-3871's
arms-only point is the exact target hostile: the arm system is a genuine
smooth local one-fold, while the omitted conserved rows are nonzero.

The required sextic sidecar is therefore the filtered/Rees kernel, transition
maps, and lower conserved equations—not merely the top scheme.

THM-3404's CRT principal-part tower is a direct model only after a literal
Danielewski suspension `XY=F(v)` or a specified DVR jet module is identified.
For the present sextic strip it is still a **proposed sidecar**: no boundary
chart and no Hamiltonian intertwiner have yet been proved.  This distinction
prevents an attractive principal-part analogy from becoming a false theorem
transfer.

There is a second exact geometry hostile for THM-3862's common-infinity
contract.  The graph curves `y=x^d`, `d=2,3,4`, all normalize to `A1` and
meet the projective line at the same point with the same tangent ray.  Their
pairwise contacts at infinity are nevertheless different:

```text
I(Gamma_2,Gamma_3)=3,                 I(Gamma_2,Gamma_4)=4. (26a)
```

Thus “all branches share one infinity address” is a support statement, not
a contact classifier.  A future completion passport must retain branchwise
Puiseux/valuation arms and their intersection matrix.

## 9. Other problems as signal, with exact stopping points

### Mahler `3/2`

THM-3848 proves that every finite strict safe prefix is compatible with the
weak closed inverse-limit shift, while a countable equality boundary is
invisible to all finite prefixes.  This does not map Mahler words to JC
coefficients.  It supplies a procedural test: whenever every finite jet is
solvable, compute the inverse-limit boundary and the native strict/effective
cocycle.  In the cubic seed that cocycle is the determinant and divisor; in
Mahler it is the unbounded suffix match plus terminal residue clock.

### The degree-six Young-subgroup gap

THM-3112 first needs cycle-weight refinement at degree six: support size
collapses the cycle types `(4,1,1)` and `(2,2,1,1)` even though the target
operator distinguishes them.  There is no map from those symmetric-group
operators to sextic coefficient buckets.  The useful common warning is typed:
degree alone is a quotient, and the first new operation may require an
orbit/factor sidecar.  The shared number six is not evidence of a theorem.

### Nonmonogenic cubic completions

THM-3868 closes a fixed seed because one monic hidden control is recovered.
THM-3801 proves that a constant-unit cubic finite completion must instead be
nonmonogenic, and THM-3811 exhibits the relevant rank-two binary-cubic/index
grammar on a different surface.  The cubic scout strengthens the task choice:
adding another coefficient to one hidden primitive merely follows the same
formal source orbit.  The next recovery experiment should use a trace-zero
rank-two module with its binary index form, different, ramification classes,
and visible/deleted companion sheets.  It must not assume the THM-3811
surface is the Russell surface.

### Incoming cusp-ideal work

Promoted THM-3872 validates the operation “quotient cosmetic additions by
the entire vanishing ideal before selecting a representative.”  The cubic
calculation applies that discipline to a different source-addition module
and finds transitivity.  THM-3881 now proves that a different full cusp-ideal
lift contracts to a rank-two pair and one norm equation, and closes its
`T=0` square lane.  That is strong signal for the rank-two hidden-control task
below, but its surface, residual, and matrix factorization have not been
identified with the Russell seed.  The cheapest lawful test is an actual
intertwiner from the Russell source-addition module to its two-coordinate
norm pair; without one, no conclusion transfers.

### Incoming place and sign packets

THM-3879 supplies a rational sextic discriminant curve whose connected `C3`
layer has a unique best two-place line and no one-place line.  THM-3882
explains the obstruction intrinsically for every everywhere-immersed rational
dual: a line divisor is twice a projection base fibre plus its ramification,
so one-place support is impossible in primal degree at least three.  The
sextic normal-strip packet is a coefficient-space boundary problem, not a
proved rational dual curve.  A decisive future test is to eliminate and
normalize the surviving `(6,4)` boundary, then determine whether a dual
projection model exists.  Until that source--target map is constructed,
THM-3882 transfers no sextic-strip obstruction.

THM-3880 gives an exact matching-sign carrier criterion on nodal--cuspidal
curves.  The quadratic extension in `(6,4)` has its own involution
`(eta,y)->(-eta,-y)`.  Sending its two sheets to signs preserves the first
descent character but forgets coefficient valuations, higher jets, and the
Hamiltonian buckets.  The next strict-transform computation should therefore
record the involution eigencharacter of `u_2,v_4`, the regular arm, and every
conserved value.  Only if its eliminated boundary is actually nodal--cuspidal
can THM-3880 be invoked; otherwise this is a parity checklist, not a theorem
transfer.

## 10. Generated task portfolio

### Anchor: resolve the `(6,4)` strict transform

Continue after the verified first dividend `(17a)`: insert `u_2,v_4`, the
regular target-arm coefficient, `E,G,A0`, and both conserved constants;
divide the `u^4` exceptional multiplicity and compute the first nonzero
constant bracket.  Carry the quadratic-extension involution through every step.
Positive signal: a finite new jet scheme with a nonzero bracket candidate.
Stopping signal: a unit ideal or a forced common-power factor after descent.

### Niche: finish all `(6,5)` valuation channels

The fourteen-point calculation closes one balanced face only.  Enumerate the
Newton fan and audit pole, zero, unit, identically-zero, and constant-scale
channels.  Every face must retain both conserved polynomials, the nonexact
`x^1` row, both arms, and the constant bracket.

### Seed lane: quotient by the legal source group

Replace the full cosmetic quotient by the formal determinant-one/Hamiltonian
source action, or retain the determinant cocycle explicitly.  Compute the
first cohomology/residual class of `(20)` under this restricted action.  Then
test rational coefficients decaying at `u=infinity`; affine polynomial
corrections on the known pole chart are exhausted.

### Geometry lane: abandon a single primitive element

Search for a finite rank-three subalgebra

```text
k[A,C] direct-sum M,                     rank(M)=2,        (27)
```

inside the Russell function field/ring, with binary-cubic multiplication and
no global index-form unit.  Audit normality, constant units, branch companions,
and all height-one valuations before asking for a plane atlas.

### Wildcard scheduler: arithmetic edge coloring

Maintain the task graph `(3)` with the sidecar

```text
(native problem tag, gcd/filter, operation word, evaluator, proof status).
```

Use ordinal rank only to choose the next edge.  The exact rank-eight hostile
must be a regression test for any cross-problem task generator.

### Incoming-signal audit

Consume THM-3877--3882 only at their promoted scopes.  Use them to generate a
global sign-kernel audit for any finite cover, a descent-character table for
`(6,4)`, a normalization-place/dual-projection test for any eliminated sextic
boundary curve, and a literal-intertwiner search from the Russell seed to
THM-3881's rank-two norm grammar.  The last search is necessary because
theorem status on one surface does not identify modules on another.

Incoming THM-3884 is only a **RESERVED / UNPROVED EMPTY STUB** for a
total-degree filtration of that norm equation.  Its useful procedural signal
is already testable: when an actual filtration is proposed, compare
`gr(ker)` with `ker(gr)` and retain the lower transport equations, using the
THM-3412 deficit as the hostile template.  The reservation itself proves no
filtration, lift, or residual obstruction.

THM-3885 is likewise only a **RESERVED / UNPROVED EMPTY STUB**, now for the
complementary `f=0` boundary.  The cheapest generated probe is nevertheless
precise: substitute `f=0` into THM-3881's two proved sidecars
`L^2 f` and `P f^2-T^2`, split off the already closed `T=0` corner, and apply
the proposed Mason/degree gate to `T!=0`.  This is a task specification, not
the reserved theorem.

## 11. Root companion reproduction

Run

```bash
python3 -B 04-computation/jc2_normal_strip_ordinal_gcd_sidecar_probe_root_20260823.py
python3 -B -O 04-computation/jc2_normal_strip_ordinal_gcd_sidecar_probe_root_20260823.py
```

Both streams byte-match
`05-knowledge/results/jc2_normal_strip_ordinal_gcd_sidecar_probe_root_20260823.out`.
The companion checks the signed triangular identity, odd-square peel,
triangular round trips, complete-graph/conductor count, the rank-eight
cross-predicate hostile, `25,914,570` valuation-ray cells, the all-row shear
criterion through degree sixty, and exact polynomial top buckets.

```text
script SHA-256   c80fa6db0e8b66c2e96d1634c779b879e5b754f5d6d5ade3cdbd90b0f500f3de
output SHA-256   6d3a7d9cbfa48cb10083f4712cf6317a08a5e7b1bdba0702feb4c3e410330765
semantic SHA-256 dc1fb2f76a1eafd90a4b7ca95ca93579b6d06d7f54130e2cd171c8418c44341a
```

## 12. Stopping certificate

This session did not prove the degree-six normal-strip theorem and did not
prove or disprove `JC(2)`.  It reduced the top sextic frontier to two rows,
closed the principal balanced `(6,5)` face, found the exact nonreduced cusp
blocking `(6,4)`, and refuted the supposed cubic quotient invariant.  It also
gave a lossless natural-number scheduler for every bounded top cell and the
minimal cross-problem hostile showing why that scheduler needs a native
arithmetic sidecar.

The next honest work is no longer “try more cubic coefficients” or “solve
only the sextic arms.”  It is the `(6,4)` strict transform, the remaining
`(6,5)` valuation fan, and a determinant-aware or intrinsically nonmonogenic
hidden-coordinate model.
