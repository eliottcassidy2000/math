---
id: THM-797
title: Odd-exception divisor grids, the q=13 signed-wall gate, and the global-component correction
status: PROVED (general odd-grid residue shell, q=13 signed-wall rigidity, and metric pin) + VERIFIED (gate-sharp rows, exact deep components, and global erosion witnesses); SCOPE CORRECTED by THM-803's quarter-grid escapes
source: codex-2026-07-14-S10 global-erosion selection analysis
depends_on:
  - THM-772   # two-sheet packet arithmetic and mandatory deep exception
  - THM-774   # folded-diamond equivalence
  - THM-789   # global erosion target
related: [THM-769, THM-775, THM-803, HYP-6800, HYP-6820]
verification:
  - 04-computation/lrc13_odd_exception_divisor_grid_codex_S10.py
  - 05-knowledge/results/lrc13_odd_exception_divisor_grid_codex_S10.out
  - 04-computation/lrc13_antigrid_all_component_selector_codex_S10.py
  - 05-knowledge/results/lrc13_antigrid_all_component_selector_codex_S10.out
---

# THM-797 — Odd-exception divisor grids and the q=13 signed-wall gate

Use THM-772's two-sheet notation

```text
A=2U union {x,y},       |U|=10,       x,y odd,
E_U={t:phi_U(t)>=1/11},
R_U={d:max_(u in U)||ud||<2/143}.
```

Put

```text
a=(x+y)/2,              b=|x-y|/2,
H_(x,y)={t:||at||+||bt||>=11/13}.
```

Thus THM-774 identifies `H_(x,y)` with the two-exception eligibility and
opposite-sheet-colour predicate, and THM-789's residual target is

```text
E_U not subset H_(x,y) minus R_U.                         (0)
```

## 1. The general odd-grid theorem

Let `q>=3` be odd.  Write `[z]_q` for the unique balanced residue in

```text
{-(q-1)/2,...,(q-1)/2}.
```

For a unit numerator `p mod q`, considered modulo sign, define

```text
D_q(U)={p mod +/- :
        |[up]_q|>=ceil(q/11) for every u in U},           (1)

A_q(x,y)={p mod +/- :
          |[xp]_q|,|[yp]_q|<=floor(2q/13),
          [xp]_q != [yp]_q (mod 2)}.                      (2)
```

Both predicates are invariant under `p->-p`.

### Theorem 1 — deep grids and exception shells

For every unit `p mod q`,

```text
p/q in E_U       iff p in D_q(U),                         (3)
p/q in H_(x,y)   iff p in A_q(x,y).                       (4)
```

Consequently,

```text
D_q(U) not subset A_q(x,y)
  implies E_U not subset H_(x,y) minus R_U.               (5)
```

If `q|x`, then (2) simplifies to the explicit shell

```text
A_q(x,y)={p mod +/- :
          q does not divide y,
          1<=|[yp]_q|<=floor(2q/13),
          [yp]_q is odd}.                                 (6)
```

In particular, if `q` divides both exceptions, this shell is empty.

### Proof

For (3), simply observe that

```text
||up/q||=|[up]_q|/q.
```

Because the numerator is integral, the inequality `>=1/11` is exactly the
threshold in (1).

For (4), write

```text
xp=qN_x+r_x,       yp=qN_y+r_y,
r_x=[xp]_q,        r_y=[yp]_q.                            (7)
```

Since `q` is odd, the balanced residues have magnitude strictly less than
`q/2`, so `N_x,N_y` are the unique nearest integers.  Closed exception
eligibility is

```text
||xp/q||,||yp/q||<=2/13,
```

equivalently

```text
|r_x|,|r_y|<=floor(2q/13).                               (8)
```

Modulo two, `q,x,y` are all one.  Formula (7) therefore gives

```text
N_x=p-r_x,          N_y=p-r_y          (mod 2).           (9)
```

The two sheet colours are opposite exactly when `N_x,N_y` have opposite
parity, which by (9) is exactly the last condition in (2).  THM-774 identifies
this lifted predicate with the folded diamond, proving (4).

Now `0 in R_U`, so morphological erosion satisfies

```text
H_(x,y) minus R_U subset H_(x,y).                         (10)
```

A class in `D_q(U)\A_q(x,y)` supplies through (3)--(4) a point of
`E_U\H_(x,y)`, and (10) proves (5).

Finally suppose `q|x`.  Then `r_x=0`.  Also `x/q` is odd, and the exact
nearest integer `N_x=(x/q)p` has parity `p`.  The two colours are opposite
precisely when `r_y=[yp]_q` is odd.  Eligibility supplies the magnitude bound;
an odd residue is nonzero, so `q` cannot divide `y`.  This proves (6). ∎

## 2. The mandatory q=13 support gate

Now assume the arithmetic conclusions of THM-772 for a hypothetical primitive
tight two-sheet packet.  In particular, no member of `U` is divisible by 13,
and at least one deep-branch exception is.

Let

```text
C=(Z/13Z)^*/{+/-1},
S(U)={+/-u mod 13:u in U} subset C.                       (11)
```

Inversion acts on these six folded classes.  Since

```text
ceil(13/11)=2,
```

condition (1) says that no `up` is `+/-1`.  Hence

```text
D_13(U)=C minus S(U)^(-1).                                (12)
```

Also `floor(26/13)=2`.  If `13|x` but `13` does not divide `y`, the only odd
nonzero balanced residues in the shell (6) are `+/-1`, so

```text
A_13(x,y)={+/-y^(-1)}.                                   (13)
```

Therefore hypothetical containment in (0), or even the weaker containment
`E_U subset H_(x,y)`, forces

```text
C minus S(U)^(-1) subset {+/-y^(-1)},
```

and inversion gives the exact support gate

```text
C minus S(U) subset {+/-y}.                               (14)
```

Thus

```text
S(U)=C
or
S(U)=C minus {+/-y}.                                     (15)
```

If both `x` and `y` are divisible by 13, (6) is empty.  Containment forces
`D_13(U)` empty, hence

```text
S(U)=C.                                                   (16)
```

This proves the following uniform eliminations before any height search:

1. every core supported on at most four folded mod-13 classes;
2. every five-class core whose missing class is not the folded class of the
   off-divisor exception; and
3. every non-full support core when both exceptions are 13-divisible.

The only q=13-invisible cases are full six-class support and the **aligned
five-class** pattern in (15).

## 3. The signed-wall upgrade and the aligned metric pin

The two support patterns left by (15) are not equally viable.  The signs
discarded by folded support supply a second, local obstruction.

### Theorem 2 — signed-wall rigidity

Assume that the two-sheet packet is hypothetically tight and label the
exceptions so that `13|x`.  Put `B=max(U)`.  Then:

1. `13` does not divide `y`;
2. full six-class support is impossible;
3. the support is aligned five-class and, as a multiset of signed residues,

   ```text
   {u mod 13:u in U}
     =(Z/13Z)^* minus {+y,-y};                            (17)
   ```

   in particular every one of the ten displayed residues occurs exactly
   once;
4. at `p=y^(-1) mod 13`,

   ```text
   phi_U(p/13)=2/13,       M(U)>=2/13;                    (18)
   ```

5. the exceptions obey the sharpened metric bounds

   ```text
   x<=2B-1,                y<=B-1,                        (19)
   13B+2xy<=2B(x+y).                                      (20)
   ```

Thus the full-support and double-13 branches are empty.  The only q=13
survivor is not merely folded-aligned: it is the exact signed complement of
the exception pair `+/-y`.

### Proof

First suppose `13` divides both exceptions.  For any folded class `c` in
`S(U)`, choose `p` with folded class `c^(-1)`.  Then

```text
phi_U(p/13)=1/13,
```

and the minimum owners are exactly the core runners congruent to `+c` or
`-c`.  The double-divisor shell is empty, so `p/13` lies in the open
complement of `H_(x,y)`.  If all minimum owners have the same signed residue
after multiplication by `p`, move a sufficiently small distance in the
direction that increases their absolute residues.  Every minimum owner then
has clearance greater than `1/13`; all other core runners start at clearance
at least `2/13` and remain above `1/13`.  This produces a one-sided interval
of `G_U` in the open complement of `H_(x,y)`, contradicting tightness.

Consequently every present folded class must be represented by both signs.
The support gate (16) requires all six classes, hence at least twelve core
runners, impossible.  This proves `13` does not divide `y`.

Now the q=13 shell contains only the folded numerator `y^(-1)`.  Repeat the
same one-sided argument for every present class `c` other than the folded
class of `y`: its inverse grid point lies outside the shell, so both signed
residues `+c,-c` must occur in `U`.  Full support would require ten runners
on the five rejected classes and at least one more on the accepted class, so
it is impossible.  In the aligned pattern the five present classes require
exactly ten runners.  Hence each sign occurs once and (17) follows.

Choose the sign of `p` so that `yp=1 mod 13`.  Multiplication by `p` turns
(17) into the ten signed residues

```text
+/-2,+/-3,+/-4,+/-5,+/-6.
```

This proves (18).  The centered interval of radius

```text
(2/13-1/13)/B=1/(13B)                                   (21)
```

about `p/13` lies in `G_U`, hence in `H_(x,y)`.  The connected `x`-eligible
tooth there has radius `2/(13x)` because `13|x`.  The `y`-error is `+1/13`,
so its eligible tooth extends only `1/(13y)` in the narrow direction.
Containment of (21) in both teeth gives

```text
x<=2B,                   y<=B.                            (22)
```

The first inequality is strict because `x` is odd and `2B` is even.  Equality
in the second would put the core speed `B=y` in the folded class omitted by
(17).  Thus (19) follows.

There is also an endpoint-owner version of this last strictness.  If one
temporarily assumes `y=B`, the narrow endpoint is

```text
p/13+1/(13B).
```

Among the normalized residues `+/-2,...,+/-6`, a core runner can have
clearance at most `1/13` there only when it is the speed `B` with normalized
residue `-2`.  Tightness would therefore force

```text
[Bp]_13=-2,             equivalently B=-2y mod 13,        (23)
```

which is also incompatible with `B=y` and `13` not dividing `y`.  This is
the signed endpoint owner hidden by the folded support quotient.

Finally (18) gives

```text
rho=(M(U)-1/13)/B >= 1/(13B).
```

Substitution in THM-772's inequality `(A*)`, followed by multiplication by
`13Bxy`, gives (20). ∎

## 4. Gate sharpness and the global-component correction

> **Scope correction (THM-803).**  The two rows in this section are sharp for
> the stated folded/signed `q=13` and exception-divisor gates, not for all
> mandatory grids.  For both rows, the quarter-grid point `t=11/52` has
> `phi_U(t)=5/52>1/11` and `Q_(9,4)(t)=1/4<11/13`.  Thus neither survives the
> combined `26,52,78` anti-grid ladder.  THM-803 supplies a different exact
> signed-wall-compatible row showing that even the combined gates still need
> an all-component selector.

The aligned pattern really can trap the entire prime grid.  Take

```text
U_0=(1,2,3,4,7,9,10,11,12,16),
(x,y)=(13,5).                                             (24)
```

This core is primitive, contains a multiple of every `m=2,...,12`, has no
13-multiple, and has `B=16`.  Exact breakpoint enumeration gives

```text
M(U_0)=2/13,
argmax phi_(U_0)={5/13,8/13}.                             (25)
```

The breakpoint certificate is finite and lossless: `phi_U` is piecewise
linear, and all its possible maxima occur at a runner cusp `k/(2u)` or an
intersection of two signed linear pieces, whose denominators divide `u+v` or
`|u-v|`.  The verifier enumerates exactly these rationals and evaluates them
with `Fraction` arithmetic.

The row also satisfies every earlier two-sheet arithmetic tax:

```text
x,y<=11B,             min(x,y)<=11B-36,
rho=(M(U_0)-1/13)/B=1/208,
1/(xy)+2rho=1/40 <= 36/845
             =2/(13x)+2/(13y).                           (26)
```

The folded support multiplicities on the classes `1,...,6` are

```text
(2,2,3,2,0,1).                                           (27)
```

Thus `S(U_0)=C\{5}` and the missing class is exactly the class of `y`.
The only deep thirteenth numerators are `p=5,8`; at both,

```text
phi_(U_0)(p/13)=2/13,
Q_(9,4)(p/13)=12/13.                                     (28)
```

So every point visible to the q=13 gate lies inside the folded diamond.

Another global escape is visible at a different component endpoint:

```text
t=7/33,
phi_(U_0)(t)=1/11,
Q_(9,4)(t)=8/33<11/13.                                   (29)
```

Thus `t in E_(U_0)\H_(13,5)`, and since `0 in R_(U_0)`, (29) proves the full
erosion failure (0).

This folded-gate-sharp row also lies in the exact **sporadic max-peel regime**.
Its full packet and maximum deletion are

```text
A_0=(2,4,5,6,8,13,14,18,20,22,24,32),
M(A_0)=2/19 at t in {2/19,8/19,11/19,17/19},
M(A_0\{32})=1/8 at t in {1/16,7/16,9/16,15/16}.
```

Hence the deletion is strictly super-lonely at the `1/12` threshold, while
the full packet is loose rather than tight because `2/19>1/13`.  Its deepest
`q=13` choices remain trapped by (28), but the different component endpoint
`7/33` supplies the escape in (29).

The complete closed `1/11`-superlevel set has eight components:

```text
[12/77,7/44],       [23/110,7/33],
[67/176,43/110],    {9/22},
{13/22},             [67/110,109/176],
[26/33,87/110],      [37/44,65/77].                       (30)
```

Their endpoint-owner pairs are

```text
7->12, 10->9, 16->10, {10,12},
{10,12}, 10->16, 9->10, 12->7.                           (31)
```

They form four reflection orbits.  The orbit containing the global maxima
`5/13,8/13` is entirely trapped: its minimum `Q` is `159/176>11/13`.
The second orbit contains (29) and escapes with a large margin.  Therefore
**selecting only global maximizers is insufficient**.  The faithful global
state must retain every closed `1/11`-superlevel component, including
threshold endpoints and endpoint owners.

### Sharpness after the signed-wall upgrade, within exception-divisor grids

The signed-wall theorem is still not a closure theorem.  Replace `U_0` by

```text
U_1=(1,2,4,6,7,9,10,11,12,16),
(x,y)=(13,5).                                             (32)
```

This row is primitive, divisor-complete through 12, has no 13-multiple, and
has `B=16`.  Its signed residues are exactly

```text
U_1 mod 13=(Z/13Z)^* minus {5,8};                         (33)
```

so it satisfies the full conclusion (17), not only folded alignment.  Exact
breakpoint enumeration gives

```text
M(U_1)=2/13,        argmax phi_(U_1)={5/13,8/13}.         (34)
```

It obeys (19)--(20) and THM-772's original size and `(A*)` taxes.  Moreover
**every nontrivial odd divisor grid of either exception is silent**:

```text
q=5:       D_5(U_1)=empty       because 5|10,
q=13:      D_13(U_1)={+/-5} subset A_13(13,5).            (35)
```

Nevertheless the non-divisor grid `q=17` supplies

```text
t=6/17,       phi_(U_1)(t)=2/17,
Q_(9,4)(t)=10/17<11/13.                                  (36)
```

Thus (36) proves the full erosion failure (0).  This row is an exact no-go
for any strategy using only exception-divisor grids, q=13 signed walls, and
the scalar arithmetic taxes: the escaping denominator is generated by a
different global component.

## 5. What remains after the signed-wall gate

THM-789 showed why refining one trapped deep anchor cannot prove global
noncontainment.  Theorem 1 supplies a different global selector whenever an
odd exception divisor has a deep unit class outside the other exception's
explicit shell.  The q=13 signed-wall specialization leaves only the exact
residue complement (17), while (32)--(36) prove that this gate remains sharp
relative to exception-divisor grids.  THM-803's quarter-grid catches both
rows, then gives a signed-wall-compatible row on which all three universal
anti-grids are silent and an owner-labelled deep component is still required.

The remaining two-sheet program should therefore branch on two axes:

1. **arithmetic grids:** use exception divisors through (1)--(6), with q=13
   always first, then THM-803's mandatory `26,52,78` anti-grids, while
   retaining other denominators generated by the core; (35)--(36) prove that
   exception divisors alone are incomplete; and
2. **geometric deep components:** when a grid is silent, enumerate or bound
   the owner-labelled components of `E_U` and compare their minimum folded-
   diamond margins, rather than refining only an argmax component.

For the mandatory q=13 branch, only the signed aligned-five pattern (17)
reaches this stage; THM-803 further requires full parity-twisted support and
silence on the anti-grid ladder, then decides the erosion predicate by a
finite all-component endpoint/cusp selector.  No claim of complete two-sheet
or uniform `n=12` closure is made here.

## 6. Tournament Analysis and challenged assumptions

There are two useful but deliberately decorated tournament views.

### Folded residue obligations

Take the six folded classes in `C` as vertices.  For (24), orient by core
multiplicity, breaking ties by the numerical class order; use inversion

```text
(2 6)(3 4)
```

as the switch/gauge and the sorted order as the tie Hamiltonian path.  Before
and after the switch the paths are

```text
(5,6,1,2,4,3),          (5,2,1,3,6,4).
```

The gauges differ on five edges.  Both tournaments are transitive: score
histogram `0,...,5`, zero directed triangles, singleton SCCs, and one
Hamiltonian path.  The bare fingerprints do not expose the missing class or
the exception shell; the decorated incidence `(class,signed occupancy,
multiplicity,inverse,accepted-by-y)` does.  The signed-wall proof sharpens
this information audit: folded support alone regards both `U_0` and `U_1` as
aligned survivors, while the twelve signed residue slots immediately reject
`U_0` and retain only the exact complement (33).

### Deep-component alternatives

Take the four reflection orbits in (30) as vertices.  The first pairwise
observable is escape margin

```text
11/13-min_(t in C) Q_(9,4)(t),
```

and the switch/gauge replaces it by component width.  With index tie-breaking,
the two Hamiltonian paths are

```text
(2,0,3,1),             (3,1,0,2),
```

again with five edge flips.  Both tournaments are transitive with scores
`0,...,3`, zero cycles, singleton SCCs, and one Hamiltonian path.  Yet only
the decorated first gauge remembers the **sign** of escape and the endpoint
owners.

Runner vertices lose the divisor grid; residue vertices lose component
geometry; a single global-max vertex loses alternative components; bare
tournaments lose all metric signs.  The minimal carrier suggested by this
theorem is a bipartite incidence object between folded divisor obligations and
owner-labelled deep components, with signed wall occupancy and exact escape
margins as sidecars.  The all-divisor-grid-silent row (32) adds a further
warning: restricting arithmetic vertices to divisors of the exceptions loses
the detecting `q=17` obligation.  This is the assumption challenge: the
underlying object is a family of global alternatives and core-generated
denominators, not one privileged deepest point or one fixed divisor list.

## Exact replay

The verifier exhausts every odd `q=3,...,101`, every unit numerator, and every
odd pair of speed classes modulo `2q`: `10,971,770` exact tests of the direct
lifted predicate, the balanced-residue shell, and the folded diamond, with no
failures.  It separately checks `143,346` divisor specializations, all `832`
q=13 support/exception profiles, all `352,716` signed ten-residue multisets
(`340,170` off-divisor support-gate pairs and `34,332` double-divisor full-
support profiles), and `201,000` endpoint-owner cases.  The signed-wall census
leaves exactly six off-divisor profiles and no double-divisor profile, with
digest

```text
6d1065da21735e471b238e0c1222cdb1063946c55f8a74d0646dc4fcb1d3efc9.
```

The replay also checks the breakpoint certificates (25) and (34), every
component in (30), witnesses (29) and (36), all exception-divisor grids of
(32), and both tournament fingerprints.  The signed-wall survivor row has
digest

```text
2b5e7a6e4593f998fd28eeb124fb5eaaba8b00747782084ee11004b457ebb546.
```

No floating point or sampled-circle verdict enters the proof.
