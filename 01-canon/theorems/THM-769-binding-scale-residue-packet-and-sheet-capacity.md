---
id: THM-769
title: Binding-scale residue packets, sheet capacity, and the exact folded two-tightener criterion
status: PROVED (elementary, conditional only on the repository's settled lower-dimensional LRC input; does not prove sporadic-branch emptiness)
source: codex-2026-07-14-S4 (n=12 residue-forcing audit)
depends_on:
  - THM-593   # unit-residue pinning, explicitly conditional on no multiple of 13
  - THM-668   # pair-sum-ruler description of rational maximizers
  - LRCUpTo13 # lower-dimensional margin for a proper on-sheet core
related:
  - THM-617   # shift-pigeonhole precursor
  - THM-761   # multi-exception sheet-cover / gcd-descent interface
  - THM-765   # hereditary primitivity / tooth decks
  - THM-768   # prime-multiple maximum perturbation
  - HYP-6820  # n=12 sporadic-branch audit
---

# THM-769 — Binding-scale residue packets and sheet capacity

Put

```text
phi_A(t) = min_(v in A) ||vt||,       delta = 1/13.
```

Let `A` be a primitive set of twelve positive integer speeds and suppose
`M(A)=delta`.  This theorem separates the familiar full-nonzero-residue
picture from the branch that it does not see.

## 1. The residue packet at an arbitrary binding scale

Choose a global maximizer in lowest terms,

```text
t* = p/Q,       gcd(p,Q)=1,       phi_A(t*)=delta.
```

Then `Q=13s` for an integer `s>=1`.  If `r_v` is the representative of
`pv (mod Q)` in `{0,...,Q-1}`, then

```text
s <= r_v <= 12s                                             (1)
```

for every `v`, and both endpoints of (1) occur.  In particular there are
active speeds `x,y` with

```text
px = s (mod 13s),       py = -s (mod 13s),
13s | x+y,              s | x, y.                           (2)
```

Thus the active pair sum is divisible by the reduced denominator and the
associated pair-sum ruler carries a factor `13`.

Define the on-sheet and off-sheet packets

```text
E = {v in A : s|v} = sU,       F = A\E.
```

At `t*`, every `w in F` lies strictly inside the band in (1).  The on-sheet
residues satisfy

```text
pu (mod 13) in {1,...,12}       for u in U,
```

and both `+1` and `-1` occur.  Moreover, on a neighbourhood of `t*`,

```text
phi_A(t) = phi_E(t) = phi_U(st),
```

so `p/13` is a local maximum of `phi_U` of height `1/13`.

### Proof

All values `||vp/Q||` are integer multiples of `1/Q`.  Since at least one is
exactly `1/13`, `13|Q`; write `Q=13s`.  The assertion that every clearance is
at least `1/13` is exactly (1).

At a local maximum of height below `1/2`, active tent functions must occur on
both slope sides: otherwise a sufficiently small displacement in one
direction raises every active constraint.  The rising active tents have
residue `s` and the falling active tents residue `12s`, proving endpoint
occurrence.  Adding the two endpoint congruences and using `gcd(p,13s)=1`
gives `13s|x+y`.  Each endpoint congruence also implies `s|px` or `s|py`;
because `gcd(p,s)=1`, this proves (2).

The same argument shows that any endpoint owner is in `E`, hence every member
of `F` is strictly interior.  Its clearance therefore stays strictly above
`1/13` near `t*`.  The two oppositely sloped endpoint owners in `E` keep
`phi_E` at most `1/13` on both sides.  This proves the local identity and the
claimed local maximum after the change of variable `tau=st`.  Notice that it
does **not** prove `M(U)=1/13`; usually `U` is looser.

## 2. The shallow branch is exactly the full-residue branch

For a tight twelve-set `A`, the following are equivalent:

1. no speed in `A` is divisible by `13`;
2. reduction modulo `13` maps `A` bijectively onto `{1,...,12}`;
3. every `a/13`, `a=1,...,12`, is a global maximizer;
4. some global maximizer has reduced denominator `13`.

Indeed, (1) makes every runner at least `1/13`-safe at every `a/13`; tightness
makes each such point a global maximizer.  The two-sided slope condition there
forces residues `+1` and `-1` after multiplication by `a`.  Given any nonzero
class `u`, take `a=u^(-1)` and read the `+1` owner: every class occurs.  There
are exactly twelve speeds, so the residue map is bijective.  This is precisely
the prime-`13` specialization of THM-593 (and the hypothesis excluding a
multiple of `13` is essential).  The remaining implications are immediate:
a full residue system has no zero, while a multiple of `13` has clearance zero
at every fraction of reduced denominator `13`.

Consequently, if `A` contains a multiple of `13`, all shallow points are
blocked and **every** global maximizer has `s>=2`.  This is the deep branch.
The common endpoint-splice target on twelve nonzero residue classes applies
only after the shallow branch has been justified; it is not a consequence of
tightness alone.

### Exact one-point criterion

For a unit `p (mod 13)`, the point `p/13` is a local maximum of height exactly
`1/13` if and only if the multiplied residue multiset

```text
{pv (mod 13) : v in A}
```

avoids `0` and contains both `1` and `12`.  The other ten entries may repeat
arbitrarily.  This is the precise reason that one binding point cannot force a
complete nonzero residue system; completeness comes from applying the slope
condition at **all** twelve shallow maximizers.

## 3. Exact sheet capacity in the deep branch

Assume that the chosen maximizer has `s>1`.  Primitivity makes `F` nonempty.
Put `e=|E|`, so `e<=11`.  The settled lower-dimensional LRC bound gives

```text
M(U) >= 1/(e+1) >= 1/12 > 1/13.                             (3)
```

Choose one global maximizer `tau_0` of `U` and consider its `s` lifts

```text
t_j = (tau_0+j)/s,       j in Z/sZ.
```

Every on-sheet speed has clearance `M(U)` at every lift.  Since `A` is tight,
each lift must therefore lie in the **closed** `1/13`-danger tooth of at least
one `w in F`.

For `w in F`, put

```text
g_w = gcd(w,s),       D_w = s/g_w >= 2.
```

As `j` varies, the phases of `w` form a translate of a `D_w`-grid, each grid
point repeated `g_w` times.  A closed circular arc of length `2/13` contains
at most `floor(2D_w/13)+1` points of a `D_w`-grid.  Hence `w` covers at most

```text
g_w (floor(2D_w/13)+1)
```

of the `s` sheets.  The necessary sheet-capacity inequality is therefore

```text
sum_(w in F) (floor(2D_w/13)+1)/D_w >= 1.                   (4)
```

For every integer `D>=2`,

```text
c(D) = (floor(2D/13)+1)/D <= 1/2,
```

with equality if and only if `D=2`.  (At `D=2` this is immediate; for
`D>=3`, `floor(2D/13)+1 <= 2D/13+1 < D/2`.)  It follows from (4) that

```text
|F| >= 2.                                                   (5)
```

There is also a useful ramification corollary.  If `r=|F|<=6`, then
`c(D)<=2/13+1/D` in (4) gives

```text
sum_(w in F) 1/D_w >= 1-2r/13,
min_(w in F) D_w <= 13r/(13-2r).                            (5a)
```

Thus for `r=2,3,4,5,6` some off-sheet runner has respectively
`D_w<=2,5,10,21,78`.  Equivalently, every small-exception deep packet has a
runner with a large gcd against `s`; this is the direct interface with gcd
descent and multi-exception sheet analysis.

If `|F|=2`, both summands in (4) must equal `1/2`, so both tighteners have
`D_w=2`.  Thus `s` is even and both off-sheet speeds are divisible by `s/2`;
the on-sheet speeds are divisible by `s`.  Primitivity forces `s/2=1`.
Therefore

> **Every primitive deep packet with exactly two off-sheet speeds has
> `s=2` and is exactly ten even speeds plus two odd speeds.**

This is sharper at the tight threshold than the generic shift-pigeonhole
bound of THM-617.  It also eliminates the one-tightener deep branch outright.

For three tighteners there is a second exact edge.  If no `D_w` equals `2`,
then `D_w>=3` and

```text
c(D_w) <= 1/3,
```

with equality only at `D_w=3`.  (For `4<=D<=6` this is immediate, and for
`D>=7` use `c(D)<=2/13+1/D<1/3`.)  Capacity forces all three `D_w=3`.
All speeds are then divisible by `s/3`, so primitivity forces `s=3`.
Consequently:

> **Every primitive deep packet with exactly three off-sheet speeds either
> has a half-sheet tightener (`D_w=2`), or has `s=3` and is exactly nine
> multiples of `3` plus three nonmultiples of `3`.**

## 4. Exact folded criteria at `s=2` and the `s=3` equality edge

Write the two-tightener residual as

```text
A = 2U union {x,y},       |U|=10,       x,y odd,
G_U = {tau in R/Z : phi_U(tau)>1/13}.
```

Represent `tau in R/Z` by its unique element of `[0,1)`.  For `w` odd and
`||w tau||<=2/13`, let `N_w(tau)` be the unique nearest integer to `w tau` and
let

```text
epsilon_w(tau) = N_w(tau) (mod 2).
```

Changing the representative by one flips both odd-runner parities, so the
inequality between the two colours below is intrinsic even without this
convention.

Then `A` is tight if and only if, for every `tau in G_U`,

```text
||x tau|| <= 2/13,       ||y tau|| <= 2/13,
epsilon_x(tau) != epsilon_y(tau).                           (6)
```

### Proof

The two lifts of `tau` under doubling are

```text
t_0=tau/2,       t_1=(tau+1)/2.
```

At both lifts the even core has clearance `phi_U(tau)`.  For odd `w`, the two
phases differ by `1/2`.  The runner `w` is closed-dangerous at one of the two
lifts exactly when `||w tau||<=2/13`.  Since `2/13<1/2`, the nearest integer
is unique; even parity kills `t_0` and odd parity kills `t_1`.  As
`1/13<1/4`, one odd runner cannot kill both lifts.

If `tau` is outside `G_U`, the even core already blocks both lifts.  If it is
inside `G_U`, both lifts are blocked exactly when the two odd runners are
eligible and have opposite parities, which is (6).  Thus (6) is equivalent to
`phi_A(t)<=1/13` for all `t`.  The lower-dimensional LRC theorem gives the
reverse inequality `M(A)>=1/13`, completing the equivalence.  The strict sign
in `G_U` and the closed signs in (6) are essential at endpoints.

This turns the last two-tightener packet into a two-colour interval-cover
problem on the loose set `G_U`: each odd runner must cover all of `G_U` under
the doubled danger map, and their nearest-integer colours must be opposite at
every point.

### The three-sheet equality packet

The `s=3` edge above has an exactly parallel description.  Write

```text
A = 3U union {x,y,z},       |U|=9,       3 does not divide xyz,
G_U = {tau in R/Z : phi_U(tau)>1/13}.
```

Use the representative `tau in [0,1)`.  Whenever
`||w tau||<=3/13`, let `N_w(tau)` be the unique nearest integer and define the
sheet colour

```text
j_w(tau) = -w^(-1) N_w(tau) (mod 3).                       (7)
```

Then `A` is tight if and only if, for every `tau in G_U`,

```text
||x tau||, ||y tau||, ||z tau|| <= 3/13,
{j_x(tau),j_y(tau),j_z(tau)} = Z/3Z.                       (8)
```

Indeed, the three lifts are `t_j=(tau+j)/3`.  An exceptional runner `w`
kills a lift precisely when `w tau` is within `3/13` of an integer `N` and
`N+wj=0 (mod 3)`, which is exactly (7).  Since `1/13<1/6`, it kills at most
one lift.  Thus all three core-safe lifts are blocked exactly when the three
eligible exceptions own all three colours.  Outside `G_U` the core blocks
every lift, and the LRC lower bound again converts global blocking into
tightness.  Changing the representative of `tau` translates all three sheet
colours together, so (8) is intrinsic.

Equations (6) and (8) expose the recursive object: a tight deep packet is a
persistent `s`-colour ownership of the lifts over the entire loose set of its
quotient core, not merely a residue condition at one maximizer.

## 5. Guardrails against overstrong residue claims

All of the following are exact elementary checks.

1. A local shallow binding point does not pin all residues.  For

   ```text
   B={1,2,3,4,5,6,7,8,9,10,12,14},
   ```

   `t=1/13` is a local maximum of height `1/13`, bound by residues `1` and
   `12`, but residue `1` is repeated and residue `11` is missing.  In fact
   `M(B)=1/11`: the core `{1,...,10}` is an upper bound and `t=1/11` attains
   it.

2. Full residues do not imply an integer arithmetic progression or tightness.
   The non-AP set

   ```text
   C={1,2,...,11,25}
   ```

   is a complete nonzero residue system modulo `13`, but `M(C)=1/12`, by the
   core `{1,...,11}` and the witness `t=1/12`.

3. Covering, hereditary primitivity, and a deep local binding point still do
   not force shallow residues.  The set

   ```text
   D={2,3,4,5,6,7,8,9,10,11,13,24}
   ```

   has a multiple of every `q=2,...,13`; every leave-one-out set is primitive;
   and `t=1/26` is a local `1/13` maximum bound by `2+24=26`.  Its residues
   modulo `13` include `0`, duplicate `11`, and miss `1,12`.  It is not tight:
   `t=1/15` already has clearance `2/15>1/13`.

4. Primitivity cannot be omitted from a no-multiple statement.  The dilation
   `13*{1,...,12}` is tight by dilation invariance, yet all of its speeds are
   zero modulo `13`.  Thus the no-multiple conjecture must be stated only
   after primitive normalization.

## 6. Tournament / information-preservation audit

The natural vertices in Section 3 are the sheet indices `j in Z/sZ`, and the
deciding data are the tightener-to-sheet closed-danger incidences.  This is a
bipartite cover (or, after projecting out tighteners, a sheet hypergraph), not
primarily a tournament on runners.  A runner-pair tournament discards the
simultaneous-cover requirement and the grid multiplicities `g_w`.

At `s=2`, the quotient becomes two sheet vertices with perfect two-colour
ownership; nearest-integer parity is the gauge.  A tournament can record which
odd runner owns which sheet, and its tie Hamiltonian path is trivial, but the
actual invariant is the persistent opposite-colour condition over every
component of `G_U`.  Its useful fingerprints are colour flips at doubled-tooth
endpoints and failures of persistent ownership, rather than directed cycles
among speeds.

This explicitly challenges the assumption that runners or binding arcs must
be the tournament vertices.  The sheet-incidence quotient preserves the LRC
predicate in (4) and (6); a runner tournament destroys exactly the joint
coverage information needed for the proof.

## 7. Honest frontier

THM-769 does not prove that the primitive tight twelve-speed sporadic branch is
empty.  It splits the missing theorem into two structurally different tasks:

- **shallow/full-residue branch:** prove lift/splice coherence strongly enough
  to turn the full residue system into the primitive AP;
- **deep branch:** for two tighteners, rule out the exact `s=2` folded parity
  cover (6); for at least three tighteners, classify the sheet covers allowed
  by (4), together with their persistence over the loose set of `U`.

The earlier single-scale residue picture was therefore missing a coordinate:
the reduced binding scale `s`, together with off-sheet incidence.  Full
nonzero residues are not the underlying object; they are the height-one fibre
of this scale-residue packet.
