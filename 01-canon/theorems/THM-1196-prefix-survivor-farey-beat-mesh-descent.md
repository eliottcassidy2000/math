---
id: THM-1196
title: Prefix-survivor Farey descent and the recursive component-span flood tail
status: PROVED exact all-cardinality beat-mesh implication and three unconditional six-comb tail bounds.  A covered prefix survivor avoids every sum/difference beat lattice joining the last comb to the prefix.  Composing THM-1198's universal dual survivor with recursive component spans proves d_6/d_5<77, d_5/d_4<189/8, and d_4/d_3<343/15 in every putative cover; the phase-aware beat refinements remain exact.  This does not bound the first three prefix ratios, prove universal six-comb noncoverage, or prove LRC(14)
source: codex-2026-07-18 slow-gap descent session
depends_on: [THM-1176, THM-1192, THM-1198]
related: [THM-1178, THM-1179, THM-1193, HYP-7715]
script: 04-computation/lrc14_prefix_survivor_farey_beat_mesh_descent_codex_20260718.py
output: 05-knowledge/results/lrc14_prefix_survivor_farey_beat_mesh_descent_codex_20260718.out
---

# THM-1196 -- prefix-survivor Farey beat-mesh descent

At radius `1/14`, put

```text
D_s={t in R: ||st||<1/14}
```

and let

```text
G=G_k(c)=[(14k+1)/(14c),(14k+13)/(14c)]              (1)
```

be a complete closed safe gap of the carrier `c`.  Suppose

```text
c<d_1<...<d_r,                    2<=r<=6,             (2)
G subset union_(i=1)^r D_(d_i).                         (3)
```

The last comb can only repair the survivor of the prefix.  The point of this
theorem is that every sum and difference beat with that last comb punctures
the putative repair set.  This is an exact, phase-aware finite descent of the
unbounded last-speed coordinate.

## 1. The prefix survivor avoids a Farey beat mesh

Define the closed prefix survivor

```text
E=G minus union_(i=1)^(r-1) D_(d_i).                   (4)
```

For every `i<r`, put

```text
q_i^+=d_r+d_i,                 q_i^-=d_r-d_i,          (5)
B_i^+={p/q_i^+ : p in Z, p/q_i^+ in G},
B_i^-={p/q_i^- : p in Z, p/q_i^- in G},                (6)
B=union_(i<r)(B_i^+ union B_i^-).                      (7)
```

All denominators in (5)--(6) are positive by (2).

> **Lemma 1 (beat-mesh puncture).**  Under (2)--(3),
>
> ```text
> E intersect B is empty.                              (8)
> ```

**Proof.**  Take `t=p/q_i^+`.  Since `(d_r+d_i)t=p`,

```text
||d_r t||=||p-d_i t||=||d_i t||.                      (9)
```

For `t=p/q_i^-`, the identity `(d_r-d_i)t=p` gives the same conclusion with
`d_r t=p+d_i t`.  If `t` also belonged to `E`, then in particular
`t notin D_(d_i)`, so (9) would imply `t notin D_(d_r)`.  Definition (4)
already excludes every earlier comb.  Thus no comb in (3) would cover `t`,
a contradiction.  Strict danger teeth are important here: equality at
distance `1/14` is safe for both defining speeds.  This proves (8).  ∎

This is the terminal, one-comb specialization of THM-1192's all-pairs
puncture law.  It gains strength because every complementary fast comb is
already deleted in (4): a safe beat is an immediate lonely point, rather
than a four- or five-comb obligation.

## 2. Exact phase-local mesh and component count

Sort the finite rational set

```text
{left endpoint of G, right endpoint of G} union B
```

and define `M_k(c;d_1,...,d_r)` to be its largest consecutive spacing.
This is an exact phase-local Farey mesh; it retains the truncated block
address `k` that a denominator density discards.  Lemma 1 implies that every
connected component of `E` has length strictly less than `M_k`.

For an exact component bound, let `n_i` be the number of open `d_i`-teeth
which meet `G`, for `i<r`.  Equivalently,

```text
n_i=#{m in Z:
      d_i(14k+1)/(14c)-1/14 < m
       <d_i(14k+13)/(14c)+1/14}.                       (10)
```

Removing `sum_i n_i` relatively open intervals from one closed interval
leaves at most

```text
C=1+sum_(i=1)^(r-1)n_i                                (11)
```

connected components.  Hence (8) gives the exact phase-aware upper bound

```text
|E|<C M_k(c;d_1,...,d_r).                             (12)
```

The interval in (10) has length

```text
6d_i/(7c)+1/7=(6d_i+c)/(7c),                          (13)
```

so the phase-free count

```text
C<=C_bar:=1+sum_(i=1)^(r-1) ceil((6d_i+c)/(7c))       (14)
```

is immediate.  Also the lattice of denominator

```text
q_max=d_r+d_(r-1)                                     (15)
```

is one of the lattices in `B`.  Consecutive points of that lattice are
`1/q_max` apart, and the two boundary remainders are shorter.  Therefore

```text
M_k<=1/(d_r+d_(r-1)).                                 (16)
```

Equations (12) and (16) are different levels of the result: (12) is the
stronger phase-aware Farey screen, while (16) is the closed scalar tail.

## 3. Harmonic-to-Farey descent

Apply THM-1176's all-cardinality survivor bound to the first `r-1` combs:

```text
|E| >=(6/49)((8-r)/c-sum_(i=1)^(r-1)1/d_i).           (17)
```

Write the prefix pressure

```text
H_(r-1)=c sum_(i=1)^(r-1)1/d_i.                       (18)
```

> **Theorem 2 (last-tooth Farey descent).**  Every cover (3) with
> `H_(r-1)<8-r` obeys
>
> ```text
> (6/49c)((8-r)-H_(r-1))
>      < C M_k(c;d_1,...,d_r)                         (19)
> ```
>
> and consequently
>
> ```text
> (d_r+d_(r-1))/c
>      < 49 C_bar/[6((8-r)-H_(r-1))].                 (20)
> ```

**Proof.**  The hypothesis makes the right side of (17) positive.  Combine
(12) with (17) to obtain (19), then use (14) and (16) and clear the positive
denominators.  ∎

For the open six-comb cone, (20) reads

```text
if H_5=c sum_(i=1)^5 1/d_i<2, then

(d_6+d_5)/c
 <49[1+sum_(i=1)^5 ceil((6d_i+c)/(7c))]/[6(2-H_5)].   (21)
```

Thus an unbounded last tooth is impossible over a fixed subcritical prefix.
More uniformly, if

```text
d_i<=Kc (i<r),            H_(r-1)<=8-r-epsilon,       (22)
```

then

```text
d_r/c
 <49[1+(r-1)ceil((6K+1)/7)]/(6epsilon).               (23)
```

This converts the infinite `d_r` tail into a finite search uniformly on
every compact prefix cone separated from its critical face.  Therefore any
sequence escaping (23) must do at least one of the following:

```text
some prefix ratio d_i/c escapes every compact set,
or H_(r-1) approaches/crosses 8-r.                    (24)
```

For `r=6`, the adjusted target left by (24) is the five-prefix critical cone
`H_5>=2`; for `r=5` it is `H_4>=3`.  For `r=4`, the three-term prefix always
has `H_3<3<4`, so every fixed four-comb prefix has a finite exact last-speed
range.  These are genuine reductions, not claims that the critical cones
are empty.  The next section uses the later five-comb dual theorem to remove
the `H_5` alternative specifically in the six-comb case.

## 4. Unconditional six-comb last-ratio bound

THM-1198 supplies exactly the missing pressure-free lower bound.  In the
unit-gap normalization, its nonnegative probability density `f` satisfies

```text
f<=7/6,                                                (25)
integral_(one closed faster danger comb) f<7/36.       (26)
```

The inequality is strict because THM-1198's equality case has normalized
slope `6/7`, whereas `d_i>c` gives slope `6d_i/(7c)>6/7`.  Subtracting five
loads from total `f`-mass one shows that the normalized five-comb survivor
has `f`-mass strictly greater than `1/36`.  Its Lebesgue measure is therefore
strictly greater than

```text
(1/36)/(7/6)=1/42.                                    (27)
```

Rescaling the unit interval to the slow gap of length `6/(7c)`, and enlarging
from the closed-comb survivor to the open-comb survivor, gives

```text
|E_5|>1/(49c),
E_5=G minus union_(i=1)^5 D_(d_i).                    (28)
```

Let `C_5` be the component bound (11) for the first five combs.  Under a
six-comb cover, `E_5 subset D_(d_6)`.  Every connected component of the
compact set `E_5` lies strictly inside one open `d_6` tooth, of width
`1/(7d_6)`.  Hence

```text
|E_5|<C_5/(7d_6).                                    (28a)
```

Combining (28) and (28a) gives the stronger uniform conclusion.

> **Corollary 3 (last-tooth component span).**  Every six-comb open cover
> of a complete `c`-slow gap obeys
>
> ```text
> d_6/c<7C_5
>  <=7[1+sum_(i=1)^5 ceil((6d_i+c)/(7c))],            (28b)
> d_6/d_5<77.                                         (28c)
> ```

For the ratio claim put `x=d_5/c>1`.  Every ceiling in (28b) is at most
`ceil(x)<x+1<2x`, while `1/x<1`, so

```text
d_6/d_5<7[1+5ceil(x)]/x<7(1+10)=77.
```

The Farey mesh supplies an independent phase-sensitive supplement.  Combine
(12), (14), (16), and the weaker non-strict version of (28):

```text
(d_6+d_5)/c
 <49[1+sum_(i=1)^5 ceil((6d_i+c)/(7c))].              (29)
```

Replacing the phase-free `1/(d_6+d_5)` by the exact `M_k` in (12) can improve
(29) on a fixed clock.  The component-span bound (28b)--(28c) is uniformly
stronger as a scalar ratio law.  Neither result yet bounds all four earlier
prefix ratios.

### The `j=4` component-span flood tail

The same dual density controls one rung deeper.  Put

```text
E_4=G minus union_(i=1)^4 D_(d_i),       q=d_5+d_6.   (35)
```

For the four **closed** prefix combs, THM-1198 gives

```text
integral_(E_4) f >=1-4(7/36)=2/9.                    (36)
```

After division by `max f=7/6`, the normalized Lebesgue survivor has length
at least `4/21`.  Rescaling by `6/(7c)`, and then enlarging from the closed-
comb survivor to the open-comb survivor in (35), gives

```text
|E_4|>=8/(49c).                                       (37)
```

The strict load improvement in (26) actually makes (36)--(37) strict, but
the displayed weak form is enough below.  Under the six-comb cover,

```text
E_4 subset D_(d_5) union D_(d_6).                    (37a)
```

A `d_6` tooth has width `1/(7d_6)`, shorter than the safe gap
`6/(7d_5)` between adjacent `d_5` teeth.  It therefore cannot meet two
distinct `d_5` teeth.  Hence a connected component of the open union in
(37a) contains at most one `d_5` tooth.  All incident `d_6` teeth together
can extend by at most one `d_6`-tooth width on either side, so that union
component has span at most

```text
1/(7d_5)+2/(7d_6)<3/(7d_5).                          (37b)
```

Let `C_4` be the component bound (11) for the first four combs.  Each compact
component of `E_4` lies strictly inside one open union component, and thus

```text
|E_4|<3C_4/(7d_5).                                   (37c)
```

> **Corollary 4 (`j=4` component-span flood tail).**  Every six-comb open
> cover of a complete `c`-slow gap obeys
>
> ```text
> d_5/c<(21/8)C_4
>  <=(21/8)[1+sum_(i=1)^4 ceil((6d_i+c)/(7c))],       (37d)
> d_5/d_4<189/8.                                      (37e)
> ```

For the ratio claim put `x=d_4/c>1`.  Every ceiling in (37d) is at most
`ceil(x)<x+1`, so

```text
d_5/d_4<(21/8)[1+4ceil(x)]/x
             <(21/8)(4+5/x)<189/8.
```

This is the functional form of the `j=4` flood tail: the final two combs can
flood a four-prefix survivor only inside short two-scale union components.
It keeps component chronology instead of replacing the two combs by their
combined density `2/7`.

### Independent safe-beat refinement

At the sum-beat points `t=p/q`, one has

```text
||d_5t||=||d_6t||.                                    (38)
```

Consequently `E_4` avoids every **safe** beat point

```text
B_safe={p/q in G: ||d_5p/q||>=1/14}.                 (39)
```

Indeed, the first four combs are safe by membership in `E_4`; (38) makes
both remaining combs safe, including equality at `1/14`, contradicting the
open cover.

It remains to bound the holes in this safe-beat mesh.  Write

```text
g=gcd(d_5,d_6),      Q=q/g,      u=d_5/g,
alpha=u/Q=d_5/q<1/2.                                (40)
```

Suppose `R` consecutive numerators are dangerous for `d_5`.  Their orbit
under addition by `alpha` stays in the open circle arc
`(-1/14,1/14)`, whose length is `1/7`.  If `R>=2`, the difference of two
consecutive points first forces `alpha<1/7`; lifting the run to that one arc
then gives

```text
(R-1)alpha<1/7.
```

The same final bound is trivial for `R=1`, and hence

```text
R<=ceil(q/(7d_5)).                                    (41)
```

Thus consecutive safe beat points, including the two boundary remainders
inside `G`, are separated by at most

```text
(R+1)/q
 <=[ceil(q/(7d_5))+1]/q
 <1/(7d_5)+2/q
 <8/(7d_5),                                           (42)
```

where the last inequality uses `q=d_5+d_6>2d_5`.

Let `C_4` be the component bound (11) for the first four combs.  The same
component-versus-mesh proof as (12), now using (39)--(42), gives

```text
|E_4|<8C_4/(7d_5).                                    (43)
```

Together with (37), the safe-beat mesh gives the independent but scalarly
weaker bounds

```text
d_5/c<7C_4,                  d_5/d_4<63.              (44)
```

The value of this route is arithmetic rather than its uniform constant: the
exact phase-local safe-beat mesh can be much smaller than the envelope (42).
It retains the equality (38), numerator positions, and truncated gap phase,
all of which the stronger component-span constant (37e) intentionally
forgets.

### The three-tooth component recurrence

One more rung closes without any new density.  Put

```text
E_3=G minus union_(i=1)^3 D_(d_i),
w_i=1/(7d_i).                                         (44a)
```

The union of `D_(d_5)` and `D_(d_6)` has component span at most

```text
w_5+2w_6<3w_4<6w_4,                                  (44b)
```

where `6w_4` is the gap between adjacent `d_4` teeth.  Hence no such
two-comb component meets two distinct `d_4` teeth.  A connected component
of `D_(d_4) union D_(d_5) union D_(d_6)` therefore contains at most one
`d_4` tooth.  The attached two-comb components can extend by at most their
full span on each side, so the three-comb component span is at most

```text
w_4+2(w_5+2w_6)<7w_4=1/d_4.                          (44c)
```

For the first three closed prefix combs, THM-1198 gives survivor `f`-mass
strictly greater than

```text
1-3(7/36)=5/12.                                      (44d)
```

Division by `max f=7/6` gives normalized length greater than `5/14`, and
rescaling gives the physical lower bound

```text
|E_3|>15/(49c).                                       (44e)
```

Let `C_3` be the component bound (11) for the first three combs.  Under the
six-comb cover, every compact component of `E_3` lies strictly inside a
three-comb open union component of span less than `1/d_4`.  Thus
`|E_3|<C_3/d_4`, and (44e) proves the next recurrence.

> **Corollary 5 (third tail rung).**  Every six-comb open cover of a complete
> `c`-slow gap obeys
>
> ```text
> d_4/c<(49/15)C_3
>  <=(49/15)[1+sum_(i=1)^3 ceil((6d_i+c)/(7c))],      (44f)
> d_4/d_3<343/15.                                     (44g)
> ```

For the ratio claim put `x=d_3/c>1`; every ceiling in (44f) is less than
`x+1<2x`, so

```text
d_4/d_3<(49/15)[1+3ceil(x)]/x
             <(49/15)(1+6)=343/15.
```

The recurrence stops naturally here.  Repeating the coarse replacement
`w_i<w_3` one rung earlier would give a seven-width attachment against only
a six-width `d_3` gap, so it would no longer forbid bridging.  Any fourth
tail rung must use phase, the already proved ratio bounds, or finer owner
chronology rather than the same black-box span estimate.

## 5. Toothpick self-similarity and the limit of quotienting

Let

```text
g=gcd(c,d_1,...,d_r).                                  (46)
```

Under the change of variable `u=gt`, (3) descends exactly to the speeds

```text
c/g<d_1/g<...<d_r/g                                   (47)
```

on the same gap address `k`.  Prefix pressure, the tooth counts `n_i`, and
`C` are unchanged, while every length and the Farey mesh scale by `g`:

```text
M_k(c;d_1,...,d_r)
 =(1/g)M_k(c/g;d_1/g,...,d_r/g).                      (48)
```

Thus (19)--(20) are covariant under the ordinary gcd descent.  A primitive
packet has `g=1`, and no further quotient of the ratios is legitimate in
general: THM-1176's normalized clock still contains the phases
`k(d_i mod c)/c`.  The Farey mesh is the correct smaller object because it
quotients each beat denominator by its own arithmetic while retaining the
actual points in the addressed gap.

## 6. Carrier and tournament audit

The speed-order tournament is again transitive and sees none of the proof.
If one insists on a pairwise observable, orient the prefix labels by

```text
O(i,j)=(d_r+d_i)-(d_r+d_j)=d_i-d_j,                   (49)
```

with increasing label as the tie gauge.  Its score histogram is
`0,1,...,r-2`, it has no directed cycle, singleton SCCs, and one Hamiltonian
path.  It preserves only denominator order.  It destroys the numerators
`p`, the gap address `k`, coincidence between sum/difference lattices, and
the mesh cells in (12).

The challenged vertex sets were runners, danger arcs, slow gaps, residues,
endpoint events, beat denominators, beat points, and proof obligations.  The
faithful vertices are the **phase-local Farey mesh cells**, with beat points
as separating walls and prefix-survivor components as occupants.  Pair
labels are sidecars recording why a wall is forbidden.  A bare tournament
cannot preserve the predicate `E intersect B=empty`.

## 7. Exact replay and honest frontier

The companion script independently checks all rational identities, the
component count, the mesh bound, common-gcd covariance, and the implication
"a beat in the prefix survivor is an exact lonely witness."  As finite
telemetry, it scans every phase-labelled packet

```text
2<=c<=7,   c<d_1<...<d_r<=3c,   r in {4,5,6}.          (50)
```

All `11,354`, `20,268`, and `27,730` respective rows contain such a beat
witness; there are no survivors of the beat screen in this bounded bank.
Those zero counts are not used to prove the all-range theorem.

Normal and optimized replays are byte-identical.  Frozen SHA-256 hashes are

```text
source  adfa80ed69f51f9980f9513f01f3d45f68d5d5e327c93ef1297f67fa48d4cefa
output  7d3bad66711432de9b2dbc5ac25d299995b9604172dd7109d71d3a31ffe8b0cb
```

What is proved is the exact exclusion (8), the phase-local inequality (19),
the subcritical bounds (20)--(23), and the unconditional six-comb bounds
(28b)--(28c), (37d)--(37e), and (44f)--(44g).  THM-1198 independently makes every
active-cardinality branch `r<=5` empty; THM-1196 shows that the remaining
`r=6` branch cannot escape through arbitrarily large fourth-to-third,
fifth-to-fourth, or sixth-to-fifth ratios.  It does not bound the first three
prefix ratios,
extend the finite-carrier BV certificates to all `c`, prove universal
six-comb noncoverage, or prove LRC(14).  The remaining escape is now a
phase-aware prefix-ratio/phase-orbit problem rather than an undifferentiated
six-ratio tail.
