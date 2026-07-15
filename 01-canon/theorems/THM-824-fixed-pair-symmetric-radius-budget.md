---
id: THM-824
title: Fixed-ratio folded target and the exact symmetric Minkowski radius budget
status: PROVED (all-size compact-set factorization and common-dilate corollary) + VERIFIED (exact affine, synthetic-packet, and 214-core replays)
source: codex-2026-07-15-S10 deep sum-arc continuation
depends_on:
  - THM-774  # folded-diamond eligibility/colour equivalence
  - THM-789  # symmetric return erosion
  - THM-821  # atomic signed-cell factorization and finite bank
related: [THM-803, THM-807, THM-817, HYP-6820]
verification:
  - 04-computation/lrc14_fixed_pair_symmetric_radius_budget_codex_S10.py
  - 05-knowledge/results/lrc14_fixed_pair_symmetric_radius_budget_codex_S10.out
---

# THM-824 — The fixed-ratio symmetric radius budget

Write `T=R/Z`, use `||t||` for circular distance to zero, and put

```text
Q(t)=||9t||+||4t||,                 theta=11/13,
C={c_-,c_+}={5/13,8/13},            h=2/169,
d_C(t)=min_(c in C)||t-c||.                                  (1)
```

For nonempty compact `E subset T` and compact `R subset T`, define

```text
rho_C(E)=max_(t in E)d_C(t),
rho_0(R)=max_(r in R)||r||.                                  (2)
```

The fixed folded target is only two short cells, and symmetry turns its full
Minkowski-containment predicate into one additive radius budget.

## 1. Exact identification of the folded target

### Theorem 1 — two cells in the `13`-cover

One has exactly

```text
H={t:Q(t)>=11/13}
 =B_T(5/13,2/169) union B_T(8/13,2/169)
 =[63/169,67/169] union [102/169,106/169].                (3)
```

The two displayed intervals are ordinary lifts in `[0,1]`; they are disjoint.

### Proof

THM-774 identifies `Q(t)>=11/13` with the conjunction

```text
||13t||<=2/13,        ||5t||<=2/13,                       (4)
```

and opposite parities for the unique nearest integers `p` to `13t` and `q`
to `5t`.  The nearest integers are unique because `2/13<1/2`.  The two
eligible teeth about `p/13` and `q/5` can meet only if

```text
|p/13-q/5| <= 2/169+2/65 = 36/845,
```

or

```text
|5p-13q| <= 36/13 < 3.                                   (5)
```

The integer in (5) is odd because

```text
5p-13q = p-q mod 2,
```

so it is `+/-1`.  Modulo `13`, using `5^(-1)=8`, this forces

```text
(p,q)=(5,2) with 5p-13q=-1,
(p,q)=(8,3) with 5p-13q=+1                               (6)
```

in one period.  Conversely both pairs have opposite parity.  Their centre
separation is `1/65=13/845`.  The farthest point of the radius-`2/169`
`13`-tooth is therefore only

```text
13/845+2/169=23/845 < 26/845=2/65
```

from the corresponding `5`-tooth centre.  Thus each intersection in (6) is
the whole smaller `13`-tooth.  These are precisely the intervals in (3). ∎

## 2. The no-switch lemma

Let

```text
I_-=B_T(c_-,h),              I_+=B_T(c_+,h).              (7)
```

### Lemma 2 — a symmetric translate cannot change target cell

If

```text
t, t+r, t-r in H,                                             (8)
```

then all three points lie in the same one of `I_-` and `I_+`.  If `c` is
that cell's centre, then

```text
d_C(t)+||r||<=h.                                             (9)
```

### Proof

The cells are disjoint because their centre distance is

```text
||c_+-c_-||=3/13=39/169>2h=4/169.                           (10)
```

Therefore the three points in (8) have unique cell labels `i,j,k` in
`{-,+}`.  The circular triangle inequality gives

```text
||2c_i-c_j-c_k||
 <=2||t-c_i||+||t+r-c_j||+||t-r-c_k||
 <=4h=8/169.                                                (11)
```

Writing the centres as `5/13,8/13`, the complete label table for the left
side of (11) is

```text
j=k=i:       0,
one switched label: 3/13=39/169,
both switched:      6/13=78/169.                            (12)
```

Only `j=k=i` is compatible with (11).

It remains to make the lift in (9) explicit.  Choose the unique real lifts
of `t,t+r,t-r` in `[c-h,c+h]`, with `c` represented by `5/13` or `8/13`.
The sum of the latter two lifts and twice the first lift are congruent modulo
one.  Their difference has magnitude at most `4h=8/169<1`, so the lifts obey
the exact real identity

```text
tilde(t+r)+tilde(t-r)=2 tilde(t).                            (13)
```

Consequently

```text
s=tilde(t)-c,
u=tilde(t+r)-tilde(t)=tilde(t)-tilde(t-r)                  (14)
```

are real representatives with `|s+u|,|s-u|<=h`.  Also
`|u|<=2h<1/2`, so `u` is the unique balanced representative of `r` and
`|u|=||r||`; similarly `|s|=d_C(t)`.  The elementary identity

```text
max(|s+u|,|s-u|)=|s|+|u|                                  (15)
```

now proves (9).  This argument includes all endpoint cases and contains no
implicit choice across the circle cut. ∎

## 3. Exact all-size factorization

### Theorem 3 — symmetric Minkowski radius budget

Let `E subset T` be nonempty and compact.  Let `R subset T` be compact with

```text
R=-R,                    0 in R.                            (16)
```

Then

```text
E+R subset H
iff rho_C(E)+rho_0(R)<=2/169.                              (17)
```

Thus the **global conjunction** of all component-by-cell stalk predicates,
equivalently the containment for their full Minkowski union `E+R`, is
equivalent at every height to two extremal radii.  An individual signed
satellite `R_k` is not symmetric and need not satisfy this factorization by
itself.

### Proof

Suppose first that `E+R subset H`.  Since `0 in R`, every `t in E` lies in
`H`.  For each `r in R`, symmetry also puts `-r in R`, so Lemma 2 applies to
`t,t+r,t-r` and gives

```text
d_C(t)+||r||<=h                                             (18)
```

for every pair `(t,r) in E x R`.  Compactness makes both maxima in (2)
attained.  Substituting a maximizing `t` and a maximizing `r` in (18) proves
the forward inequality in (17).

Conversely suppose the radius inequality holds.  For `t in E`, choose a
centre `c in C` attaining `d_C(t)`.  For every `r in R`, the circular
triangle inequality gives

```text
||t+r-c|| <= ||t-c||+||r||
            <=rho_C(E)+rho_0(R)<=h.                        (19)
```

Hence `t+r in B_T(c,h) subset H`, proving the reverse implication.  Notice
that this direction does not require symmetry; symmetry and the zero element
are what make the converse necessity true. ∎

### Monotonicity

The signed budget

```text
kappa(E,R)=2/169-rho_C(E)-rho_0(R)                         (20)
```

is decreasing under either inclusion `E subset E'` or `R subset R'`.
It is therefore a genuine monotone obstruction: once negative, no refinement
that enlarges either exact set can repair it.  This is different from ranking
individual stalk margins, whose signed target cells can move against one
another.

The assumptions in (16) are essential for the reverse implication.  For
example,

```text
E={5/13},             R={0,3/13}
```

has `E+R={5/13,8/13} subset H`, but
`rho_C(E)+rho_0(R)=3/13>2/169`; the one-sided `R` is not symmetric.

## 4. The common-dilate branch

For a positive integer `d`, let `m_d(t)=dt` on `T` and put

```text
H_d={t:||9dt||+||4dt||>=11/13}=m_d^(-1)(H).               (21)
```

In the two-sheet application this is the common-ratio exception pair

```text
(x,y)=(13d,5d),                                           (22)
```

with `d` odd.  Define

```text
rho_(C,d)(E)=max_(t in E) min_(c in C)||dt-c||,
rho_d(R)=max_(r in R)||dr||.                               (23)
```

### Corollary 4 — scale-uniform fixed-ratio budget

Under the hypotheses of Theorem 3,

```text
E+R subset H_d
iff rho_(C,d)(E)+rho_d(R)<=2/169.                          (24)
```

Indeed multiplication by `d` is a circle homomorphism, so

```text
m_d(E+R)=m_d(E)+m_d(R),                                   (25)
```

and the image of `R` remains compact, symmetric, and contains zero.  Apply
Theorem 3 to the two images.

When (24) holds, every `dt` lies in one of the two cells in (3), and hence

```text
d_C(dt)=||13dt||/13.
```

Thus (24) has the exact phase form

```text
max_(t in E)||13dt|| + 13 max_(r in R)||dr|| <=2/13.       (26)
```

If these are the THM-789 deep and return sets for a hypothetical two-sheet
packet, its individual exception cap gives `13d<=11B`, where `B=max(U)`.
The closed central return interval then gives

```text
rho_d(closure(R_U))>=2d/(143B)<1/2.
```

Substitution in (26) recovers exactly THM-789's `w=13d` thickness tax

```text
max_(t in E_U)||13dt||
 <=2/13-26d/(143B)=2(11B-13d)/(143B).                    (27)
```

For the fixed branch `d=1`, take the THM-789 sets

```text
E=E_U={t:min_(u in U)||ut||>=1/11},
R=closure(R_U),
R_U={r:max_(u in U)||ur||<2/143}.                          (27a)
```

Settled `LRC(<=13)` makes `E_U` nonempty, and `closure(R_U)` is compact,
symmetric, and contains zero.  Therefore the entire THM-803/821 containment
is exactly

```text
max_(t in E_U)||13t||
 +13 max_(r in closure(R_U))||r|| <=2/13.                  (28)
```

The mandatory central return cell has radius `2/(143B)`, where `B=max(U)`.
Replacing the second maximum in (28) by this lower bound recovers exactly
THM-789's pointwise `w=13` thickness tax:

```text
max_(t in E_U)||13t||
 <=2/13-26/(143B)=2(11B-13)/(143B).                        (29)
```

Every return satellite farther from zero strictly strengthens (29).  This is
the promised scale-free use of THM-817's signed cells: their full numerical
effect is their outer circular radius, even when their count is linear in
`B`.

## 5. Linear evaluation and endpoint ancestry

Once exact component intervals are known, (17) is evaluable in linear time
in their number.  On `[0,1]`, the function `d_C` is affine between

```text
0, 5/13, 1/2, 8/13, 1.                                   (30)
```

Hence `rho_C(E)` is attained at an input endpoint or at a contained Voronoi
boundary `0` or `1/2`.  Likewise `||r||` is maximized on a return interval at
an endpoint or a contained half-integer.  For THM-817's primitive return
cells the self-antipodal cell is absent, so the cell endpoints suffice.

This reduces evaluation after component construction from the
`c_E*N_R` atomic sum arcs of THM-821 to `O(c_E+N_R)` exact comparisons.  It
does **not** make endpoint owners disposable under recursion.  For one fixed
pair `(E,R)`, the two numbers in (2) decide containment.  Under speed
replacement, gcd descent, common-dilate transport, or erosion iteration,
owners record which inequality moves the maximizing endpoint.  The
theorem-bearing dynamic carrier remains

```text
(owner-decorated deep intervals, owner-decorated signed return cells)
 -> (rho_C(E),rho_0(R)) -> kappa(E,R).                     (31)
```

The final arrow is evaluation-exact; forgetting owners before the first
arrow is not claimed to be functorial.

## 6. Tournament Analysis and challenged vertices

Tournament vertices are six candidate **proof carriers**, not runners:

```text
topology counts,
signed component/cell labels,
the radius pair,
exact sum arcs,
exact input interval families,
owner-decorated input families.                             (32)
```

The evaluation observable first orders by loss of the fixed predicate and
then by payload size.  The switch/gauge first orders by retained endpoint
ancestry; declaration order is the fixed tie Hamiltonian path.  The replay
gives the two transitive orders

```text
radius pair < exact input intervals < owner-decorated inputs
            < exact sum arcs < signed labels < topology counts,

owner-decorated inputs < exact input intervals < signed labels
            < exact sum arcs < radius pair < topology counts.  (33)
```

They have score histogram `(0,1,2,3,4,5)`, no directed cycle, six singleton
SCCs, one Hamiltonian path, and the recorded edge-flip count in the replay.
The tournament is telemetry: (17) proves evaluation sufficiency, while
THM-821's liar fibres prove that topology and labels alone do not.  The
transport order records retained ancestry, not a theorem that owners by
themselves determine a future recursion.

The challenged assumption is that disconnected return satellites force a
quadratic component-incidence object at this fixed ratio.  Symmetry of the
**full return union** pairs every satellite with its negative and prohibits
target-cell switching in the global conjunction; a lone satellite can switch
cells, exactly as the one-sided counterexample after (20) shows.  Hence only
the full union's outer radius affects the present global predicate.  Alternate
vertices such as sum arcs, cusps, endpoint owners, and proof obligations
remain useful at different stages, but none should be confused with the
minimal fixed-evaluation quotient.

## 7. Verification and scope

The exact replay independently:

1. reconstructs (3) from every affine cell of `Q` and from the parity/Bezout
   tooth intersections;
2. emits the complete no-switch centre table;
3. checks a deterministic exact bank of connected and disconnected compact
   interval packets, including translations near the forbidden `3/13`
   switch;
4. regenerates THM-821's `213` seeded cores and THM-817's disconnected row,
   compares the raw all-stalk predicate with (17) on all `214` rows, and
   independently recovers all `9,974` atomic verdicts; and
5. emits certificate and stored-output SHA-256 digests.

No floating point or sampled-circle verdict is used.

Reproduce with

```bash
python3 04-computation/lrc14_fixed_pair_symmetric_radius_budget_codex_S10.py
```

The frozen source and stored-output SHA-256 digests are

```text
source  d9c77bce064e23c9891012612cc80e3c8e1fbc4e148685d6c28063cc37f841b5,
output  7c2c4b124893f07659918f9da45903cee7fc93acfeacc0c5b2b8ec327a125e6b.
```

The canonical mathematical certificate digest emitted inside the output is

```text
8eac54e4513394b3857199e3e5528a4bd6780b8ff5f53efcbe0c9ba26c06f95c.
```

An independent byte-equal replay used `4.62` seconds of user CPU (about
`19.6` seconds wall time under the shared concurrent load).

The theorem does not prove that every arithmetic core violates (28), does not
close the global `n=12` sporadic branch, and says nothing uniform about odd
pairs outside the common ratio `(13d):(5d)`.  It does not bound the number of
deep or return components and does not replace the stronger tightness target
`G_U subset H_d`; it is an exact simplification of the necessary deep-set
erosion predicate.  The `214`-core bank is a regression check, not an
exhaustion over heights.
