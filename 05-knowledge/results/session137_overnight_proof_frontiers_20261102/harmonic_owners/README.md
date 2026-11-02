# First-missing-speed chambers and complete harmonic boundary certificates

Session 137, 2 November 2026. **PROVED** for the specified second-band
family; **FINITE-EXACT** tests supplement the proofs. **LRC(14) remains OPEN.**
No claim of literature priority or formal-kernel verification is made.

## Useful results

The inherited owner atlas extends substantially farther than its uniform
`1/(q+2)` floor when more initial second-band speeds are present. The missing
coordinate is the **first absent second-band speed**, not just the unit-owner
spectrum. At the new floor the complete missing stratum is an explicit
primitive rational grid, including every isolated point and all tied owners.

This gives an all-family harmonic continuation, complete maximizer census
at the formerly excluded boundary, and a precise test distinguishing an
attained boundary maximum from a maximum strictly below it. Neither measure
nor average is substituted for an actual witness.

## 1. Inheritance, novelty boundary, and concept board

Read-only source: session137 checkout `bb7cb2d1f61a14632389814469317131d897b640`; no git mutation was done.
The closest mechanisms are:

- THM-819, `01-canon/theorems/THM-819-primitive-harmonic-law-for-interval-cores.md`:
  primitive harmonic **measure**, with its corrected null-boundary caveat.
- THM-826, `01-canon/theorems/THM-826-farey-profile-theorem-interval-core-measure.md`:
  exact Farey tents/remnants; this already gives an algorithmic interval-core
  profile and must not be presented as new.
- Session114, `session114_continuation_symmetry_20261030/lrc/GENERAL_CORE.md`:
  complete top windows for `{1,...,q-1,q+1}` above `1/(q+2)`.
- Session116, `session116_dimension_owner_towers_20261030/lrc/SECOND_BAND.md`,
  and its credited concurrent `odd_eighth_owner_families_20261030/PROOFS.md`:
  the second-band unit-owner promotion and harmonic transfer above that
  floor are inherited. Its O5 also already proves the **scalar first-tail**
  law `M(C union {q})=1/(q+d)`, with primitive-grid maximizers. That
  `m=1,k=1` special case below is recovered, not newly claimed. Session116
  explicitly excludes the core-atlas floor and gives a counterexample there.
- Session135, `session135_fractran_tm_structure_20261101/lrc_machine/README.md`:
  the original interval-core formula retains isolated primitive-grid points.
  Its `A=empty,d=1` case is recovered here, not claimed again.

The new statements are the full core first-missing-speed chamber extension, its boundary
and all-owner classification, and the resulting deeper harmonic ceilings
and boundary census. A targeted source search found no identical statement;
that is not a literature novelty claim.

The live concepts are: missing-speed Dirichlet confinement; promoted unit
owners; primitive equality grids; harmonic first arrival; and rational
positive/negative certificates. The hostile is the nonunit attachment that
changes a boundary while leaving every old unit owner unchanged. The old
`q<=25` shortcut remains refuted by THM-762/764; all denominators here are
family/row dependent.

## 2. The universal closed-chamber theorem

Let `q>=3`, with **any** subset `A subset {1,...,q-1}`. Define

```
C = {1,...,q-1} union {q+b : b in A},
d = min({1,...,q-1} minus A), or d=q if A is full,
n = q+d,
U = {u:1<=u<q, gcd(u,q)=1},
e(u)=u+q*1_A(u),   E={e(u):u in U},   s=min(E).
```

In the requested family `1 in A`, hence `d>=2`. Allowing arbitrary `A`
exposes the interval-core special case without changing the proof.

For `r0=1/n <= r <= 1/q`, put `epsilon=1/q-r`. For each primitive `a/q`,
write `p=a^(-1) mod q` in `{1,...,q-1}` and set

```
W_a(r) = [a/q-epsilon/e(p), a/q+epsilon/e(q-p)],
W(r) = disjoint union_(a in U) W_a(r),
P_n = {b/n : 1<=b<n, gcd(b,n)=1}.
```

Then, as equalities of **closed sets**,

```
G_C(r) = W(r)                   for r0 < r <= 1/q,
G_C(r0) = W(r0) union P_n.                                  (1)
```

The core maximum is exactly `1/q`: the q-point circle bound for its AP
subcore gives the upper bound, and every primitive q-grid point survives
all attachments (none is divisible by q).

The intervals are nondegenerate when `r<1/q` and become singleton centers
at the top. At the lower boundary there are exactly

```
phi(q) nondegenerate components,
phi(n)-2*1[gcd(q,d)=1] isolated components.                  (2)
```

The only intersections `W(r0) intersect P_n` are two reflected endpoints
when `gcd(q,d)=1`; otherwise the intersection is empty. This is a sufficient
universal chamber, not a claim that every individual core changes immediately
below it.

### Proof A: missing-speed confinement and wall audit

Every speed in `{1,...,n-1}` except `q` belongs to `C`. Place the `n` indexed
circle points `0,t,...,(n-1)t` in cyclic order, retaining multiplicity. If
one gap is smaller than `r`, its endpoint indices differ by `j` with
`1<=j<=n-1` and `||jt||<r`. Safety forces `j=q`. Consequently

```
t=a/q+z,       |z|<r/q.
```

The center must be primitive: if `g=gcd(a,q)>1`, the present speed `q/g`
has distance `|qz|/g<r`, a contradiction. If `r>1/n`, a smaller gap exists
because the average gap is `1/n`. If `r=1/n` and no smaller gap exists,
every gap is `1/n`; hence the orbit is the full `n`-grid and `t=b/n` with
`gcd(b,n)=1`. This isolates exactly the equality alternative, without
throwing away any finite null set.

Consider the left side `t=a/q-z`, `z>=0`, in the strict-confinement case.
Its residue-1 speeds are `p` and, if present, `q+p`; the largest is `e=e(p)`.
The unreduced residue is `1/q-ez`. It cannot wrap through zero and become
safe on the other side, since

```
ez < e*r/q < 2r <= 1/q+r.
```

Thus safety necessarily gives `z<=epsilon/e`. Reflection gives the right
constraint. Conversely take any `0<=z<=epsilon/e`, whether or not the
confinement inequality is strict. For a speed `v in C`, let
`c=av mod q`, represented in `{1,...,q-1}`. Its raw fractional phase is
`c/q-vz`. The upper bound `<=1-r` is automatic. If `c=1`, then `v<=e`
and the lower bound is exactly the displayed owner inequality.

If `c>=2` and `p` is unpromoted, then `p>=d` and
`v<=q+cp`. Therefore

```
(v-p)*epsilon <= p*(c-1)/q,
```

because `epsilon<=d/[q(q+d)]` and
`d*(q+(c-1)p) <= (c-1)p*(q+d)` is equivalent to
`d<=p(c-1)`. This proves `c/q-v*epsilon/p>=r`.

If `p` is promoted, then `e=q+p>q`, `v<2q`, and
`epsilon<=1/(2q)`. Consequently `v*epsilon/e<1/q`, so for `c>=2`
the raw phase is strictly greater than `1/q>=r`. No hidden wrap occurs:
all phases have been bounded inside `[r,1-r]`.

Every primitive `n`-grid point is safe: `n` is absent from `C`, every
speed is positive and less than `2q<=2n`, and no speed can be a multiple
of `n`. Thus all its distances are at least `1/n`. This proves (1).

### Proof B for the boundary, by an equality-limit argument

The proved open chamber implies `W(r0)` is safe by continuity. Any boundary
point not in `W(r0)` cannot be approached by safe points at larger radii.
The circle-gap argument above then rules out a strict gap (which would put
it in the audited owner windows). All gaps must equal `1/n`, producing
exactly the primitive equality grid. This is an independent way to organize
the last step, not a second independent proof of the underlying wall audit.

### Disjointness and exact component counts

Since every `e>=d`, every side width at the floor is at most `1/(qn)`.
Distinct `q`-centers are separated by at least `1/q`, and
`2/(qn)<1/q` because `n>=4`; the windows are pairwise disjoint.
For `b/n` and `a/q`, a nonzero distance is at least `1/(qn)`. A grid point
can therefore enter a window only at equality: its directional owner must
be the unpromoted `d`, requiring `gcd(q,d)=1`, and `|qb-an|=1`.
There are precisely two such reflected incidences. The lower endpoint for
`a=d^(-1) mod q` is

```
a/q - 1/(qn) = ((an-1)/q)/n.
```

This is primitive since `qb=an-1`. Removing these two already-contained
points, or none when `d` is nonunit, proves (2).

## 3. All endpoint owners, not just one winning speed

Signs mean `{vt}=r` for `+` and `{vt}=1-r` for `-`.

- For `r0<r<1/q`, a window's left endpoint has exactly owner `e(p),+`;
  its right endpoint has exactly owner `e(q-p),-`.
- At `r=r0`, the same statement holds **except** on a side whose unit
  index is `d`, if `2d in A`. That side has the two same-sign owners
  `d` and `q+2d`. Necessarily `2d<q`. There are no other ties.
- At `r=1/q`, the singleton `a/q` has `p,+` and `q-p,-`, together with
  `q+p,+` if `p in A` and `2q-p,-` if `q-p in A`.
- At **every** primitive boundary point `b/n`, let `u=b^(-1) mod n`.
  Its complete owner sets are

```
positive: C intersect {u,u+n},
negative: C intersect {n-u,2n-u}.                           (3)
```

These formulas include overlap points, isolated points, and all ties.
For the exceptional endpoint tie, equality in the unpromoted wall audit
requires simultaneously `epsilon=d/[qn]`, `d=p(c-1)`, and `v=q+cp`.
As `p>=d,c>=2`, this forces `p=d,c=2,v=q+2d`. The promoted case is strict.
Formula (3) is the exact residue equality, with `max C<2n` bounding all
possible lifts. At each isolated point both signs occur; a missing base
owner in `{1,...,n-1}` could only be `q`, which is exactly one of the two
nonisolated overlap cases.

## 4. Harmonic continuation and the boundary decision

Append the actual harmonic bank

```
V = C union {qm,2qm,...,kqm},     m,k positive integers.
```

The inherited family is recovered when `d=2`. For the new whole range
`1<=k<=n-1`, the following is exact:

1. If `dm>=s`, then

```
M(V) = min{m/(qm+s), 1/(k+1)}.                              (4)
```

2. If `dm<s`, put `N=n/gcd(qm,n)`. Then

```
M(V)=1/n   iff k<N;
M(V)<1/n   iff k>=N.                                       (5)
```

The below-floor value in the second line is **not** classified here.
For `dm>s,k<=n-2`, (4) lies strictly inside the new chamber. If its first
bound is strictly smaller, the only two maximizers are

```
t=(am+h)/(qm+s), 1-t,
a=s^(-1) mod q, h=(as-1)/q.
```

Equality of the two bounds is handled by the complete ceiling census below,
not discarded.

**Proof.** The `(k+1)`-point circle-gap argument bounds the k-term tail AP
by `1/(k+1)`; at equality all gaps are equal and its phase is a primitive
`(k+1)`-grid point. The core bound `1/q` excludes larger radii before the
chamber argument is applied. On a side with owner `e`, write
`y=q|t-a/q|`. Core safety means
`qr+ey<=1`. The first tail satisfies `r<=||my||<=my`. Hence any safe point
in the open chamber forces `r<=m/(qm+e)<=m/(qm+s)`. The tail AP separately
bounds clearance by `1/(k+1)`. Conversely choose the side of owner `s` and
`y=r/m` at the minimum in (4). Every tail phase `jmy=jr` is in `[r,1-r]`
when `r<=1/(k+1)`. Since `dm>=s`, the proposed radius is at least `1/n`,
where (1) proves the core safe. This supplies an actual witness.

If `dm<s`, the same argument excludes every radius above `1/n`; at the
floor it excludes every point of `W(r0)` as well. On a primitive `n`-grid
point the first tail has order `N=n/gcd(qm,n)`. Its first `k` multiples
are all nonzero iff `k<N`; their distances are then at least `1/N>=1/n`.
Thus exactly all or none of `P_n` survive. Compactness of the full safe
clearance function makes absence at the floor a strict maximum inequality.
The strict-bottleneck uniqueness uses equality in `e>=s` and `my>=r`;
`E` has one least element, producing one side and its reflection.

## 5. Every harmonic ceiling, including its isolated boundary stratum

Choose `1<=h<=d`, put `n_h=q+h`, and append `k=n_h-1` harmonic tails.
Their optimum `1/n_h` is retained **if and only if**

```
mh>=s.                                                     (6)
```

For each `e in E`, let `a_e=e^(-1) mod q`. All q-window maximizers are

```
t_(e,j)=a_e/q-j/[q*m*n_h],  1-t_(e,j),
1<=j<=floor(mh/e), gcd(j,n_h)=1.                            (7)
```

These points have no duplicates. If `h<d`, (7) is the complete list.
If `h=d`, also include the **entire** primitive `n`-grid precisely when
`gcd(qm,n)=1`, taking a set union rather than double-counting its two
q-window endpoints. Consequently the exact boundary count is

```
2 sum_(e in E) #{1<=j<=floor(md/e):gcd(j,n)=1}
 + (phi(n)-2)*1[gcd(qm,n)=1].                               (8)
```

For (6) at the boundary, a grid-only exception cannot violate the criterion:
`gcd(qm,n)=1` forces `gcd(q,d)=1`, whence `s=d` and `md>=s` automatically.
In the nonpreserving range both sets are empty. Formula (7) follows by
forcing the tight tail AP phase to be primitive `j/n_h` and imposing its
literal owner inequality `ej<=mh`. This proves necessity, sufficiency,
completeness, and the equality cases of (4).

### Terminal mass, with the boundary correction retained

For `n_h=q+h`, define

```
w_(n_h)(j)=1/ell_j+1/(n_h-ell_j),  ell_j=j^(-1) mod n_h,
H_(n_h)=sum_(1<=u<n_h,gcd(u,n_h)=1) 1/u.
```

In the preserving range, let `B` be the left-hand terminal mass coefficient:
`measure G_V(1/n_h-eta)=B*eta` for all sufficiently small positive `eta`.
The exact formula is

```
B = 2/(qm) sum_(e in E)
       sum_(1<=j<=floor(mh/e),gcd(j,n_h)=1) w_(n_h)(j)

  + 1[h=d and gcd(qm,n_h)=1] * 2/(qm)
       * (H_(n_h)-w_(n_h)(m)).                             (9)
```

The second line is included only under its stated condition, so `w_(n_h)(m)`
is evaluated only when `m` is coprime to `n_h`. It is the equality-grid sidecar. It is not present in the
inherited open-chamber formula. It sums all primitive-grid weights and
subtracts the two overlap points already counted in the first line.

**Fastest-owner proof, including ties.** At any tail-AP maximum, the two
tight tail owners are `qm*ell` and `qm*(n_h-ell)` for the appropriate
primitive phase. A tight tail speed is at least `2q`, hence faster than
every core speed, except possibly when `m=1` and that tail speed is `q`.
In this exception `qt=z+sigma/n_h`, with `sigma=+1` or `-1`. A core speed
with the same signed equality would satisfy

```
v*(n_h*z+sigma) = q*(n_h*z'+sigma),
```

so `v congruent q mod n_h`. Since `v<2q<n_h+q` and `q` is absent from
`C`, this is impossible. Thus the tail always supplies the fastest tight
owner in each direction, even at the new core tie `d,q+2d`.
The local interval grows at the sum of the reciprocals of those two speeds.
Summing (7) gives the first line of (9). Multiplication by the unit `qm`
permutes primitive residues when the grid is present, so its total weight
is `2H_(n_h)/(qm)`. Its two q-window overlap points have phases `+/-m/n_h`
and combined weight `2w_(n_h)(m)/(qm)`. This proves the correction.

For a fixed `h`, put `L=n_h*lcm(E)` and `H_E=sum_(e in E)1/e`. Once in the
preserving range, the exact count `N_h` and coefficient `B_h` satisfy

```
N_h(m+L)-N_h(m) = 2Lh*phi(n_h)*H_E/n_h,
(m+L)B_h(m+L)-m B_h(m) = 4Lh*H_E*H_(n_h)/(q*n_h).
```

Indeed each owner cutoff increases by a multiple of `n_h`, and the boundary
correction is periodic in `m` modulo `n_h`. Therefore

```
N_h(m)/m -> 2h*phi(n_h)*H_E/n_h,
B_h(m) -> 4h*H_E*H_(n_h)/(q*n_h).
```

The extra isolated maxima are invisible in these leading limits. This is
why the asymptotic or measure-only statement cannot replace the exact set.

## 6. Hostiles and the quotient-loss ledger

1. **Deleting isolated points loses maximizers.** For `q=3,A={1}`,
   `d=2,n=5`, the lower boundary contains isolated `1/5,4/5` in addition
   to its two positive q-windows. With `m=1,k=4`, all four primitive fifths
   are maxima, but a q-window-only census retains just `2/5,3/5`.
2. **The unit-owner spectrum alone is not a boundary invariant.** For
   `q=4`, the cores `A={1}` and `A={1,2}` both have `E={5,3}`. At `r=1/6`,
   `t=1/6` is safe for the first and killed by speed6 in the second.
   Their first-missing indices are2 and3. Nonunit promotions are invisible
   in the old open chamber but not at adjacent-band entry.
3. **Boundary owners cannot be arbitrarily chosen.** For
   `q=5,A={1,4}`, at `r=1/7`, the overlap endpoints `3/7,4/7` have tight
   owners2 and9 with the same sign. The faster owner is required when
   continuing below the boundary; a unique-owner continuation is false.
4. **Even `(q,d,E)` does not preserve boundary ownership or continuation.**
   Within the requested `1 in A` family take `q=8`, `A={1,2}` and
   `A'={1,2,6}`. Both have `d=3,E={9,3,5,7}` and the same entire closed
   atlas down to `1/11`. At `t=4/11`, the first has sole positive owner3,
   while the second has tied positive owners3 and14. For
   `eta=1/10000`, the point `t-eta/3=119989/330000` is safe for the first
   at `r=1/11-eta=9989/110000` but unsafe for the second. Thus the
   actual attachment set, or the boundary co-owner sidecar, is necessary
   for continuation below the chamber. `check_boundary_loss.py` verifies
   both integer-residue verdicts; its normal/optimized receipts agree.
5. **No universal lower-floor claim.** For `q=3,A={1,2},m=1,k=3`, the exact
   maximum is `1/7`. The new floor is `1/6`; here `d=3,s=4,N=2,k>=N`,
   so (5) correctly predicts a strict drop. This includes the inherited
   hostile below the old `1/5` floor.

```
source: labelled physical speeds C, exact phase t and threshold r
target: (q,d,E), oriented q-windows, primitive equality grid, owner lifts
preserved: full closed safety throughout the specified chamber; all maxima
lost by E alone: nonunit attachments, first missing speed, equality components
restoration: d and, for tied owners, the actual attachment set A
cheapest decisive tests: q4 same-E cores; q3 isolated fifths; q5 owner2/9 tie
```

There is no forced tournament interpretation: the relevant carrier is a
labelled closed interval complex with a primitive-grid equality stratum.

## 7. Exact finite certificates and reproduction

`check.py` is dependency-free and uses `Fraction`, integer residues, and
explicit exceptions rather than optimization-elided assertions. It compares
(1) with direct closed safe-interval intersection, independently checks a
smaller universe by arrangement walls and open-cell midpoints, and verifies
all owner sets, harmonic ceiling maximizers and terminal mass coefficients.
A smaller bank independently measures the exact safe-set growth below the
maximum and checks (9). Equality intersections are
retained throughout.

The maximum-value checks also use a separate generic certificate fact:
for any finite positive integer row with maximum speed `H`, a maximum
clearance has reduced denominator at most `2H`. Indeed at a local maximum
of the lower envelope of triangular distance functions, two active slopes
have opposite signs (possibly the two sides of a single half-integer peak).
Solving their equalities gives `r=integer/(v+w)` with `v+w<=2H`.
Thus for a proposed `r` of denominator at most `2H`, a safe witness at `r`
and an exact empty safe-set check at `r+1/(8H^2)` prove it is the exact
maximum: two distinct rationals with denominators at most `2H` differ by
at least `1/(4H^2)`. This supplies positive and negative rational
certificates, not an appeal to a numerical optimizer or a sample average.
The negative verifier is magnitude-dependent interval intersection; no
polynomial binary-time or uniform small-denominator-bank claim is made.
It does not reduce arbitrary unbounded LRC to the fixed chamber family.

```
python3 -B check.py > normal.json
python3 -O -B check.py > optimized.json
cmp normal.json optimized.json
python3 -B check_periods.py > periods.normal.json
python3 -O -B check_periods.py > periods.optimized.json
cmp periods.normal.json periods.optimized.json
python3 -B check_boundary_loss.py > boundary_loss.normal.json
python3 -O -B check_boundary_loss.py > boundary_loss.optimized.json
cmp boundary_loss.normal.json boundary_loss.optimized.json
sha256sum -c MANIFEST.sha256
```

The primary receipt records 4,092 complete attachment sets for q=3..12,
16,368 closed-atlas tests, 184,512 endpoint-owner checks, 248 independent
wall-oracle tests, 2,226 larger prefix-family tests through q=40,
2,649 exact-or-strict-below harmonic decisions, 729 complete ceiling lists
and terminal mass checks, and 108 independent exact mass-slope checks.
The separate no-producer-import residue-block checker adds 24,522 affine-period
tests, including 12,276 boundary-correction cases, across all 4,092 attachment
sets q=3..12 and every `1<=h<=d`. Both programs have byte-identical normal
and optimized receipts.
The scope of each count is also stated in `normal.json`. The independent
agent's separate audit uses no producer imports and is packaged by the
parent session; its first report passed all attachment sets q=3..10 at
five rational radii (5,100 cases), plus 1,020 boundary classifications and
the two named isolated/tie controls. Its extended independent run also
passed 447 harmonic maximum/boundary decisions, 464 complete ceiling sets,
292 exact terminal-mass slopes, and 428 extra primitive-grid points. These
are separate audit universes, not added to the primary counts.
Its final mathematical review accepted
the chamber, all owners, harmonic maximum dichotomy, complete ceilings, mass
correction, periods and exact-maximum certificates after the factor repair
recorded below. This is a finite audit, not the source
of the universal quantifiers. No external literature audit is claimed.

### Draft correction retained

The independent proof audit caught one missing factor `p` in the displayed
unpromoted-wall inequality. The correct statement is
`(v-p)*epsilon <= p*(c-1)/q`, as already supplied by the next algebraic
line. The earlier display without `p` fails for `q=5,p=d=2,v=9,r=1/7`.
It was repaired before publication; neither the theorem nor the exact
checker changed. The same-sign owner2/9 control makes the mistake decisive.
