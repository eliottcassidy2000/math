# Cofinal pair planes, owner descent, and multiplicity windows

**Research synthesis, 2026-08-26.** The direct proved advances from this
session are [THM-4231](../01-canon/theorems/THM-4231-arbitrary-pair-cofinal-depth-six-haar-repair-and-exact-outsider-lift.md),
the strengthened cross-scale clause in
[THM-4082](../01-canon/theorems/THM-4082-mahler-renormalized-linear-chart-and-exact-bit-defect.md),
and the multiplicity-window clause in
[THM-4036](../01-canon/theorems/THM-4036-sun-2468-energy-and-support-exponent.md).
[MISTAKE-520](../01-canon/MISTAKES.md) repairs the inherited fixed-pool
coverage ledger. The independently audited planar-Jacobian wall promoted as
[THM-4232](../01-canon/theorems/THM-4232-weight-eleven-u-zero-primitive-cm-planar-jacobian-exclusion.md)
is part of the same synthesis. Incoming [THM-4227](../01-canon/theorems/THM-4227-two-outsider-scale-separated-depth-eight-haar-wedge.md)
and [THM-4228](../01-canon/theorems/THM-4228-common-gcd-two-outsider-periodic-observable-haar-ray.md)
were independently audited before use here. LRC(14), JC(2), Mahler's `3/2`
problem, and the classification of Sun's `2-4-6-8` holes remain **OPEN**.

Reserved THM-4229/4230 files are not proved inputs. Historical drafts
and scratch packets below are evidence only at the status explicitly stated.

## Inheritance pass

- **LRC closest mechanism:** THM-4170's exact component-discrepancy bound.
  **Hostile:** THM-4207 shows that marginal one-outsider decks do not compose.
  **Corrected near miss:** MISTAKE-520 observes that THM-4156's safe full pool
  already closes every fixed-pool subset by heredity; a more elaborate repair
  census is not new body coverage. **Least-used sidecar:** the component count
  of each literal safe union.
- **Planar-Jacobian closest mechanism:** THM-4222's complete lower Newton model
  and primitive `Q(zeta_11)` CM obstruction. **Hostile:** THM-4218's hidden
  `j=0` elliptic tail. **Corrected near miss:** deleting a top owner requires a
  complete replacement-face audit, not reuse of the old polygon. **Least-used
  sidecar:** the first surviving pure-`p` owner and every aggregate coefficient
  collision.
- **Mahler closest mechanism:** THM-4082's exact defect from the limiting
  isometry. **Hostile:** the denominator-19 moving-start programmer.
  **Least-used sidecar:** the same literal point through every scale.
- **Sun closest mechanism:** THM-4036's second moment. **Hostile:** the exact
  hole `896315812331399`. **Least-used sidecar:** representation multiplicity,
  as distinct from support.

## Live concept board

| lane | discrete object | operation | retained invariant | information that must not be quotiented |
|---|---|---|---|---|
| LRC anchor | labelled deletion edge / exceptional body | add two literal combs | exact Haar mass, components, transversal incidence | outsider labels, phase, direct margin |
| JC anchor | complete valued support with owner labels | kill one coefficient and descend to the next face | polygon, components, genus, CM/Hom, degree | omitted affine divisors, collisions, attachments |
| Mahler niche | one point in all renormalized charts | compare scales | exact first differing bit | integer/positive/terminal status |
| Sun niche | multiplicity fibre `a(n)` | truncate low and high fibres | tuple mass and square energy | zero fibres and height-sensitive lifts |
| tournament crossfeed | exact bad-owner mask on a cyclic word | lower/upper Boolean zeta | owner tensor and successor incidence | cyclic adjacency and realizability |

The shared useful idea is not “all five problems are the same.” It is that a
first nonzero coordinate becomes useful only when the sidecar needed by the
next operation remains attached.

## 1. LRC: three complementary cofinal geometries

[THM-4227](../01-canon/theorems/THM-4227-two-outsider-scale-separated-depth-eight-haar-wedge.md)
and THM-4231 solve different parts of the outsider plane.

THM-4227 is sequential. For an ordered pair it first intersects a depth-eight
base repair with `G_q`, retains the new component count, and then intersects
with `G_r`. It proves the directed wedge

```text
q>=3391,
10633545731 r >= 321902813232(q+130),                  (1)
```

or its transpose. Its strength is scale separation: one outsider may be much
smaller than `17548`.

THM-4231 is symmetric. For a six-deletion `R`, it applies Bonferroni directly
inside `U_R=G_(P\R)`:

```text
mu(U_R intersect G_q intersect G_r)
 >=(5/7)M_R-(6c_R/49)(1/q+1/r).                        (2)
```

The exact activation

```text
kappa_2(R)=ceil(108c_R/(7(45M_R-4)))                   (3)
```

produces a `54,566`-edge deck at `Q=17548` with no nine-cover. Two independent
paths agree on all `593,775` deletions and literally scan all
`14,307,150` nine-bodies. This closes every pair of distinct outsiders
`q,r>=17548` and, with THM-4156/4191, every one of the
`C(32,11)=129,024,480` eleven-faces in the resulting chart.

THM-4228 retains arithmetic correlation instead of discarding it. Writing
`(q,r)=g(u,v)`, it treats `G_u intersect G_v` as one periodic observable and
proves the sharp primitive-pair density floor

```text
mu(G_u intersect G_v)>=66/91,                           (3a)
```

with equality only at ratio `1:13` or `13:1`. Component discrepancy at the
single outer frequency `g` then closes every pair with `gcd(q,r)>=3467`.
This reaches arithmetic rays missed by a purely size-based picture. The three
methods retain different coordinates: directed scale order, common divisor
plus primitive shape, or symmetric marginal size.

### Pair-plane reduction

For `m=min(q,r)` and `M=max(q,r)`, THM-4231 handles `m>=17548`. If
`3391<=m<=17547`, THM-4227 handles

```text
M>=ceil(321902813232(m+130)/10633545731).               (4)
```

The right side is increasing and is `535125` at `m=17547`. Therefore every
still-uncovered outsider pair satisfies

```text
m<=3390,
or
3391<=m<=17547 and M<=535124,
and in either case gcd(q,r)<=3466.                       (5)
```

This changes the topology of the open problem. Pair entry is no longer one
undifferentiated infinite quadrant: it is one finite exact region plus a
finite list of fixed-small-label rays. The natural next analytic obligation is
one cofinal theorem for each small literal outsider; the remaining computation
then becomes finite. THM-4228 removes the gcd-large arithmetic subrays, but it
cannot touch `min(q,r)<=3390` through its `gcd>=3467` hypothesis.

### Repair-or-direct descent

The `Q=17547` activation deck has exactly one nine-cover,

```text
W={85,88,143,168,193,240,252,264,290}.                 (6)
```

That cover is not dangerous. A direct exact pool calculation gives

```text
mu(G_W) ticks=4,802,564,195,362,
c(G_W)=506,
kappa_2(W)=995.                                         (7)
```

Thus `(6)` is uniformly safe far before the repair deck loses its last cover.
At this writing `(7)` is **VERIFIED-SCRATCH**, pending its independent
composite-deck audit. It suggests the stronger certificate:

```text
every B either misses an active repair edge
or B itself has direct two-comb activation <=Q.         (8)
```

The correct next experiment is to enumerate every cover at descending `Q`,
attach each cover's exact `(M_B,c_B,kappa_2(B))`, and locate the first cutoff
where `(8)` holds. This is the LRC analogue of keeping moderate-multiplicity
fibres rather than asking only whether the zero set is empty.

## 2. Planar JC: owner descent removes whole walls

At exact weight ten, THM-4218 closes the dense `zeta!=0` row and THM-4220
closes `zeta=0`. Their union therefore closes the entire coefficient locus

```text
upsilon*xi*(upsilon+xi)!=0                              (9)
```

with `zeta` arbitrary.

At exact weight eleven, THM-4222 closes `U!=0`. THM-4232 descends across
`U=0`: if `Delta!=0`, the replacement rational face is owned by `Delta p^4`;
if `Delta=0`, it is owned by the fixed nonzero term `-1376 p^3/135`.
The genus-five primitive-`Q(zeta_11)` component is unchanged, both replacement
faces are rational, and the complete special-fibre genus remains fifteen.
Consequently THM-4222 and THM-4232 together close

```text
A*B*Z*(A+B)!=0                                          (10)
```

with `U` arbitrary. The remaining named walls in this M11 coefficient chart
are `A=0`, `B=0`, `Z=0`, and `A+B=0`, together with other cells and seam entry.

This is a reusable **owner-descent** operation:

1. enumerate the complete support before deleting an owner;
2. identify the first surviving term in every branch of the wall;
3. recompute lower planes, aggregate collisions, and compactified boundaries;
4. preserve every positive-genus component and its attachment/Hom data; and
5. compare the conserved generic degree with the complete special response.

The fixed residual support point `(2,0,1)` from `-Qs^2/2`, omitted by an old
helper but lying strictly above all relevant faces, is now included and
gap-checked. This is a harmless certificate repair, not a change to THM-4222.

### M12 hostile

The dense weight-twelve main component has an explicit degree-four route to a
`j=0` elliptic curve through

```text
v^2=(W^2-4UZ)P^6+4Z,
x=P^2.                                                   (11)
```

Therefore the M11/M13 statement “every positive-genus source component is
Hom-orthogonal to the good elliptic target” cannot be copied to M12. A viable
M12 proof must retain the quartic quotient, polarization or integral Hom
lattice, attachment images, and the actual carrier-response degree. Any
claimed divisibility mismatch here remains **OPEN** until that integral audit
is complete; reserved THM-4230 is not evidence.

## 3. Mahler: exact pairwise first defects

For every nonzero `x in Z_2` and `0<=s<t`, the strengthened THM-4082 clause is

```text
v_2(H_t(x)-H_s(x))=s+2+2v_2(x).                        (12)
```

It follows immediately because the two errors from the common limit have
unequal valuations; the shallower error survives by strong triangle equality.
Hence a fixed nonzero point is Cauchy but never stabilizes at a finite chart,
and the corresponding carry words agree for exactly the displayed number of
bits. The boundary `x=0` is exact equality.

The lawful transfer to JC is methodological: compare two strict transforms
and find the exact first nonzero normal coefficient. The transfer does not
carry Keller identities, component genus, or Hom data; those remain mandatory
sidecars. The lawful transfer to LRC is a firewall: compatible finite-depth
certificates with moving `q,r,B` do not yield one fixed object. THM-4231 avoids
that error by keeping the same literal triple and proving a uniform tail.

## 4. Sun: tuple mass in a moderate window

For fixed `0<theta<1` and `delta>0`, THM-4036 now records

```text
W_(theta,delta)(X)
 ={n<=X:theta A(X)/X<=a(n)<X^(1/24+delta)},             (13)

sum_(n in W_(theta,delta)(X)) a(n)
 >=(1-theta-o(1))A(X),                                  (14)

#{n<=X:a(n)>=X^(1/24+delta)} <<_delta X^(1-delta).      (15)
```

Low fibres carry at most `theta A(X)` and the energy bound makes the high
tail negligible. These are tuple-mass statements, fully compatible with a
density-one hole set. They neither repair the known hole nor imply positive
support density.

The best LRC transfer is to define an intrinsic repair multiplicity

```text
c_(q,r)(B)=#{lawful labelled repairs R:R intersect B=empty}, (16)
```

not the order-dependent number resolved by a greedy separator. Its first and
second moments should reveal whether most residual bodies have many robust
certificates and isolate a small exceptional support for direct treatment.

## Connection contracts

### Tournament Boolean hierarchy -> LRC component deck

```text
source:       cyclic cell failure masks F_i with lengths ell_i
target:       deletion masses and component counts
map:          lower zeta plus adjacent-union transition count
preserved:    labels, exact mass, cyclic components
destroyed:    literal q/r phase and component addresses
sidecar:      F_(i-1) union F_i and direct boundary scan
test:         independent reverse scatter on every C(30,6) deletion.
```

No tournament is forced: the intrinsic operation is cyclic adjacency. The
transfer from THM-4223/4225 is the Boolean owner tensor plus its operational
sidecar.

### Primitive-pair observable -> arithmetic LRC rays

```text
source:       two danger combs with common factor g and primitive shape (u,v)
target:       one periodic observable sampled at frequency g
map:          factor (q,r)=g(u,v) and integrate the joint safe indicator
preserved:    common divisor, primitive ratio, exact density, oscillation
destroyed:    literal phase alignment inside a coarse universal floor
sidecar:      primitive-shape overlap and component count of the base repair
test:         sharp 1:13 control plus literal joint-wall boundary scans.
```

The next refinement should tabulate shape-specific density and oscillation,
not merely replace the sharp universal `66/91` floor by an unsupported average.

### Mahler first defect -> JC owner descent

```text
source:       two compatible charts and their first nonzero difference
target:       adjacent strict transforms on a coefficient wall
map:          valuation/order of first surviving owner
preserved:    contact order and literal coefficient label
destroyed:    global component, attachment, Keller, and Hom data
sidecar:      full Newton support and regular-model ledger
test:         one complete wall with simultaneous collision deletions.
```

### Sun energy -> LRC repair multiplicity

```text
source:       multiplicity fibres with first/second moments
target:       intrinsic repair counts c_(q,r)(B)
map:          low/high truncation and collision energy
preserved:    total certificate incidence and concentration
destroyed:    exact margins, phases, and outsider labels unless retained
sidecar:      labelled repair list with exact Haar surplus
test:         full c(B) histogram on a frozen residual universe.
```

### Sun zero fibres -> JC coefficient walls

Aggregate generic mass does not classify zero fibres; likewise a dense JC
cell does not classify its coefficient walls. The lawful map is only a warning
to retain lower-tail/wall strata. The decisive JC test remains a complete
source-first resultant, component, genus, attachment, and degree audit.

## Procedurally generated next tasks

1. **LRC anchor:** finish the composite repair-or-direct descent from
   `Q=17547`; freeze all covers and direct thresholds, then independently
   audit the first lower transition.
2. **LRC entry:** build a fixed-prefix/component atlas for every outsider
   `q<3391` not in `P`; seek one cofinal `r` threshold per literal `q`.
3. **LRC finite closure:** after those rays, enumerate the finite remainder in
   `(5)` with a native joint-wall or exact direct-margin verifier, using
   primitive-pair-specific oscillation where the gcd theorem is wasteful.
4. **JC anchor:** apply owner descent to the cheapest of `A=0`, `B=0`, `Z=0`,
   `A+B=0`, selecting by replacement-face genus rather than coefficient name.
5. **JC hostile:** compute the M12 elliptic quotient's integral Hom and
   polarization lattice before using any degree-mod-four proposal.
6. **Multiplicity niche:** emit the intrinsic `c_(q,r)(B)` histogram, second
   moment, minimum margin, and failure masks for the THM-4207 residual bodies
   and the THM-4231 transition covers.
7. **Boolean wildcard:** compare lower-zeta mass, upper-zeta forced-failure
   counts, and cyclic boundary tensors without inventing a cosmetic
   tournament orientation.

The strongest session-level lesson is operational: first use ambient
monotonicity, then descend to the first surviving owner, and only then count
certificates. Scalar summaries become useful after—not before—the label,
boundary, and fixed-object sidecars are secured.
