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

Incoming THM-4229 is proved finite-chart input and THM-4230 is a proved
relative M12 squeeze; neither closes arbitrary LRC entry or the hidden M12
locus. Historical drafts and scratch packets below are evidence only at the
status explicitly stated.

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

## 1. LRC: direct cofinal quadrant and method sidecars

THM-4231 now gives the strongest size-only coverage. For each labelled
`B in binom(P,9)`, retain the literal safe set `U_B=G_B`, its exact mass `M_B`,
and its cyclic component count `c_B`. THM-4170 and Bonferroni give

```text
mu(U_B intersect G_q intersect G_r)
 >=(5/7)M_B-(6c_B/49)(1/q+1/r).                        (1)
```

All `14,307,150` bodies have `45M_B-4>0`. Define

```text
kappa_body(B)=ceil(108c_B/(7(45M_B-4))).               (2)
```

Two structurally independent full censuses agree that the maximum is exactly
`1307`, uniquely at

```text
B*={170,176,190,193,240,252,264,286,290},
M_B* ticks=4,579,301,272,924,
c_B*=618.                                               (3)
```

They agree on the complete seven-field ledger fingerprint
`5b8c08ad9d02622a`. The uniform direct cutoff is `1307`, and the unique
second-largest threshold is `1290`. Thus only `B*` needs treatment below
`1307`. For `1290<=q<=1306`, its first direct-safe larger label is exactly
`r_0(q)=2614-q`; the intervening triangle has `33+31+...+1=289` pairs. Two
independent literal geometries exhaust that triangle and find every pair
strictly safe, with minimum at `(1300,1305)`. Hence every distinct
`q,r>=1290` is safe; with THM-4156/4191 this fills all
`C(32,11)=129,024,480` faces. This is a finite literal patch to a
certificate-sharp direct theorem, not a claimed literal threshold.

### Pair-plane reduction

Every still-uncovered pair now has

```text
min(q,r)<=1289.                                         (4)
```

For genuine outsiders, the infinite remainder is exactly the `1,259` rays
indexed by `q in {1,...,1289}\P`, not the previous finite-box-plus-rays
region. A cofinal theorem on each ray would leave a finite exact remainder. It
is important not to call `(4)` a finite reduction before those tails are proved.

The first ray is now proved by two independent full censuses:

```text
q=1, r>=542,        max_B kappa_1(B)=543.               (4a)
```

Here `543` is the exact certificate maximum. Its unique extremal body is
literally safe at `r=542`, while every other body has certificate cutoff at
most `530`; this closes the ray from `542` without claiming literal
minimality.

The same `B*` is the unique extremizer, but is literally safe at `r=542`.
Thus only `1,258` unbounded outsider rays remain, plus a finite initial segment
on the `q=1` ray.

A broader exact primary also scans all `1,276` labels
`q in {1,...,1306}\P` against all bodies: `18,255,923,400` cases. Every
limiting surplus is positive and the global tail maximum is

```text
K=931 at q=1305,
B={20,170,190,193,240,252,264,286,290}.                 (4b)
```

The load-bearing ray has a separate exact replay, but the full range lacks an
independent census. Hence the exterior safety claim `max(q,r)>=931` is
**VERIFIED-SCRATCH**, not proved canon, pending that audit; its finite
remainder is the box `max(q,r)<=930`.

[THM-4227](../01-canon/theorems/THM-4227-two-outsider-scale-separated-depth-eight-haar-wedge.md)
and [THM-4228](../01-canon/theorems/THM-4228-common-gcd-two-outsider-periodic-observable-haar-ray.md)
are now coverage-subsumed, but not methodologically obsolete. THM-4227's
sequential intersection retains directed scale order and component birth.
THM-4228 factors `(q,r)=g(u,v)`, retains primitive shape, and proves the sharp
two-comb floor

```text
mu(G_u intersect G_v)>=66/91,                           (5)
```

with equality only at ratio `1:13` or its swap. These coordinates are natural
inputs for the fixed-ray problem, where symmetric Bonferroni discards too much
phase information.

### Why the repair route disappeared at the optimum

The old depth-six activation transition remains exact:

```text
Q=17547: one nine-cover
W={85,88,143,168,193,240,252,264,290};
Q=17548: no nine-cover.                                 (6)
```

The new direct census proves

```text
kappa_body(W)=995,                 rho(W)=17548.         (7)
```

The first strict depth-six repair activates only at `3077`. Consequently, for
the composite threshold

```text
theta(B)=min(kappa_body(B),rho(B)),                      (8)
```

the exact repair-or-direct certificate maximum is `1307`, with every optimizer
resolved by the direct branch. The `17547/17548` transition is a useful exact hypergraph
sidecar, but its old safety headline is superseded. This is a concrete version
of the Sun-inspired multiplicity lesson: count certificates only after testing
the target object's own mass and boundary complexity.

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
Hom-orthogonal to the good elliptic target” cannot be copied to M12.
THM-4230 now performs the needed visible integral audit on the exact gate

```text
U Z (W^2-4UZ)(U+W+Z) != 0.
```

Its visible `j=0` Hom lattice is saturated and contributes degrees divisible
by four, while the complete carrier responses are `42` and `34`. Hence every
hypothetical Keller point on this gate lies in the hidden locus

```text
H_0={kappa:Hom(A_12(kappa),E_0)!=0}.
```

This is a proved squeeze, not an M12 exclusion: `H_0` is countable and proper
but nonempty, with `W=0` giving
`A_12~E_0^2 x E_1728^2`. The next owner is therefore an extra correspondence,
not a coefficient. It must retain the polarized integral Hom lattice, all
twelve attachment evaluations, and the node-annihilator condition. The
independent support audit also restores the omitted fixed point `(2,0,1)` at
gap one; this repairs completeness without changing the theorem.

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

The best LRC transfer begins only after fixing a depth `d`, a finite labelled
deck `R_d(q,r)`, a lawfulness predicate `Lawful_(q,r)(R)`, and an exact signed
Haar margin `sigma_(q,r)(R)`. Then define the typed multiplicity

```text
c_(d;q,r)(B)=#{R in R_d(q,r):Lawful_(q,r)(R),
               R intersect B=empty, sigma_(q,r)(R)>0}.       (16)
```

The strict/non-strict margin convention must be frozen with the deck. This is
not the order-dependent number resolved by a greedy separator. Its first and
second moments can supply only Cauchy--Schwarz and Markov support/concentration
bounds; Sun's exponents and zero-fibre conclusions do not transfer. Those
bounds may still isolate a small exceptional support for direct treatment.

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
target:       typed intrinsic repair counts c_(d;q,r)(B)
map:          low/high truncation and collision energy
preserved:    total certificate incidence and concentration
destroyed:    exact margins, phases, and outsider labels unless retained
sidecar:      labelled repair list with exact Haar surplus
test:         full c(B) histogram on a frozen residual universe.
```

Only the moment inequalities transfer through this contract, not the Sun
exponents or a classification of zero fibres.

### Sun zero fibres -> JC coefficient walls

Aggregate generic mass does not classify zero fibres; likewise a dense JC
cell does not classify its coefficient walls. The lawful map is only a warning
to retain lower-tail/wall strata. The decisive JC test remains a complete
source-first resultant, component, genus, attachment, and degree audit.

## Procedurally generated next tasks

1. **LRC anchor:** independently audit the full `18,255,923,400`-case ray
   census; if `931` survives, promote the finite-box reduction with its
   orientation proof and retain `931` as certificate-only.
2. **LRC entry:** refine the remaining `q<=1289` thresholds by reversing
   orientation and retaining primitive phase rather than only symmetric loss.
3. **LRC finite closure:** after those ray tails, enumerate their finite
   remainders with a native joint-wall or exact direct-margin verifier, using
   primitive-pair-specific oscillation where the universal gcd floor is
   wasteful.
4. **JC anchor:** apply owner descent to the cheapest of `A=0`, `B=0`, `Z=0`,
   `A+B=0`, selecting by replacement-face genus rather than coefficient name.
5. **JC hostile:** on the M12 `W=0` split, compute the hidden Prym Hom lattice,
   all twelve attachment evaluations, and the node-annihilator sublattice;
   test whether its degree form represents `34` or `42`.
6. **Multiplicity niche:** freeze `(d,R_d,Lawful,sigma)` and emit the intrinsic
   `c_(d;q,r)(B)` histogram, second
   moment, minimum margin, and failure masks on the fixed-small-ray residuals;
   compare it with direct `(M_B,c_B)` thresholds before invoking repair.
7. **Boolean wildcard:** compare lower-zeta mass, upper-zeta forced-failure
   counts, and cyclic boundary tensors without inventing a cosmetic
   tournament orientation.

The strongest session-level lesson is operational: first use ambient
monotonicity, then descend to the first surviving owner, and only then count
certificates. Scalar summaries become useful after—not before—the label,
boundary, and fixed-object sidecars are secured.
