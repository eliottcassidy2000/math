---
id: THM-2531
title: "Prime-necklace guard boundary selector and pair-labelled odd lift"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  THM-2439's translation-equivariant
  lexicographic marker, applied to every nonconstant thirteen-root Boolean
  mask with the retained guard slope, canonically selects the last occupied
  root and first empty root of its maximal initial run.  This is an affine-
  covariant adjacent occupied-to-empty domain wall.  Its binary word scores
  form one length-thirteen doubling orbit modulo 8191 and orient a transitive
  no-tie tournament whose source is the marker; complement is its exact
  global converse.  Reflection is covariant only when the oriented slope is
  transported.  Exactly 126 of the 630 nonconstant rotation necklaces have
  a reflection stabilizer, so forgetting that slope creates genuine
  dihedral orientation ties.  On the full thirteen-arm live fibre, the
  selector refines every positive THM-2527 mask layer into disjoint Boolean
  atoms carrying one actual ordered predecessor boundary.  After enlarging
  THM-2528's four-arm observable to retain that full mask, every positive
  selector class still has strict positive-versus-negative path imbalance,
  pointwise.  On THM-2529's actual singleton/adjacent-pair deep-comb bank the
  ambient degree-twelve marker collapses to a quadratic local wall.  In the
  THM-2530 target-anchored path chart, its class masses are linear Gram
  coordinates and give a genuine skew selector star from an occupied root
  to the excluded target on every fibre; rational positivity forces all
  twelve primitive root colours.  The adjacent wall itself reaches the
  target only on endpoint rays, whose mass need not be positive.  This is a
  pair-labelled fibre lift, not a projection of the general selector into
  the current scalar algebra: THM-2439's homometric hostile and ambient
  degree-twelve marker cost remain.  The empty target is an excluded event,
  not semantic arrival or owner; no owner loop, row exclusion, or LRC(14)
  proof is claimed.
source: codex-2026-07-27-prime-necklace-boundary-selector
depends_on:
  - THM-2439-cyclic-marker-replica-degree-and-homometric-gram-boundary
  - THM-2527-owner-weighted-all-mode-odd-bank-and-boolean-cut-coordinate
  - THM-2528-intrinsic-four-arm-boolean-path-and-joint-autocorrelation-scalarization
  - THM-2529-deep-comb-adjacent-fibre-odd-consumer-and-zero-target-boundary
  - THM-2530-anchored-deep-gram-cone-and-lossless-skew-target-refinement
related:
  - THM-2518-perron-inverse-branch-owner-word-cospan-recovery
  - THM-2526-affine-skew-orientation-gauge-boundary
script: 04-computation/lrc14_prime_necklace_boundary_selector_thm2531.py
output: 05-knowledge/results/lrc14_prime_necklace_boundary_selector_thm2531.out
script_sha256: 386ca6f9c87feee4c899505b1c8913215e90347cee1140499eb2070ad4443693
output_sha256: 9f8eb5eef6428498e861784b0ec6a4aae79ea2b33ab2ba74f542ef22a3235276
hash_basis: working-tree bytes (LF)
---

# THM-2531 -- the prime necklace chooses one guard domain wall

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The raw lexicographic root marker is not new here.  THM-2439 already proves
its translation equivariance, exact Boolean replica degree, and failure to
descend through cyclic autocorrelation.  The new point is to put that marker
on THM-2527's positive full-mask layer *before* the autocorrelation quotient
and to retain it when THM-2528's signed operator path is opened.  This turns
the previously unmarked positive layer into disjoint atoms carrying one
actual ordered root boundary.

There is also a useful dynamical and tournament description:

```text
prime Boolean necklace
  -> length-thirteen doubling orbit modulo 8191
  -> transitive word tournament with a unique source
  -> canonical first guard-directed 1-to-0 domain wall.          (1)
```

The construction uses the full thirteen-arm mask.  It does not claim that
the old four-arm projection or the scalar collision profile reconstructs
the selected edge.

## 1. Prime freeness, the Mersenne orbit, and the word tournament

Let

```text
X=F_13,                    tau in F_13^*,
e:X->{0,1},                1<=sum_r e_r<=12.                    (2)
```

The nonzero `tau` is the oriented guard slope retained by THM-2526.  For
each root `a`, read the cyclic word and its binary value in guard order:

```text
W_(tau,a)(e)=(e_(a+j tau))_(j=0)^12,

L_(tau,a)(e)=sum_(j=0)^12 2^(12-j)e_(a+j tau).                 (3)
```

The thirteen words are distinct.  Indeed, equality of the words based at
`a` and `b` makes `e` invariant under translation by `b-a`.  If `a!=b`,
that translation generates the prime cyclic group, so `e` is constant,
contrary to (2).  Thus the rotation stabilizer of every mixed mask is
trivial.

Let

```text
alpha_tau(e)=the unique a maximizing L_(tau,a)(e).             (4)
```

This is exactly THM-2439's lexicographic marker, transported from slope one.
There is an equivalent Mersenne recurrence.  Put `M=2^13-1=8191`.  Directly
shifting (3) gives

```text
L_(tau,a+tau)=2L_(tau,a)-M e_a,
L_(tau,a+tau)=2L_(tau,a)                  mod M.                (5)
```

The number `8191` is prime and `2` has order `13` modulo `8191`.  Hence the
thirteen binary words are the ordinary representatives of one exact
length-thirteen doubling orbit.  Primality of the root cycle, rather than
Mersenne primality, is the conceptual freeness proof above; (5) is a second
exact realization of the same object.

Orient the complete graph on the thirteen roots by

```text
a ->_e b                  iff L_(tau,a)(e)>L_(tau,b)(e).        (6)
```

This is an intrinsic tournament: vertices are root predecessors, the
pairwise observable is binary-word comparison, `tau` is the orientation
gauge, and there are no ties.  It is transitive because (6) is a strict
total order.  Its unique source is `alpha_tau(e)`.  Moreover every occupied
root beats every empty root: an occupied word begins with `1`, while an
empty word begins with `0`.  The support bits and `tau` are load-bearing
sidecars; the unlabelled tournament order alone is not being asserted to
recover the physical mask.

## 2. The canonical adjacent occupied-to-empty boundary

The maximal word begins with `1`.  Since the mask is proper, define

```text
q_tau(e)=min{j in {1,...,12}:e_(alpha_tau(e)+j tau)=0},

s_tau(e)=alpha_tau(e)+(q_tau(e)-1)tau,
t_tau(e)=alpha_tau(e)+q_tau(e)tau=s_tau(e)+tau.                (7)
```

Then

```text
e_(s_tau(e))=1,                  e_(t_tau(e))=0.               (8)
```

The initial run length `q_tau(e)` is the maximum cyclic run length of ones.
Otherwise a longer occupied run would give a lexicographically larger word.
When several maximum runs have the same length, the remaining cyclic word
breaks the tie uniquely.  Thus

```text
partial_tau(e)=(s_tau(e),t_tau(e))                             (9)
```

is a canonical adjacent guard-directed domain wall.  It is also an arc of
the word tournament (6).  We retain the fuller selector packet

```text
Sigma_tau(e)=(alpha_tau(e),q_tau(e),partial_tau(e)),           (10)
```

because `(q,partial)` recovers the marker `alpha`.

The exhaustive ambient census provides a useful scale check.  There are

```text
(2^13-2)/13=630                                               (11)
```

nonconstant rotation necklaces.  Grouped by the selected maximum one-run
length `q`, their counts are

```text
q:       1   2    3    4   5  6  7  8  9 10 11 12
count:  40 172  178  114  62 32 16  8  4  2  1  1.           (12)
```

Across all `8,190` labelled mixed masks, every possible anchor occurs `630`
times and every adjacent directed guard edge occurs `630` times.  These
uniform census facts are consequences of covariance, not a uniformity claim
about the live LRC distribution.

## 3. Exact affine, reflection, and complement laws

For `u!=0` and `b in F_13`, let the affine pushforward of a mask be

```text
((u,b)_*e)_r=e_(u^(-1)(r-b)).                                  (13)
```

Transport the oriented slope at the same time, `tau -> u tau`.  Then

```text
W_(u tau,ua+b)((u,b)_*e)=W_(tau,a)(e),                         (14)

alpha_(u tau)((u,b)_*e)=u alpha_tau(e)+b,
q_(u tau)((u,b)_*e)=q_tau(e),
partial_(u tau)((u,b)_*e)
 =(u s_tau(e)+b,u t_tau(e)+b).                                (15)
```

Thus the selector uses no affine origin.  It is equivariant, not invariant:
the selected physical roots move when the chart is relabelled.

For the reflection `J_c:r->c-r`, fixed-slope and transported-slope laws
must not be conflated:

```text
alpha_tau(J_c e)=c-alpha_(-tau)(e),
partial_tau(J_c e)=c-partial_(-tau)(e) coordinatewise,         (16)

partial_(-tau)(J_c e)
 =(c-s_tau(e),c-t_tau(e)).                                    (17)
```

Equation (17) is the lawful simultaneous reflection of the roots and guard
direction.  Replacing `tau` by `-tau` at a fixed mask is generally neither
the same word tournament nor its converse.  For the singleton mask `{0}`,
the two descending root orders are

```text
(0,-1,-2,...,-12)              and (0,1,2,...,12),            (18)
```

which are neither equal nor reverse.  Therefore guard reversal, root
reflection, and semantic converse cannot be silently identified.

Complement has a different and exact law.  Since the binary weights sum to
`M`,

```text
L_(tau,a)(1-e)=M-L_(tau,a)(e).                                (19)
```

Consequently complement takes the word tournament to its global converse,
and its marker is the old unique minimum:

```text
alpha_tau(1-e)=argmin_a L_(tau,a)(e).                          (20)
```

The complement selector chooses a guard-directed zero-to-one boundary of
the original mask.  It is a second canonical wall, not in general the
reverse of (9).

## 4. Dihedral ties are real and exactly classified

Prime length removes rotations but not reflections.  Each reflection of a
thirteen-cycle has one fixed root and six transposed pairs, so it fixes
`2^7=128` masks, `126` of them nonconstant.  A nonconstant mask cannot be
fixed by two different reflections: their product would be a nontrivial
rotation.  Hence exactly

```text
13*126=1,638                                                   (21)
```

labelled mixed masks, or `126` of the `630` rotation necklaces, have a
reflection stabilizer.

The maximum forward and backward cyclic words agree exactly on those `126`
necklaces.  Indeed, equality of the two maxima based at `a,b` says
`e_(a+j tau)=e_(b-j tau)` for every `j`, which is invariance under the
reflection `r->a+b-r`; the converse is immediate.  They are genuine ties if one tries to canonicalize after
quotienting the two guard directions.  The other `504` necklaces pair under
reflection, giving

```text
126+(630-126)/2=378                                           (22)
```

dihedral mask orbits.  Retaining `tau` resolves every tie before quotient.

There is no hidden complement-dihedral tie.  If a dihedral permutation sent
`e` to `1-e`, cardinality preservation would force

```text
sum e=13-sum e,                                               (23)
```

impossible at odd length.  Thus complement and reflection are separate
operations throughout this construction.

## 5. Pair-labelled Boolean refinement of the positive mask lift

Return to THM-2527's live full root fibre.  At a future base point `z`, let

```text
e(z)=(e_r(z))_(r in F_13) in {0,1}^13,
g(z) in {0,1}                                                 (24)
```

be the complete predecessor mask and common late owner or owner--word
factor.  Define the selector cylinders

```text
Lambda_(a,q)^tau
 ={z:alpha_tau(e(z))=a, q_tau(e(z))=q},
                                  a in F_13, 1<=q<=12.         (25)
```

They are finite Boolean predicates of the thirteen root bits.  Equivalently,
`alpha=a` says that the integer `L_(tau,a)` is strictly larger than its
twelve competitors, and `q` records the first zero in that word.  The
`156` cylinders are disjoint and partition the mixed-mask locus:

```text
1_(e nonconstant)=sum_(a,q)1_(Lambda_(a,q)^tau).               (26)
```

On `Lambda_(a,q)^tau`, the selected actual predecessor pair is

```text
(s,t)=(a+(q-1)tau,a+q tau),        e_s=1, e_t=0.               (27)
```

Let `Psi_tau(e)` be THM-2527's integer cut score.  Its exact properties are

```text
Psi_tau(e) in {0,...,98},
Psi_tau(e)=0 iff e is constant.                               (28)
```

The old positive threshold lift can now be refined without an external
copy coordinate:

```text
g Psi_tau(e)
 =sum_(j=1)^98 sum_(a,q)
    1_{g=1, Psi_tau(e)>=j, Lambda_(a,q)^tau}.                  (29)
```

Every summand in (29) is a Boolean event on the same future base and carries
the ordered pair (27).  Integrating gives the exact pair-labelled form of
the positive odd coordinate:

```text
O_tau(-4tau)
 =sum_(j=1)^98 sum_(a,q)
    measure{g=1, Psi_tau(e)>=j, Lambda_(a,q)^tau}.             (30)
```

This uses no preferred root: an affine chart change permutes the summands
according to (15).  It also introduces no Bernoulli or formal copy tag.

For a rational owner weight `0<=g<=1` with common denominator `D`, use the
finite threshold identity

```text
g=(1/D)sum_(ell=1)^D 1_(Dg>=ell).                             (31)
```

Substituting (31) in (29) gives a finite positive Boolean layer expansion
with the harmless scalar `1/D`.  For an arbitrary real owner weight the
same selector gives a positive weighted partition, but not a finite Boolean
owner expansion.  The live owner--word case is Boolean.

## 6. Every positive selector class retains the four-arm imbalance

This step has a type boundary.  THM-2528's original path carrier records
the four addresses

```text
(r,r+t,r+t+v,r+t+v+h)                                        (32)
```

and the two occupied endpoint bits.  Those four observed arms do **not**
determine the thirteen-bit marker.  Before restricting the path sets, enlarge
the observable to the full thirteen-arm natural-extension fibre `e(z)` and
then append the old path coordinates `(r,v,h)`.  This is the same base `z`
and uses no independent copy, but it is strictly richer than the four-arm
projection.

Let `P_t,N_t` be THM-2528's two disjoint Boolean path unions, with uniform
counting measure in their three finite address coordinates.  Inserting any
full-mask Boolean predicate `E(z)` into the same path expansion gives, at
the fixed coordinate `t=-4tau`,

```text
13*12^2 [measure(P_(-4tau) intersection E)
         -measure(N_(-4tau) intersection E)]
 =integral g(z)1_E(z)Psi_tau(e(z))dz.                         (33)
```

Take `E=Lambda_(a,q)^tau`.  By (28), every selector class of positive
owner measure satisfies the strict inequality

```text
measure(P_(-4tau) intersection Lambda_(a,q)^tau)
 >measure(N_(-4tau) intersection Lambda_(a,q)^tau).           (34)
```

There is no cancellation between selector classes: (34) holds for **every**
non-null class, not merely for one class obtained after summation.  The same
is true after aggregating cylinders that select a fixed adjacent edge, or
after intersecting with any positive threshold layer `Psi_tau>=j`.

Expanding (33) back into the paired `h,-h` wedges of THM-2528 shows more.
Inside every positive selector class there exist `v`, a positive
THM-2528 tournament-kernel direction `h`, and an absolute root `r` for which
the corresponding two occupied endpoint intersections have the strict
signed comparison of THM-2528 equation (19), now with the canonical boundary
(27) on the same full mask.  Thus the odd path and selected domain wall
coexist on one ancestry fibre.

They are not the same pair.  Both endpoints of a `P/N` collision path are
occupied, while the target in (27) is empty.  The selector is a retained
domain-wall sidecar on the full mask, not a claim that the empty root is one
of the positive collision endpoints.

## 7. The sharp wall and the autocorrelation quotient boundary

THM-2527's sharp cut inequality has `52` nonconstant equality masks in four
rotation necklaces.  The selector gives a common structural description of
all of them:

```text
q_tau(e)=4.                                                   (35)
```

Thus on the sharp `13/42` wall, the canonical edge is the terminal edge of
a lexicographically selected maximum four-run.  In slope-one canonical
words the four necklaces are

```text
1111000011000,  1111000011100,
1111000110000,  1111001110000.                                (36)
```

The selector must be attached before autocorrelation.  THM-2439's exact
homometric pair

```text
A={0,1,3,9},                    B={1,2,5,7}                    (37)
```

has the identical complete cyclic correlation profile

```text
c(A)=c(B)=(4,1,1,...,1),                                      (38)
```

and therefore the same `Psi` and every downstream collision-power scalar.
Nevertheless, at `tau=1`,

```text
partial_1(A)=(1,2),                 partial_1(B)=(2,3).        (39)
```

The masks are not dihedrally equivalent.  Hence no function of the complete
self-correlation, including THM-2527's odd scalar coordinate, can recover
the selected boundary universally.

There is a matching algebraic invoice.  THM-2439 proves that the marker
indicator has intrinsic independent-replica degree `12` on all nonempty
proper thirteen-bit masks (`10` on its smaller cardinality-at-most-ten
domain).  The full packet (10) recovers that marker, so projecting this
specific selector into a same-gauge scalar/current algebra inherits the
degree-twelve obstruction.  Existing four-arm products do not pay it.

Thus (29)--(34) are lawful Boolean statements on the retained full fibre,
but they are not a theorem that the present scalar current algebra emits the
marker or selected-edge current.  Conflating those two levels would repeat
exactly the quotient error isolated by THM-2439.

## 8. On the actual deep comb the marker collapses to degree two

The ambient degree-twelve invoice is not the right invoice on every live
subfamily.  THM-2529 proves that the actual deepest-comb mask, in its natural
cyclic root chart, is always

```text
{r}                 or {r,r+1}.                              (40)
```

There is only one occupied run.  At the aligned positive deep slope, the
marker and terminal-wall indicators reduce pointwise to

```text
1_(alpha_1(e)=a)=e_a(1-e_(a-1)),

1_(partial_1(e)=(s,s+1))=e_s(1-e_(s+1)).                    (41)
```

Thus the exact lex selector is quadratic on this `26`-mask bank.  No
degree-twelve replica service is needed to select its wall.  The ambient
lower bound in Section 7 remains correct for arbitrary mixed masks and for
any argument that has already forgotten the one-run restriction.

The same collapse holds for every oriented slope, not only `tau=1`.  Let
`k_tau in {1,...,12}` represent `tau^(-1)`, and put

```text
epsilon_tau=+1,                     1<=k_tau<=6,
epsilon_tau=-1,                     7<=k_tau<=12.             (42)
```

For a two-root mask `{r,r+1}`, the second occupied bit occurs at position
`k_tau` when the word starts at `r` and at position `13-k_tau` when it
starts at `r+1`.  Therefore

```text
1_(alpha_tau(e)=a)=e_a(1-e_(a-epsilon_tau))                  (43)
```

on the entire deep-comb bank.  This gives a quadratic marker in the actual
physical guard direction, whatever nonzero unit that direction is.

The target-anchored THM-2530 chart now supplies a stronger ordered pair than
the adjacent wall alone.  Relative root `0` is forbidden, so the only masks
are

```text
{j},             1<=j<=12,
{j,j+1},         1<=j<=11.                                  (44)
```

For every one of these masks and every retained `tau`, the marker root is
occupied and root `0` is empty.  Hence

```text
(alpha_tau(e),0)                                                (45)
```

is a canonical ordered occupied-to-excluded-target pair on **every** fibre.
It is an arc of the word tournament because every occupied word beats every
empty word.  Usually it is a chord rather than an adjacent guard edge.  In
absolute target coordinates it is `(t+alpha_tau(e),t)`, and affine chart
transport moves both endpoints covariantly.

Let `F` be the nonnegative weight in one lawful target cell and `K` be
THM-2530's anchored Gram matrix.  Define the marker-class masses

```text
gamma_a^tau
 =integral F e_a(1-e_(a-epsilon_tau))
 =K(a,a)-K(a-epsilon_tau,a).                                (46)
```

The `gamma_a^tau` are nonnegative, their pointwise events are disjoint, and

```text
gamma_0^tau=0,
sum_a gamma_a^tau=integral F.                               (47)
```

Thus every positive cell has at least one positive canonical target-arc
class.  More importantly, the selector masses are linear coordinates of the
already-lawful Gram matrix.  THM-2530's lossless inverse from its skew matrix
back to `K` therefore recovers the entire vector `gamma^tau`; the ambient
degree-twelve current obstruction has genuinely disappeared on this bank.

The selected arcs themselves form a clean skew star current.  Pointwise put

```text
S_tau(e)=E_(alpha_tau(e),0)-E_(0,alpha_tau(e)).               (48)
```

After weighted integration,

```text
S_tau(a,0)= gamma_a^tau,
S_tau(0,a)=-gamma_a^tau,
S_tau(a,b)=0                         if a,b!=0.                (49)
```

Every nonzero positive entry is the measure of a literal selected
occupied-to-excluded-target arc.  This is an oriented star, not a complete
tournament, and it is lossless for `gamma^tau`.

There is also no root-colour loss.  Suppose the target cell and its weight
are rational and positive.  If one nontrivial Fourier coefficient of
`gamma^tau` vanished, its rational degree-at-most-twelve polynomial would
vanish at a primitive thirteenth root.  Irreducibility of `Phi_13` would
force all thirteen `gamma_a^tau` equal.  Equations (47) make that impossible.
Consequently

```text
gamma_hat^tau(k)!=0                      for every k!=0.       (50)
```

This is a genuine nonnegative selector current with all twelve root colours;
it avoids the autocorrelation phase loss on the actual deep-comb bank.  It
does not assert a nonzero target-character mode across the `(s,t)` cell
indices.

If one insists that the selected pair also be an adjacent domain wall, the
endpoint issue from the earlier construction remains exact.  In the positive
deep direction, the adjacent wall ends at the target only for

```text
{12}, {11,12},              partial_1(e)=(12,0),              (51)
```

while in the negative direction it does so only for

```text
{1}, {1,2},                 partial_(-1)(e)=(1,0).             (52)
```

Writing `alpha_j,beta_j` for THM-2530's singleton/adjacent-pair masses, the
two adjacent target-wall masses are

```text
K(12,12)=alpha_12+beta_11,
K(1,1)  =alpha_1+beta_1.                                    (53)
```

Neither endpoint diagonal is forced: an interior ray has both zero.  This
does not affect the uniform target chord (45); it only blocks upgrading that
chord to the first adjacent zero in a prescribed deep direction.

THM-2529's exact zero-target hostile has relative singleton mask `{1}`.
Its canonical target chord is `(1,0)` for every slope, and its negative-deep
selector makes that chord the literal adjacent wall

```text
(1,0)=occupied root -> excluded target root.                 (54)
```

This sharpens the hostile's geometry but does not repair its semantics: root
`0` is still empty and has zero target charge.  Moreover `-1` must be
retained as the oriented deep chart used in (52).  It may not be silently
identified with THM-2526's physical guard slope when that live slope is a
different unit.  The nonadjacent target chord (45) is valid at the physical
guard slope; only the adjacent-wall identification in (52) spends `-1`.

## 9. Exact semantic boundary

The advance is precise:

1. every live mixed mask has one affine-covariant, guard-directed actual
   predecessor pair `(occupied,empty)`;
2. every positive THM-2527 threshold atom can carry that pair without an
   external origin or copy coordinate; and
3. after retaining the full mask, every non-null pair class preserves the
   strict THM-2528 four-arm sign; and
4. on the actual target-anchored deep comb, a quadratic selector produces a
   genuine occupied-to-excluded-target skew star whose nonnegative class
   vector has all twelve primitive root colours.

In the ambient construction, the selected target is empty because the old
collision predicate `F` is absent there.  In the anchored deep construction,
root `0` is more strongly typed as the excluded target/safe-complement root,
but it still has no positive target charge.  Absence is not a semantic
arrival event, a new owner, or a future packet.  Both selected roots lie at
the same predecessor horizon, while the common owner factor occurs later.
No temporal old-to-new map, source-to-arrival cospan, blocker-to-owner loop,
signed scalar-cover current, row exclusion, or proof of LRC(14) follows.

For arbitrary mixed masks, the next lawful use must either give the empty
side a positive semantic meaning or construct the degree-twelve full-mask-
to-current intertwiner.  Calling the selected empty root an arrival without
one of those maps is forbidden.

On the actual one-run deep-comb bank, Section 8 already removes the ambient
marker cost and gives the target-ending chord uniformly.  The remaining
obligation there is to send the excluded/safe target type through a lawful
semantic exit or couple `gamma^tau` to a nonzero target-character/scalar-
cover mode.  Endpoint positivity in (53) is needed only if that exit insists
on an **adjacent** target wall.

## 10. Exact dependency-free referee

Run

```bash
python3 04-computation/lrc14_prime_necklace_boundary_selector_thm2531.py
python3 -O 04-computation/lrc14_prime_necklace_boundary_selector_thm2531.py
```

Both executions reproduce

```text
05-knowledge/results/lrc14_prime_necklace_boundary_selector_thm2531.out
```

byte-for-byte.  The referee uses only integer arithmetic and checks:

- all `8,190` mixed masks, their trivial cyclic stabilizers, Mersenne
  recurrences, unique anchors, maximum runs, selected boundaries, and word
  tournament orders;
- `106,470` translation and `98,280` multiplicative-slope covariance
  instances, the fixed-slope reflection law, and exact complement converse;
- the `630/126/378` rotation/reflection/dihedral counts and absence of every
  complement-dihedral tie;
- every cut score and threshold lift, all `52` sharp masks, the four sharp
  words, and selected run length four;
- the uncollapsed `72+72` four-arm kernels and pointwise strict imbalance on
  every mixed mask and all `156` selector classes; and
- THM-2439's homometric correlation hostile with its two distinct selected
  boundaries;
- all `26` singleton/adjacent deep masks, the quadratic marker/wall formulas,
  all twelve slope transports, all `23` target-anchored path masks, the four
  endpoint cases, and the hostile boundary `(1,0)`; and
- linear Gram recovery of the marker/wall masses, their positive target-arc
  partition, and the nonconstant cyclotomic selector vector.

The finite script audits the complete ambient Boolean bank.  Equations
(14)--(20), (26)--(34), (43)--(50), and the semantic carrier distinction are
the exact symbolic proofs and typing argument above.  The all-colour claim
in (50) is the `Phi_13` irreducibility proof, with the script checking its
nonconstant positive coefficient hypotheses exactly.

**QED.**

The independent audit rederived prime-cycle freeness, the doubling-orbit and
dihedral counts, the affine/reflection/complement laws, the selector-class
path imbalance, and the homometric quotient hostile.  On the actual deep
comb it separately checked the all-slope quadratic marker, the uniform
occupied-to-excluded-target chord, the Gram formula for `gamma^tau`, and the
`Phi_13` all-colour argument.  Normal and optimized executions reproduced
the stored output byte-for-byte.
