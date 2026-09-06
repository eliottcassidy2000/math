# Virtual pair walls enter an actual primitive ten-pack component

**Status: PROVED elementary signed-pair and event lemmas + PROVED
unbounded constructive LRC(14) subclass + FINITE-EXACT; independently
audited.** LRC(14) remains open. The family below has a positive interval
of strict physical witnesses. No actual THM-3818 decoder entry is claimed.
[Independent audit: PASS](overnight9_20260906_lrc_virtual_pair_wall_audit.md).
The [repaired sharp margin](overnight9_20260906_lrc_virtual_pair_wall_margin.md)
strengthens the one-sided band below. Incoming THM-4448 also guarantees
existence for the large-scale family; this report retains an explicit
virtual-wall witness and owner word, not a new first existence claim.

## 1. Inheritance and the interface recovered

Current frontier and guardrails were read first. The closest mechanisms are
**THM-4335 / owner permutations, component addresses, and minority renewal**
(`01-canon/theorems/THM-4335-lrc14-owner-permutation-component-address-and-minority-renewal.md`)
and the independently audited fourth-round event interface
`05-knowledge/results/overnight4_20260906_lrc_body_event.md`.

Incoming **THM-4446 / primitive ten-pack descent and dilation rays**
(`01-canon/theorems/THM-4446-lrc14-primitive-ten-pack-descent-and-dilation-rays.md`)
removes common pack dilation. Incoming **THM-4447 / composite-clock capacity
and small-clock reduction**
(`01-canon/theorems/THM-4447-lrc14-composite-clock-capacity-and-small-clock-reduction.md`)
was read directly from origin/main: it packages the exact gcd-sensitive
capacity and recovers the prior
`05-knowledge/results/lrc14_effective_clock_empty_core_sep06.md`.
Their remaining clock-three signature needs component placement or owner
overlap, not a smaller sum of three individual capacities.

The canonical hostile is the body from **THM-530, Section A**
(`01-canon/theorems/THM-530-lrc-gp-intersection-global-witness-floor.md`): its
measure 14249/252252 is a finite bounded-body minimum and lies below 6/77.
It is not a universal ten-body Haar floor. A second hostile, already in the
fourth report, blocks all ordinary tail event cosets while retaining a safe
row. The corrected near miss here is a tempting deletion inference: no
d-owned surviving endpoint does not imply the d constraint is redundant.
It can delete whole components. Section 4 gives the minimal two-frequency
example and an actual primitive ten-pack embedding.

The least-used sidecar is a **signed pair frequency of the tails, divided
by three**, together with its own body-wall coset. It is usually not a tail
frequency and is not the gcd of the body. The arithmetic sum/difference
clocks in THM-4106 concern a different pair observer and do not themselves
give the owner-band statement below.

Concept board: tail carries; relative owner labels; signed pair frequency;
virtual body walls; positive component addresses; hereditary gcd filters.
Source is an actual tail owner word above a body phase. Target is the
distance of the signed pair frequency from an integer. The map preserves
the necessity of distinct owners under full spoil, but loses individual
activity and tooth addresses. The full wall code and a retained body point
restore those coordinates. The cheap hostile at y=1/4 below prevents an
iff claim after taking this quotient.

Portfolio: this is the actual LRC(14) Anchor; the simultaneous Laurent/jet
work is the Niche; the independently audited no-three-line density theorem
is the Wildcard. Their shared method retains an actual collision coordinate
before using a scalar bound; no mathematical reduction between them is
being asserted.

## 2. A sharp signed-pair band for any ternary-unit tail

For finite positive C and three distinct positive ternary units T, write

    G_C={y mod1: ||cy||>=1/14 for all c in C},
    K_w(y)={j in Z/3Z: ||w(y+j)/3||<1/14},
    F_T={y: union_(w in T) K_w(y)=Z/3Z}.

Every strict mask has at most one owner. A tail w is active exactly when
||wy||<3/14; then, for the unique nearest integer n_w to wy,

    K_w(y)={-n_w w^(-1) mod3}.                         (1)

Equality at clearance 3/14 is inactive for the strict mask. A full row
3C union T is unsafe precisely when G_C is contained in F_T, but the
sufficient gates below will not be promoted to equivalent conditions.

For a pair u,v of distinct ternary units, choose epsilon in {+1,-1} so
u+epsilon v is divisible by three, and put

    d=|u+epsilon v|/3 > 0.

**Signed-pair lemma.** At every phase y,

    y in F_T  implies  ||dy||>4/21.                   (2)

More locally, if ||dy||<=4/21, the two tails are either not both active,
or they have the same bad-sheet owner. The third tail then cannot fill
all three sheets. This holds for every pair in T simultaneously.

**Proof.** Suppose the pair is active and write wy=n_w+e_w with
|e_w|<3/14. Since u=-epsilon v modulo three, equality of the two owners
in (1) is equivalent to n_u+epsilon n_v=0 modulo three. Distinct owners
therefore give

    ||(n_u+epsilon n_v)/3||=1/3,
    |(e_u+epsilon e_v)/3|<1/7.

The reverse triangle inequality for circle distance proves (2), strictly.
The strict inequality includes the weak-safe boundary in its contrapositive.
Changing the representative y by an integer changes the combined tooth
integer by a multiple of three, so this is a relative-owner statement
independent of sheet gauge. Taking the absolute value in d loses no distance.

The constant is sharp as a uniform statement, even for primitive additive
tails. For any positive even N, set

    T=(9N-1,33N-1,42N-2),  y=1/(42N-2),  d=8N.

These tails are primitive: gcd(9N-1,24N)=gcd(9N-1,24)=1. Their nearest
tooth integers are (0,1,1), all are active, and the owners are (0,1,2).
Indeed (9N-1)/(42N-2)<3/14. Yet

    ||dy||=4N/(21N-1)=4/21+4/[21(21N-1)] -> 4/21.     (3)

No larger universal constant can replace 4/21 in (2). Conversely the three
pair bands do not characterize F_T: for T=(1,4,5), its virtual pair
frequencies are 1,2,3, and at y=1/4 all three distances exceed 4/21, while
tail one is inactive. Thus y is not fully spoiled.

## 3. Virtual wall codes: a pointwise body-to-owner bridge

For any virtual frequency d from Section 2, use the numerator-one coset

    V_d^+={y_k=(14k+1)/(14d): k in Z/dZ}.             (4)

Every point satisfies ||dy_k||=1/14<4/21. Consequently

    V_d^+ intersects G_C  implies  G_(3C union T) nonempty.   (5)

A witness is obtained by evaluating the three literal masks and choosing
a missing sheet j, then x=(y_k+j)/3. If d belongs to C this is an actual
body wall; (5) does not require that membership. Reflection handles the
negative coset. These are not the ordinary effective tail endpoints, whose
numerator is three and whose frequencies belong to T.

For body frequency c, its exact virtual blocker code is

    B_d(c)={k: dist_(14d)(c(14k+1),0)<d}.            (6)

Put g=gcd(c,d), q=d/g and h=c/g. Multiplication by h permutes the q reduced
classes; every class has multiplicity g. The admissible integer distances
are -q<n<q and n=h modulo fourteen. Therefore

    |B_d(c)|=g[floor((q-1-h)/14)+floor((q-1+h)/14)+1],
    |B_d(c)|<=g ceil(q/7),
    B_d(c)=Z/dZ iff 14d divides c.                  (7)

This is the inherited affine cyclic count, with a different numerator; no
novelty is claimed for counting residues in an open interval. Its use at
a signed-pair virtual wall is the new body-to-owner application.

For example, suppose C has ten members, at least four are divisible by d,
none is divisible by 14d, and all the other members are coprime to d.
The multiples block no virtual point. At most six remaining frequencies
each block ceil(d/7), so

    6ceil(d/7)<d                                     (8)

forces an actual full-row witness by (5). Condition (8) holds for all
integers d>=37 by checking residues 37,...,43 and advancing by seven.
It also holds at d=31; every integer from 32 through 36 has a factor in
common with 42. Hence it holds whenever d>=31 and gcd(d,42)=1.

Unlike a common-clock descent, the whole body need not be divisible by d.
Unlike the ordinary tail event criterion, the frequency d can be absent
from T. Both distinctions are used in the following family.

## 4. Why an owner-free deletion need not be redundant

Take C=(1,14) and d=1. The safe set G_(14) consists of fourteen components

    [(14k+1)/196,(14k+13)/196],  k=0,...,13.

Imposing the frequency-one condition removes exactly the first and last
components, leaving the other twelve unchanged. Neither surviving endpoint
is a frequency-one wall. Both such walls, 1/14 and 13/14, are instead
strictly blocked by frequency fourteen. Thus the d constraint is not
redundant even though it owns no surviving component endpoint.

This mechanism occurs in a primitive ten-pack too: take C={1,...,9,14}.
The phase 1/28 is safe for C without 1 and unsafe for 1, whereas 1/11 is
safe for all of C. Frequency fourteen still blocks both possible 1-walls.
No component-level deletion dichotomy may ignore the fully removed
components. The positive theorem below exhibits an actual wall and a
positive component cell explicitly, avoiding that missing implication.

## 5. A constructive unbounded primitive family

Let d>=31 satisfy gcd(d,42)=1. Define

    a=1,  b=3d+1,  c=3d+2,
    L=42dbc,  h=1+mL,  m>=1,
    C={d,2d,3d,4d,14,14b,14c,h,2h,1},
    T={1,b,c},  S=3C union T.                        (9)

All ten body frequencies are distinct; m>=1 ensures h exceeds every
base body frequency. The tail is primitive, additive and ternary-unit.
The body and the thirteen-speed full row are primitive because C and T
contain 1. No lower-dimensional LRC theorem is needed for this construction.

**Constructive theorem.** Every row S in (9) has M(S)>1/14. Put

    k=floor((35d+18)/42),
    y=(14k+1)/(14d),  x=y/3.                         (10)

Then x and (y+1)/3 are weak-safe phases. The point y is the left endpoint
of a positive body component, with frequency d its unique endpoint owner.
On an explicitly positive initial cell all three tails are active and
have the same owner label 2. Both free sheets yield strict physical witnesses
at interior points of that cell.

**Proof.** Since d is 1 or 5 modulo six, the address is exact:

    y=5/6+5/(21d)   if d=1 mod6,
    y=5/6-2/(21d)   if d=5 mod6.                    (11)

Write delta=y-5/6, so |delta|<=5/(21d). On the virtual coset,

    dy=k+1/14,
    jd y=j/14 mod1 for j=1,2,3,4,
    14b y=14y mod1,  14c y=28y mod1,
    hy=y mod1,  2hy=2y mod1.                        (12)

The last line uses 14d dividing L. Thus the remaining body conditions
reduce exactly to the four frequencies 1,2,14,28 at y. At y=5/6 their
clearances are 1/6,1/3,1/3,1/3. For d>=31 their changes are bounded by
|delta|, 2|delta|, 14|delta| and 28|delta|. The strongest required estimate
is 20/(3d)<11/42, valid already for d>280/11. All four are therefore
strictly above 1/14. In (12), only d itself is at equality. This proves
body safety and the unique endpoint owner.

The three tail phases have the following nearest integers:

    ay=y                 nearest integer 1;
    by=3k+y+3/14         nearest integer 3k+1;
    cy=3k+2y+3/14        nearest integer 3k+2.

At delta=0 the signed errors are -1/6, 1/21 and -5/42; the errors change
by delta, delta and 2delta respectively. They stay strictly inside
(-3/14,3/14) for d>=31. Since a,b are 1 modulo three and c is 2 modulo
three, (1) gives the same owner 2 in all three cases. The two labels 0,1
are free, proving the stated physical endpoint witnesses.

For a fully explicit positive addressed cell take

    eta=min( 1/(14d),
        min_(z in C, z!=d) (||zy||-1/14)/(2z),
        min_(w in T) (3/14-||wy||)/(2w) ) > 0.        (13)

Circle distance is 1-Lipschitz. Thus throughout [y,y+eta] the other body
frequencies remain strictly safe and every tail remains active with its
same nearest integer and owner. The d phase moves from 1/14 into its safe
interval. For every interior point y+t, 0<t<=eta, sheets 0 and 1 are
strictly safe for every physical speed. This proves M(S)>1/14 and exhibits
the actual component address, rather than only an ambient phase statistic.

## 6. What this family defeats, and which entry conditions it does not claim

**Every ordinary tail endpoint coset is completely blocked.** For each
w in T, C contains 14w. At every positive or negative ordinary event phase
(14k+-3)/(14w), the corresponding body phase is integral and hence strictly
unsafe. All three ordinary event deficits are zero. The virtual d endpoint
in (10) succeeds by owner collision even though every tail is active.

**The virtual count also works without the special h residue.** Four
members of C are multiples d with ratios 1,2,3,4 and block no virtual
point. The other six are coprime to d. Thus the independent criterion
(8) proves existence as well. More generally h can be any sufficiently
large integer coprime to the eight base body frequencies; then the same
virtual count closes the row, although formula (10) need not remain its
witness. The special choice h=1+mL restores an exact phase and owner word.
Composite d is allowed: every gcd and its branch multiplicity is retained.

**The necessary numerical filters of THM-4446 can all hold.** Require
h>Q=91^6 and distinguish P={h,2h}. Every base frequency z divides L, so
gcd(h,z)=1 and gcd(2h,z)<=2. Consequently

    kappa(C;P)
      =min_(u in P,z in C\P) max(u/gcd(u,z),z/gcd(u,z))
      =h>Q.                                         (14)

The minimum is attained at (h,1). The primitive ratio of the distinguished
pair is (1,2), whose sum is three, and max T<11 max C. Increasing m makes
these cross heights arbitrarily large while the pack stays gcd one.

These are **necessary numerical filters only**. They do not establish
actual THM-3818 entry, the identity W=V_dec, or the absence of all short
multi-support crossing relations. No such inference is made. The rows
are explicitly safe and may be closed by other inherited mechanisms too.

**There is no alternative large ten-pack clock hidden here.** In S, exactly
ten speeds are divisible by three and their gcd is three: they are 3C.
For a prime p!=3, at most seven speeds are divisible by p. More precisely,
there are seven even speeds; a prime divisor of d divides four speeds; 7
divides at most four; a prime divisor of h divides only 3h and 6h; and any
other prime can occur in at most one of b,c and its paired body multiple,
so divides at most two. Coprimality of h with every base frequency makes
these cases exhaustive. Hence 3C is the only ten-subset with nontrivial
gcd, and every eleven- and twelve-subset is primitive. The family occupies
the exact clock-three signature left by THM-4447, while still being closed
by the explicit owner collision.

**All denominator clocks through fourteen can also be blocked.** Restrict

    d=715(11+168r),  r>=0.                            (15)

Then d is divisible by 5,11,13, is coprime to 42, and is 1 modulo eight.
The physical speeds 12d,3d,42,42b,9d respectively block all grid phases
with denominators in the groups

    {2,3,4,6,10,12}, {5,11,13}, {7,14}, {8}, {9}.

For q=8 use b=3d+1=4 modulo eight, so 8 divides 42b. Blocking means a
speed becomes integral, strictly below the weak target; no endpoint
convention can rescue those clocks. The virtual denominator 42d succeeds.

## 7. Exact controls and what remains open

The standard-library verifier checks all 120 ternary-unit triples with
maximum at most fourteen, at every activity/band boundary and intervening
cell; both signs of the pair map; the sharp primitive sequence at fifty
even N; and exact numerator-one wall counts for d<=80 and c<=100, with
additional full-cover boundaries. It uses literal physical masks, not a
reused aggregate tail measure.

It then checks 84 family controls through d=180 at two m values each,
including large h>91^6, plus three controls in (15). The full virtual codes
and ordinary tail endpoint blockers are replayed on selected controls.
For the three strongest controls every ten-, eleven- and twelve-subset gcd
is checked independently. Rational physical inequalities verify the
endpoint, midpoint and end of each explicit component cell.

Two examples are:

| d | h | T | y | physical x |
|---:|---:|---|---|---|
| 31 | 567879096121 | (1,94,95) | 365/434 | 365/1302 |
| 7865 | 183926030247961 | (1,23596,23597) | 91757/110110 | 91757/330330 |

Both have kappa=h, zero ordinary event deficits for every tail, unique
body endpoint owner d, and literal tail masks ({2},{2},{2}). The second
blocks every denominator from two through fourteen. Its positive cell
length from (13) is 2622/10126047595301492855; the very small size is retained
rather than replaced by a phase-independent body measure.

Reproduction:

    python3 -B 04-computation/overnight9_20260906_lrc_virtual_pair_wall.py
    python3 -B -O 04-computation/overnight9_20260906_lrc_virtual_pair_wall.py

All checks use integers/Fractions and explicit exceptions that remain
active under optimization. They verify the all-parameter proofs above,
not supply an extrapolation from a tail or height census.

The remaining universal entry question is whether an actual primitive
ten-pack produced by the decoder must meet some signed-pair band complement,
virtual wall, inactive-owner cell, or owner-collision cell. The present
family supplies one unbounded actual component construction and a sharp
necessary band; it does not force one of those events for every live row.

**Frozen reproducibility:** 65,423 active gates; normal and optimized
outputs are byte-identical LF. Source SHA256:
`69cbf5285d3c65272e1c4200e5c095e9ea9bb687c1cd926ce95eefb6cc375d48`.
Output SHA256:
`50fd7538753f1495d91e4bf517fda84b9d13b2a54b7752a5c86428b6f82f8e51`.
Semantic payload SHA256:
`0fb1cb24e4cbe69af6c48cb2cac53fe94b95a111c2949ea1dabae211520d2d53`.

**Filing:** root integrated these audited artifacts in the ninth checkpoint;
reproduction paths are relative to the repository root. Earlier outside-worktree
notes describe author provenance, not the present file location.
