# Where to move the Collatz proof frontier: inverse shields, lifted joins, and exact search limits

**Status: PROVED elementary mechanisms; FINITE-EXACT bounded experiments;
OPEN universal certificate coverage and the first-coalescence lift question.**
2026-10-04, `difference-families-oct04`. No literature-priority claim.

This continuation changes the allocation of proof work. Holding the forward
window at six odd steps, increasing inverse depth from 16 to 28 leaves the
same 16 of the inherited 239 seed obligations resolved. A separate complete
audit of every smaller source proves that **only 17** of those obligations
can be resolved at that forward window by *any* inverse depth. Advancing
the forward window to eight, with inverse depth 16, resolves 93; advancing
to sixteen resolves 198. These are exact finite comparisons, not coverage
of arbitrary integers or reductions in a newly rerun compiler's seed count.

The new all-height content is a reset-two inverse shield, a sharp delayed
escape from that shield, an exact lift of a common-future diagram to an
arithmetic family, and a sign-sensitive sibling/inverse-run descent rule.
The search has a proved completeness bound on its infinite inverse fan.

## 1. Inheritance, portfolio, and the changed target

Closest mechanisms:

- [Checked reset switches](checked_switch_phase19_20261004.md), especially
  the unbounded first-reset rule and the unresolved reset-two debt.
- [The concurrent grounded observation union](adaptive_observation_union_20261004.md):
  keep unfinished actual edges, close backward from 1, and reuse every
  available proof. Its 802-edge initial graph grounds 509 of 512 requested
  sources; precisely 32 additional literal edges complete that instance.
- [The sibling grammar](creative_sibling_20260925.md): retain intermediate
  sibling heights and distinguish a smaller obligation from a completed
  certificate. Its forward-closure theorem already explains why pure
  inverse ports alone cannot enlarge that particular generated family.
- [The boundary compiler](collatz_boundary_compiler_20261004.md) and
  [THM-4512, coefficient-descent cylinders](../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md):
  exact words need the final oddness bit; coefficient inequalities have
  their own scope and do not force an arbitrary orbit to enter their domain.

Canonical hostile: 27 has no smaller-source common future until odd step
37. Corrected near miss: a real join to a smaller source need not lift to
a decreasing family; 231 -> 233 supplies an explicit counterexample below.
Least-used sidecars: the first exposure to the union of smaller completed
routes, the two separate word lengths, and the exact mixed binary/ternary
period of a common-future family.

The live board is **grounded graph / original source / inverse depth /
ordered carries / family period / reset debt / sign**. Anchor: extend
checked common-future certificates. Niche: explain inverse-search failure
at the reset-two join. Wildcard: the three-adic sibling resonance and its
sharp rank boundary. META-PATTERNS used: search the statement before the
method; correct the object before sharpening the technique; keep the
missing second coordinate; retain selected guards. No new meta-pattern is
promoted from this one thread.

## 2. Every checked join has an exact arithmetic-family lift

Write `U(n)=oddpart(3n+1)` on positive odd integers. A positive valuation
word w of length r and total A has carry B and affine map

    U_w(n)=(3^r n+B)/2^A.

Suppose actual words w,v satisfy

    U_w(n0)=U_v(m0),   0<m0<n0,
    lengths r,s,      total valuations A,D.

Empty words are allowed and have length, cost, and carry zero. Define

    P=2^(A+1) 3^max(s-r,0),
    Q=2^(D+1) 3^max(r-s,0).                            (1)

**PROVED family-lift theorem.** For every integer t>=0, the same exact
words give

    U_w(n0+Pt)=U_v(m0+Qt).

If `P>=Q`, then `0<m0+Qt<n0+Pt` for all t>=0, so the diagram transports
a supplied certificate for a strictly smaller source at every height.
It preserves that supplied basin on either signed sheet. It does not
provide the child certificate by itself.

**Proof.** The two exact word guards are fixed by increments divisible
by `2^(A+1)` and `2^(D+1)`. Their endpoint increments are both
`2*3^max(r,s)*t`. Positivity and the rank comparison follow directly.
Conversely, requiring increments to preserve both exact words gives
`3^r k=3^s l`; its primitive positive solution produces P,Q. Thus (1)
is the least simultaneous positive period for these exact words. QED.

The ratio and intercept of the dependency are

    m = lambda*n+beta,
    lambda=2^(D-A)3^(r-s),   beta=m0-lambda*n0.          (2)

Different word lengths naturally introduce powers of 3 into the source
period. Purely dyadic family bookkeeping would lose these joins.

Two concrete families, independently replayed:

    111+8748t  --(1)--> J <--(1,1,2,1,1,1,2,3)-- 103+8192t;
    283+1296t --(1,2)--> K <--(1,1,1,1,3,2)---- 223+1024t.

Their slopes are `2048/2187` and `64/81`. The denominator gaps are
`139=3^7-2^11` and `17=3^4-2^6`. These are actual coefficient differences
of the displayed proof transports. The first is also the inherited -17
cycle's coefficient gap, but the carry and the word are different; the
gap alone does not identify cycles.

**Hostile.** The actual 17-step word from 231 to 233 is

    (1,1,2,1,1,2,6,1,1,1,1,2,2,1,2,1,1), total27.

As a diagram certifying source 233 from smaller source 231, its inverse
slope is `2^27/3^17>1`. The t=0 diagram is valid; its least-period lift
eventually has child larger than source. Retain the slope test even after
checking one successful join.

## 3. Why the first reset-two join resists inverse switching

For sign `epsilon=+1` or `-1`, let
`U_epsilon(n)=oddpart(3n+epsilon)`. Suppose q>=2, t is positive odd,

    n=2^q t-epsilon,
    3^q t-epsilon=2 mod4.

The actual word `(1^(q-1),2)` ends at

    J=(3^q t-epsilon)/2.                              (3)

**PROVED inverse shield.** For every 1<=d<=q, the unique smallest positive
odd d-step predecessor of J is

    x_d=2^d 3^(q-d)t-epsilon.                         (4)

In particular it is at least n, with equality only at d=q. Therefore no
inverse word of depth at most q gives a smaller dependency at this join.
There is no hidden bound on the individual inverse exponents.

**Proof.** At J the least permissible inverse exponent is 2. The resulting
integer is x1. Each of the next q-1 steps permits exponent 1, yielding
(4). Every inverse map `(2^a x-epsilon)/3` increases with both a and x.
All branches use at least that first exponent and at least one afterward,
so none can lie below the displayed branch. Strictness gives uniqueness.
Finally `x_d+epsilon=(3/2)^(q-d)(n+epsilon)`. QED.

There is a stronger, useful boundary.

**PROVED delayed escape.** Any smaller predecessor of J at depth at most
q+3 must first retrace the entire canonical reverse path back to n.
Equivalently, its forward path passes through n before J. Such a witness
is an ordinary ancestor of n; it introduces no new branch around n at
this early join. A genuinely different branch can first beat n at q+4,
and this bound is attained.

**Proof.** Suppose a reverse path first differs at position j<=q. The
allowed exponents have fixed parity, so its exponent is at least two
larger than the canonical one. Since

    E_(a+2,epsilon)(x)=4 E_(a,epsilon)(x)+epsilon,

after q+ell total steps its endpoint z obeys

    z+epsilon >= 4(2/3)^ell(n+epsilon)
                 -2 epsilon (2/3)^(q+ell-j).           (5)

For the plus sheet and ell<=3 the worst comparison is ell=3,j=q;
the difference from n+1 is at least `(5(n+1)-16)/27>0` for n>=3.
The minus sheet has a positive correction instead and n-1>0. Shallower
depths are covered by (4). Thus the changed branch cannot be smaller.
At ell=4 the leading factor is `64/81<1`, so this obstruction ends.
The example n=283,q=2,J=319, with predecessor 223 and reverse depth 6,
attains it: its path does not pass through 283. QED.

This shield is sign-blind and applies at the minus basins too. It explains
a search obstruction, not why the plus sheet has only one positive basin.
It also does not forbid an ordinary ancestor of n or a later forward join.

## 4. Climb a sibling ladder, then undo exponent-one steps

Let `S(n)=4n+1`, so `U(S(n))=U(n)`. Fix k>=0 and ell>=1. If

    3^ell divides S^k(n)+1,
    m=2^ell (S^k(n)+1)/3^ell-1,                        (6)

then m is positive odd and has ell actual exponent-one steps to S^k(n).
It therefore shares the future U(n) with n.

**PROVED sharp rank criterion.** For every positive integral instance (6),

    m<n iff 3^ell>2^(ell+2k).                         (7)

To see the sharpness, put `t=(S^k(n)+1)/3^ell`, a positive even integer,
and `D=3^ell-2^(ell+2k)`. Then

    n-m=[D*t+2(4^k-1)/3]/4^k.

If D>0 the difference is positive. If D<0 it is less than 1 and is an
integer, so it cannot be positive. Equality of the two powers is impossible.
This uses integrality, not an asymptotic slope approximation.

Let ell(k) be the least positive exponent satisfying (7). The source
guard is one complete ternary residue class:

    4^k(3n+1)+2=0 mod3^(ell(k)+1).                    (8)

| k | ell(k) | source class | least positive odd source | child |
|---:|---:|---|---:|---:|
| 0 | 1 | 2 mod3 | 5 | 3 |
| 1 | 4 | 40 mod81 | 121 | 95 |
| 2 | 7 | 273 mod2187 | 273 | 255 |
| 3 | 11 | 94109 mod177147 | 94109 | 69631 |

For k=1 the source class is already inside the pure inverse rule n=4 mod9;
the new diagram supplies another dependency, not new local applicability
against that rule. For k=2, every source is divisible by3 and hence has
**no odd predecessor at all**, while (6) still supplies a smaller common-
future source. This proves an infinite local extension beyond *all pure
inverse-word rules*. It does not claim an extension beyond every earlier
forward/common-future grammar.

A subfamily lies entirely in the hard reset-two branch:

    n=4647+69984h,   m=4351+65536h,   h>=0.

Here n has q=3 and first reset2. The paths use words `(1)` and
`(1,1,1,1,1,1,1,5)` to the same future. The coefficient gap is139.
This ties the delayed inverse escape to the same exact `3/2` threshold
as [the virtual contraction ladders](virtual_contraction_ladders_20261004.md).
Their direction of certificate transport and source guards remain distinct.

**Strict extension of a specified short bank.** At h=2, n=144615 and
the sibling rule supplies child135423. Complete search through the first
six forward positions and all inverse words of length at most6 finds no
smaller child with a nonexpanding lift. Already at the first forward
position, inverse length7 supplies the alternative child101567, with word
`(1,1,1,1,1,2,3)` and join216923. Its entire family is

    n=1731+2916t,   m=1215+2048t,   t>=0.

The slope is512/729 and the coefficient gap is217=7*31. The source144615
is t=49. This is a strict local extension of the complete F=6,D=6 bank,
not a claim that the children are automatically grounded. Retain both
children when a supplied certificate may make either route preferable.

**Signed hostile.** On `3n-1`, the inverse of `(1,2)` is `(8n+5)/9`.
At n=5 it equals5 despite coefficient8/9<1. The seven-step -17 pattern
likewise has inverse coefficient2048/2187 but fixes17. The positive carry
sign in (6)--(7) must not be discarded when discussing basin uniqueness.

## 5. Search every possible inverse exponent at a fixed depth

For positive odd z, every d-step predecessor satisfies

    predecessor >= (2/3)^d(z+1)-1.                   (9)

This follows by iterating `(2^a z-1)/3 >= (2z-1)/3`. A search for a
predecessor below the **original threshold N** can prune the entire
remaining depth d whenever

    2^d(z+1)>=3^d(N+1).                              (10)

For a candidate next inverse exponent a, it can omit that and all larger
exponents when

    2^(d-1)(2^a z+2)>=3^d(N+1).                      (11)

Thus a depth cap gives a complete finite search with no arbitrary exponent
cap. Divisibility chooses a's parity; a multiple of3 has no predecessors.
The implementation includes the empty inverse word when the forward
endpoint is already smaller.

Repeated positive vertices can also be removed when searching for a
nonexpanding lift. A positive 3n+1 cycle with length r and cost A satisfies
`(2^A-3^r)n=B>0`; its inverse coefficient is therefore greater than1.
Deleting such a loop shortens the witness and improves its inverse
coefficient. This reasoning does not assume that the only positive cycle
is1. It would not justify the same pruning under 3n-1.

The complete algorithm was independently checked on every smaller candidate
source for 10,560 combinations of threshold, endpoint, depth and coefficient
allowance. The independent path simply walks those candidate sources
forward; it does not use (9)--(11).

## 6. Frozen benchmark: which direction is productive?

Input: the exact 239 seeds saved by the inherited compiler on all odd
sources <=10000. Their list is hashed in the new JSON. Each experiment
allows at most F actual forward steps and D reverse steps; it retains the
nonexpanding family criterion. All valuations are included by (11).

| Forward cap F | Reverse D=6 | D=8 | D=12 | D=16 |
|---:|---:|---:|---:|---:|
| 6 | 0 | 7 | 12 | 16 |
| 8 | 78 | 84 | 89 | 93 |
| 12 | 159 | 164 | 166 | 169 |
| 16 | 192 | 195 | 196 | 198 |

Entries count original seed obligations with a newly found smaller-source
join, including direct descent as an empty reverse word. They do not say
the source has a short actual route to1: the smaller source's certificate
is still needed. These counts also do not represent a rerun of the old
online learner, whose subsequent learned bank could change.

At F=6, D=20,24,28 still resolve16. The inverse-node visits rise from
91,321 at D=16 to1,568,023 at D=28. At F=8,D=6 there are78 resolutions
and9,109 inverse-node visits. This comparison shows where the present
search is spending its work; it does not equate reverse-node visits with
CPU time or charge the unequal forward budgets as equal.

**Complete all-depth audit.** Independently replay every odd source through
10000 to1 and record which smaller sources visit each integer, with their
ordered costs. For each seed and its first six forward positions, inspect
all those smaller trajectories. There are only finitely many smaller
sources, and their entire completed paths are known, so this audit covers
*every possible inverse depth*, including depths not searched above.

Exactly17 seeds have any such smaller ancestor, and all17 admit a
nonexpanding family. The other222 have no smaller ancestor at any of those
positions, even if the family-slope restriction is dropped. Full routes
are used only for this saturation audit, not as inputs to the bounded search.

The seventeenth source, missed even at inverse depth28, is1567. Its sixth
forward endpoint is4465; child1055 reaches that endpoint in48 odd steps.
If child1055 already has a completed certificate, the shared graph can
reuse its suffix directly. It need not rediscover48 inverse edges. The
17-source ceiling concerns witnesses starting below the original source;
an independently certified larger source is also legitimate for graph
grounding, but does not supply this particular smaller-child induction.

The exact forward-window ceiling is:

| F | 1 | 4 | 6 | 8 | 12 | 16 | 24 | 32 | 48 | 64 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| seeds with a smaller-source common future | 6 | 13 | 17 | 93 | 169 | 198 | 220 | 228 | 238 | 239 |

The last first coalescence is source703 at odd step51. Source27 first
coalesces at step37, also its first actual descent. These are finite facts
about the supplied universe, not uniform bounds.

## 7. A further signal, with the correct quantifiers

For each odd n through100000, let G_<n be the union of the completed
trajectories of smaller positive odd sources. Define kappa as the first
j>=0 for which U^j(n) lies in G_<n. In this finite universe all paths are
explicitly completed. Among the49,999 nonroot inputs:

- 23,400 already belong to G_<n, so kappa=0;
- 26,599 are new sources; 2,419 of them coalesce strictly before their
  first actual descent;
- every first coalescence has some smaller witness whose diagram admits
  a nonexpanding all-height lift;
- the largest saving among genuinely new sources is70 odd steps:
  source75007 coalesces at8 instead of descending at78, using child37503.
  This is an instance of the inherited reset switch.

**OPEN secondary question.** Assuming smaller sources have completed
certificates and a first coalescence exists, must that first coalescence
always admit a nonexpanding lift? The finite signal supports this
statement, but the arbitrary-pair hostile231 ->233 shows why it needs the
first-coalescence condition and a choice among all smaller witnesses.

Even a proof would not establish that kappa is finite for an arbitrary
source. That is the remaining coverage issue. It should not be quietly
substituted for the finite signal or for the distinct stopping-time equality
question in THM-4512.

## 8. Reformulated proof program

The object worth retaining is a **grounded graph of parametric proof
obligations**. A node carries the original integer N, an actual current x,
the exact word from N to x, and candidate common-future diagrams with their
two costs, periods, slopes and child proofs. A finite phase such as the
[joined mod36 observer](triplet_crt_223_233_20261004.md) annotates this
object; it does not replace N or the ordered carries.

The operational changes justified by this session are:

1. Preserve unfinished observations and all available certified children,
   as in the concurrent grounded-union work. An unresolved cycle of proof
   obligations cannot certify itself.
2. Lift successful joins by (1), retaining mixed binary/ternary source
   periods. Store the exact family theorem alongside a ground instance.
3. At a reset-two join, apply the shield before deepening inverse search.
   Before depth q+4, new branches cannot beat the original source; ordinary
   ancestors can be queried separately.
4. Use the complete inverse bound (11), then advance the forward frontier
   when the reverse search is saturated. Keep the original N when comparing
   costs or declaring descent.
5. Replace recurrent long forward blocks by the existing unbounded
   rational-shadow macros, with their actual guards and source-relative
   payment inequalities. A macro count is not a bound on orbit length.

**Precise remaining theorem sought:** a well-founded symbolic rule for
the unresolved forward frontier, covering every positive source while
retaining exact guards and child identity. The reset shield identifies
where a fixed inverse bank cannot supply it. The family lifts and sibling
rules enlarge the available transformations, but neither proves that every
frontier enters one of their domains. No global decreasing meta-rank has
been proved in this session.

This reformulation keeps the target on universal grounded coverage, with
finite modular geometry serving as a coordinate system and with observed
plateaus giving rigorous stopping reasons for particular search directions.

## Reproduction and audit boundaries

The continuation [Binary debt and ternary sibling addresses](collatz_binary_ternary_guard_fusion_20261004.md)
integrates incoming commits9b320b7701 and7c15ea9559. Their guarded half-source
rewrites complement the joins here. The continuation extends their eight-row
debt bank to sixteen, proves that sibling height is a ternary isometric
address, and uses positivity to establish a natural-density product law
for the binary and ternary necessary tests. The combined direct tests
still leave positive density; the source-relative forward frontier remains
the target. This does not change the frozen239-source counts above.

Run `python3 -B 04-computation/experiments/collatz_join_shields_and_lifts_20261004.py`
and repeat with `python3 -O -B`. The [script](../../04-computation/experiments/collatz_join_shields_and_lifts_20261004.py),
[JSON](collatz_join_shields_and_lifts_20261004.json), and
[stdout](collatz_join_shields_and_lifts_20261004.out) specify the frozen
seeds, every search cap, successful lifted witnesses, exact finite-universe
saturation and first-coalescence data. The checks use integer arithmetic
and Fractions, with no optimization-sensitive assertions.

Positive controls include every returned word's literal replay, the
independent affine identity, very large positive family parameters, and
the sharp q+4 escape. Hostiles include the expanding231/233 join, the
minus-cycle minima5 and17, and the inverse-depth plateau. The all-depth
finite audit and the first-coalescence probe use fully completed routes;
they are clearly separated from the bounded search and supply no hidden
oracle to it. No external literature theorem is needed for the new proofs.
