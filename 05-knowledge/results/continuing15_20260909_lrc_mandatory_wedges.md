# Mandatory actual wedges: three more zero-budget clocks and a synchronized stopping tree

**Status: PROVED ANALYTICALLY + FINITE-EXACT; INDEPENDENTLY AUDITED.**
The complete inherited zero-budget slice is surveyed, and a bounded actual-ratio
decoder excludes the clocks **5880, 6804, and 7056**. The necessary clock array
therefore decreases from 6,889 to **6,886**, with maximum 11934. Its zero-budget
slice decreases from 356 to **353**, with maximum 6552. These are consequences
only in the inherited **primitive thirteen-speed, selected-six, actual
connected-complement** setting. General LRC(14) remains open.

The new structural coordinate is the ratio on a **forced bridge** of the
complete possible zero graph. Two such bridges meeting at a vertex must be
realized simultaneously. Their complete ratio banks can conflict with global
gcd depths, duplicate an actual speed, or force a positive overlap between
their endpoints. The endpoints need not themselves be a strict-atlas edge.

The declared extension stops at 5376. Its one remaining word has an explicit
primitive seven-tail realization whose *entire actual strict graph* is a tree
with every edge zero at one common phase. Nevertheless 1,600 grid points remain
safe. This is a precise hostile to mistaking synchronized native edges for an
unsafe covering; it is not a counterexample to LRC.

The [independent referee](continuing15_20260909_lrc_independent_audit.md) accepts this result;
this filed copy changes only audit status and adds this link.

## 1. Inheritance, scope, and retained coordinates

The closest proved mechanism is the complete native-zero classification in
[continuing14](continuing14_20260908_lrc_zero_classification.md), with its
[independent referee](continuing14_20260908_lrc_independent_audit.md). That result
retains all 126 proper nonempty subset profiles, the complete strict atlas,
actual quotient clocks, strict walls, marginal-attainment chambers, all seven
position vertices, and the distinct-speed slot bound. Its necessary array has
raw supplier certificate SHA

    863d33accf19a5515b3ecd2048aac35d77f30f4e004371c73f5760d345bbb240

and semantic array SHA

    b603ca0186c1ebaa2318782348a7446a16f529543147eb53062f1ba360e8e57b.

The ratio/depth sidecar is
[THM-4106, pair-owner reciprocity and valuation-tree decoder](../../01-canon/theorems/THM-4106-lrc-pair-owner-reciprocity-and-valuation-tree-decoder.md),
as developed in
[continuing9's ratio-tree report](continuing9_20260907_lrc_ratio_tree.md).
The local-depth hostile has clock 2 and requested margins `(1,2,1)`:
ratios `u_1/u_2=1/4` and `u_3/u_2=3/2` propagate to primitive
`(1,4,6)`, whose actual margins are `(1,2,2)`. Independently compatible
edges need not have simultaneously compatible absolute depths.

The corrected near miss is the 7200 native tree: native zeros, including a
common phase, do not imply coverage because non-atlas intersections were
discarded. The present 5376 object pays the same missing coordinate on a new
surviving full-profile word. The least-used sidecar is the **bridge cut** of
the complete possible graph: it turns a possible edge into a logically
required actual edge. We never infer actual adjacency merely from membership
in the possible graph.

The concept board for this bounded session was: marginal saturation; actual
bridge cuts; oriented ratio banks; primitive depth synchronization; forced
non-atlas endpoint intersections; and common-phase covering. The map from an
actual tail tuple to its possible graph preserves each actual strict edge but
forgets ratios and phases. The mandatory-wedge decoder restores two ratios
and their exact common primitive scale. The final physical hostile shows what
still disappears after that restoration: intersections on the other actual
pairs and the union of all seven danger sets.

No arbitrary new clock census, partial minimum-spanning-tree pruning,
scale analogy, or weakening of actual connectedness is used.

## 2. Native definitions and why zero budget forces chamber-zero edges

For a clock `t` and an actual positive tail speed `u_i`, write
`d_i=gcd(t,u_i)`. The inherited domain is the complete divisor alphabet
`D_t={d in V:d|t}`, where V is the pinned 42-element six-profile alphabet.
For each sorted seven-word S over D_t, require `gcd(S)=1` and, for every
proper nonempty position subset I, set

    c = gcd(S_i : i in I),
    w = sort(gcd(c,S_j) : j not in I),

and require `(c,w)` in the full inherited profile bank of size `7-|I|`.
All 126 subsets are tested. The source enumerates complete seven-multisets
without rejecting partial words by a heuristic bound.

Every surviving zero-budget clock is divisible by seven. No member of V has
a factor seven. Consequently `7d|t` for every d in D_t. On a t-point grid,
each individual strict danger set `||u_i x||<1/14` has at most `t/7` points.
An unsafe covering by seven tails would attain every marginal bound and make
all pair intersections zero. At an individual endpoint wall the corresponding
danger count loses d points, so such an unsafe phase is in an attainment
chamber for each actual tail. This justifies using **chamber-zero** native
input banks. The endpoint contradictions below use the stronger **all-phase**
minimum, so they require no new endpoint-attainment assumption.

An actual strict edge has coprime reduced ratio p:q, `p<q`, with `p+q<=356`
and each prime divisor of p+q congruent to 2 modulo 3 with exponent at most
two. The complete atlas contains 5,855 pairs. If endpoint margins are a,b,
then `e=gcd(a,b)` is the actual pair sheet, the quotient clock is `n=t/e`,
and the two possible orientations satisfy

    (e*gcd(n,p), e*gcd(n,q)) = (a,b) or (b,a).

For each ordered pair of margins define `B_t(a,b)` to be every oriented
ratio `u_a/u_b` from this full atlas whose exact intersection minimum over
attainment chambers is zero. When a=b, the two reciprocal orientations are
distinct and both retained. Their counts reproduce the pinned continuing14
oriented pair table. No bank is inferred from a selected minimizing owner.

The possible zero graph G has the seven **position** vertices, with edge ij
exactly when `B_t(d_i,d_j)` is nonempty. The actual connected strict graph
at an unsafe phase is a connected spanning subgraph of G. Every bridge of G
therefore belongs to the actual graph. If two bridges meet, their three
endpoints form a **mandatory actual wedge**. This implication fails for a
merely possible nonbridge edge and is not used there.

The inherited distinct-speed slot filter is also retained: if margin a has
multiplicity m_a and `B_t(a,a)` is empty, then

    m_a <= sum_(b != a) m_b * |B_t(a,b)|.

All original seven-word residuals are recomputed from the complete graph and
this filter. This is a necessary filter, not a realization theorem.

## 3. Exact primitive-wedge decoder

Consider a mandatory wedge with margins `(a,b,c)` and center b. Put

    g=gcd(a,b,c),  T=t/g,  (A,B,C)=(a/g,b/g,c/g).

Enumerate the **entire Cartesian product**

    r in B_T(A,B),  s in B_T(C,B).

The quotient clocks and reduced margins in these two banks are exactly those
in the unnormalized banks: dividing the three margins and t by g loses no
ratio. Propagate `(r,1,s)` to its primitive positive integer triple v.

**Exact local scale criterion.** There exist positive integer values on the
three vertices with ratios r,s and original margins `(a,b,c)` if and only if

    (gcd(T,v_1),gcd(T,v_2),gcd(T,v_3))=(A,B,C).             (1)

Indeed every realizing triple is H v. Because `gcd(v)=1`, the gcd of its
three original clock margins is `gcd(t,H)=g`. Write H=gH'. Then
`gcd(T,H')=1`, and its margins are `g*gcd(T,v_i)`, proving necessity.
Conversely H=g gives a realization whenever (1) holds. This proof permits
an arbitrary common gcd of the original three tails; it does not silently
assume those tails are primitive on the original clock.

If two entries of v coincide, distinct actual speeds are impossible.
If (1) fails, no global rescaling repairs that local depth mismatch.
Otherwise the endpoint ratio is exactly r/s. Reduce it to coprime p:q;
there is no condition that p:q belong to the strict atlas. The endpoint
sheet is `gcd(a,c)` and its actual quotient clock is `t/gcd(a,c)`.
Let m be the exact all-phase intersection minimum for p:q on this quotient.
If m>0, those actual endpoints overlap at every phase by at least
`gcd(a,c)*m`, contradicting zero-budget unsafe coverage.

Thus a mandatory-wedge signature closes if every ratio pair is a duplicate,
a depth conflict, or a positive all-phase endpoint intersection. A pair with
correct depths, distinct entries, and zero endpoint minimum is retained as
an honest residual. Even all such local tests passing does not imply the
seven-tail tuple is realizable, nor that its pairwise phases agree.

The capacity engine constructs the strict circle intersection at denominator
`L=14pq`, joins the interval through the origin, and sweeps **all** start/end
events and exact strict-wall holes on the n-point grid. On each interval
`(a,b)` its count is

    ceil((n*b-phase)/L) - floor((n*a-phase)/L) - 1.

The finite event sweep computes both the minimum over all phases and the
minimum over open chambers. Independent direct interval evaluations at every
wall and one point per chamber, followed by literal native grid counts,
are retained as controls. No expected-overlap approximation enters (1) or
the endpoint contradiction.

## 4. Declared finite universe and stopping protocol

The input is **all 356 remaining zero-budget clocks**, not a sampled subset.
Their full alphabets use 121 of the 203 continuing14 domains. The native
companion enumerates **7,105,596** seven-multisets and obtains **291,811**
accepted distinct domain words. Every regenerated whole bank has exactly
the inherited raw SHA. Reapplying graph connectedness and the complete
oriented-slot filter gives **47,568** surviving clock-labelled words.

First survey every mandatory wedge with **at least one singleton oriented
ratio bank**. This produces 94 clock-labelled signatures and 53 distinct
normalized signatures. Their complete arithmetic universe contains 1,980
ratio products. Exactly 33 normalized signatures close, excluding 1,087
original residual words and completely closing 5880 and 7056.

Next make one cheap **topological** boundary survey, on that same entire
356-clock universe: retain mandatory wedges whose product of bank sizes is
at most 256. This yields 561 clock-labelled signatures, 304 normalized
signatures, and 5,424 eligible words on 77 clocks. The catalogue is complete;
its arithmetic products are **not** all evaluated. In particular no claim
is made to have examined the 26,729 products of the full 304-key catalogue.

Among clocks not already closed by the singleton stage, exactly 6804, 5376,
and 5124 have every residual word eligible for this bounded survey. Declare
the descending order

    [6804, 5376, 5124],

evaluate every bounded signature at each successive clock, and stop at the
first clock with a complete arithmetic residual. Clock 6804 closes. Clock
5376 retains one word, so the procedure stops there. **5124 is untested
arithmetically**; its participation in the cheap topology catalogue is not
a claim that its ratio products have been consumed.

The union actually evaluated is exactly **64 normalized signatures and
3,055 complete ratio pairs**:

| Exact classification | Ratio pairs |
|---|---:|
| Positive all-phase endpoint minimum | 2,158 |
| Impossible simultaneous primitive depths | 630 |
| Duplicate actual speed values | 9 |
| Retained arithmetic residual | 258 |

Every one of the 2,158 positive cases has positive **all-phase** minimum;
there is no extra chamber-only composite positive in this universe.
Exactly 37 signatures close. Applying these proved signatures to every
eligible occurrence in the already complete topology catalogue removes
**1,119** full-profile residual words. Exactly the three clocks
`[5880,6804,7056]` have no word left. Other clocks retain their explicit
remaining word counts; no partial deletion is promoted to a clock theorem.

## 5. Three distinct arithmetic explanations

### 7056: a forced off-atlas endpoint credit

The three continuing14 residual words all have the unique margin-4 leaf
joined only to the unique margin24, and the unique crossing from the
`{9,9,18}` group to its complement is the edge18--24. Therefore the
actual wedge4--24--18 is mandatory in all three words. Its complete banks are

    u4/u24 = 289/30,
    u18/u24 in {111/244,117/236,135/212,219/100,
                 225/92,243/68,291/4}.

The seven endpoint ratios u4/u18 are

    35258/1665, 34102/1755, 30634/2025, 2890/657,
    13294/3375, 9826/3645, 578/4365.

Each has exact all-phase minimum **55** on the actual quotient clock3528.
The endpoint sheet is2, so every product forces **110** overlapping points
on the7056 grid. E=0 makes that uniformly impossible. The proof uses the
entire seven-product bank, not a convenient selected path ratio.

### 5880: a common-depth impossibility

The mandatory normalized signature `(5880;4,20,5)` has exactly two complete
ratio products. Both violate (1). This is a scale synchronization closure,
independent of estimating an endpoint intersection. The certificate retains
both rational edge labels, the propagated primitive triples, and every
actual gcd margin so the obstruction can be checked prime by prime.

### 6804: bounded banks without a singleton shortcut

There are two distinct normalized signatures at the declared first boundary
clock: `(756;1,3,2)` has114 products, including43 residuals, whereas
`(2268;4,6,9)` has171 products, all with positive all-phase endpoint minima.
The latter mandatory wedge accounts for the full clock closure. Retaining
the former's residuals is necessary for honest local-versus-global scope.

## 6. Exact stopping state and a full physical hostile

At5376 the complete bounded catalogue has nine normalized signatures and
790 ratio pairs. After applying all37 signatures proved anywhere in the
declared arithmetic universe, exactly one original residual word remains:

    (3,3,6,8,12,16,16).                                  (2)

Its two margin-3 leaves must both meet the unique margin6. Their complete
oriented bank is `{289/58,295/22}`. Reusing a ratio duplicates a speed;
the two opposite choices propagate to `(3179,638,8555)` and its reversed
endpoints, with correct1792-clock margins `(1,2,1)` and zero all-phase
endpoint minimum. The complete four-product bank consists of two duplicates
and two retained residuals. In particular the prior singleton-slot mechanism
cannot exclude two distinct leaves from these two distinct actual slots.

A further **literal tail-only hostile**, without any additional clock
enumeration, is

    U=(157007631,422522895,31510182,129461800,
       36729660,8763568,5333680).                          (3)

These are seven distinct positive integers with gcd1 and actual5376-clock
margins exactly (2). Computing all21 actual reduced pair ratios shows that
the **entire** actual strict atlas graph consists of the following six edges;
indices are zero-based in (3):

| Endpoints | Coprime unordered ratio |
|---|---|
| 0,2 | 58:289 |
| 1,2 | 22:295 |
| 2,4 | 163:190 |
| 3,5 | 22:325 |
| 4,5 | 68:285 |
| 4,6 | 44:303 |

It is connected and acyclic. Thus all cycle conditions hold. At the single
original grid phase `alpha=1/2`, meaning points `(j+1/2)/5376`, **each**
individual danger set has768 points and **every** one of these six actual
edge intersections is empty. The test computes actual huge integer residues
on all5376 points; it does not separately minimize each edge and infer a
common phase. Yet their union leaves exactly **1600 safe points**.

The certificate retains all21 actual pair ratios and intersections, the
full actual strict graph, and the hash of each of the seven literal danger
sets. Positive intersections on actual non-atlas pairs account for the
missing union information. This is a safe tail-grid object satisfying every
listed graph, ratio, common-depth, distinctness, cycle and native-phase
condition. It is **not** an unsafe row, a complete thirteen-speed instance,
or an exclusion of5376. No selected-six phase realization is asserted for
this standalone tail-only control.

The bounded computation stops here. An extension must retain a further
coordinate, such as forced longer-path endpoint intersections on the
remaining tree banks, or a certified lower bound for overlap away from the
strict atlas coupled to actual global ratios. Merely completing more
independent zero-edge or cycle tests does not resolve the demonstrated
coverage loss.

## 7. Replay and audit boundary

The [Python source](../../04-computation/continuing15_20260909_lrc_mandatory_wedges.py)
and [native profile companion](../../04-computation/continuing15_20260909_lrc_mandatory_wedges.cpp)
are self-contained consumers of the two pinned repository inputs. They import
no earlier mathematical producer and do not depend on private exploratory
word files. Run from the repository root:

    python 04-computation/continuing15_20260909_lrc_mandatory_wedges.py
    python -O 04-computation/continuing15_20260909_lrc_mandatory_wedges.py

Outside the repository, run the equivalent absolute source path from the
repository working directory. A private temporary native build regenerates
every requested complete word bank. On relocation under `04-computation`,
the certificate is written to `05-knowledge/results`; outside it is written
beside the source.

The [certificate](continuing15_20260909_lrc_mandatory_wedges_certificate.json)
contains the full declared clock array, all121 domain pins, the complete
singleton and bounded topology catalogues, every evaluated oriented native
bank, every one of the3055 ratio products, all37 proved signatures,
1119 removed-word witnesses with **position** wedges, every remaining clock
count, the exact final array, the declared stopping order, and the full
physical hostile. The final array semantic SHA is

    acbcbd2e63f01572e3af4789b058ce71a557bcc6b8317be8af8b21b1ff45f2cd.

The
[transcript](continuing15_20260909_lrc_mandatory_wedges.out)
records the universe and gates. All checks raise explicitly and remain active
under `-O`; normal and optimized stdout/certificate byte comparisons are
required before the packet is frozen. Gate counts and final array hash are
reported by the frozen transcript and certificate, not inferred from a
prototype. Each mode passes **77,404 always-active gates**: 67,007 Python
gates and 10,397 native gates. The 91 complete requested oriented banks
examine 57,438 compatible native atlas pairs with no positive-lower-bound
skip. These native inputs and all 2,408 distinct computed endpoint
capacities are separate from the full 3,055 ratio-product records.

This session promotes a reusable sufficient obstruction, not a necessary
and sufficient decoder for actual tails. The all356-clock **topology** census,
the64-key **arithmetic** universe, and the single physical **phase** control
are three separately declared finite scopes.
