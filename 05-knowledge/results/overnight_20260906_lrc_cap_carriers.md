# A sharp gap for LRC carrier dictionaries with no three collinear points

**Status: PROVED ANALYTICALLY + FINITE-EXACT finite base + INDEPENDENTLY
AUDITED.** If the **full** live carrier dictionary of a
primitive sorted odd ternary-unit triple is noncollinear and contains no
three collinear points, then

```text
min_i S_i <= 204/5957 < 6/77.                            (1)
```

The first constant is sharp, uniquely attained at `(23,29,37)`. This is an
all-height closure of a geometric class of genuinely multi-direction
carriers; it does not repeat the signed norm-four one-ray closures.
LRC(14), physical entry, and the general projection inequality remain
**OPEN**. THM-4423 was an unproved reservation at the inherited commit and
is not a dependency.

## Inheritance and the actual transfer from no-three-in-line

The closest proved mechanism is
[THM-4422 / projection-deficit-and-beatty-row-reduction](../../01-canon/theorems/THM-4422-lrc14-projection-deficit-and-beatty-row-reduction.md):
it retains the full carrier support and turns the target into a maximum
of boundary deficits. The canonical hostile to global pair selection is
`(1,19,79)`; the first count-dense multi-direction dictionary in its finite
atlas is the six-point `A_2` configuration at `(19,23,29)`. The corrected
near miss is a fixed convex average of the three projections, explicitly
refuted in THM-4422. The least-used sidecar is the pair of **owner colors**
on the integral rank-two lattice; parity alone does not preserve liveness.

The primary seed was Guy and Kelly's
[*The No-Three-In-Line Problem*, 1968](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/B126DA7E4957722BAC70AC7B7F6E1FA2/S0008439500056770a.pdf/the-no-three-in-line-problem.pdf).
Their paper counts grid collinear triples by primitive directions and then
uses a stated independence assumption for a conjectural extremal heuristic.
The full primary PDF was read. Neither the heuristic nor a random-model
conclusion is imported here. The transferred operation is to inspect the
**complete collinearity hypergraph**, rather than replace the configuration
by its number of directions or a selected cap.

This transfer changes the ambient problem: an arbitrary no-three-in-line
subset of a square can have many points, whereas our dictionary contains
**every** lattice point of the two admitted owner colors in a convex open
region. That completeness makes integral midpoints compulsory.

The concept board and its completed update are:

| Concept | Question | Result / surviving sidecar |
|---|---|---|
| Raw projection roofs | Can one global pair be bounded? | A complete cap dictionary has a sharp strict gap |
| Owner-color lattice | Does a midpoint stay live? | Yes, for equal owner color and equal parity |
| Primitive directions / collinearity | Are directions enough? | No; the exact midpoint triple is the obstruction |
| Boundary deficits | Can the new triple count pay the deficit? | Still open; position inside the weighted strip matters |
| Convex completeness | Can a convenient subset be used? | No; a selected cap can erase a mandatory carrier |

**Completed continuation:** [the midpoint-deficit analysis](overnight_20260906_lrc_midpoint_deficit.md)
now proves an optimal endpoint-disjoint curvature payment and its exact
affine residual. A complete infinite two-line family has arbitrarily many
compulsory progressions but zero curvature. Thus the weighted question in
the board requires a boundary offset; the hard residual problem remains open.

The anchor was the residual projection selector; the niche was owner-color
and parity incidence; the wildcard was the complete collinearity hypergraph.
The cheap hostile was the four collinear carriers of `(1,11,23)`: selecting
two endpoints produces a cap but deletes precisely the relevant information.
The META-PATTERNS cards used were “Search the statement before the method,”
“Correct the object before sharpening the technique,” and the quotient-loss
discipline. No new meta-pattern is promoted from this one lane.

## 1. The full dictionary and its two colors

Let `w=(a,b,c)` satisfy

```text
1<=a<b<c, gcd(a,b,c)=1,
a,b,c odd and nonzero modulo 3.
```

Define the integral rank-two lattice and its full open carrier dictionary

```text
L_w = {C in Z^3 : C dot w=0},
p_i(C)=3(w_j+w_k)-14|C_i|,
Omega_w={C in L_w : p_i(C)>0 and C_i!=0 mod 3 for every i}.               (2)
```

The region cut out by `p_i>0` is open, convex, bounded, and centrally
symmetric in the real kernel plane. These strict inequalities exclude
tangencies consistently with THM-4414/4422.

For a live carrier, the three nonzero residues `w_i C_i` have sum zero
over `F_3`, so they are all equal. Set

```text
tau(C)=aC_1=bC_2=cC_3 mod 3 in {1,2}.                                  (3)
```

The sets `tau=1` and `tau=2` are individual affine cosets of `3L_w`.
Negation swaps them, so they each contain `N/2` members, where
`N=|Omega_w|`.

One can retain the intermediate index-three lattice

```text
Gamma_w={C in L_w : aC_1=bC_2=cC_3 mod 3}.
```

Then `Omega_w` is exactly the convex section of `Gamma_w` with its
zero-color sublattice `3L_w` deleted. Both indices
`[L_w:Gamma_w]` and `[Gamma_w:3L_w]` are three. This is an arithmetic
description of the true admissible point set, not a replacement by a
uniform random subset.

## 2. Midpoint forcing and a quantitative collinear-triple tariff

**Lemma (PROVED).** If two distinct carriers `C,D` have the same owner
color and the same class in `L_w/2L_w`, then

```text
M=(C+D)/2
```

is a third distinct live carrier, of the same owner color. Therefore the
three points `C,M,D` are collinear.

**Proof.** Equal parity makes `M` integral. Since `C,D` lie in the real
kernel, so does `M`. The lattice is saturated: if an integral vector lies
in that plane, it belongs to `L_w`; equivalently the coordinatewise parity
kernel on `L_w` is exactly `2L_w`. Convexity preserves all three strict
roof inequalities. Finally, modulo three,

```text
w_i M_i = (w_i C_i+w_i D_i)/2 = (tau+tau)/2 = tau !=0.
```

Thus the owner gate is preserved. The midpoint is distinct from both
endpoints. QED.

The rank-two free lattice `L_w` has exactly four parity classes. More
concretely, because all three speeds are odd, a carrier's parity vector
has even coordinate sum and is one of
`000,011,101,110`; all four occur in the lattice. Consequently:

**Corollary (PROVED).** If the full dictionary has no three collinear
points, then

```text
N<=8.                                                                  (4)
```

Indeed, each of the two owner colors has at most one point in each of the
four parity classes. This is a sufficient universal bound; no assertion
that cardinality eight is attainable by an LRC dictionary is made.

There is also a useful quantitative survivor when (4) fails. Write
`N=2n`, `n=4k+r`, `0<=r<4`. Let `A(Omega_w)` count unordered triples
which are equally spaced and lie entirely in one owner color. Every
unordered pair of equal-color, equal-parity endpoints forces one such
triple, and distinct endpoint pairs have distinct triples. Balancing the
four parity buckets in each color gives

```text
A(Omega_w) >= 2[4 binom(k,2)+r k].                                      (5)
```

In particular `N>=10` forces at least two such triples, exchanged by
negation. This is deterministic incidence, not independent-event counting.
It supplies no positive numerical deficit by itself: all three points can
lie where a particular deficit function is affine or zero. Endpoint
distance from the actual weighted boundary strip remains a needed sidecar.

The midpoint bound extends verbatim to a complete convex section of a
rank-`d` lattice restricted to `q` cosets of an **odd** index dilation:
a cap has at most `q*2^d` points. Oddness is what makes division by two
preserve an individual admitted coset. This generalization is algebraic;
it supplies no claim about arbitrary selected subsets of a square grid.

## 3. All-height projection gap for multi-direction caps

Use the exact projection capacities of THM-4414/4422:

```text
q=3/(7c),
S_i=sum_(C in Omega_w) min(q,p_i(C)/(14w_jw_k)).                         (6)
```

Suppose the full dictionary is not contained in a line and has no three
collinear points. By (4), each projection satisfies the trivial but
consumer-correct upper bound

```text
S_i <= N q <= 24/(7c).                                                   (7)
```

For every integer `c>=101`,

```text
24/(7c) <=24/707 <204/5957<6/77.                                         (8)
```

The complementary base is exactly `1<=a<b<c<=100` with the assumptions in
Section 1. This includes every allowed smaller height; `c=100` itself is
excluded by the stated oddness condition, rather than silently omitted.
The exact script independently enumerates all `5,409` such primitive
triples. There are `1,723` full cap dictionaries and `1,230` noncollinear
full cap dictionaries. Their cardinalities are

```text
N=0: 36 rows; N=2: 457 rows; N=4: 723 rows; N=6: 507 rows; N=8: none.
```

All `N=4` and `N=6` cap rows are noncollinear. The exact best values by
cardinality in this finite base are

```text
N=4: 12/413   at (1,7,59),
N=6: 204/5957 at (23,29,37).                                              (9)
```

The overall maximum in (9) occurs at exactly one triple. Together with
the strict tail (8), this proves the universal sharp bound (1) and its
equality classification for this entire geometric class.

The full sharp dictionary is

```text
w=(23,29,37),
Omega_w=+/-{(1,-11,8),(10,1,-7),(11,-10,1)},
(11,-10,1)=(1,-11,8)+(10,1,-7),
(S_1,S_2,S_3)=(276/7511,204/5957,7776/172753).                            (10)
```

Its exact physical mass is only `3100/172753`. Thus the sharp claim is
about the specified network/projection consumer, not the physical comb
maximum. The strict gap to the target is

```text
6/77 - 204/5957 = 2862/65527.                                             (11)
```

The older equality comb `(1,5,11)` has only two carriers on one line and
is outside the new multi-direction hypothesis. The first dense `A_2`
example `(19,23,29)` is included and strictly closed by (1).

## 4. The full-dictionary firewall

At `w=(1,11,23)` the full dictionary is

```text
{+/- (1,2,-1), +/- (2,4,-2)}.
```

It has four collinear points. Keeping only the two outer endpoints makes
a cap, but does not make the original dictionary a cap. Likewise a chosen
six-point cap inside a larger carrier support inherits neither (4) nor the
capacity bound (7). The source in every theorem above is the **complete**
set in (2), with the original owner and roof filters. Cardinality,
direction support, or a visually attractive cap subset is not a substitute.

The first implication that fails for a selected subset is midpoint
completeness: the midpoint is guaranteed to be an admitted carrier, not a
member of the chosen subset. The repaired object is the full incidence
hypergraph equipped with owner color and parity. The new remaining question
is how its forced arithmetic-progression edges interact with the three
weighted boundary-deficit functions.

## 5. Verification and connection contracts

Source:
[overnight_20260906_lrc_cap_carriers.py](../../04-computation/overnight_20260906_lrc_cap_carriers.py).
Output:
[overnight_20260906_lrc_cap_carriers.out](overnight_20260906_lrc_cap_carriers.out).

```powershell
python -B 04-computation/overnight_20260906_lrc_cap_carriers.py
python -B -O 04-computation/overnight_20260906_lrc_cap_carriers.py
```

The exact integer box enumerator and an independent Bezout-progression
enumerator produce identical full dictionaries on every base row. Direct
three-point determinants and a separately implemented direction-at-each-
vertex test agree on every cap verdict, including rows rejected by the
midpoint theorem. The proof's forced midpoint triples are checked directly
for membership, color, and injectivity of their endpoint address. The
sharp row also agrees with literal physical interval construction and exact
augmenting max-flow. All comparisons use integers or rational fractions;
there are no optimizable Python assertions or floating-point decisions.

The normal replay executes `115,344` explicit checks and audits `10,610`
forced same-owner midpoint triples. Normal and optimized outputs agree
byte for byte. SHA-256 values of the raw LF source and output are:

```text
source: 9a717d42f3cf1e2d29b2c9d5f5f6a492e3c42e5c9e599d2e009b6e816a35c166
output: 0c1ccbdcb61acaa823c7bb85137ef752affc944a6b5ce5db2887ee6177478324
```

The root agent independently accepted the universal midpoint/coset proof,
the injective endpoint-pair count, and the strict tail arithmetic. Its
separate no-import
[middle-coordinate audit](../../04-computation/overnight_20260906_root_cap_audit.py)
enumerates carriers by solving for the middle coordinate and computes all
projection roofs using the integer denominator `14abc`. It independently
reproduces all `5,409` base rows, `1,723` caps, `1,230` noncollinear caps,
and the unique sharp row. Its normal and optimized streams match the saved
[independent output](overnight_20260906_root_cap_audit.out).

| Source -> target | Map / preserved predicate | Loss / required sidecar | Decisive test |
|---|---|---|---|
| Guy–Kelly line counting -> LRC incidence | Use exact collinear triples on the full lattice carrier | Random-grid distribution and its heuristic do not transfer | A same-color equal-parity pair forces its midpoint |
| Convex owner dictionary -> eight parity buckets | `(tau(C),C mod 2L_w)` | Locations and boundary deficits are forgotten | A repeated bucket produces a certified live triple |
| Cap carrier -> projection ceiling | `N<=8`, each actual summand `<=3/(7c)` | No entry or other-speed synchronization is supplied | Complete `c<=100` base and strict `c>=101` tail |
| Large dictionary -> arithmetic-progression hypergraph | Endpoint pair maps to its compulsory midpoint triple | Counts alone lose weighted strip location | Apply all three exact convex deficit functions to each triple |
| Full dictionary -> selected cap | Inclusion only | Midpoint closure is destroyed | `(1,11,23)` endpoint-subset hostile |

**Remaining frontier.** A genuinely multi-direction dictionary exceeding
the bound (1) must contain a collinear triple. A dictionary with at least
ten carriers must contain the stronger same-owner midpoint structure in
(5). Neither condition alone is an exclusion. Collinear multi-direction
dictionaries, their location-dependent boundary deficits, arbitrary local
nonresonance, physical entry, and LRC(14) remain open.
