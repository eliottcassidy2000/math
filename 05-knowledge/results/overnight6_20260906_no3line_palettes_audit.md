# Independent audit of the palette-conditioned no-three-line expansion

**INDEPENDENT PROOF AUDIT + EXACT CONTROLS PASS.** This accepts the full
scope of [the palette report](overnight6_20260906_no3line_palettes.md),
including local injection normalization, actual conditional cumulants,
outer palette corrections and the positive cylinder algebra. It is an
audit of that synthesis, not a separate discovery claim. No producer,
repository navigation or Git files were edited by this audit.

## 1. Probability and automorphism normalization

The object is a fixed simple bipartite 2-regular skeleton with independent
uniform full bijections on its two shores. Naming the components is
essential when cycle lengths repeat. Conditioning on each component's
actual row and column palettes leaves a Cartesian product of uniform
internal bijections. Each named palette pair has precisely
`product_i (k_i!)^2` full labelings, so all palette pairs of the specified
sizes have the same probability `(product_i k_i!/n!)^2`.

For the complete edge union F of a family of grid events, a row/column
palette mismatch makes its conditional probability zero. Otherwise, in
component i the inverse images of the incident grid rows and columns are
uniform injections. Their denominator is exactly
`(k_i)_(r_i)*(k_i)_(s_i)`. The number that preserve the required edges is
`inj_sh(F_i,C_(2k_i))`. Extra edges are allowed: this is a non-induced
embedding problem. Dividing both this numerator and its possible-host
denominator by the same shore-preserving automorphism group gives the
stated copies/N formula. There is no extra shore swap or endpoint flip.

This calculation uses the complete union of the events. Repeated events
do not change joint presence, but the list of occurrences must be restored
when cumulant partition terms are formed. Local degree three, a closed
C4 inside a C6, and palette-crossing edges have zero probability. In
contrast, two specified disjoint edges have probability one in C4=K22;
an induced-subgraph interpretation would incorrectly reject this control.

For one nonaxis line triple, each local component receives a matching,
because the three rows and three columns are distinct. A specified
s-edge matching has `s!*m_s(C_(2k))` shore-preserving injections: only the
s edges are permuted, and their row/column orientation is fixed. This is
why s!, rather than `2^s*s!` or one, occurs in the matching kernel.

For completeness, on a cycle with L edges a nonempty s-matching may be
counted by choosing a marked selected edge and the positive cyclic gaps
of unselected edges. This gives
`m_s=(L/s)*binom(L-s-1,s-1)` for `1<=s<=floor(L/2)`. Substituting L=2k
gives exactly the three matching counts and q formulas in the producer.
The k=2, s=3 case is excluded by s>k before its singular displayed formula
is used. The independent check includes `q_(3,3)=1/3` and `q_(2,2)=1`.

## 2. Genuine conditional connectedness and outer corrections

Once palettes are fixed, different skeleton components have independent
internal permutation variables. Four-cycles have no residual randomness:
their full 2x2 palette rectangles are already determined. Removing those
components from an event's dependency set is therefore valid. If the graph
on event occurrences formed by overlap of the remaining component sets is
disconnected, its two parts are functions of independent variable blocks.
The joint moment-generating function factors, so the mixed coefficient in
its logarithm is zero. The statement remains true for repeated events,
and for occurrences that are conditionally constant. Connectedness is
only a necessary condition, as stated.

Multilinearity expresses the ordinary conditional cumulant of X as the
sum over ordered event lists, with repetitions. At order three the stated
five-term kernel is the ordinary joint cumulant and uses the correct
complete-union probabilities. No factor six is silently introduced or
discarded: ordered occurrences belong to the cumulant expansion itself.

The third total-cumulance formula is proved directly by conditioning the
cube of `(X-mu_P)+(mu_P-E[X])`. The term with one centered conditional
factor vanishes, while the term with its square gives
`3 Cov(mu_P,Var(X|P))`. Thus both outer terms in the report are required.

The all-order statement can also be checked formally without an extra
probabilistic assumption. Write the conditional log moment-generating
series as K_P(t). The unconditional series is
`log E[exp K_P(t)]`; expanding its outer joint cumulants and extracting
a multilinear coefficient partitions the event occurrences into blocks.
Each block contributes its conditional joint cumulant, and the outer
cumulant combines those random block quantities. All variables and
palettes are finite here, so the formal coefficient operations are valid.
The report appropriately attributes this classical identity to
[Brillinger's conditioning paper](https://www.ism.ac.jp/editsec/aism/pdf/021_1_0215.pdf);
the local graph normalization and the geometric controls are separate
repository work.

The ordinary/factorial warning is essential. For conditionally deterministic
X=x, ordinary cumulants of order>=2 vanish, but the third factorial
cumulant is 2x. At n=4 the mean conditional third factorial cumulant is
4, while the unconditional third factorial cumulant is 4/3. Reusing the
ordinary conditional graph rule for factorial cumulants would be invalid.

## 3. Independent board universe and finite values

The referee uses a third generation route. It chooses exactly two column
neighbors for each row, rejects column degrees other than two, and finds
components by graph BFS. It uses neither the producer's palette-board
generator nor its complete shore-permutation census. This enumerates each
labelled binary board exactly once, with the following complete counts:

| n | component type | distinct boards |
|---|---|---:|
| 4 | 2C4 | 18 |
| 4 | C8 | 72 |
| 5 | C4+C6 | 600 |
| 5 | C10 | 1,440 |

For n=4, either actual C4 may represent the first named component. Thus
the 18 boards give 36 named palettes, each with one conditional board.
For n=5, the C4 component is uniquely identifiable; there are 100 palettes
with six conditional boards apiece. The full-label multiplicities 32 and
24 satisfy `18*32=4!^2` and `600*24=5!^2`, respectively. The report's
named palette law, conditional board law and unconditional board law
therefore agree with the independent universe.

The n=4 named-palette histogram is exactly
`{0:18,2:8,4:4,6:4,8:2}`. Direct central moments reproduce
`(mean,var,kappa3)=(2,56/9,16)`, with total-cumulance terms `(0,16,0)`.

At n=5 the independent census reproduces every producer histogram entry
and gives

```
(mean,var,kappa3)=(13/3,1769/225,26159/675),
(E conditional kappa3, kappa3 conditional mean, 3 covariance)
 =(-871/675,6467/360,7949/360).
```

Every singleton event probability and every palette mean is checked.
At two prescribed n=5 palettes the referee also checks every list of one,
two or three grid events up to permutation of occurrences, including
repetitions, against complete-union injection probabilities. There are
52 nonaxis grid events; these tests retain their actual Euclidean labels.

For n=6, row-neighbor enumeration inside the two prescribed 3x3 rectangles
gives exactly 36 conditional boards. It independently recovers

```
(pA,pB,pC,pAB,pAC,pBC,pABC)
 =(1/3,1/3,1/3,1/9,2/9,1/6,1/9),
kappa(A,A,B|P)=0,
kappa(A,B,C|P)=1/54.
```

The minimality claims are correctly restricted. Two disconnected simple
bipartite cycle components require n>=4. Two components with residual
internal randomness after conditioning require both cycle lengths at
least six, hence n>=6. Neither statement claims global minimality outside
this model.

## 4. All-four-cycle representation and positive algebra

For n=2m, the row labels and column labels are each partitioned into
unordered pairs, and a bijection associates the two collections of pairs.
Each board made of m disjoint K22 components uniquely determines these
three objects. This proves the count
`(n!)^2/(2^(2m)*m!)` and uniformity. Named palettes count each board m!
times; each palette has `2^(2m)` internal labelings, giving the same
automorphism normalization.

The referee enumerates all 1,350 such n=6 boards and checks the entire
geometric occupancy decomposition on each one. A triple using two points
of a block must use one of its two opposite-corner diagonals. Its majority
block and remaining block are uniquely determined, so the ordered block
sum has no unwanted factor two. A triple using three blocks is counted
once by an unordered choice of blocks and one point from each. A point
at the intersection of another block's diagonals yields two different
triples; the separate diagonal sum preserves both. No event wholly in
one K22 can have three distinct row and column labels.

The partial-assignment cylinders are actual indicator functions. Their
product rule retains both source conflicts and repeated target labels on
a shore. Extension relations resolve each partial cylinder as a sum of
full-assignment atoms. These atoms are orthogonal idempotents summing to
one, so the stated state is positive on every real square, not merely
on the sample squares in the producer.

Independently, at n=3 the referee constructs all 34 one-shore partial
injections and their indicator vectors on the six full permutations.
Its entire 34x34 Gram matrix equals `A A^t/6` and agrees entrywise with
the collision-sensitive product/expectation rule. The two-shore version
is the tensor product, proving the same positivity. Conditioning on
named palettes gives exactly the product state on internal permutations;
averaging over palettes gives the full state.

This repairs the multiplication only after actual labels are retained.
It does not make the old unlabelled disjoint-union weighting map positive,
nor does it assign an old formal defect a new random-variable meaning.
The report preserves that distinction throughout.

## 5. Reproduction, audit limits and freeze

```
python -B 04-computation/overnight6_20260906_no3line_palettes_audit.py
python -B -O 04-computation/overnight6_20260906_no3line_palettes_audit.py
```

Both LF outputs agree: **67,111 optimization-live gates PASS**. The
auditor imports only Python's standard library and no producer or prior
event-census module. Its input hash gate freezes the reviewed producer.
No new n=8 third-event census was performed. The all-size proof is the
conditional-injection and cumulant argument, not a finite extrapolation.

```
audit source 9cff9a74ef208c54bc5964a315367d088761f47dc61148762f7210db8c61c926
audit output 885415823263c2fbfdd041931e332ef3e0ac5973a496f25801bac3124d2179fd
producer source cc0fc2774d8475766da8f90ab4a5bb12d070170bfc5ba2b4d3bbf0c74e05f80f
```

The strongest conclusion is the exact conditional local kernel and genuine
connected expansion with its outer mixture terms. No unconditional
component-additive formula, all-size sign claim, tail estimate, or positive
probability of no collinear triples follows from those identities alone.
No mathematical correction to the producer report was needed.

**Filing:** root integrated this independently audited report after
`f5f0f7f75`. The proof and transcript are frozen in the sixth manifest.
