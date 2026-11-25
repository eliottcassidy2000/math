# Overnight Chan–Pak / degree-57 Moore tandem

Research span: 2026-11-24 evening to 2026-11-25 morning (UTC−7). The session alternated between constructive tournament height calculations and necessary conditions on a hypothetical degree-57 Moore graph. It does **not** settle either open problem.

## Results at a glance

| Side | Verified result | Scope |
|---|---|---|
| Chan–Pak | A new eight-site core and three **fixed** physical \(+2\) channels with an integral inverse at every modulus | Joint readings encode the payload; a prescribed single absolute odd target still need not occur |
| Chan–Pak | \(\Omega(X^{2/3}/\log X)\) **distinct** literal one-arc \(+2\) pairs below \(X\), even in any fixed admissible lower-height class | Holds separately in the six- and eight-site variable-tail families; not positive-density coverage |
| Chan–Pak | \(E_{31}\pm1\) absent from both full variable-tail guard families by exact search; \(E_{61}\pm1\) absent from the offset product-spaced subfamilies by a prime-divisor argument | Full variable-tail exponent 61 remains undecided; all compatible fixed moduli pass |
| Moore | Both native complement-spine endpoints exceed \(2^{18927}\), uniformly; the existing upper bounds stay \(2^{28860}\), \(2^{28863}\) | Conditional on existence; only native endpoint charts |
| Moore | Two distinct valid normalized decks differ in at least 1410 literal arc positions in a common displayed order | Conditional on both valid decks; frozen-cell spectral tests can certify uniqueness from partial data |
| Moore | Point-preserving Gewirtz transports excluded; full action has no sharp set with 7/49 swaps or containing a complete order-seven coset; no branch-sign gauge symmetry of prime order at least 13 | Full Gewirtz transport model is not excluded at remaining parameters; this says nothing about arbitrary automorphisms of a Moore graph |

A final pass [sharpens the distinct-pair denominator](ROUGH-DIVISORS-AND-CONSTANT-SEEDS.md) to \((\log X)^{0.086071\ldots}(\log\log X)^{3/2}\) without a prescribed residue, constructs a 45-site literal seam \(E_{23}+1\to E_{23}+3\), and reduces the five-cycle order-11 sign case to eight orbit-constant seeds. The improved bound and the independently counted example are not an absolute selector.

## The transfer between the two problems

The starting model was a native tournament written across a balanced cut. Its **common** paths survive the same-side arc edit; the exact signed difference is instead a contraction count read by endpoint type. That distinction suggested building a guarded transfer whose four large count channels stay fixed while one exact reading changes by 2.

Going back to Moore, an edit in a real normalized graph deck has the opposite character: every changed row forces distinct deletions at frozen probes, and the difference must cancel in each branch and label. In a degree-57 chart that prevents tiny edits. The rectangular completion test, circulation test and localized Gram spaces make those cancellations usable without a full graph witness. The arithmetic guard and the Moore graph payload are therefore separate objects; a literal one-arc count seam does not imply a one-arc trade of the graph.

## Tournament proofs and replays

- [Cut layers and Petersen controls](CUT-CANCELLATION.md): ten literal graph-preserving gap-2 swaps, five scalar pairs, and grading by same-side steps. Large common counts and small differences can coexist.
- [Complete single-apex inverse and gap-window counts](FIBRE-TRANSFER.md), [endpoint determinant](ENDPOINT-AND-SIGNED-DECK.md), and [minimum bare-gap proof](CIRCULATIONS-AND-HEIGHT-LAW.md). Any prescribed even height difference \(2d\) is possible; the one-apex class optimally begins at \(h(d-1)\). That property by itself does not prescribe either endpoint.
- [Six-site guard census](GUARDS-AND-DECK-DISTANCE.md) and [two short-guard charts](GEWIRTZ-AND-DUAL-GUARDS.md). A directly counted 14-site example has twin-prime heights \(15269\to15271\). The two charts can share a scalar height without inducing a common physical middle tournament.
- [Variable tails and exact scalar inverse](VARIABLE-TAIL-SELECTOR.md), [new eight-site guard and three fixed linear readings](EIGHT-SITE-FIXED-CHANNELS.md). The exact exponent-31 checks are in the second note.
- [Complete fixed-modulus laws](TAIL-LOCAL-LAWS.md) and [distinct product-spaced pairs](PRODUCT-SPACED-DISTINCT-SEAMS.md). The product quotient makes the lower bound on distinct outputs possible, and simultaneously explains the Mersenne-prime obstruction for that particular subfamily. A fuller absolute selector has to keep more flexibility than this quotient, or use a different chart.

## Moore proofs and replays

- [Weighted balanced-cut augmentation](WEIGHTED-AUGMENTATION.md): the uniform endpoint lower exponent 18927. [Original native-window proof](DEGREE57-HEIGHT.md): upper exponents.
- [Signed support and 1410-arc bound](ENDPOINT-AND-SIGNED-DECK.md), [frozen rows and columns](TWO-COORDINATE-RIGIDITY.md), and [arbitrary frozen-cell circulations and projector Gram minors](CIRCULATIONS-AND-HEIGHT-LAW.md). For example freezing 22 grid rows and 19 columns pins any valid completion if all incidences of those lines are prescribed.
- [Gewirtz point-incidence obstruction](GEWIRTZ-AND-DUAL-GUARDS.md) and [exhaustive seven/49 swap and coset check](GEWIRTZ-SEVEN-EXCLUSION.md). The exact-cover phase is necessary: the local Fano geometry alone passes at two order-seven subgroup fixtures. The full-group cases at 10,13,…,46 swaps are not all resolved.
- [No order 19 or larger-prime sign symmetry](NO-LARGE-PRIME-SIGN-SYMMETRY.md) and [no order 13](NO-ORDER13-SIGN-SYMMETRY.md). These constrain the type matrix before fitting actual transport permutations.

The proof notes point to their executable replays. None of the numerical or modular tests substitutes for the asserted scope of its accompanying theorem.
