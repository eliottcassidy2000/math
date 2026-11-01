# From an executable counter to a spectrum theorem

## The distinction that drives the session

Three logically different tasks must not be conflated:

1. Count classes at a fixed order, such as A_n, S_n or M_n.
2. Evaluate H(T) for a supplied tournament.
3. Decide or prove whether a target h is realized by some tournament of any order.

Actual finite FRACTRAN programs now implement the first two tasks. Their
existence alone does not solve the third. A counter must be paired with a
proved reduction or a constructive family whose range is controlled.

## Four exact bridges

**Certified atoms -> infinitely many counting terms.** Busch's bound turns
H<=21 into a finite inventory of strong SCCs. Unique ordered SCC decomposition
then yields an all-order positive linear recurrence for eleven multiplicity
sequences. The397-fraction stream,67-fraction tail and43-fraction alias tail
are actual implementations. The latter share registers only where the entire
sequences are proved equal. They do not identify the underlying graph classes.

**Counter -> bounded core, with the right quantifier.** The core-floor g(N)
has literal20/14-fraction implementations. Its exact inverse K(B) bounds the
vertices remaining after deleting all singleton SCCs. K(21)=7 is attained
with H15, not H21. This is a path-count budget theorem, not a witness for
an impossible exact count. All-order induced H-preserving deletion is canonical;
arbitrary contraction of directed twins does not preserve H.

**Budget -> local obstruction certificate.** Minimal induced obstructions
to H<=B have at most3r vertices, where3^(r-1)<=B<3^r, and the bound is sharp.
Except for the r-triangle ordered sum, the stronger boundK(B)+1 applies.
The complete B21 basis consists of56 explicit adjacency classes. Exact basis
plateaus at5/7 and19/21 give an alternative structural certificate for the
already-established exclusions7 and21. Equal basis cardinalities alone would
not suffice: the actual forbidden isomorphism-class sets are compared.

**One-apex word -> arithmetic -> converse-fixed fibre.** Keeping a marked
apex and the unique order on its transitive deletion gives an exact binary
word, not a guessed graph encoding. The32-fraction program computes H and the
odd-cycle count. Its inverse fibres and the square-shell theorem retain this
mark. The H5 example shows that forgetting the apex changes parity; inherited
classification of the exceptional double-apex tournaments supplies the
specific correction. The general adjacency evaluator keeps all arcs and has
no such family restriction, but it does not prove range completeness.

## What remains open

The statement that every positive odd integer except7 and21 occurs as H(T)
is not proved here. The exclusions were already proved and attributed to prior
literature; the new aim is stronger certified structure, not a priority claim.
The one-apex family misses13, which unrestricted tournaments realize, so it
cannot alone prove the universal statement. Adding a finite bank of seeds or
observing a long interval of outputs is not enough without a uniform reduction
for every remaining integer.

The positive finite-state multiplicity streams are possible because H is
bounded. They do not make the unrestricted A000568 sequence C-finite. Likewise,
small fraction counts do not imply efficient binary computation: output
valuations themselves impose execution costs. Exact microstep results are
stated only for the specified machines and valid initializations.

The closest prior converse-count sidecar is
`feedback_word_hamilton_spectrum_20261020/BINARY_PRODUCT_DENSITY.md`:
it already proves W(X)=X log(2)+O(sqrt(X)) for reduced words and bounds the
converse-fixed population by O(sqrt(X)). The square-shell result sharpens
that fixed-population sidecar to an exact count and resolves parity of each
pointed H-fibre. It does not claim the prior word-density asymptotic as new,
or improve the leading asymptotic density bound for distinct H-values.
