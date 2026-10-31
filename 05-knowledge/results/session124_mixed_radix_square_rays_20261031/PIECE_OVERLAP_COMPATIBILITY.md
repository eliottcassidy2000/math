# Exact classification of positional-parity compatibility by piece overlaps

**PROVED.** This classifies the compatibility sidecar, not the existence of
Hamiltonian connector words. Ownership, square edges, endpoints, and connected
traversal remain separate obligations.

## Two outputs

Fix one affine offset and two selected input paths with a common reference
position parity ε(x) for every label x common to them. Cut each into contiguous
pieces. Every oriented piece A placed in an output word has a well-defined
shift σ_A in F₂ such that every label x of that piece has output position
parity ε(x)+σ_A. Reversal changes the shift, not its constancy on the piece.

Build a bipartite overlap graph:

- left nodes: pieces used in output0;
- right nodes: pieces used in output1;
- edge A–B precisely when their source-label sets have a common label.

**Theorem.** The two outputs agree in position parity on all common labels
of this offset if and only if σ_A=σ_B along every overlap edge. Equivalently,
σ is constant on every connected component of the overlap graph.

Proof. A common label x belongs to one piece A in output0 and one piece B in
output1. Its output parities agree exactly when
ε(x)+σ_A=ε(x)+σ_B, or σ_A=σ_B. Taking every common label gives exactly the
edge equations. Equality propagates along paths; conversely a componentwise
constant assignment satisfies every edge equation. ∎

Thus compatibility alone has one free bit per connected component, before
adding prescribed-prefix, endpoint, connector-color, and realizability
conditions. This is a classification of solutions of the compatibility
equations, not a claim that all such assignments have a connector realization.
A piece containing only a new unshared label produces no overlap constraint;
actual output lengths and the parities of all shared labels still constrain
that new label in a real word.

## Affine ownership and different offsets

For radix R and |c|<R/2, each source label has its own disjoint affine band.
Within a band, compatible source parities identify the same offset owner in
both outputs. Distinct offsets therefore do not share output labels. Apply
the graph theorem separately to each offset, and compare each low singleton's
position parity directly. The exceptional top label appearing in only one
output is omitted from the overlap graph by definition.

## Why separating a known prefix works

For the cut after input1, the common label sets are {1} and all common labels
other than1. The graph has two components: head-to-head, and tail-to-tail.
Their shifts may be different. In actual-head-color notation H and tail-shift
notation S, compatibility requires H0=H1 and S0=S1, not H=1+S within a word.
The successful relaxed Q25 certificate at even input/residue18 explicitly
uses this freedom: head19 has even output parity while its tail has shift0.

For a cut at a nonprefix edge2–7, the labels on each side can differ between
the two selected paths. Cross-overlap edges may connect the nominal left and
right pieces into one component. Then the graph theorem requires a uniform
shift. A certified membership-overlap graph can justify more freedom; without
it, assigning independent shifts is not safe. A mere named cut edge does not
supply that missing membership information.

## More than two outputs

Use one node for each output piece, and a hyperedge for each source label that
occurs in at least two outputs, incident to all pieces containing that label.
The same proof gives equality of all shift bits on each hyperedge, hence one
bit per connected component of the incidence graph. This is the natural
multi-output compatibility sidecar; fixed or prescribed colors add further
linear equations over F₂.

## Correction lineage

Earlier frozen tables imposed one uniform shift on every piece of each
source block. They therefore satisfy these equations and remain valid.
That stronger sufficient architecture was unnecessarily restrictive for the
known-prefix head/tail partition. The new result does not retroactively
invalidate any old certificate; it identifies exactly which equality can be
removed and which overlap information must be retained when cuts vary.
