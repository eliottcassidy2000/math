# Zero-heavy cores, exact ports, and product symmetry

**2026-10-20. Current scope:** unrestricted A000568 compression and all-order
black-core connectivity remain **OPEN**. This packet separates exact
operation-preserving encodings from counts that forget their ports.

| Lane | Result | Status and source |
|---|---|---|
|Zero-heavy boundary|One-feedback-vertex class count A_1(n)=2^(n-1)-(n-1)-binom(n-2,2) for n>=2; exactly two pair-free classes at each n>=4|PROVED, [anchor quotient](zero_heavy.md)|
|Fixed anchor number|Canonical re-anchoring by bounded-delay port-word transduction gives rational ordinary and pair-free class series for every fixed feedback-vertex bound|PROVED, independently audited, [construction](zero_heavy.md)|
|Separating anchors|The ordered join of m directed triangles needs exactly floor(3m/2) anchors to distinguish all remaining vertex signatures|PROVED, [sharp hostile](zero_heavy.md#4-why-a-uniformly-small-separating-anchor-cannot-suffice)|
|Closed blue endpoints (n-3,2)|Entire layer has three black moves to the literal partner for every n>=6 after integrating the concurrent whole-block repairs|PROVED with FINITE-EXACT head6–8, [proof integration](endpoint.md#9-incoming-completion-of-all-three-residuals)|
|One-vertex numeric growth|Every directed path-cover count has an exact quadratic attachment response, with specified component endpoints; its labelled average depends only on the old profile|PROVED, [port tensors](profile_growth.md)|
|False smaller growth state|Pair-free six-cores with equal full profile, scores, feedback and automorphism cycle index have maximum child H values159 and147|FINITE-EXACT counterexample, [witness](profile_growth.md#2-smallest-reinforced-profile-hostile-in-the-audited-universe)|
|Square/diagonal microcosm|A product of pointed Boolean sets has literal commuting involutions whose orbit counts give the primitive Pythagorean triple; actual deletion fibres have sizes4,2,1|PROVED, [symmetry and growth](microcosm.md)|
|Capped weighted core|Same-order eleven-vertex family members have identical cap and scores but blue counts17 and0; actual chain-length equality is lost|PROVED example, [predicate boundary](microcosm.md#6-the-same-five-vertex-core-reaches-the-zero-heavy-boundary)|

The board comprises feedback anchors, canonical word ownership, homogeneous
chains, two-component endpoint tensors, blue positional involutions, and
the exceptional state in the square product. The anchor is the unresolved
zero-heavy/strong-blue boundary; the niche is exact one-vertex response;
the wildcard is the commuting-involution interpretation of (3,4,5).

The closest inherited mechanisms and hostile cases are named at the start
of each note. The inherited run-word clone dichotomy is checked again but
not claimed as new. THM-1900's incorrect deletion-of-path premise is repaired
in place, retaining the true neutrality theorem and its orphan terms;
THM-012b's three-cycle-only range is corrected from n<=5 to n<=4.

The exact checks are ordinary Python with explicit failures. Their declared
universes, independent paths, and limits are recorded beside each proof.
No finite census is used as an all-order connectivity argument.

For the endpoint result write the marked path as `(a,b,W,c,z)`, with
`W=(u,...,v)`, and let x be the extra middle vertex pointing to a.
The initial large-order cases requiring another ledger were: open
middle with x=v; closed middle with b not pointing to v and x=c; closed
middle with b pointing to v and x=u. Two natural endpoint-path repairs
fail on explicit infinite families. These are failures of the specified
repairs, not non-leaf bridges or disconnected black components. The incoming
[complete endpoint-two theorem](../endpoint_two_repairs_20261020g/README.md)
now repairs all three cases by whole-block words, independently reviewed
here. The complete layer is closed; the next direct layer is (n-4,3).
Exact clone tests retain the two outer weights in `C3[TT2,W,TT2]`.

## Incoming work and the next precise questions

The incoming [feedback-word spectrum packet](../feedback_word_hamilton_spectrum_20261020/README.md)
already counts H with fixed feedback anchors. Our additional operation is
canonical re-anchoring: it removes isomorphism duplicates while retaining
the complete word. Both are finite-state for a fixed anchor bound, but they
count different objects. Its [binary-product density bound](../feedback_word_hamilton_spectrum_20261020/BINARY_PRODUCT_DENSITY.md)
concerns the number of attainable numerical H-values; it is not a bound on
the number of canonical tournament words. The six-vertex response twins
both have feedback vertex number two, so this is a concrete next state
space on which to combine the two machines.

The incoming [score-fibre return Gram matrix](../regular_return_gram_icosahedra_20261020/gram_return_addendum.md)
keeps actual intermediate tournaments as columns and has raw H on its
diagonal. Its five-vertex icosahedral faces are actual incidence fibres.
This is a useful counterpart to the square product: in both cases an
arithmetic or geometric structure follows from a specified map of objects.
Neither map alone preserves the black-connectivity predicate. A useful next
test is to keep only certified black incidences in that return carrier and
compare its components with the next endpoint-defect layer, retaining labels
and endpoint pairs throughout.

For unrestricted A000568 the new finite-state theorem changes the obstacle:
one need not enumerate every core at a **fixed** feedback-vertex bound.
What remains is a controlled relation between the owner automata as that
bound increases. The ordered joins of triangles show why simply declaring
the anchor alphabet uniformly bounded cannot work. No efficient uniform
transition or closed unrestricted scalar recurrence has been established.

See [AUDIT.md](AUDIT.md) for the proof reviews, computational universes,
correction lineage, source pins, and reproduction command.
