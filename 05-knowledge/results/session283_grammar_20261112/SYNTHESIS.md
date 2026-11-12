# Session283: exact families, metagraph fracture, and physical state codes

**Chan–Pak remains OPEN.** This session improves several exact recursive
interfaces without pretending that one of them already classifies every
unbounded tournament core or selects every desired path count.

## 1. An actual missing link in exact A000568 correction

Session279 had tight stratified moment certificates but still generated their
rich packets by listing all odd partitions. Its efficient low-order jets did
not preserve the required minimum-defect strata. The new
[conditional-family compiler](enumeration/PROOF.md) reverses the elimination:
retain the long cycle core carrying the defects, then integrate all short
cycles through their exact gcd-response field. This produces the **same full
stratified packet**, not an approximation to it.

For example the order64 packet uses266 long cores instead of16,444 full odd
partitions. Larger exact backbones now make the degree-four certificate
collapse to the exact integer at orders150,200 and250 without supplying
low residue bits. The independently inherited full counter and exact cache
are read only after recovery as controls; they do not choose the recovered
integer. The full decimal integers and operation/timing receipts are saved in
[the completion receipt](enumeration/completion.json).

This is a better certificate producer, not a faster full counter claim.
The optimized full C++ counter remains faster than the Python certificate
producer on the matched tested sizes. Exact values at these orders already
existed in the repository; the new content is the family compilation and
explanation of the certificate's failure boundary.

### An infinite classification from53 small patterns

The [threshold theorem](enumeration/THRESHOLD_PROOF.md) classifies the first
rational ambiguity of the degree-four ground-defect certificate for every
initial odd backbone ending at b>=27. Put a=b+2. Then

    first ambiguity = 4a+t(a),

where t is determined by53 four-length offset patterns and is periodic on
odd a with period6930. Its eight possible values are

    12,16,18,20,22,24,26,28.

The proof first bounds the obstruction's total mass, then excludes every
five-cycle core and every repeated-length four-core. What remains is a
finite difference-gcd grammar, valid at arbitrarily high levels. For a
backbone through35 the certificate is rationally exact through167 and first
becomes ambiguous at168. A backbone through49 is exact through215.
Rational ambiguity is not the same as failure of integer rounding. The theorem
also gives an all-order backbone scheduler: within b>=27 the smallest sufficient
backbone tops for n=150,200,250 are33,45,57. The proof, not a numerical trend,
justifies using the rule beyond the tested orders.

## 2. A complete local edit grammar on infinite weighted families

Every tournament has the inherited representation as a pair-free core Q
whose vertices are inflated to transitive chains. The new
[arc-edit grammar](grammar/PROOF.md) treats arbitrary chain lengths without
expanding their vertices:

- one arc reversal changes the reduced core rank by at most four, sharply;
- the complete one-step rank distribution uses at most9*binom(b,2) templates
  on at most b+4 vertices, where b is the old core rank;
- retaining the actual split lengths also gives the exact child tournament;
- if every chain has length at least three, the entire rank-transition row
  depends only on total size n and b, not the core shape or allocation.

If chains of length two are allowed, one extra count suffices. Singleton
boundaries require the endpoint templates. This is a full classification of
the specified one-edit response of every weighted family, not just finite
experimental agreement. A three-million-vertex example is computed from
small templates, without constructing those vertices.

Uniformity fails after arbitrary repeated edits unless its missing sidecars
are retained: one flip may leave the no-singleton sector. Explicit hostile
pairs also show that neither rank nor the whole separator histogram is a
Markov state for the actual rank distribution.

## 3. A closed all-time expectation grammar

The [separator heat law](metagraph/PROOF.md) nevertheless finds an exact
manageable evolution. Let P_T(z) count unordered vertex pairs by their number
of distinguishing exterior vertices. Under a uniformly random arc reversal,
the expected polynomial obeys a closed linear recurrence. Classical
Krawtchouk coordinates diagonalize this particular observable space.

It gives exact expected core rank at every time without enumerating the
metagraph. Under independent rate-one arc flips, every bounded-core family
has the macroscopic expectation law

    E[core rank at time tau/n]/n -> 1-exp(-2tau).

This is explicitly an expectation theorem, not a concentration claim or a
uniform-isomorphism-class law. It explains how small-core descriptions
fracture under editing while distinguishing expectation closure from a full
stochastic grammar.

## 4. One real path count can encode the whole substitution profile

The [physical state code](paths/PROOF.md) turns the complete ordered path-cover
profile of an n-vertex tournament into the Hamiltonian count of one fixed
strong triangle-cap construction. Positive cap coefficients become
superincreasing places. Greedy exact integer division recovers every profile
coordinate. A universal cap uses O(n² log(n+1)) vertices, and a certified
finite search supplies much shorter sufficient caps in small orders.

This improves the inherited many-query tomography into a single physical
observation. It does not identify the original isomorphism class. Two
five-vertex H9 first creators have the same full path-cover profile and hence
the same code, but different one-vertex extension menus. Uniform module
substitution and unrestricted growth therefore require different state
spaces even in this tiny example.

A constructive [contraction-deck addendum](paths/ADDENDUM.md) gives the exact
count polynomial for independently replacing transitive three-vertex modules
by directed triangles. Modulo2^r, only contraction subsets of size less than r
are needed. It classifies precisely when such a family covers all odd residues
modulo8, but proves that every cube with at least two independent switches
misses an odd number between its extrema. Compensated interactions or overlap
between different absolute bases are therefore necessary to turn this
particular recursive switch language into uninterrupted coverage.

## 5. The remaining Chan–Pak obligation

[The bridge note](BRIDGE.md) separates exact class census, exact operation
response, and prescribed-target existence. The new compilers can evaluate
and compare recursive families exactly. They do not yet construct a member
with an arbitrary required odd H. The previous renewable nine-bank remains a
useful base, but a well-founded all-target inverse descent is still missing.

The unifying practical principle is to eliminate an interior only after
fixing the interface read by the next operation. Here those interfaces are
respectively a gcd field, an endpoint split template, and a full path-cover
profile. Their honest loss boundaries are as important as their compression.

[Independent review](audit/REVIEW.md) and [reproduction/release receipt](REPRODUCE.md)
record the final source and controls. Full global core classification and
Chan–Pak completeness are not claimed.
