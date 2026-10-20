# The golden21 bank separates junction labels but does not force junction support

**FINITE-EXACT:** the specified21 Fano junction probes fail the proposed
zero-support certificate for every one of the105 target junctions. They
nevertheless distinguish all105 targets by their complete measurement
vectors. **PROVED implication:** the incoming physical Kempe-uniformity
theorem cannot be applied through this failed support-cover argument.
No physical counterexample, full21-measurement nonuniqueness theorem,
or rooted-minor obstruction is asserted.

## 1. Inheritance and the exact proposed map

The new [physical target certificate](../deficiency_three_20261020/boundary.md),
section2, proves that14 suitable zero closure measurements force the
boundary support of a physical cubic seven-pole into one junction
kernel J. An independent Kempe argument then forces its tensor to be
`mJ`: a two-colour path starting at one target wire end must finish at
the other end of that same wire, since any other endpoint would produce
a proper colouring outside J after switching. These switches biject
the nine wire-palette assignments; global colour permutations handle
the six ordered junction palettes. The self-closure count is `54m`,
or `9m` on colour-orbit coordinates.

That mechanism was read and audited before this test. The checker below
also independently transports and verifies the incoming14-probe support
cover for every target. The closest earlier mechanism is the
[marked golden-character bank](arithmetic.md), section5: mark zero in
the character cube F2^3, choose a nonzero direction d, delete the pair
containing zero from its parallel matching, and turn one of the other
three pairs together with d into a junction. The seven directions give
21 closures. The corrected near miss is treating rank deficiency as
either an automatic obstruction or an automatic physical certificate.
The retained sidecar here is the complete set of nonnegative boundary
support columns, not just a Gram eigenvalue.

Number the seven nonzero character addresses1,...,7. A label
`abc;de,fg` has an all-distinct triple abc and equal-colour pairs de,fg.
The21 bank consists exactly of labels whose triple satisfies
`a xor b xor c=0`, so its junction triples are Fano lines. This is checked
against the parallel-direction construction, not merely inferred from
the count21.

For a target J, take **all** bank rows L with `<L,J>=0`. Since the
rows are0/1 and J is nonnegative, these are precisely the bank supports
disjoint from J. Their zero measurements exclude the union of their
positive columns. The decisive question is whether that union is the
entire complement of the nine-orbit target support.

## 2. Complete negative support-cover test

The universe is all `3^7=2187` labelled colour words. Parity retains546
words; global colour symmetry gives91 orbits of size six. Every junction
supports54 words, or nine orbits. All105 target labels are tested.

| Target type | Targets | Disjoint bank probes | Surviving orbits | Surviving orbits outside J |
|---|---:|---:|---:|---:|
|Target itself belongs to the21 bank|21|10|33|24|
|Target lies outside the bank|84|8|39|30|

Thus the bank zeros leave198 or234 labelled words possible, respectively,
instead of the desired54. Every available zero row was used; taking
fewer of them cannot repair the cover. The two uniform cases are actual
GL(3,2) orbits of target labels: the group of168 invertible linear address
maps preserves the bank and acts transitively on its21 members and on
the84 other junctions. Both orbit claims are independently enumerated.

Here are lexicographically least surviving orbit representatives for
one target in each orbit, with colours written in port order1,...,7:

| Target | Surviving word outside target |
|---|---|
|`123;45,67`|`0010211`|
|`124;35,67`|`0001020`|

Each displayed word has odd multiplicities for all three colours. It
fails its target junction but is accepted by none of the bank rows
whose target measurement vanishes. These are minimal words in the
declared lexicographic orbit order, not claims about smallest physical
graph realizations.

For a purely tensor-level hostile, add the indicator of the entire
six-word colour orbit of such a witness to J. The resulting nonnegative
invariant parity tensor keeps all those zero measurements and keeps the
target self-closure count54, but has support outside J. Its physical
realizability has not been established. This example only shows why
the support inference and subsequent Kempe-uniformity inference cannot
be performed from these zero constraints and normalization alone.

## 3. The positive control and the precise remaining question

The complete21-component measurement vectors of the105 junction
targets are all distinct. Entries in the receipt are orbit inner
products; multiply each by six for labelled-colour closure counts.
Therefore the bank identifies the target label **if the unknown tensor
is already known to be one of these105 normalized junction kernels**.

The failed support cover does not settle whether all21 measurement
equalities, together with physical realizability or additional Kempe
constraints, could certify a junction target. Those nonzero equalities
were not replaced by the weaker set of zero equations in a claimed
equivalence. Conversely, distinguishing105 candidate labels does not
certify that an arbitrary unknown tensor belongs to that finite family.

The incoming14-probe certificate remains valid. Its successful support
cover uses target-adapted rows from the full105 bank; the fact that21 is
larger than14 does not make this particular21-row subset adequate.
Neither result supplies internal disjoint paths or a Petersen-minor
certificate without an additional graph-side argument.

## 4. Reproduction and source fingerprints

```text
python3 05-knowledge/results/scc_basins_20261020d/arithmetic_junction_certificate.py
python3 -O 05-knowledge/results/scc_basins_20261020d/arithmetic_junction_certificate.py
```

The checker uses only the Python standard library and explicit failures.
It compares normalized-orbit computations with direct labelled-word
enumeration, verifies all105 transported incoming covers, reconstructs
the21 bank in two ways, and checks both GL(3,2) orbits. Its full
[JSON receipt](arithmetic_junction_certificate.json) retains every
target, its zero-probe labels, all surviving outside orbits, and its
complete bank measurement vector.

Normal and optimized runs have identical stdout and byte-identical JSON.
The JSON SHA256 is
`67f49854009884f6aeda3a0fc90a91217b75cdd2f803ced03813958b6445ac87`.
The normalized-LF source hash of the incoming boundary theorem read here is
`638670dd9128cf712a9bb617c501c89cab0c3364c581c42caa46432abab86e99`.
The checker records its own normalized-LF source hash. The incoming
packet first arrived in commit `a852aca417`; its later integration does
not change the declared test universe.
