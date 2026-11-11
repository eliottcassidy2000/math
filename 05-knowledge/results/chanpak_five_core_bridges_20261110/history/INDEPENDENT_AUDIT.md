# Independent audit of the positive compiler and binary selector obstruction

**PASS, with one boundary clarification requested.** This audit is by a
different agent from the two proof producers. It covers the written
mechanisms in [POSITIVE_COMPILER](../algebra/POSITIVE_COMPILER.md) and
[SELECTOR](../cover/SELECTOR.md), with independent exact controls. It does
not audit the five-core all-height covering conjecture or promote it.

## Positive compiler

Grouping the inherited THM-016 sum by vertex support gives the signed
endpoint identities with the stated endpoint reversal: an end-marked
negative sum is minus the full-core start count, and a start-marked one
is minus the full-core end count. Summing over the distinguished endpoint
gives k_-=-H. The apex identity h=d+j_++j_- then gives equation(3) with
the correct signs. All its terms are nonnegative; h>=d follows from
insertion monotonicity or the even-path OCF sum, and k_+>=d from spanning
paths. The parity constraints and divided-by-two expression check.

The seven physical samples invert the coefficient map exactly. In
particular p20 contains 2d+2i, with no parity term; this is the sample
which separates d from the linear and endpoint terms. The largest listed
context has only two flank vertices plus the apex. The rank claim is
properly restricted to scalar linear measurements of the rational span;
it makes no unsupported claim about nonlinear encodings or arbitrary
future operations.

The independent [audit script](audit_compiler_selector.py) imports no
producer compiler, algebra or selector code. It generates the complete
20 ordinary tournament isomorphism classes through order5, retains all
470 cuts, enumerates directed paths on every subset by permutations,
and computes the endpoint weights directly. On 6,110 actual flank/apex
constructions through order12 it checks the positive formula against
this lane's subset-DP counter. The seven physical reconstructions,
endpoint signs, parity constraints and rational rank7 all pass. Uniform
cuts, empty flanks, and non-strong cores are included as boundary controls.

## Selector obstruction

The proof's essential positive term is d(D_A+D_B)m. The AB term is a
finite sum of exponentials under a pumped three-tape loop, even for
unequal input loop lengths; it contains no linear drift. Parities become
constant after at least one traversal. A nonzero exponential polynomial
with distinct positive bases, plus an affine term of positive slope,
has only finitely many zeros at nonnegative integer m: the largest
surviving base eventually dominates, or the affine remainder does.
A finite union of bank rows cannot cover all pumped m. This correctly
handles an unlabelled selector whose successful row changes with m.

Consequently every successful loop is monochromatic or empty on each
input tape. In a trimmed finite automaton, two different input digits
within one strongly connected component would create a mixed closed
walk. An accepting path crosses only finitely many components, giving a
uniform run bound. If transitions carry words rather than one digit or
epsilon, first subdivide those finitely many words into digit transitions.
The proof should explicitly specify ordinary most-significant-bit-first
representations; that is the convention used by its pumping formula.

The degenerate uniform-cut axis is treated correctly. There the count
is d*g(v); deleting trailing binary ones preserves g, never increases
the run count, and reduces to v=2q with g(v)>=q+1. Elsewhere a positive
linear or bilinear coefficient bounds both addresses by O(X). Thus the
normalization maps witnesses to bounded-run pairs of O(log X) bit length
without needing the selector itself to implement normalization. The
polylogarithmic output bound follows. This is a restriction on rational
ordinary-binary witness selection, not on output-support recognition or
unbounded-register algorithms.

**Boundary clarification:** periodic formula(2) is valid for m>=1 (or
eventually), not necessarily at m=0. For prefix1, loop0 and empty suffix,
g(2^m)=2^m+1 for m>=1 but g(1)=1 at m=0. The written pumping proof already
uses fixed parities for m>=1, so this does not affect the theorem.

Independent finite controls use the direct g recurrence, not the
producer's digit sum or virtual-prefix drift formula. For all510 words
of lengths1..8 and nine fixed prefix/suffix contexts, 4,590 exact checks
verify the exponential-annihilator identity, positive drift exactly for
mixed words, and independence of prefix and suffix. These are controls
of the written infinite proof, not a census of finite-state selectors.

## Reproduction

```text
python -B 05-knowledge/results/chanpak_five_core_bridges_20261110/history/audit_compiler_selector.py
```

[compiler_selector_audit.json](compiler_selector_audit.json) records the
universes, counts, rank, hostile m=0 case and script hash. The separate
[five-core census](receipt.json) concerns the new historical bridges;
it is not part of the producer proof dependency.

## Seven-persistent-vertex reset

**VERIFIED independently for the literal construction; PROVED algebraic
family boundary.** This audit covers sections1–2 of
[RESET](../alternatives/RESET.md). The broader replacement censuses and
neutral edge component in its sections3–4 were not independently rerun here.

The audit constructs the old source from the historical rooted word
`11^k,01,11,10`, using the history lane's constructor. It imports no producer
constructor, attachment, isomorphism routine or path-count implementation.
Independent subset DP checks every intermediate graph for k1..6. Every
arc between the persistent old vertices agrees at every stage, and all
six stages are strong. At k2, the exact count trace is

    141 -> 33 -> 69 -> 141 -> 215 -> 425.

The original vertices `[0,1,2,3,4,5,7]` correspond to the new source U
vertices `[1,2,3,6,7,8,0]`, respectively. The two reset vertices correspond
to U4,U5. This checks the full isomorphism, not just H. The receipt retains
all adjacency rows at every stage, every attachment mask and the label map.

For the all-k calculation, put D=2^k and t=2D+1. The binary word
`1 0^k 1 0^r` is the integer2^r*t. Since g(t)=g(D)=D+1, the first doubling
adds t+1 and all later doublings add2^j*t. Therefore, for k,r>=1,

    g(2^r*t) = D+1+(2^r-1)*t+1
              = (2^(r+1)-1)*D+2^r+1.

The deletion and the two reset insertions have r2,3,4, respectively;
their counts are7D+5,15D+9,31D+17. The last two bits remain00, so the
converse1-ready continuation from ALTERNATE141 gives3(31D+17)+2.
Equality with the original source count34D+5 requires D4, hence exactly
k2. This independently confirms both the positive141 reset and its
first failed implication as a uniform source-return construction.

Reproduce with

```text
python -B 05-knowledge/results/chanpak_five_core_bridges_20261110/history/audit_reset.py
```

The receipt is [reset_audit.json](reset_audit.json). These six finite
controls support the explicit all-k calculation; they are not an
all-order search for other resets.
