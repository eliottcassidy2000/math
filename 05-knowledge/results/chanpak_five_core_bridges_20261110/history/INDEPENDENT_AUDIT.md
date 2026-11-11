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
