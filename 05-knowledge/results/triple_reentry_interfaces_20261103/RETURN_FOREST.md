# The six-root ternary forest and the missing lift

**PROVED:** P1–P3 below are elementary statements about integer addresses and
certified interfaces. **FINITE-EXACT:** literal tournaments realize the six
root counts. **OPEN:** a source class and legal return maps satisfying P3.
Raw Hamiltonian-path counts are used throughout, never `H/Aut`.

## P1. Deleting the two forbidden values leaves exactly six roots

On positive odd integers put

    f_j(h)=3h+2j,       j=0,1,2.

This directed graph is the rooted ternary tree with root 1. To see this,
write `u=(h-1)/2`. The three maps become

    u -> 3u+1, 3u+2, 3u+3.

For `h>1` there is a unique parent: choose `j` so `h-2j` is divisible by
three, then `h_parent=(h-2j)/3` is a smaller positive odd integer. Equivalently
`j=((h-3)/2) mod3`. This is bijective ternary addressing, with digits 1,2,3
on the nonnegative `u` coordinate; it is not an ordinary base-three digit
shift on `h`.

Delete the vertices 7 and 21 and their incident edges. The surviving forest
has precisely the roots

    R = {1,23,25,63,65,67}.

Indeed 7 has children 21,23,25 and 21 has children 63,65,67. The deleted
vertices are adjacent. Every other surviving non-root retains its unique
parent, and inverse descent terminates at one of the six printed roots.
The roots are forced: none has a surviving incoming edge in this grammar.

Thus a completely closed ternary system starting at 1 would be wrong: its
`j=2` branch creates 7. A faithful system must exclude that transition and
start new components at the five nontrivial surviving children of the
removed nodes. If the singleton source is kept outside the main carrier,
the corresponding seven roots are `{3,5,23,25,63,65,67}` instead.

This is an exact arithmetic forest. It does not identify tournaments with
integers or prove a single new tournament realization by itself.

## P2. The six scalar roots have literal tournament witnesses

`scalar_seeds.json` stores outgoing row bitmasks, with original vertex names.
Independent permutation enumeration gives:

| Raw H | Vertices | Provenance |
|---:|---:|---|
|1|1|singleton|
|23|6|one birth on a labelled transitive five-vertex parent|
|25|6|one birth on parent rows `(4,1,2,7,15)`|
|63|8|inherited S124 constructor, then removal of a universal source|
|65|7|one birth on rows `(8,1,3,6,15,31)`|
|67|7|same six-vertex parent, a different named cut|

For 63 the inherited recipe has seed `(1,0,0)` and blocks
`1,0,0,1,1,0,0`. Its nine-vertex output has a universal source; deleting
that forced first vertex preserves H. The remaining strong tournament has
rows `(126,252,248,112,96,192,128,25)`. This recovers existing realizability,
not a new discovery of count 63. The older THM-4094
[prime-lane reduction](../../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md)
already uses63 as a strong carry atom; the new forest is a different,
sufficient return-lifting organization, not a replacement for that stronger
existing multiplicative reduction. Credit:
[the million-count constructor](../session124_million_path_spectrum_20261031/construct_witness.py).
The new verifier does not depend on running that constructor.

These witnesses do **not** certify the observer/guard predicate of the
current ternary source family. That qualification is the point of P3.

## P3. A precise sufficient return-lifting theorem

Let `C` be a class of actual marked tournaments, with raw-count map `H` and
explicitly recorded shared ports and any other coordinates used by later
operations. Suppose:

1. There is a certified member of C at each count in R.
2. For every `T in C` and each `j` with `3H(T)+2j notin {7,21}`, a terminating
   legal construction produces a member `Return_j(T) in C` with
   `H(Return_j(T))=3H(T)+2j`.
3. The construction transports the declared marks, ports and source data;
   its count and membership certificate apply to that actual output.

Then every positive odd integer other than 7 and 21 is a raw tournament
Hamiltonian-path count. Proof: descend its unique P1 address to R, start
with the certified root source, and replay the finite address using (2).
The already proved all-order exclusion supplies the other direction.

Condition (3) is stronger than bare existence needs, but matches the
source-preserving grammar under investigation. If a different replacement
source is used, its membership and count must be authenticated; numerical
equality alone does not transfer a source certificate. No uniform runtime
bound is required for this existence theorem. Any claimed efficient
construction needs its own cost bound.

This is a sufficient theorem schema, not a newly proved return mechanism.
One could also prove completeness by a different grammar without satisfying
it. The six-root reduction is useful because it cleanly separates a finite
seed problem from a universal closure problem.

## What the current return construction already proves

The inherited [moving-guard theorem](../session163_ternary_triples_20261103/strong_returns/MOVING_GUARD.md)
gives genuine all-depth `3H` returns with a moving guard. It is stronger
than a bounded number of successful zero branches.

The [first-step analysis](../session163_ternary_triples_20261103/strong_returns/FIRST_STEP.md)
explains the nonzero obstruction. Its general scalar count is `3H+2G`.
For the marked observer triple the correction is

    (2G_0, G_1-G_0, G_2-G_1),       G_0 <= G_1 <= G_2.

The desired nonzero offset is `j(2,-1,0)`. At a fixed observer port its
negative second coordinate is incompatible with this gadget's nonnegative
second correction. Moreover, its automatically propagated guard has zero
charge, so a one-off scalar `+2` or `+4` does not give a repeatable nonzero
branch. These are source-coordinate obstructions, not failures of ternary
arithmetic.

The incoming S164 complementary construction now repeats scalar offsets and
provides universal entry to a canonical source bank, with a count-changing
price. The old exclusion applies to its specified ear, not to every changed
attachment. [INCOMING](INCOMING.md) integrates those stronger results and
measures the price on all six literal seeds.

The next decisive target is therefore a return with either (a) an explicitly
transported new observer chart that changes this sign constraint, or (b) a
larger interface whose joint path-cover data permits the required correction.
It must survive a second nonzero return. First test the same literal source
and named ports for words `11,12,21,22`; scalar-only success is insufficient.

## Why checking every modulus still does not fill the gap

The inherited [S157 universal construction](../session157_gaps_20261103/universal/NOTE.md)
already realizes every eligible residue modulo every modulus by strong
one-apex tournaments (odd residues when the modulus is even). Therefore
both 7 and 21 are congruent to attainable strong counts modulo **each**
modulus, although neither is attainable exactly. Taking any finite list of
moduli does not help: apply the theorem to their least common multiple.

This is a direct, problem-internal local/global counterexample. The
realizing source is allowed to change with the modulus. No fixed integer
witness or compatible finite source is thereby produced.

For a *fixed finite* interface bank the situation is different. Its output
set is finite, so an absent target is excluded by some sufficiently large
modulus (take a modulus exceeding every nonzero difference from that target).
The prime-support compiler in README handles its unit predicate over all
moduli exactly. Increasing the bank changes this finite certificate; it
must not be silently treated as a uniform source-return theorem.
