> **CORRECTED HISTORICAL SYNTHESIS (2026-08-25).**  The former version
> conflated cubic edge colorings with `Z/3Z` flows, primality with Eisenstein
> splitting type, and the geometry of lines on cubics with an incorrect
> expected dimension.  MISTAKE-507 records the correction.  The analogies
> below are not theorem transfers.

# Fano Schemes, Eisenstein Splitting, and Snark Obstructions

## Three correctly typed objects

For a cubic hypersurface of dimension `n` in projective `(n+1)`-space, the
Grassmannian of lines has dimension `2n`, and restricting a cubic to a line
gives four coefficient conditions.  Thus the expected dimension of its Fano
scheme of lines is `2n-4`.  A smooth cubic surface (`n=2`) has a zero-
dimensional scheme of 27 lines; it is not an expected-dimension `-1`
exception.

For a **rational prime** `p`, its splitting in the Eisenstein integers is
determined by `p mod 3`: `p=3` ramifies, primes congruent to `1 mod 3` split,
and primes congruent to `2 mod 3` remain inert.  This classification assumes
primality first.  A residue class does not test primality: `7` and `25` are
both `1 mod 3`.

For a cubic graph, a proper three-edge-coloring labels the three incident
edges at every vertex by the three nonzero elements of `F_2^2`; their sum is
zero.  Conversely, a nowhere-zero `F_2^2` flow gives those three colors.
Hence the exact standard equivalence is with a nowhere-zero 4-flow, not a
nowhere-zero `Z/3Z` flow.

## What survives the comparison

All three subjects involve an extension or realization locus, but their
interfaces differ:

| Object | Local data | Global question | Faithful carrier |
|---|---|---|---|
| line on a cubic | restriction coefficients | does the line lie on the cubic? | section of `Sym^3(S^*)` |
| rational prime in `Z[omega]` | residue after primality is known | split, inert, or ramified? | ideal factorization |
| cubic graph piece | colors on a cut | does the coloring extend and glue? | boundary extension relation |

No map here preserves all three predicates.  The proof-grade snark mechanism
is THM-4116: a piece is replaced by its ordered boundary extension-count
vector, and gluing is an exact dot product.  This retains extendability while
discarding interior geometry; multiplicities, cut order, and color gauge are
the required sidecars.

The old “local -> algebraic -> global” hierarchy remains a useful question-
asking mnemonic, but it carries no implication between the three theories.
