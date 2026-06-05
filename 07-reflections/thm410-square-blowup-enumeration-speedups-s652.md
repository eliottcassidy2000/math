# THM-410 Square-Blowup Enumeration Speedups S652

The user's "square the node count; each node becomes a tournament of the
original size" is the tournament substitution product in disguise.  That was the
right tangent because it separates two enumeration questions that raw search
usually conflates:

```text
How many prime/indecomposable tournament classes exist?
How many imprimitive block templates exist once the module system is retained?
```

A000568 answers the first question only after quotienting every label
symmetry.  The square-blowup carrier asks a smaller companion question: choose
an `n`-vertex base tournament and assign an `n`-vertex tournament class to each
base vertex, modulo the automorphisms of the base.  Burnside gives this without
ever building an `n^2`-vertex labelled tournament:

```text
SqTempl(n) = sum_B (1/|Aut(B)|) sum_g A(n)^{cycles(g)}.
```

The first values from S652 are:

```text
1, 1, 12, 704, 2127984, ...
```

This is not A000568 at `n^2`.  It is a decomposable shadow of it, and that is
the point.  It is the modular-decomposition layer that should be peeled before
canonicalizing a giant blowup.

The Hamiltonian-path side initially tempts a false product formula:

```text
H(B[F_i]) ?= H(B) * product_i H(F_i).
```

S71o already warned that lexicographic product is not generally multiplicative.
The S652 replacement is the exact run formula:

```text
H(B[F_i]) =
  sum_r MacroWords_B(r) * product_i OrderedPathCovers(F_i, r_i).
```

The direct 9-vertex checks are the important sanity test.  A cyclic 3-base with
three transitive 3-blocks has `H=2721`, while the naive product predicts `3`.
The path-cover/macro-word formula returns `2721` exactly.  That is the missing
interleaving mass: Hamiltonian paths are not just paths through blocks; they can
return to a block many times when the macro tournament is cyclic.

The THM-410 half of the session is a different kind of speedup.  It says that a
matching of reversed arcs off a transitive order has exact triangle count

```text
c3 = sum interval-interiors.
```

This is a theorem-grade pruning key.  For the H=21 window, where `H >= 1+2c3`
forces `c3<=10`, interval-matching states can be counted by a small weighted
matching DP.  At `n=14`, all matchings are `2,390,480`, but only `188,273` have
THM-410 weight at most `10`.  Again, this is not all low-`c3` tournaments.  It
is a cheap exact carrier that can seed, bound, or sanity-check broader searches.

The S9 high-gap unlock result is the cautionary backdrop.  All 157 high
permanent-through-`n=9` H gaps unlocked at `n=10` under biased sampling.  That
means high-end absence often reflects search geometry, not obstruction.  S652's
working rule should be:

```text
sample aggressively for witnesses;
prune only with exact carriers;
record which quotient was destroyed.
```

For Tournament Analysis, the challenged assumption is especially sharp here.
The vertices of the speedup tournament are not runners or arcs by default.  They
can be reversed intervals, modules, block-run words, path-cover segments,
automorphism cycles, H-values, or proof obligations.  Each vertex choice keeps a
different invariant and throws away a different side channel.  The enumeration
tool should name that quotient explicitly before trusting the speedup.

The next implementation target is a real substitution enumerator: compute
path-cover polynomials for all classes up to `n=6`, cache base cycle indices,
and generate square-blowup H spectra by template rather than by labelled
Held-Karp.  That would make the new sequence useful instead of merely pretty.
