# THM-410 Square-Blowup Enumeration Speedups S652

The user said: square the node count; each node becomes a tournament of the
original size.  That is lexicographic substitution, and it is exactly the kind
of object where a raw enumerator can make itself unnecessarily miserable.

Flatten it, and the result is an `n^2`-vertex tournament.  Held-Karp sees a
huge subset lattice, and isomorphism search sees a huge symmetric group.  Keep
the carrier, and the same object is:

```text
base tournament + tournament blocks + automorphism cycles + macro-word runs.
```

That distinction separates two questions that raw search usually conflates:

```text
How many prime/indecomposable tournament classes exist?
How many imprimitive block templates exist once the module system is retained?
```

A000568 answers the first question only after quotienting every label symmetry.
The square-blowup companion asks a smaller question: choose an `n`-vertex base
tournament and assign an `n`-vertex tournament class to each base vertex, modulo
the automorphisms of the base.  Burnside gives this without building an
`n^2`-vertex labelled tournament:

```text
SqTempl(n) = sum_B (1/|Aut(B)|) sum_g A(n)^{cycles(g)}.
```

The first values from the S652 Burnside scout are:

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

S710 already warned that lexicographic product is not generally
multiplicative.  S652 makes the replacement explicit:

```text
H(B[F_i]) =
  sum_r MacroWords_B(r) * product_i OrderedPathCovers(F_i, r_i).
```

A Hamiltonian path in a strong base can leave a block and return to it.  So the
unit is not a block path; it is a block run.  The direct checks make this hard
to ignore:

```text
cycle3[chain3,chain3,chain3]: H=2721, naive product=3.
C3[C3]:                         H=3159, naive product=81.
C3[T2,T2,T2]:                    H=45,   naive product=3.
```

The supplemental scout also gives the scale smell test.  For a 25-vertex square
example, a raw Held-Karp proxy is about `8.4e8` states, while the module-state
proxy is about `4.0e4`.  By `n=7`, the proxy gap is `2.758e16` raw states
versus `1.470e7` module states.  That is not magic; it is just what happens
when the computation does not throw away the module carrier.

THM-410 contributes the other half of the picture.  It says that if we start
from a transitive tournament and reverse a matching of long edges, `c3` is an
additive interval sum:

```text
c3 = sum interval-interiors.
```

That is a pruning oracle with theorem backing.  Instead of saying "try edge
masks with low `c3`," say:

```text
try interval ledgers whose interior sums stay under budget.
```

At `n=10`, there are only `9496` interval matchings, compared with `2^45`
labelled tournaments and `9733056` isomorphism classes.  The Burnside scout
pushes the same weighted-matching idea further: at `n=14`, only `188273` of
`2390480` matchings have THM-410 weight at most `10`.  This is not an exhaustive
low-`c3` enumerator.  It is a certified near-transitive witness generator.

The arbitrary-upset formula widens the lane.  Relative to a fixed linear order,
for `x<y<z`, the triple is cyclic exactly when `r_xy=r_yz!=r_xz`.  That gives a
bitset side channel for fixed-order deformations.  It does not solve
isomorphism, but it makes the `c3` observable cheaper and more typed.

The S9 high-gap unlock result is the cautionary backdrop.  All 157 high
permanent-through-`n=9` H gaps unlocked at `n=10` under biased sampling.  That
means high-end absence often reflects search geometry, not obstruction.  The
next move should be to replace "biased sample" with "typed witness menu":

```text
THM-410 = additive witness cache.
Upset formula = fast c3 side channel for a fixed order.
Square blowup = module-preserving H automaton.
Burnside templates = global orbit side channel.
```

For A000568, this does not replace orderly generation.  It says that a good
enumerator should first ask whether the candidate has a modular quotient,
whether its `c3` comes from a cheap interval/upset certificate, and whether a
decorated fiber canonical form will separate it before full `S_N`
canonicalization.

For Tournament Analysis, the challenged assumption is especially sharp here.
The vertices of the speedup tournament are not runners or arcs by default.  They
can be reversed intervals, modules, block-run words, path-cover segments,
automorphism cycles, H-values, or proof obligations.  Each vertex choice keeps a
different invariant and throws away a different side channel.  The enumeration
tool should name that quotient explicitly before trusting the speedup.

The next implementation target is a production substitution enumerator: compute
path-cover polynomials for all classes up to `n=6`, cache base cycle indices,
and generate square-blowup H spectra by template rather than by labelled
Held-Karp.  The rule of thumb is short enough to keep:

```text
speedups are side channels.
```
