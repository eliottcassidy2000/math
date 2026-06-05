# HYP-2228: THM-410 and Square-Node Blowups Give Enumeration Speedup Carriers

**Status:** OPEN method hypothesis with finite computation evidence
(codex-2026-06-05-S652).

## Claim

THM-410 interval-reversal matchings and square-node tournament substitution are
enumeration carriers.  They do not replace full tournament isomorphism
enumeration, but they expose two compressed slices that raw labelled search sees
only after paying for too many labels.

The two carriers are:

1. **Interval matching kernel.**  Start from the transitive tournament on a
   linear order and reverse a matching `M` of arcs.  By THM-410,

   ```text
   c3(T_M) = sum_{(a,b) in M} (b-a-1).
   ```

   Thus a low-`c3` near-transitive branch can carry an additive interval weight,
   instead of recomputing all triples at every node.

2. **Square substitution kernel.**  Given an `n`-vertex base tournament `B`,
   replace each vertex by an `n`-vertex tournament block.  The resulting
   tournament has `n^2` vertices, but the imprimitive template count is

   ```text
   SqTempl(n) =
     sum_[B in A000568(n)] (1/|Aut(B)|) sum_{g in Aut(B)} A000568(n)^{cycles(g)}.
   ```

   For Hamiltonian paths, naive multiplicativity fails for cyclic bases, but an
   exact formula survives:

   ```text
   H(B[F_i]) =
     sum_r MacroWords_B(r) * product_i OrderedPathCovers(F_i, r_i).
   ```

   Here `r_i` is the number of runs of block `i`; `MacroWords_B(r)` counts
   directed block-word paths in the base; and `OrderedPathCovers(F_i,r_i)` counts
   ordered covers of the fiber by `r_i` directed paths.

## S652 Evidence

`04-computation/thm410_square_blowup_enumeration_s652.py` validates the THM-410
matching formula exhaustively through `n=8` matchings, with zero failures.

The low-`c3` matching DP gives the following cap sizes for the H=21-relevant
bound `c3<=10`:

```text
n   all matchings   c3<=10
8   764             722
10  9,496           5,538
12  140,152         34,284
14  2,390,480       188,273
```

These are not all low-`c3` tournaments; they are a cheap theorem-certified
near-transitive carrier.

The square-substitution Burnside companion sequence through `n=5` is:

```text
n  A000568(n)  A(n)^2 uniform pairs  raw A(n)^(n+1)  square templates
1  1           1                     1               1
2  1           1                     1               1
3  2           4                     16              12
4  4           16                    1,024           704
5  12          144                   2,985,984       2,127,984
```

So the new related sequence is:

```text
1, 1, 12, 704, 2127984, ...
```

It counts imprimitive square templates, not all isomorphism classes on `n^2`
vertices.

The block path-cover formula was checked against direct Held-Karp on 9-vertex
lexicographic products:

```text
chain3[chain3,chain3,chain3]       H=1     naive product=1
cycle3[chain3,chain3,chain3]       H=2721  naive product=3
cycle3[cycle3,chain3,cycle3]       H=2961  naive product=27
chain3[cycle3,chain3,cycle3]       H=9     naive product=9
```

Thus the exact block formula catches the large run-interleaving contribution
that simple product heuristics miss.

## Enumeration Use

The practical speedup thesis is:

- use THM-410's interval weight as a branch-and-bound key for near-transitive
  `c3`-capped searches, including H-spectrum gap work;
- use modular decomposition and base automorphism cycle-index data before
  canonicalizing large blowups;
- compute `H` of substitution products by path-cover polynomials and a
  macro-word DP, not by Held-Karp on all `n^2` vertices;
- treat A000568 as one layer in a family of companion counts: rooted,
  self-converse, q-Burnside, and now square-substitution templates.

## Assumption Challenge

The vertices of the analysis need not be the original tournament vertices.
Useful vertex sets here include reversed intervals, matching arcs, block
modules, macro-runs, ordered path-cover segments, automorphism cycles, H-values,
and proof obligations.

The interval quotient preserves exact cyclic-triangle witnesses for matching
flips, but destroys interactions among nonmatching flip sets.  The square
substitution quotient preserves modular decomposition and block `H` data, but
destroys the prime/indecomposable tournaments outside the imprimitive slice.

## Limits

This is a carrier, not a full enumeration theorem.  THM-410 is exact for
matching reversals; arbitrary flip sets need higher-order corrections.  The
Burnside formula counts square templates modulo base automorphisms; it does not
claim that these are all distinct `n^2`-vertex isomorphism classes.  The block
path-cover formula is exact for substitution products, but its path-cover
polynomials still need their own optimized implementation for larger fibers.

**See:** `04-computation/thm410_square_blowup_enumeration_s652.py`;
`05-knowledge/results/thm410_square_blowup_enumeration_s652.out`;
`07-reflections/thm410-square-blowup-enumeration-speedups-s652.md`; THM-410,
HYP-2227, HYP-2226, HYP-2209, HYP-2200, THM-115.
