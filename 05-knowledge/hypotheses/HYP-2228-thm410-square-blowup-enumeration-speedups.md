# HYP-2228: THM-410 and Square-Node Blowups Give Enumeration Speedup Carriers

**Status:** OPEN method hypothesis with S652 computational evidence.

## Claim

Tournament enumeration gets faster when the quotient remembers the witness
carrier.  THM-410 remembers reversed intervals.  Arbitrary upsets around a
linear order remember a bitset deformation carrier.  Square-node substitution
remembers modules, automorphism cycles, and block-run words.  If those carriers
are flattened into raw edge masks or raw `S_N` canonicalization too early, the
computation loses the structure that made it small.

S652 records two compatible speedup lanes.

1. **Interval and upset carrier.**  Start from the transitive tournament on a
   linear order.  If a matching `M` of arcs is reversed, THM-410 gives

   ```text
   c3(T_M) = sum_{(a,b) in M} (b-a-1).
   ```

   For a general upset set relative to the same order, each `x<y<z` is cyclic
   exactly when

   ```text
   r_xy = r_yz != r_xz.
   ```

   Thus low-`c3` and near-transitive searches can keep interval ledgers or
   bitset side channels before scanning all triples.

2. **Square substitution carrier.**  Given an `n`-vertex base tournament `A`,
   replace each vertex `i` by an `n`-vertex tournament block `B_i`.  Between
   blocks, all arcs follow `A`.  The result has `n^2` vertices, but still has a
   visible module carrier.

   For the uniform square `Sq(T)=T[T]`, S652 verifies

   ```text
   score((i,a)) = n * score_T(i) + score_T(a)
   c3(Sq(T)) = (n^3 + n) * c3(T).
   ```

   For general substitution,

   ```text
   c3(A[B_i]) =
     sum_i c3(B_i)
     + sum_{directed 3-cycles abc in A} |B_a| |B_b| |B_c|.
   ```

   Hamiltonian paths require an additional run carrier:

   ```text
   H(A[B_i]) =
     sum_r MacroWords_A(r) * product_i OrderedPathCovers(B_i, r_i).
   ```

   Here `r_i` is the number of runs of block `i`; `MacroWords_A(r)` counts
   directed block-word paths in the base; and `OrderedPathCovers(B_i,r_i)`
   counts ordered covers of the fiber by `r_i` directed paths.

## S652 Evidence

`04-computation/thm410_square_blowup_enumeration_s652.py` validates the
THM-410 matching formula exhaustively through `n=8` and uses a weighted matching
DP to count the H=21-relevant cap `c3<=10`:

```text
n   all matchings   c3<=10
8   764             722
10  9,496           5,538
12  140,152         34,284
14  2,390,480       188,273
```

These are not all low-`c3` tournaments.  They are a cheap theorem-certified
near-transitive carrier.

The same script adds a Burnside companion sequence for square-substitution
templates.  It chooses an `n`-vertex base tournament and assigns an `n`-vertex
tournament class to each base vertex, modulo base automorphisms:

```text
SqTempl(n) =
  sum_[B in A000568(n)] (1/|Aut(B)|) sum_{g in Aut(B)} A000568(n)^{cycles(g)}.
```

The first values through `n=5` are:

```text
1, 1, 12, 704, 2127984, ...
```

This counts imprimitive square templates, not all isomorphism classes on `n^2`
vertices.

`04-computation/thm410_square_blowup_speedups_s652.py` adds a supplemental
bounded scout.  It verifies the THM-410 matching cache through `n=8`, verifies
the arbitrary-upset rule over all labelled tournaments through `n=6`, and
checks the uniform-square score and `c3` formulas over all labelled tournaments
through `n=5`.

The two Hamiltonian-path scouts agree on the key warning: naive multiplication
fails for strong outer quotients, and the missing mass is block-run
interleaving.

```text
cycle3[chain3,chain3,chain3]: H=2721, naive product=3
C3[C3]:                         H=3159, naive product=81
C3[T2,T2,T2]:                    H=45,   naive product=3
```

The supplemental scout also records the square-module state proxy:

```text
n  total vertices  raw Held-Karp states  module-state proxy  raw/module
3               9                   4608                 273       16.88
4              16              1.049e+06                2824      371.31
5              25              8.389e+08               40095    20921.83
6              36              2.474e+12              710268  3483053.10
7              49              2.758e+16           1.470e+07 1877090681.38
```

The proxy is not a generic speedup for random unmarked tournaments.  It is the
payoff for retaining the substitution carrier.

## Enumeration Program

For H-spectrum gap hunting:

- use S9's lesson that sampling proves achievability but not permanence;
- build certified witness menus: THM-410 intervals, general upset bitsets,
  square/module substitutions, and other modular decompositions;
- compute `H` of substitution products by path-cover polynomials and a
  macro-word DP, not by Held-Karp on all `n^2` vertices;
- record which quotient proves the count, rather than recording only the final
  `H`.

For A000568-style isomorphism-class enumeration:

- run modular decomposition before full canonicalization;
- canonicalize the prime quotient and decorated fiber types before falling back
  to `S_N`;
- cache base automorphism cycle indices and decorated-fiber Burnside counts;
- use score, `c3`, block path-cover polynomials, macro-word automata, and
  upset bitsets as fast reject/cache keys.

## Tournament Analysis

The vertices of the analysis need not be the original tournament vertices.
Useful vertex sets here include reversed intervals, matching arcs, upset
patterns, block modules, macro-runs, ordered path-cover segments, automorphism
cycles, H-values, and proof obligations.

The interval quotient preserves exact cyclic-triangle witnesses for matching
flips, but destroys interactions among nonmatching flip sets.  The upset
bitset quotient preserves a fixed-order deformation ledger, but destroys
order-free canonical information.  The square substitution quotient preserves
modular decomposition and block `H` data, but destroys prime/indecomposable
tournaments outside the imprimitive slice.

## Limits

This is a carrier, not a full enumeration theorem.  THM-410 is exact for
matching reversals; arbitrary flip sets need the broader upset rule or
higher-order bookkeeping.  The Burnside formula counts square templates modulo
base automorphisms; it does not claim that these are all distinct `n^2`-vertex
isomorphism classes.  The block path-cover formula is exact for substitution
products, but its path-cover polynomials still need optimized production
implementations for larger fibers.

**See:** `04-computation/thm410_square_blowup_enumeration_s652.py`;
`05-knowledge/results/thm410_square_blowup_enumeration_s652.out`;
`04-computation/thm410_square_blowup_speedups_s652.py`;
`05-knowledge/results/thm410_square_blowup_speedups_s652.out`;
`07-reflections/thm410-square-blowup-enumeration-speedups-s652.md`; THM-410,
HYP-2227, HYP-2226, HYP-2209, HYP-2200, THM-115, OPEN-Q-055.
