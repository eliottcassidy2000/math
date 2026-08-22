# Dyadic fibres: a ternary carry inside the binary tree

## Status

**SYNTHESIS AROUND PROVED THM-3435.**  The fibre theorem and exact companion
have now passed an independent derivation and hostile audit.  This reflection
adds interpretations, not further proved dependencies, and gives no LRC(14)
decrement.

## 1. Inheritance board

The closest proved mechanisms are:

- THM-3432, which splits a fixed-zero block on `2R` into even zero and odd
  half charts;
- THM-3429/3434, which distinguish inactive pullbacks from active complete
  prime grids and retain an affine lift sidecar;
- THM-3416, whose two maximal half-block formulas already warn that normalized
  coefficient parity survives quotient order;
- the exact `Q=8` rank-four atom, which was known as a cover but not isolated
  as a four-branch Boolean partition.

The corrected near miss was “an odd coefficient gives one section of every
dyadic fibre.”  It is true only through active depths one and two.  At depth
three, `(Q,r)=(8,1)` already has two points on its sole eight-fibre.

The least-used sidecar was the affine **start** of the cyclic interval on each
fibre.  Counts and extra masks retain its length but not its location.

The live board became

```text
inactive depth | active grid | cyclic interval | mod-7 carry |
Boolean quartet | eight-branch overlap | affine start
```

## 2. The exact decomposition

Write `Q=2^aN`, `N` odd, and `r=2^bs` with `b=min(v_2(r),a)`.  If `b=a`, the
mask is a full `2^a`-fold pullback from `N`.  Otherwise remove those `b`
inactive levels and put

```text
M=2^(a-b)=7q+c,                    c in {1,2,4}.
```

On the reduced fibre `y+Nt`, the change `n=st mod M` turns the mask into one
open cyclic interval of length `M/7=q+c/7`.  It therefore contains `q` or
`q+1` integers.  The extra-point bases are exactly

```text
E=H_N^(c)(s+qN),
```

where the superscript `c` means radius `c/14`.  The original fibre count is
`2^b(q+1_E(y))`.  This is both more rigid and
more limited than a descent theorem: it retains every count and the cyclic
order, but a cover still depends on how the intervals start relative to one
another.

The strict-boundary split is clean.  In an active branch, an endpoint would
equate quantities of two-adic valuations `b<a`, impossible.  Inactive
pullbacks retain odd-base endpoints; `(Q,r)=(14,2)` is the first even-modulus
example.

## 3. Why a ternary recurrence appears

The binary depth `k` is read modulo seven:

```text
2^k mod 7:       2 -> 4 -> 1 -> 2,
carry to q:      0 -> 1 -> 0 -> 0.
```

Thus the state machine has three states because `2` has multiplicative order
three modulo `7`.  Every third level, `q` is odd and the extra coefficient
switches from `s` to `s+N`; at the other two levels it remains `s` while the
radius is `2/14` or `4/14`.

This is the precise version of the speculative “ternary shadow of a binary
branch tree.”  The tree itself is binary.  The **carry automaton** on fibre
mass is ternary.  Conflating those two statements would erase the operation
that creates the period three.

## 4. The size-four Boolean realization

At `Q=8`, four transverse coefficient branches give

```text
coefficient 1  -> {0,7}
coefficient 5  -> {1,6}
coefficient 9  -> {3,4}
coefficient 13 -> {2,5}.
```

They partition the eight sheets.  Hence for every

```text
J subset {1,5,9,13},
```

the map

```text
J |-> union_(r in J) H_8(r)
```

is an injective Boolean-algebra homomorphism: union is OR and symmetric
difference is literal sheetwise XOR.  This is an exact size-four object, but
it is **not** intrinsically a tournament.  Its intersection graph is empty;
an orientation would be imported decoration.

Periodically lift the sheet support to a subset `S_J` of the natural numbers.
Each selected owner contributes two residue classes modulo eight, so

```text
density(S_J)=|J|/4,
sum_(n<=x,n in S_J) 1/n=(|J|/4)log x+C_J+O(1/x).         (1)
```

This realizes a four-bit cylinder algebra inside the harmonic series.  It
does not say every subset of the naturals has a density or a harmonic
coefficient.  An arbitrary subset is always an indicator/XOR object, and its
harmonic subseries is well-defined, but the asymptotic in `(1)` needs the
periodicity present here.

## 5. Why depth three is not another tournament

Eight coefficient preimages are spaced by eighths in phase.  Their arcs have
length `1/7>1/8`, so overlap begins exactly here.  Nevertheless the grid
covering radius is `1/16<1/14`, so all eight branches cover every sheet, with
multiplicity one or two.

On a fixed sheet, multiplication by the odd unit `2ell+1 mod 8` orders the
labels in a local `C_8`; a double hit uses adjacent labels in that order.
The order changes with the sheet.  After forgetting it, any odd label
difference can become adjacent, while an even difference never can.  The
global conflict graph is therefore a subgraph of `K_(4,4)`, not a fixed
cycle and not a tournament.

The first two transverse controls go in opposite directions:

```text
Q=8, seed 1:  every sheet has multiplicity 2; any seven branches cover;
Q=16,seed 1:  deleting branch 0 leaves only 14 of 16 sheets.
```

So “eight branches minus one” is not a cover certificate.  The missing
coordinate is again affine position, now encoded by the sheet-dependent odd
unit and interval start.

## 6. How this changes the frontier

- Even moduli no longer need to be attacked as an opaque new census.  Every
  owner has an exact inactive depth, baseline fibre mass, extra mask, and
  affine-start sidecar.
- For target-free cap-seven work, moduli divisible by eight already lie over
  the rank-four atom.  The genuinely new all-rank lane is concentrated at
  two-adic depth one or two, where active owners are honest partial sections
  and branch fusion is XOR-exact.
- The `1,2,4` remainder cycle is a lawful recurrence connection to ternary
  state machines.  It is not evidence for a C3 tower until a map preserving a
  target predicate is supplied.
- The size-four quartet gives the requested Boolean/harmonic realization,
  while its empty intersection graph gives the corresponding tournament
  no-go.

The next decisive test is not another unrestricted even-modulus search.  It
is a typed cap-seven solver on `Q=2N` and `Q=4N` that stores, for each active
owner, both the expanded-radius base support and its affine section bit.
That is the minimum sufficient state suggested by the decomposition.
