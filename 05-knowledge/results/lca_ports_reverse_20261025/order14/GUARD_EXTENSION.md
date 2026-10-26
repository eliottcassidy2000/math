# One source bit selects an order-fourteen three-black ledger

2026-10-25. **PROVED, root proof audit and independent FINITE-EXACT
controls passed:** on each of
the four full source cells in `README.md`, the new relocation ledger
works whenever T(3,5)=0. When T(3,5)=1, the inherited ledger works.
Thus the decision can use one original arc; it need not test the whole
failed-word anti symmetry. This concerns only the four stated order14
cells and preserves literal vertices, ports, and all clone adjacencies.

## Scope and legality

Use exactly the four28-dimensional source cells and words from
`README.md`. They all have incoming ports containing10 and omitting
2,4,8,11, directed closed blue P=(0,...,13), and the closed middle12->1.
The two assigned chord orientations remain part of their source data.

Assume T(3,5)=0. The new Q1 requires2->13 and13->3, which follow by
P-anti from0->11 and10->0. Q2 then requires4->0 and3->5, which are the
complements of omitted source pairs0->4 and5->3. Q3 restores2->3 and
12->13 after two complements. All other joins are retained. The same
edge rectangle therefore reaches the exact B=K_P(T). None of these
steps uses F-blue once the one-bit guard is supplied.

## Three uniform pair witnesses

1. At T, Q1's actual reversal sends the ordered pair(0,11) to(1,12).
   Their values are1 and0:11 is not a port, and the middle is closed.

2. At S=K_Q1(T), Q2's actual reversal sends(0,8) to(13,11). Both
   unordered pairs are omitted by Q1. Since8 is not a port,
   T(0,8)=1. P-anti gives T(13,11)=T(2,0)=0, because2 is not a port.
   Consequently S(0,8)=0 and S(13,11)=1.

3. At B, Q3's actual reversal sends(0,11) to(1,12). These are both
   nonpath pairs of P, so B(0,11)=0 and B(1,12)=1. Its blackness at B
   is equivalent to blackness at the actual preceding source, because
   positional reversal preserves its own path-edge set.

Each mismatch disproves the relevant anti condition, proving all three
words black. The pair witnesses do not use the exceptional tournament's
degree sequence. They hold uniformly on every guarded completion.

## The one-bit decision is exhaustive

The original symmetry classification proves F-blue implies T(3,5)=0.
Its contrapositive shows that T(3,5)=1 makes F black. The inherited first
two words are already uniformly black and legal on each full source cell;
its literal target is B. Therefore select at the original source T:

    T(3,5)=0: use the new relocation rectangle;
    T(3,5)=1: use the inherited branch-specific ledger.

This partitions each28-dimensional cell into two nonempty27-dimensional
affine cells. It is a sufficient selector, not an equivalence between
T(3,5)=0 and F-blue: each zero-guard half-cell has many completions,
whereas exactly one completion per full cell has blue F. The selector
does not splice prefixes of different ledgers.

The original missing-symmetry hostile flips{3,5} and its reflected pair
{8,10}; it crosses from the zero guard to the one guard, where the new
ledger becomes illegal and the inherited one succeeds. Thus the guard
is substantive.

The independent parity-equation code embedded in
[INDEPENDENT_REVIEW.md](INDEPENDENT_REVIEW.md) certifies all completions
of all eight guard half-cells, their dimensions, exact legality and target,
and actual-stage blackness. It imports no producer and passes normally
and under `-O`. The original producer files are unchanged.
