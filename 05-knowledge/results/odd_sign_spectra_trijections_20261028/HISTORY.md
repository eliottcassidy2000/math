# Inheritance, repairs, and typed connections

The session started from snapshot 3eb49fa63d7f in a clean isolated worktree.
It recovered current canon before historical notes. No chronological log was
used as a truth source. A concurrent fetch then brought in Session90's
[even-graph / odd-descent synthesis](../session90_even_graphs_odd_descent_20261029/README.md),
based on the same starting snapshot. Its claims are credited as inherited
concurrent work, not presented as discoveries of this packet.

Before close-out, incoming commit3c3861df0d supplied
[M1-M3 of the minimal-blue packet](../minimal_blue_selectors_20261028/THEORY.md).
The previous singleton minimality question is now closed for its declared
language: a four-singleton endpoint decoder proves necessity, and the unique
signed boundary offset reduces the extra interface to4n-8 states. We read
the decoder and necessity/sufficiency argument and updated the current routes;
we do not claim this concurrent result or independently rerun its census here.

## What was already present

- THM-4721 [Walsh duality](../../../01-canon/theorems/THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md)
  and the [October19 even-graph note](../even_graphs_20261019_sylow_refinement_components_walsh_duality.md)
  already explain ordinary equinumerosity through signed orbit lines.
  THM-4531 [twisted Mallows-Sloane](../../../01-canon/theorems/THM-4531-twisted-mallows-sloane-principle-and-no-natural-bijection-for-switching-classes.md)
  already supplies the ordinary converse trace.
- THM-4722 [Sylow Hall obstruction](../../../01-canon/theorems/THM-4722-no-equivariant-bijection-tournaments-even-graphs-for-every-n-sylow-hall-violator.md)
  blocks a natural equivariant class-bijective map for n>=5.
  This packet does not circumvent it by renaming a linear basis correspondence.
- THM-4835 [marked carriers](../../../01-canon/theorems/THM-4835-euler-characteristic-criterion-for-marked-excess-minus-one-carriers-in-tournaments.md)
  proves the monomial Euler criterion, arbitrary-tree end-edge map, functional
  root-swap map, bigon weak/nonlinear separation, and theta nonweak family.
  Its repaired tier law is essential: an odd subgroup defect alone is not
  sufficient to refute weak trijectivity.
- THM-4836 [coherent carriers](../../../01-canon/theorems/THM-4836-coherent-carrier-sets-of-tournaments-and-the-pair-or-triple-law.md)
  concerns which coherent carrier kinds occur in a tournament. It is not the
  claim that two or three unmarked sparse graphs encode every tournament.
- [Session81 graph carriers](../session81_counterexample_families_20261028/GRAPH_CARRIERS.md)
  separates a marked completion chart from a larger faithful sparse graph
  encoding. Session90 adds exact-symmetry codecs and sharper padding bounds.
- The [October20 trijection reset](../trijection_reset_20261020/README.md)
  and [graph trijections](../tournament_graph_trijections_20261020/README.md)
  retain order and host marks. Forgetting those marks is not an inverse map.
- Session89's odd self-converse lift and stronger low-bit formulas precede
  P2. Session90's joint-decoration Fourier theorem precedes this integration.
  Our new outputs are the explicit two-layer sign cancellation census, free
  flip actions, and their operation-specific spectral/geometric applications.

## Older inspirations and their repairs

The July reflection [One word, two even graphs](../../../07-reflections/one-word-two-even-graphs.md)
already warns that degree-even and automorphism-sign-even graphs are different.
The present carrier polynomial keeps the sign relative to the joint stabilizer,
which is a further distinction from global graph evenness.

The June3 [parity-sector reflection](../../../07-reflections/lrc-tournament-parity-sectors-worry-is-the-odd-sector-s587.md)
suggested that the LRC worry-set was the self-converse sector.
That identification was not supplied with a witness-preserving map.
Its heading and maintained warning now distinguish tournament parity facts
from the proposed LRC transfer. L1-L2 give an exact, limited comparison;
they do not purport to refute every possible future definition of a worry-set.

The June3 [unit-distance synthesis](../../../07-reflections/unit-distance-n22-and-its-deep-tie-to-tournaments-and-lrc-s599.md)
combined a triangular-lattice construction with unrestricted language.
The degree-six step fails without separation; arbitrary oriented triangles
need not be cyclic; its pair correlation is not standard additive energy.
A warning preserves its provenance and marks the missing implication.
The old 22-point maximum is not adopted as a current unrestricted result.

## A common algebraic mechanism, with different objects

For a binary matrix M and a selected support x, the zero-syndrome indicator is

    1_(Mx=0) = 2^(-r) sum_(y in F2^r) (-1)^(y dot Mx),

where r is the number of rows. This follows by pairing characters when Mx
is nonzero. It applies to arithmetic square classes: each prime valuation
modulo two gives one row, and square products are kernel supports.
It is an orthogonal projection, just as averaging a group action projects
onto invariant functions in S1. Both retain relations by testing all required
characters; they need not have the same matroid or the same object set.

A single sign can lose most of that information. For example Liouville's
lambda(n)=(-1)^(sum_p v_p(n)) is one linear character of the valuation parity
vector. It gives +1 at 6 although 6 is not a square. Mobius additionally
tests whether any valuation exceeds one, information absent from parity:
3 and 12 have the same square class, but mu(3)=-1 and mu(12)=0.
Thus the full parity vector, its total parity, and squarefree support are
three different representations. These elementary definitions explain why
a sign cancellation is not automatically an arithmetic reconstruction.

The earlier
[consecutive flags / signed incidence analysis](../consecutive_flags_order10_20261028/THEORY.md)
uses actual equations for graceful differences and square sums. Those
equations retain endpoint coordinates and signs; the current Walsh count
does not solve their distinct-label or prescribed-value constraints.
No graceful-tree theorem follows merely from a shared incidence notation.

## Map ledger

| Connection | Source -> target and map | Preserved predicate / loss / sidecar | Cheapest decisive test |
|---|---|---|---|
| Even graphs | Completion invariant functions -> signed support lines, Walsh transform | Counts and specified linear operators; loses a positive object matching; retain signs and A | K2 endpoint swap |
| Sparse trijection | Completion cube -> pair cube on E, if the Euler criterion holds | Individual flips and free complement; retain D and the action | Bigon k=2 degree-one mismatch |
| Paired oddness | Aligned classes -> odd self-converse classes, inherited Eplett codec | Class and transported flips; full Aut is replaced by a centralizer at a mark | n=2 negative-times-negative support |
| LRC reflection | Good-time components -> permutation representation of reflection | Component trace; discards positions, lengths, and existence at trace zero | Three (1,2) thresholds |
| Geometric theta | Abstract path carrier -> algebraic planar point configuration | All unit/nonunit pairs for length3; forgets oriented automorphism action | Length2 K_(2,3) and rigid eight-vertex carrier |
| Binary arithmetic | Selected factors -> valuation parity syndrome | Square-product predicate; discards exponents and squarefreeness | 6 versus a square; 3 versus 12 |
| Compression | Tournament classes -> abstract sparse or faithful geometric graphs | Exact reconstruction only if enough output states exist | Quadratic input entropy versus m log m capacity |

These distinctions are the cohesive theory supported here: specify the group,
its signed representation, the retained operators, and the output contract.
Counting, geometric realization, and lossless reconstruction impose different
conditions on that same data. The new examples show exactly where those
conditions agree or separate.

## Stopping reasons and next questions

- Equal signed spectra need not classify objects or prove an equivariant
  coordinate map. Search for a smallest isospectral/non-monomial pair only
  after fixing the acting group and coordinate set.
- Free four-element orbits explain a lower divisibility bound, but the
  sharper Session89 valuations need additional structure.
- LRC component/reflection parity is saturated by L1-L2. Further work needs
  a different object, with the empty/isolated/interval controls retained.
- Geometric encoding lower bounds do not constrain exact real-coordinate
  side channels or a representation for just one observer.
- The 81/27 hostile and the singleton carrier's seven-vertex collision remain
  distinct tests. No marginal spectrum replaces joint named-port data.
