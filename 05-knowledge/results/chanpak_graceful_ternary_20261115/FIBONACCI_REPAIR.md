# A Fibonacci repair contract on a path of conflicting switches

The [graceful matching lemma](GRACEFUL.md) says that raw consecutive
cut weights already permit an interval of repairs when neighboring
cuts conflict. Those raw weights yield only quadratic width, which is
not enough on a logarithmic-order strong tournament. A useful target
is to preserve the same conflict graph, but change its weights.

## Exact weighted matching decoder

Use m switches arranged on a path. A legal choice is any set with
no consecutive indices (equivalently a matching on the corresponding
m-edge path). Give switch i weight F_(i+1), where F_1=F_2=1. Then
the sums of legal choices are **bijective** with

    {0, 1, ..., F_(m+2)−1}.

Indeed the choices omitting m give [0,F_(m+1)−1] by induction.
Those using m give F_(m+1)+[0,F_m−1]. The intervals are disjoint
and consecutive. The inverse reads the top weight: choose m iff the
target is at least F_(m+1), subtract it, and skip the conflicting
neighbor. This is the finite Zeckendorf proof, not a claim of a new
integer representation theorem. Grouping neighboring switches in
pairs gives the three symbols 00,10,01, but block boundaries carry
one compatibility bit. It is not the unrestricted ternary tree.

These are the only positive nondecreasing weights giving a bijection
onto a full initial interval for this fixed conflict path. Once all
earlier weights are forced, smaller next weight would collide with
an earlier subset (the new singleton is legal); larger next weight
would skip the next integer permanently because all later weights
are at least as large. The induction starts at weight one. This
pinpoints what must be designed and what cannot be traded for
any arbitrary list of graceful differences.

## Physical criterion for an exact repair

For a target N preceded by a literal congruent baseline N₀≤N,
let M be the placement modulus. It suffices to build m physical
switches for which every legal set I has count

    N₀ + M Σ_(i in I) F_(i+1),

with N₀, the marked continuation information, and both +2 options
retained. If (N−N₀)/M≤F_(m+2)−1, the constructive inverse above
selects an exact witness. Negative residuals can use the opposite
orientation of an already-proved switch array, where available;
the theorem does not assume a negative array exists.

The requirement is only additivity on legal matchings, not on all
subsets. But it is an integer identity. A modularly vanishing mixed
difference still moves the residual quotient and cannot be discarded.
The one-apex obstruction in [GRACEFUL.md](GRACEFUL.md) says this
contract cannot be supplied by mere separated mark intervals in that
reader; an interleaved or different physical reader is needed.

## A vertex-budget obstruction for any such strong repair array

More generally, even arbitrary weights on m switches with this
conflict path give at most F_(m+2) outcomes. Suppose a single physical
array is required to handle an entire quotient interval of length
Ω(N^β), with its witnesses in an absolute band [cN,CN] for fixed
positive c,C and all of them strong. If each switch consumes at
least s private vertices, then N^β≲F_(m+2)≲φ^m whereas
s m≤n≤1+3 log_5(CN) by
[Busch's lower bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3).
Necessarily

    β ≤ (3 log φ)/(s log 5).

For s=3 this is about 0.299. In particular three private vertices
per switch are too expensive even for width Ω(N^(1/3)), and also
for the full O(N^(3/8−1/160)) residual envelope in the prior
[polynomial-precision result](../chanpak_bott_cayley_periods_20261115/POLYNOMIAL_PLACEMENT.md).
This does not assert that the worst-case residual actually takes
that size; it rules out a specified uniform array as the way to
cover such an interval. Shared vertices, changing the reader with
the target, or a smaller residual are not excluded.
