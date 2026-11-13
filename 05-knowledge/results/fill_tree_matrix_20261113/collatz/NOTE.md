# Guarded Collatz trees, a literal relaxation, and unavoidable spatial excursions

**PROVED elementary reductions and obstruction / FINITE-EXACT controls.**
Ordinary Collatz remains OPEN. Universal reachability for the user's literal
relaxation is not proved here. A genuine spanning-tree interpretation is
available, but making that tree span all positive integers is the substantive
surjectivity obligation, not a consequence of its branching shape.

## Inheritance and live concepts

The closest proved mechanism is the guarded inverse closure in
[the old inverse-tree note, Theorem 1](../../collatz_procgen_20260924_inverse_tree_mod192.md).
The canonical hostile example is the root cycle: unfolding backward paths
without removing it repeatedly creates the same integer vertices. The corrected
near miss is independent subdivision of odd edges: the paths13→5 and53→5
share the native tail40→20→10→5. Their odd LCA is5, but their first native
meeting is40. A spatial embedding must identify that shared tail.
The least-used sidecar is the peak of a macro edge, retained below.

The board is: guarded inverse branches / shortest-path rank / odd macro edges
with native peak / dyadic classes / basin surjectivity. The new forced family
tests every lane: extra doubling does not change an odd macro edge, the inverse
guard can remain unavailable for arbitrarily long runs, an arborescence can
require unbounded spatial excursions, and finite box holes do not establish
global non-reachability.

Other recovered current work is
[THM-4514, descent-tree and spine blocks](../../../../01-canon/theorems/THM-4514-collatz-value-time-poset-spine-blocks-and-two-place-descent-tree.md)
and
[THM-4523, backward-tree rigidity](../../../../01-canon/theorems/THM-4523-collatz-functional-graph-is-rigid-backward-trees-identify-vertices.md).
The first does not supply a total rank; the second explicitly does not prove
connectivity. The mod192 note distinguishes a model-tree density from the
density of one integer basin. These boundaries are retained.

## 1. What the binary inverse tree says exactly

Let C(n)=n/2 for even n and C(n)=3n+1 for odd n, on positive integers.
The exact preimages of n are

    2n, and (n-1)/3 precisely when n=4 mod6.

The parity guard is essential: the second predecessor must be odd.
Each native inverse node therefore has one or two children before the root
cycle is removed. Let B be the set of positive integers whose forward C-orbit
reaches1. Delete the edge1→4 from C restricted to B. The result is an
in-arborescence toward1: every nonroot has one parent and a finite path to1,
with no cycle. Reversing these edges gives the actual rooted inverse tree.
It spans all positive integers if and only if Collatz holds.

This is the correct replacement for a Hamiltonian-path picture. Different
integer orbits merge; they need not form one path visiting every integer once.
The inverse tree is a tree on the known basin B, whose equality with the full
vertex set is still the question. The standard guarded closure formulation is
also described in [Lagarias's overview, section2](https://arxiv.org/html/2111.02635v1).

If every original edge is made traversable both ways, root reachability is
equivalent to the original Collatz conjecture, not a weaker substitute.
Indeed define x~y when their forward orbits intersect. For a functional graph
this is an equivalence relation: common forward tails make transitivity hold.
The endpoints of every original edge lie in one equivalence class. Therefore
an undirected path to1 already places its starting point in the ordinary
basin of the1→4→2→1 cycle. The converse is immediate.

On odd vertices the accelerated map is A(u)=oddpart(3u+1). Its inverse fibre
is instead infinitely branching:

    v=(2^a u-1)/3 >0, with a>=1 and v odd.

For u=1 mod3 the valid a are even, for u=2 mod3 they are odd, and there
are no children for3|u. Thus “binary” describes the native step clock, not the
odd accelerated tree. The old
[reverse-tree pieces, sections1 and4](../../collatz_mod6_20260917_reverse_tree_pieces.md)
also retain the residue phase of these fibres. An odd macro-tree drawing cannot
be subdivided independently, because its native even tails can overlap.

## 2. The literal relaxation and its exact reduction

Read the user's rules literally, retaining all ordinary edges and adding

    n→2n when n is even;
    n→(n-1)/3 when n is odd, n>1 and n=1 mod3.

Call this graph L. In the second rule n=6r+1 maps to2r. This is not a reverse
ordinary Collatz edge: its even endpoint would ordinarily halve, not map back
by3x+1. For instance7→2 is new. The genuine odd-predecessor rule in section1
has an even input4 mod6; these are different guards.

**Doubling-elimination lemma.** Added even doubling does not change which
vertices reach1, even when paths are required to remain inside a fixed
interval[1,X].

Proof. Between successive odd vertices every available move is a doubling or
halving on one dyadic chain. A finite walk that eventually reaches the chain's
odd endpoint can be replaced by successive halvings. This replacement has the
same endpoints and never increases the maximum visited integer. Apply it to
the initial even segment and to each segment following an odd arithmetic
move. Conversely all halving paths were already allowed. No inverse3 move is
available on the interior even vertices under these literal rules. QED.

Consequently all odd-to-odd reachability is captured exactly by

    A(u)=oddpart(3u+1),                         native edge peak 3u+1;
    B(u)=oddpart((u-1)/3), u>1 and u=1 mod6,   native edge peak u.

For a general starting integer, first halve to its odd part, retaining the
starting value as a possible peak. This is a reachability reduction with a
peak sidecar; it is not an assertion that all vertices of a dyadic fibre form
one strongly connected component.

Ordinary Collatz implies that every vertex of L reaches1, since L contains
the original edges. A converse or a proof of universal L-reachability is not
established. It is therefore a relaxation, but strict logical weakness is not
claimed. The earlier E-graph relaxation in
[reverse-tree pieces, Proposition4.2](../../collatz_mod6_20260917_reverse_tree_pieces.md)
adds even forward3n+1 edges; it is yet another graph. Its Q1/Q2 results must
not be assigned to L by exchanging arrows without retaining the guards.

## 3. An all-order obstruction to bounded spatial fill

**Forced-excursion theorem.** For every even K>=2 put n_K=2^K-1. Every
L-path from n_K to1 must visit an integer at least

    P_K=2*3^K-2.

In particular no constant C can guarantee a root-reaching L-path inside
[1,Cn] for every starting integer n. This also holds for ordinary Collatz.

Proof. For j=0,...,K-1 define

    x_j=3^j*2^(K-j)-1.

All x_j are odd and greater than1. Since K is even, x_0=0 mod3; for j>=1,
x_j=2 mod3. Therefore none admits the added inverse3 edge. For j<K-1,

    3x_j+1=2x_(j+1),

and x_(j+1) is odd. The only available odd transition is consequently
x_j→x_(j+1). Even doubling can only delay arrival at that same odd endpoint,
by the doubling-elimination lemma. Thus every successful path must traverse
all x_j in order. At the final x_(K-1)=2*3^(K-1)-1 the next native move is
forced and reaches3x_(K-1)+1=P_K. Finally
P_K/n_K is asymptotic to2*(3/2)^K and is unbounded. QED.

The underlying Mersenne start identity is classical; see
[Lagarias's overview, section3.2](https://arxiv.org/html/2111.02635v1).
The added conclusion here is that the user's additional choices cannot avoid
this excursion because their guards fail throughout the forced odd run.
This is a lower bound, not a claim that every n_K reaches1. The finite
certificates below establish reachability and sharpness for20 named members.

The connection contract is now concrete. Source: the literal directed graph
L. Target: paths or a spanning in-tree laid out at integer positions. Map:
retain each vertex's integer coordinate and each macro edge's largest native
coordinate. Preserved predicate: root reachability. Lost by an unweighted
macro tree: excursion height and shared even tails. Needed sidecar: peak and
native identifications. Decisive test: the forced family above. It rules out
uniform constant-factor buffers, while leaving a construction with an
unbounded buffer or a different spatial coordinate possible.

## 4. The precise spanning-tree proof obligation

For any countable directed graph G with distinguished root1, the following
are equivalent:

1. Every vertex has a finite directed path to1.
2. There is a function d:V→N, with d(1)=0 and d(v)>0 for v≠1, such that
   each v≠1 has an outgoing edge v→w with d(w)<d(v).
3. G contains a spanning in-arborescence toward1.

For1⇒2 choose shortest-path distance. For2⇒3 choose one decreasing edge
at each nonroot; strict descent in N rules out both cycles and infinite
forward paths, so every chosen path ends at1. For3⇒1 use the tree paths.
On integer-labelled graphs choose the least eligible neighbour to make the
selection explicit. This construction applies unconditionally to the known
basin, but supplying d for every integer is exactly what is still missing.

Neither local branching counts, a tree on every currently reached vertex,
nor a density tending to1 replaces that requirement. A density-one subset
can omit an infinite sparse set. For original functional Collatz there is
not even a choice of forward parent to optimize; the arithmetic orbit is
already fixed. L adds a real choice at B-vertices, which could support a
different rank, but section3 proves it cannot be bounded by a fixed spatial
buffer around the starting value.

## 5. Where5,17 and7,21 really meet

For ordinary positive3n+1, the exact odd paths include

    7→11→17→13→5→1,       21→1.

For positive3n-1, the odd accelerated map is A_-(u)=oddpart(3u-1). It has
the explicit cycles and feeder

    5→7→5;
    17→25→37→55→41→61→91→17;
    21→31→23→17.

Thus5 and17 are the minima of these two explicit nontrivial minus-map
cycles;7 lies on the first, and21 feeds the second. The native first cycle
is5→14→7→20→10→5. These are exact dynamical relations, not merely matching
small integers. Negation transports positive3n-1 to negative3n+1, as used
in the recovered
[clone-clock COLLATZ note](../../clone_clocks_negative_collatz_20261113/COLLATZ.md).
No statement that these exhaust all positive minus-map behaviour is needed.

The tournament facts attach different predicates to the same numbers:
5 and17 occur in the square-plus-one intersection of the strong unit-ready
minimum family, whereas7 and21 are excluded ordinary Hamilton counts.
There is no established map taking the minus-map cycles to tournament
non-realizability. The cycle membership supplies a precisely typed comparison
and possible arithmetic guide; it does not explain those exclusions.

## 6. Exact controls and reproducibility

`check.py` compares three independent finite descriptions: native reverse
reachability, the graph with added even doubling removed, and guarded odd
inverse fibres with macro-edge peaks retained. All agree at the following
box bounds; ordinary Collatz is a separate positive control.

| Native cap X | Ordinary nodes reaching1 within box | L nodes reaching1 within box | First L hole |
|---:|---:|---:|---:|
|100|49|64|15|
|1000|340|643|63|
|10000|4235|6407|255|
|100000|39706|64000|1023|

Every hole here is a failure to stay in that box, not a global counterexample.
The deterministic policy “halve evens; use B when legal; otherwise3n+1”
is directly certified for all starts1..100000, with a10000-step cap per
unresolved chain. It certifies161649 distinct vertices and visits no value
above86093440. This is a finite control, not a theorem of convergence.

For K=2,4,...,40 the saved literal paths attain exactly the proved lower
bound P_K. Hence their minimum possible native peaks are established exactly,
including starts as large as2^40-1. The guard/forced-path identities are also
checked for every even K through200. All named plus/minus paths are checked.

Run from the repository root:

```text
python -B 05-knowledge/results/fill_tree_matrix_20261113/collatz/check.py
python -O -B 05-knowledge/results/fill_tree_matrix_20261113/collatz/check.py
```

The two JSON receipts agree. Their source pins use LF-normalized bytes for
the five recovered current documents. No theorem ID was reserved, no global
Collatz claim was promoted, and no density observation is used as coverage.
