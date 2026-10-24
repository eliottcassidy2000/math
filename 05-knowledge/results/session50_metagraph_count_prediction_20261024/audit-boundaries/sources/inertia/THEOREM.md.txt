# Native observer holonomy: finite saturation certificates and a window obstruction

Status: **PROVED** certificate, counting sandwich, inertia/cycle-index identities,
adaptive rescue, and all-size window obstructions. **FINITE-EXACT** complete
native relocation census through order seven. **OPEN** all-size connectivity
of the single-vertex-relocation Hamiltonian-path graph. No new Burnside formula,
compressed cross-order recurrence, complexity advantage, or priority claim.

## 1. Inheritance and the precise new target

The fixed-path fibre identity `h(T)=H(T)/a(T)`, where `a=|Aut|`, is inherited
from THM-781, *Hamiltonian-path inverse metagraph fibre*, and Session48's
observer-moments theorem. THM-796, *Three-sorted recursive tiling-line-node
incidence*, Section 12 already describes the marked tournament groupoid.
THM-1430, *The tiling/class/metagraph dictionary*, already gives exact class
BFS and full Hamiltonian-path enumeration engines. None is renamed here.

The October20 metagraph-growth packet supplies the hostile unmarked order-three
metagraph and the inflation `C3[TT_k]`, including `H(C3[TT_2])=45`, `a=3`,
`h=15`. THM-4731, *Tiling line metagraph self-lines by twisted Burnside,
sigma twins, depth, and the tribonacci partner*, already identifies the
backward-nonpath partner as rigid for `n>=4` and studies its path counts.
Those source objects are inherited, not discoveries of this packet.
All quoted source documents are byte-pinned inert `.txt` archives in `sources/`.

The additional contribution is an owner-permutation certificate saying when a
collection of native local reobservations has recovered BOTH a whole word
fibre and its source isotropy. It gives a finite exact test without presuming
all-size local connectivity. Two explicit families show why bounded contiguous
window moves cannot substitute for that test.

Concept board: normalized word fibre; free marked observers; source isotropy;
owner-labelled closed walks; exact path-count saturation; window size versus
number of moved vertices. The canonical hostile is a frozen observer with
other global observers, not a failure of the underlying tournament to have
Hamiltonian paths.

## 2. Native records and reobservation arrows

Fix `n>=1`, `P0=(0,1,...,n-1)`, and

    X_n = {tournaments on [n] containing P0 as a directed path},
    |X_n| = 2^m,     m=(n-1)(n-2)/2.

The nonpath pairs are listed lexicographically; bit one directs the smaller
owner to the larger owner. Path arcs are always forward. Write `T_x` for the
entire source matrix of word `x`. An actual Hamiltonian order `p` in `T_x`
produces `y=normalize(T_x,p)`. The arrow

    x --p--> y

retains `p` as a bijection FROM the owner indices of `y` TO those of `x`.
It is an isomorphism of the entire source tournament, not only of a score,
path fragment, or unaligned source/partner pair. Its inverse is also a lawful
reobservation. Normalizing means literally reading every source arc in order
`p`; it does not use an isomorphism canonicalizer.

Our tested move family relocates one vertex of the current observer to another
position, retaining the relative order of all other vertices, and accepts the
move exactly when the resulting full order is Hamiltonian. Its inverse is a
single-vertex relocation. The number of MOVED VERTICES and the SPAN of this
move differ: moving the first vertex to the end spans the whole observer.

This graph is NOT the original blue/black complement-line metagraph, nor the
single-arc-flip graph. No source arc is changed. It is an enriched, fixed-source
observer correspondence on the same word carrier. A claimed extraction from
only the unmarked blue/black graph does not acquire these labels for free.

## 3. The holonomy saturation certificate

Take any connected finite reverse-closed graph `B` of authenticated arrows,
with distinct word vertices, all reached from a root `x0`. Choose a spanning
tree. Its composed edge labels give an isomorphism

    r_x : T_x -> T_x0

for each vertex, with `r_x` also an actual Hamiltonian order in the root source.
For each arrow `x --p--> y`, form the closed owner permutation

    g_e = r_x o p o inverse(r_y).

Every `g_e` is an actual element of `Aut(T_x0)`. Let `G_B` be the subgroup they
generate, and put `b=|B|`, `g=|G_B|`, `H=H(T_x0)`.

**Theorem 1 (finite saturation).**

    b g <= H.                                                (1)

Moreover,

    b g = H
      iff B is the ENTIRE normalized word fibre of [T_x0]
          and G_B = Aut(T_x0).                              (2)

The criterion is valid for every size and every authenticated arrow family.
It does not assume the corresponding local path graph is connected.

**Proof.** Source automorphisms act freely on actual ordered Hamiltonian paths:
fixing an entire order fixes every owner. If `h=H/a` is the full word-fibre
size, then `b<=h` and `g<=a`. This proves (1). Equality in the product of these
two positive inequalities occurs exactly when both are equalities, proving
(2). More directly, the `b g` orders `gamma o r_x`, for `gamma in G_B`, are
all distinct. Equality of two such orders implies equality of their normalized
words `x`, then equality of the owner permutations. All are paths of the root
source. This proves the bound without treating any marked observer as having
nontrivial stabilizer. QED.

There is also an exact lift interpretation. A walk in B lifts from P0 to an
actual root-source observer. If it ends at x, append the inverse chosen tree
path: the resulting loop lies in G_B, so its endpoint order is gamma o r_x.
Conversely, every product of fundamental loops followed by the tree path
realizes gamma o r_x. Thus the connected component of P0 in the actual-path
lift has exactly b g paths. Saturation also certifies that this lifted
component is the whole actual observer set. This argument uses the actual
edge permutations and their inverses, not unlabelled graph connectivity.

Changing the spanning tree with the same root and root gauge changes the
generators but not the holonomy subgroup. Changing the root or its owner gauge
conjugates the subgroup; its order and the completeness result are unchanged. It is enough to retain any
verified subgroup of source automorphisms: saturation still proves completeness.
Using all fundamental loops is a constructive way to obtain that subgroup
from the actual reobservation records, without an automorphism oracle.

The certificate consists of word vertices, a source matrix, authenticated tree
arrows, loop owner permutations, a verified finite group closure/order, and an
independent exact value of `H` (for example the standard subset recurrence).
A verifier checks permutations against all source arcs and checks `b g=H`.
No candidate class label or known A000568 value enters that verification.

### A count sandwich before saturation

If an authenticated reverse-closed move graph covers every vertex of `X_n`,
and its components are `B`, then

    sum_B b_B g_B/H_B <= A000568(n) <= number of components B. (3)

For the lower bound, every word `x` carries a certified subgroup of its actual
source automorphism group, of order `g_B`. Sum `g_B/H_B <= a(T_x)/H(T_x)` over
all words. A true class has `h=H/a` words, so the right sum contributes one
per class. For the upper bound, no component crosses source isomorphism
classes, and every class contains a path and therefore at least one word.

If every component saturates, both bounds equal the number of components.
Saturation is SUFFICIENT, not necessary for equality in the lower bound:
a rigid fibre can be disconnected while every component already has the full
trivial isotropy. The exact lower-bound deficit is

    sum_x (a(T_x)-g_B(x))/H(T_x).                            (4)

It vanishes precisely when every word's certified subgroup is its whole
source automorphism group. It is not governed merely by connectivity.

### Exact rescue if a local graph fails

Suppose `b g < H`. There is a Hamiltonian order `q` of the root source outside
`{gamma o r_x}`. Enumerate actual paths until finding one; membership is tested
by its normalized word `y` and, if `y in B`, by `q o inverse(r_y) in G_B`.

- If `y` is absent, add its full reobservation arrow and the word; `b` grows.
- If `y` is present, add that arrow's new owner holonomy; `g` grows.

Close under any desired local moves and repeat. Each rescue strictly increases
the integer `b g<=H`, so termination at saturation is proved. For a global
census, merge any existing components met by the added arrows. This is an
exact fallback, not a bounded-window or efficient all-size theorem. In the
worst case it approaches the inherited full path enumeration. The checker
records explicit terminating rescue histories for the frozen sources at
orders three through six.

## 4. What the closed loops count, and what they do not

For a fixed word `x`, full closed reobservation arrows are precisely the
permutations `p` with `normalize(T_x,p)=x`; these are `Aut(T_x)`. A nonidentity
one changes the ACTUAL observer to `p(P0)` while returning its normalized word.
It does not stabilize the marked record `(T_x,P0)`.

Therefore three distinct quantities must stay distinct:

    number of actual observers                 H(T),
    number of normalized observer words         h(T)=H(T)/a(T),
    marked-record stabilizer order              1.

The source isotropy `a(T)` is retained by owner-labelled closed reobservations,
not by the stabilizer of the marked observer. The directed triangle already
has `(H,h,a)=(3,1,3)` while every marked-record stabilizer is trivial.

The exact loop-weighted identity is

    A000568(n) = sum_(x in X_n) a(T_x)/H(T_x).               (5)

Its identity-loop part is

    sum_(x in X_n) 1/H(T_x) = 2^binom(n,2)/n!,             (6)

and its nonidentity-loop part is the missing symmetry correction

    Delta_n = sum_x (a(T_x)-1)/H(T_x).                     (7)

Equations (5)-(7) are an observer-coordinate version of standard groupoid
mass/Burnside, not a new enumerative formula. Raw loop counts do not give
A000568: `sum_x a(T_x)=sum_[T] H(T)`, which is `35,180`, not `456`, at `n=7`.
In groupoid language, the source inertia groupoid has cardinality A000568;
its number of isomorphism classes instead sums conjugacy-class counts and must
not be substituted. The marked-observer groupoid is free and has cardinality
`2^m`; taking its inertia does not restore the forgotten source isotropy.

### Full cycle-type correction

For indeterminates `p_1,p_2,...`, let `p_lambda(g)` be the product of one `p_l`
for each length-`l` cycle of `g`. Then

    Z_n(p) = sum_x (1/H(T_x)) sum_(g in Aut(T_x)) p_lambda(g)
           = sum_[T] (1/a(T)) sum_(g in Aut(T)) p_lambda(g). (8)

The proof is the same exact `h=H/a` change of measure. At all `p_l=1`, this is
(5). At all `p_l=q`, it counts tournaments whose vertices have `q` available
colours, up to colour-preserving isomorphism. A full labelled observer mark
kills the nonidentity terms, so it cannot be substituted into (8).

For comparison only, the inherited classical Burnside coefficient for an
odd partition `lambda=(l_1,...,l_s)` is

    [p_lambda] Z_n = 2^e(lambda)/z_lambda,
    e = sum_i (l_i-1)/2 + sum_(i<j) gcd(l_i,l_j),
    z_lambda = product_l l^(m_l) m_l!;                    (9)

partitions with an even part contribute zero. An even permutation cycle has
a power exchanging two of its owners, which cannot preserve their directed
arc. For odd cycles, internal unordered-pair orbits number `(l-1)/2`, and
between two cycles they number `gcd(l_i,l_j)`; each orbit freely chooses one
orientation. Grouping the usual permutation fixed-point sum by cycle type
gives (9). This recovered degree-local engine is only an independent checksum.
The native certificate generation never calls it. Matching (9) coefficient
by coefficient checks more than matching its specialization A000568.

## 5. All-size bounded-window obstructions

A width-`w` move changes the observer only inside ONE contiguous interval of
at most `w` positions, retaining its outside owners and their order. It may
arbitrarily reroute the inside, provided the full resulting order is directed.
This is broader than the tested width-restricted single-vertex relocations.

**Theorem 2 (fibre-completeness obstruction).** Let `T*_n` have the forward
path `0->1->...->n-1` and every nonpath arc directed backward. For `n>=3`,
its displayed observer admits no nonidentity reobservation supported in any
proper contiguous interval. Nevertheless it has at least `n` actual
Hamiltonian observers. For `n>=4` it is rigid and thus has at least `n`
normalized observer words.

**Proof.** For a rerouted interval `[a,b]`, if `a>0`, the preceding owner `a-1`
can point to only `a` among the interval's owners. Thus the first owner is
forced to be `a`. After the ascending prefix `a,...,i` is used, owner `i`
points to only `i+1` among the remaining owners. The whole interval is forced.
If `a=0`, properness gives `b<n-1`. The following owner `b+1` receives an arc
from only `b` in the interval, forcing the last owner; backward induction
forces the rest. Meanwhile `n-1->0`, so the `n` cyclic rotations of the
displayed path are distinct valid paths. For `n>=4`, the scores are
`d_0=d_1=1`, `d_i=i` for `2<=i<=n-3`, and
`d_(n-2)=d_(n-1)=n-2`; each score class has one or two owners. A score
class of size two cannot be swapped by an automorphism, because its internal
arc would reverse. Therefore every owner is fixed. QED.

This source is the inherited tribonacci-partner object. The theorem adds its
explicit window-freezing mechanism, not a new count for its Hamiltonian paths.
It shows why component counting and saturation from any fixed window size
fail beyond that size. It does not contradict the unrestricted-span relocation
census: first-to-last and last-to-first rotations span the whole observer.

**Theorem 3 (isotropy obstruction, not just disconnection).** For every fixed
window bound `w>=1`, choose `k>=max(1,w-1)` and let

    T = C3[TT_k]

in its concatenated three-chain observer. Every window of length at most `w`
induces a transitive tournament, hence has a unique Hamiltonian order. The
observer word is isolated under ALL width-`w` reroutings, with trivial local
holonomy. But cyclically rotating the three equal chains is a genuine source
automorphism of order three. Consequently local closed-walk isotropy does
not recover all source isotropy, and the lower bound in (3) is strict for
this order's full width-`w` move graph.

**Proof.** An interval meeting all three consecutive chains must include the
entire middle chain and at least one owner from each outer chain, requiring
at least `k+2>w` positions. Thus an allowed window lies in one chain or in
two consecutive chains. In either case the induced source is transitive,
with the displayed increasing order its only directed Hamiltonian order.
No nonidentity move leaves the root word. Rotation of the three chains,
preserving their internal chain order, preserves every source arc. Its order
is three. The frozen word alone has a positive missing weight at least `2/H`
in (4). QED.

The exact block-inflation construction is inherited from the growth packet;
the new point is the all-size failure of bounded contiguous-window isotropy
recovery. At `k=2`, the explicit source rows are `(14,12,56,48,35,3)`,
`H=45`, `a=3`, `h=15`. Every width-three rerouting is frozen at the displayed
observer, losing `2/45` from its loop-weighted contribution. The checker also
covers `k=1,3,4` (orders `3,9,12`), with independent block-run and subset-DP
path counts `3,2721,421425` respectively.

These theorems concern CONTIGUOUS SPAN, not a bound on the number of moved
vertices, source arcs, or arbitrary noncontiguous support. They do not refute
all enriched metagraph formulas, cross-size constructions, or nonlocal moves.

## 6. A small exact global law for the tested move graph

**Proposition 4 (relocation-arrow total).** For `n>=3`, the number of directed
accepted single-vertex relocation arrows over all words of `X_n`, including
loops and retaining the relocation permutation, is

    E_n = 2^m (n^2+n-4)/8.                                (10)

Adjacent swaps can never be valid because they reverse a mandatory path arc.
For an interval `[i,j]` of length at least three there are two relocations,
left endpoint to right end and conversely. Each needs `1+epsilon_L+epsilon_R`
distinct free bits, where the epsilons indicate an outside boundary on that
side. Thus both orientations together contribute a fraction
`2^(-epsilon_L-epsilon_R)` of all words. There is one full interval, `2(n-3)`
one-boundary intervals, and `(n-3)(n-4)/2` two-boundary intervals. Their sum is
`1+(n-3)+(n-3)(n-4)/8=(n^2+n-4)/8`. QED.

This global edge total is a consistency check, not a class-count law. It does
not certify fibre connectivity, holonomy saturation, or cross-order closure.

## 7. Native finite-exact census and independent controls

`check.py` generates each whole fixed-path cube, accepts actual relocations,
constructs components and spanning-tree transports, and closes the groups
of actual fundamental-loop owner permutations. Only then does it compare
against independent path counts, automorphisms, and inherited class counts.

| n | Words | Components/classes | Directed arrows | Saturated |
|---:|---:|---:|---:|---:|
| 1 | 1 | 1 | 0 | 1/1 |
| 2 | 1 | 1 | 0 | 1/1 |
| 3 | 2 | 2 | 2 | 2/2 |
| 4 | 8 | 4 | 16 | 4/4 |
| 5 | 64 | 12 | 208 | 12/12 |
| 6 | 1,024 | 56 | 4,864 | 56/56 |
| 7 | 32,768 | 456 | 212,992 | 456/456 |

The component test proves the resulting bounded census exactly. It is not
merely numerical agreement with A000568. At order seven the source-automorphism
histogram is `1:399, 3:47, 5:4, 7:1, 9:4, 21:1`.

The identity-loop masses at `n=3..7` are
`4/3, 8/3, 128/15, 2048/45, 131072/315`; the symmetry corrections are
`2/3, 4/3, 52/15, 472/45, 12568/315`. Raw closed-loop counts are
`4,12,92,1224,35180`; these do not equal the class counts.

Restricting relocations to span three gives the rigorous partial-graph
sandwiches `2<=2<=2`, `4<=4<=6`, `176/15<=12<=35`, and
`80758/1485<=56<=420` at orders three through six. The order-four lower bound
is exact despite three unsaturated components, illustrating the distinction
between isotropy completeness and connected-fibre completeness.

Independent controls include:

1. Literal Hamiltonian permutations and direct automorphism checks on every
   discovered representative through `n=7`, separate from subgroup generation.
2. The complete labelled tournament cube through `n=5`, direct full-permutation
   canonicalization, marked mass, and coefficientwise fixed-point cycle indices.
3. Coefficientwise comparison with inherited odd-partition Burnside through
   `n=7`, not used in generation; (10) checked for every applicable order.
4. All nonidentity proper-window permutations of the explicit `T*_n` observers
   through `n=8`, global rotations, and direct source automorphisms.
5. Every window permutation for the four stated cyclic-chain inflations,
   with independent block-word/run-multinomial Hamiltonian counts.
6. Triangle hostile separating a trivial marked stabilizer from order-three
   source isotropy, partial-graph sandwiches, and explicit global rescue paths.

Every check is exception-based and remains active under `python -O`. Normal
and optimized reports are byte-identical. The manifest-only verifier copies
only pinned packet files to an empty temporary directory and reruns both.
There are no imports from another checker, network calls, Git dependencies,
third-party packages, or hidden parent-workspace reads.

## 8. Remaining target and boundary

The specific next mathematical question is whether unrestricted-span single-
vertex relocations always connect all Hamiltonian paths of a tournament.
The finite census supports it through seven but does not prove it. Generic
undirected backbite/reconfiguration theorems are inapplicable and are not
imported. Even a proof would give an exact fixed-size extraction engine, not
a compressed recurrence predicting the next A000568 term or a demonstrated
improvement over the inherited odd-partition engine.

What is closed here is the finite-certificate alternative: explicit owner
holonomy plus an exact `H` detects missing fibre/isotropy information, and a
strictly increasing authenticated global-reobservation rescue finishes if
local moves do not. The all-size window families show why the scope matters.
