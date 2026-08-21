# Tournament subset decks, higher-order faces, and odd-middle Kneser duality

**Session date:** 2026-08-21

**Status:** combinatorial identities `PROVED`; catalogue and local-rule searches
`FINITE-EXACT`; reconstruction and extremal background `CITED`; the
`c_5 >= 9/64` face-potential lift `DERIVED + LITERATURE-AUDIT-NEGATIVE`;
oriented-matroid characterization and higher lifts `OPEN`.

Exact replay:
[`tournament_subset_deck_hierarchy_scout_codex_20260821.py`](../04-computation/tournament_subset_deck_hierarchy_scout_codex_20260821.py),
with stored output
[`tournament_subset_deck_hierarchy_scout_codex_20260821.out`](../05-knowledge/results/tournament_subset_deck_hierarchy_scout_codex_20260821.out).

That replay is exhaustive through `n=6`, checks the odd-middle transform
through `r=6`, and contains the targeted higher-rule/GP controls.  The full
`n=7` deck census in Section 3 is a separate exact canonicalization of
McKay's cited `456`-representative catalogue; the four displayed bitstrings
make its affine hostile independently reproducible, but the catalogue sweep
is not bundled into the main replay.

## 1. Outcome first

There are two mathematically different ways to replace the edges of a
tournament by triangles or squares.

1. **Induced-card hierarchy.**  An ordinary edge tournament `T` on `V` gives
   one ordinary tournament `T[S]` on every subset `S`.  For `|V|=7` this gives

   ```text
   layer size k:       0  1  2   3   4   5  6  7
   number of cards:    1  7 21  35  35  21  7  1.
   ```

2. **Higher tournament.**  One may instead assign an independent orientation
   to every triangle, square, or general `d`-set.  These are the higher-order
   tournaments of Leader--Tan.  They have many more degrees of freedom than
   the corresponding induced cards.

The distinction is load-bearing.  Within the induced hierarchy, the main
exact results are these.

- Every higher labeled layer contains the lower layers by a binomial
  incidence transform.  Unlabeled histograms retain its marginal identities
  but lose the overlap maps.
- When `n=2r+1`, the two middle layers have the same size and their additive
  incidence transform is exactly invertible.  After complementing the
  `(r+1)`-sets, the matrix is the odd/Kneser graph `KG(2r+1,r)`, with

  ```text
  (A^-1)_(S,R)=(-1)^|S intersect R| /
               ((r+1) binom(r,|S intersect R|)).        (1.1)
  ```

  Thus the user's `35` triangle and `35` square layers at `n=7` are a
  lossless dual pair for **addressed additive face fields**.
- The four-vertex tournament type splits exactly into an additive triangle
  boundary, a square interaction, and a dual-odd chirality.  The square
  interaction is the first pair-correlation of cyclic triangles.
- Layerwise powers and Boolean-lattice inversion recover the complete overlap
  spectrum of local motifs.  At `n=7`, quadratic triangle data on layers
  `3,4,5,6` separates identical, edge-sharing, vertex-sharing, and disjoint
  pairs of cyclic triangles.
- Multiplicity-aware unlabelled `6`-decks reconstruct all `456` tournaments
  on seven vertices, but this is discrete and nonlinear.  There is no affine
  functional of that deck recovering even the Hamilton-path count.

For genuinely higher tournaments, a second group of exact results emerged.

- Edge signs produce triangle orientations by their switching-invariant
  triangle flux.  The image consists exactly of the flat oriented two-graphs:
  it has dimension `binom(n-1,2)`, codimension `binom(n-1,3)`, and each fibre
  contains `2^(n-1)` edge tournaments.
- An exhaustive search over all isomorphism-equivariant rules that orient a
  `5`-set from its induced ordinary tournament gives

  ```text
  max Pr(random T_6 produces a directed 6-set)
      =4608/32768=9/64.                                (1.2)
  ```

  Hence

  ```text
  c_5 >= 9/64,                                         (1.3)
  ```

  improving Leader--Tan's general `d=5` construction `1/16`; their upper
  bound remains `1/6`.  The successful `4608` edge tournaments are exactly
  the switching/relabeling orbit of Gunderson--Semeraro's special six-vertex
  oriented two-graph.  Their numerical `9/64` Turan bound is cited; the
  explicit face-potential lift `(1.2)--(1.3)` was not found in the audited
  literature and is not asserted as a formal priority claim.

## 2. Inheritance pass and concept board

The closest proved canon mechanisms are:

- `THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration`:
  a coarse subset action erases the three perfect-matching channels on a
  four-set;
- `THM-3380-hamiltonian-deletion-layer-monoid-semiring-and-small-order-boundaries`:
  deletion cards form a useful hierarchy only after their operation and
  multiplicity types are fixed;
- `THM-1965-tournament-invariant-refinement-lattice-through-order-six` and
  `THM-1966-signed-redei-information-is-independent-at-order-seven`:
  scores, cycle profiles, spectra, and signed path data already separate into
  incomparable information channels.

The canonical hostile is Stockmeyer's family of nonisomorphic tournaments
with the same point-deleted deck.  Kocay later supplied a corrected
construction after finding a flaw in the original proof.  By Kelly's
identity, such a pair has the same unlabelled profile at every proper order:
adding more marginal layers does not restore the missing gluing.

The corrected near miss is the claim that equal binomial layer sizes imply
dual tournament structures.  Complementation pairs their **indices**, not
their induced tournament values.  The least-used required sidecar is the
address/overlap map telling which cards share which vertices.

The live concept board was:

| Object | Representation | Preserved coordinate | Missing sidecar |
|---|---|---|---|
| ordinary tournament | induced subset cards | all local edge orientations when addressed | overlap labels after quotienting |
| cyclic-triangle field | scalar on `3`-sets | additive incidence and overlap moments | source/sink chirality |
| odd middle layers | `KG(2r+1,r)` transform | every additive face harmonic in characteristic zero | card alignment; bad characteristics |
| higher tournament | alternating sign on `d`-sets | directed-simplex compatibility | lower-face origin and gauge |
| oriented two-graph | triangle flux of edge signs | switching class | a chosen edge lift |
| determinant support | rank-two chirotope/Pluecker faces | signs of determinant circuits | magnitudes and zero minors |

## 3. Addressed cards, unlabelled profiles, and Kelly incidence

Let `T` be a tournament on an `n`-set `V`.  Its addressed labeled `k`-deck is

```text
D_k^lab(T)=(T[S])_(S in binom(V,k)),                     (3.1)
```

including the subset address and the restriction maps on overlaps.  Its
unlabelled multiplicity deck is the vector

```text
D_k(T)=(N_F(T))_(F an isomorphism type on k vertices).   (3.2)
```

Every edge occurs in exactly

```text
r_k=binom(n-2,k-2)                                      (3.3)
```

addressed `k`-cards.  Thus a proposed family of addressed cards comes from a
single tournament exactly when all occurrences of each labeled edge agree.
Over `F_2` this realization space has

```text
dimension  binom(n,2),
codimension binom(n,2)(r_k-1)                           (3.4)
```

inside the raw collection of card-edge bits.  Any one addressed layer with
`k>=2` therefore reconstructs the original edge tournament.

For an `r`-vertex tournament type `F` and `k>=r`, double counting pairs
`(copy of F, containing k-set)` gives Kelly's identity

```text
sum_(|S|=k) N_F(T[S])
    =binom(n-r,k-r) N_F(T).                             (3.5)
```

Consequently a full unlabelled `k`-profile linearly determines every lower
unlabelled profile.  It does **not** determine their joint placement.  For
fixed `k`, the number of possible profile vectors grows only polynomially in
`n`, while the number of tournament types grows exponentially in `n^2`; so
fixed-order profiles must eventually collide.  More strongly, the
Stockmeyer--Kocay point-deck counterexamples have all proper profile vectors
equal by `(3.5)`.

The first small hostile is already on four vertices.  A universal source over
a directed triangle and a directed triangle over a universal sink have the
same unlabelled `3`-deck: one cyclic and three transitive cards.  They are
nonisomorphic and dual to one another.

### Exact seven-vertex census

On McKay's complete catalogue of `456` seven-vertex tournament types, the
number of types separated by one multiplicity-aware `k`-deck is

```text
k                         0   1   2   3    4    5    6
distinct deck signatures  1   1   1  15  164  424  456
largest collision fibre 456 456 456  79   18    3    1. (3.6)
```

Thus the `6`-deck reconstructs all seven-vertex tournaments.  The `5`-deck
still has `31` collision buckets containing `63` types.  If multiplicities
are discarded and only the **set** of occurring `6`-card types is retained,
there are `453` signatures and three collision pairs, agreeing with McKay's
2022 reconstruction audit.

A sharp `5`-deck hostile is given by catalogue lines `163` and `401`:

```text
001100000101000110101
101000100101100110101.                                  (3.7)
```

They have the same `5`-deck, score sequence `(2,2,3,3,3,3,5)`, cyclic-triangle
count `11`, and six unordered disjoint cyclic-triangle pairs, but their
Hamilton-path and directed Hamiltonian-cycle counts are respectively

```text
(H,t_7)=(109,8) and (111,9).                            (3.8)
```

Thus layer `5` loses a genuinely full-support coordinate even after several
natural lower moments are supplied.

Reconstruction here is nonlinear.  In zero-based McKay catalogue indices,

```text
D_6(T_19)+D_6(T_64)=D_6(T_23)+D_6(T_57),                (3.9)
```

coordinatewise, while their Hamilton-path totals are

```text
31+113=144 != 29+117=146,                              (3.10)
```

and their directed Hamiltonian `7`-cycle totals are

```text
0+11=11 != 0+12=12.                                   (3.11)
```

The representative edge bitstrings are

```text
19: 000000001001000010101
64: 000001100001010100101
23: 000000010000100010101
57: 000001001000010000101.                             (3.12)
```

Therefore neither invariant is an affine functional of the unrooted
`6`-deck, even though the deck separates individual types.

Converse duality provides a second useful grading.  The `456 x 56` six-deck
incidence matrix has rank `56`, splitting as ranks `34` and `22` in its
converse-even and converse-odd sectors.  At `n=7`, a self-dual `6`-deck comes
exactly from a self-converse tournament.  This fails at layer `5`: `38`
non-self-converse tournaments, in `19` converse pairs, already have
self-dual `5`-decks.

## 4. The odd-middle Kneser inverse

Let `|V|=2r+1`, let `x_S` be a scalar on each `r`-set, and aggregate it on
the next layer:

```text
y_K=sum_(S subset K, |S|=r) x_S,       |K|=r+1.         (4.1)
```

Index `K` by its complementary `r`-set `R=V\K`.  Then `(4.1)` is

```text
y_R=sum_(S intersect R=empty) x_S=A x,                  (4.2)
```

where `A` is the adjacency matrix of `KG(2r+1,r)`, the odd graph
`O_(r+1)`.  Define

```text
B_(S,R)=(-1)^t / ((r+1)binom(r,t)),
t=|S intersect R|.                                     (4.3)
```

For fixed `S,U`, put `q=|S intersect U|`.  The `r+1` sets `R` disjoint from
`U` are obtained by omitting one point from the `(r+1)`-set `V\U`.
Exactly `r-q` omissions give `|S intersect R|=r-q-1`, and `q+1` give
`|S intersect R|=r-q`.  Their two contributions in `(AB)_(U,S)` cancel when
`q<r`; when `q=r`, they sum to one.  Hence `B=A^-1`, proving `(1.1)`.

The spectrum is

```text
lambda_j=(-1)^j(r+1-j),
multiplicity binom(2r+1,j)-binom(2r+1,j-1), 0<=j<=r.    (4.4)
```

Thus the transform is singular over `F_p` exactly for primes `p<=r+1`, and
invertible over `Z/mZ` exactly when `gcd(m,(r+1)!)=1`.  This is an explicit
model of how an “exceptional prime” should arise from a finite relation
matrix: it is visible in the determinant or Smith data, not merely in a
numerical coincidence.

For the user's `n=7`, `r=3` case, the inverse is

```text
x_S =  1/4  sum_(|S intersect R|=0) y_R
      -1/12 sum_(|S intersect R|=1) y_R
      +1/12 sum_(|S intersect R|=2) y_R
      -1/4 y_S.                                         (4.5)
```

The spectrum is `4^1,(-3)^6,2^14,(-1)^14`; only characteristics `2` and `3`
are bad.  Equivalently,

```text
A^-1=(-A^3+2A^2+13A-14I)/24.                           (4.6)
```

The inverse requires addressed cards.  An unlabelled multiset has forgotten
which square complement is disjoint from which triangle, so `(4.5)` is then
ill-typed.

## 5. Layers as an overlap-moment hierarchy

The Boolean-lattice transform does more than propagate first moments.  Let a
weight `w_A` live on `r`-sets and put

```text
F_k(K)=sum_(A subset K, |A|=r) w_A,
M_k^(m)=sum_(|K|=k) F_k(K)^m.                            (5.1)
```

Expanding the power groups ordered `m`-tuples of local objects by the size of
their union.  If every union has size at most `q<n`, binomial inversion gives

```text
(sum_A w_A)^m
 =sum_(k=r)^q (-1)^(q-k) binom(n-k-1,q-k) M_k^(m).       (5.2)
```

This is a finite cluster/moment expansion: each successive layer resolves a
new overlap scale.

For a seven-vertex tournament, let `c_A` be the cyclic-triangle indicator,
let

```text
C_3=sum_A c_A,
q_K=sum_(A subset K, |A|=3)c_A,
M_k=sum_(|K|=k)q_K^2.                                  (5.3)
```

An ordered pair of triangles has union size `3,4,5`, or `6`, corresponding
respectively to the same triangle, two triangles sharing an edge, sharing
only one vertex, or being disjoint.  If `N_u` counts ordered pairs with union
size `u`, then

```text
N_3=M_3,
N_4=M_4-4M_3,
N_5=M_5-3M_4+6M_3,
N_6=M_6-2M_5+3M_4-4M_3,                                (5.4)

C_3^2=M_6-M_5+M_4-M_3.                                 (5.5)
```

In particular `N_4=2C_4`, where `C_4` is the number of four-sets inducing
the strong four-vertex tournament, equivalently the number of unoriented
directed Hamiltonian `4`-cycles.  So a square is literally the first
two-body interaction of cyclic triangle faces.

This gives a systematic recursive program.  Replace a single profile count
by its overlap spectrum, then replace each overlap count by its signed,
rooted, or flag-refined version.  The hierarchy is analogous to moments and
cumulants: unflagged layers are marginals, while overlap flags carry the
connected information that reconstruction needs.

## 6. The exact four-card packet

For a labeled four-set `K`, define

```text
q_K = number of cyclic triangle faces of T[K],
h_K = number of universal sources in T[K]
      -number of universal sinks in T[K].                (6.1)
```

Then `q_K` is only `0,1,2`, and `(q_K,h_K)` classifies all four tournament
types:

| type | local score multiset | `q_K` | `h_K` |
|---|---:|---:|---:|
| transitive | `0,1,2,3` | 0 | 0 |
| source over `C_3` | `1,1,1,3` | 1 | +1 |
| `C_3` over sink | `0,2,2,2` | 1 | -1 |
| strong `C_4` type | `1,1,2,2` | 2 | 0 |

The first coordinate is an additive boundary:

```text
q_K=sum_(A subset K, |A|=3)c_A.                         (6.2)
```

At `n=7`, every labeled `c_A` is recovered from all `q_K` by `(4.5)`.  The
second coordinate is the missing dual-odd coorientation.  Globally it is a
score moment:

```text
D_4=sum_K h_K
   =sum_v [binom(d^+(v),3)-binom(d^-(v),3)].             (6.3)
```

Let `R_4` count the strong four-sets/directed `4`-cycles, and let
`N_0,N_+,N_-,N_2` count the four rows of the table.  Then

```text
N_2=R_4,
N_+ =((n-3)C_3-2R_4+D_4)/2,
N_- =((n-3)C_3-2R_4-D_4)/2,
N_0 =binom(n,4)-(n-3)C_3+R_4.                           (6.4)
```

Thus the complete unlabelled square profile consists of:

1. the inherited triangle count;
2. the new even square interaction `R_4`;
3. the new dual-odd chirality `D_4`;
4. the total number of four-sets.

Neither “triangles determine squares” nor “squares are independent of
triangles” is correct.  The exact statement is the decomposition `(6.4)`.

## 7. Edge tournaments as flat triangle tournaments

Fix antisymmetric edge signs

```text
x_ij=-x_ji in {+1,-1},
```

and define the alternating triangle flux

```text
h(i,j,k)=x_ij x_jk x_ki.                                (7.1)
```

Vertex switching `x_ij -> s_i x_ij s_j` fixes every `h`.  On each ordered
four-set `i<j<k<l`, the fluxes satisfy the Bianchi/flatness equation

```text
h_ijk h_ijl h_ikl h_jkl=1.                              (7.2)
```

Conversely, every alternating triangle field satisfying `(7.2)` comes from
an edge tournament.  Choose an apex `o`, switch so `x_oi=1`, and set

```text
x_ij=-h(o,i,j).                                         (7.3)
```

Equation `(7.2)` supplies all remaining faces.  Therefore

```text
dimension of edge-derived triangle fields =binom(n-1,2),
codimension inside all triangle fields     =binom(n-1,3),
size of every edge-lift fibre               =2^(n-1).    (7.4)
```

Only the fraction

```text
2^(-binom(n-1,3))                                      (7.5)
```

of independent triangle tournaments is edge-derived.  The first hostile is
a four-set whose four triangle signs have product `-1`.

Flatness is not directedness.  On one tetrahedron, `(7.2)` permits eight
face-sign patterns, whereas only the two alternating boundary patterns are a
directed tetrahedron.  Leader--Tan's rule—keep the cyclic orientation on a
directed triangle and follow the minority edge on a transitive triangle—is
exactly `(7.1)`.  Under this lift, a four-set is a directed tetrahedron if and
only if its ordinary tournament is a universal source or sink over a cyclic
triangle.  Hence its directed-tetrahedron count is

```text
N_+ + N_-=(n-3)C_3-2R_4.                               (7.6)
```

A random edge tournament makes this event with probability
`(8+8)/64=1/4`, giving Leader--Tan's sharp `c_3=1/4` lower construction.

More generally, after division by a coherent reference orientation, an
orientation on every `r`-set is an `(r-1)`-cochain on the full simplex.  It
comes from signs on `(r-1)`-sets exactly when it is flat on every `(r+1)`-set.
Since the simplex has no positive-dimensional cohomology,

```text
rank(delta_p)=binom(n-1,p+1),                            (7.7)
```

so the lower-face-derived `r`-tournaments have dimension
`binom(n-1,r-1)` and codimension `binom(n-1,r)`.  Literal iteration cannot
produce new higher signs because `delta^2=0`; useful edge-to-higher-face rules
must be nonlinear.

## 8. Square holonomy is a different square object

The four-card type above must not be conflated with Hamiltonian square
holonomy.  Fix a reference tournament `rho` and put

```text
a_ij=x_ij rho_ij=a_ji.
```

For each unoriented Hamiltonian four-cycle `C`, define

```text
Q_C=product_(e in C)a_e.                                (8.1)
```

There are `3binom(n,4)` raw square signs.  Their realizable span is the
even-cycle subspace of `K_n`, of dimension

```text
binom(n-1,2)-1=n(n-3)/2.                                (8.2)
```

The edge-to-square map therefore has an `n`-dimensional kernel: `n-1`
vertex-switch cuts plus global edge reversal.  A square field is realizable
exactly when the product of its values is `+1` for every family of four-cycles
whose symmetric difference is empty as an edge set.

On one `K_4`, its three Hamiltonian-square values multiply to `+1`, but these
local equations do not suffice from `n=5` onward.  A minimal hostile sets

```text
Q_(12341)=Q_(13241)=-1
```

and every other square value to `+1`.  Every individual `K_4` triple then has
product `+1`, but the three cycle vectors `13241`, `13251`, and `14251` xor
to zero while their assigned values multiply to `-1`.  The replay supplies
the explicit edge lists and checks the local/global distinction.

This is the same type boundary encountered in gain networks.  Triangle
holonomies generate the Abelian cycle field on a complete graph; a square
holonomy is their product under either triangulation, and the tetrahedral
Bianchi identity makes the result triangulation-independent.  On a bipartite
graph, however, triangles are absent, so squares become the first observable
phase carrier.  That is why the dephasing expansion sees triangle flux first
in a general graph but square flux first in a bipartite Jacobian-response
dilation.

## 9. A finite-exact `c_5 >= 9/64` lift

Leader--Tan call an alternating orientation of every `d`-set a
`d`-tournament.  A `(d+1)`-set is directed when its `d` face orientations are
the alternating boundary of one of its two orientations.

An isomorphism-equivariant local rule can orient each `d`-set using only its
ordinary induced edge tournament.  Such a rule is especially finite.  Every
automorphism group of a tournament has odd order: an even-order group would
contain an involution, but a nontrivial involution cannot preserve the edge
between two swapped vertices.  Hence every tournament automorphism is an
even permutation.

Choose one labeled representative `R_i` of every ordinary `d`-tournament
type and one sign `epsilon_i`.  For an arbitrary labeled tournament `U`, let
`pi:U->R_i` be an isomorphism and orient `U` by

```text
epsilon_i sign(pi).                                     (9.1)
```

This is well-defined because any two choices differ by an even automorphism,
and every equivariant alternating local rule has this form.

For `d=4`, there are four types and only `16` rules.  Exhaustion over all
`2^10` edge tournaments on five vertices gives a unique rule up to global
negation, with `144/1024=9/64` directed five-sets.  This reproduces the
Leader--Tan lower bound and proves it optimal within this local-rule class.

For `d=5`, there are twelve types and `4096` rules.  Use ascending edge-bit
order

```text
01,02,03,04,12,13,14,23,24,34.                          (9.2)
```

The twelve least-mask representatives are

```text
0,2,4,5,8,9,10,11,12,40,41,76,                         (9.3)
```

and one maximizing sign vector has positive bits

```text
1,1,1,0,0,0,0,1,1,0,1,0,                              (9.4)
```

namely rule integer `1415`.  Its global negative is the only other maximum.
For every one of the `2^15=32768` edge tournaments on six vertices, orient
its six five-faces by `(9.1)` and test the alternating boundary condition.
Exactly `4608` pass, proving `(1.2)--(1.3)`.

The table has a compact invariant formula.  For edge signs on a labeled
five-set, put

```text
h(v,a,b)=x_va x_ab x_bv,
H^(v)=(h(v,a,b))_(a,b != v),                             (9.5)
```

where the four remaining indices are in ascending order.  Then

```text
f_5(T)=1/4 [sum_(v=0)^4 (-1)^v Pf(H^(v))
            -product_(i<j)x_ij] in {+1,-1}.             (9.6)
```

Each Pfaffian has three terms, so the first sum consists of the `15`
bowties—pairs of triangle fluxes sharing one vertex.  The last term is the
full ten-edge parity, itself a switching-invariant Eulerian cycle product.
The Fourier expansion of `f_5` has exactly these `15` degree-six monomials,
each with coefficient `+/-1/4`, and the degree-ten monomial with coefficient
`-1/4`.  Exact reduction on the Boolean edge cube gives

```text
[sum_v (-1)^v Pf(H^(v))-product x_ij]^2=16,              (9.7)
```

and agrees with rule `1415` on all `2^10` labeled five-tournaments.

There is an even cleaner recursion.  Let `f_4` be the optimal four-face rule
(representative sign bits `-,-,-,+` on the replay's canonical four types).
Given a five-tournament `T` and any apex `v`, switch it by

```text
s_v=1,       s_u=x_vu                                   (9.8)
```

so that `v` dominates every other vertex, and call the normalized deletion
`U_v`.  Then, independently of the chosen apex,

```text
f_5(T)=(-1)^(v+1) f_4(U_v).                              (9.9)
```

This identity holds for every tournament, not only the favorable ones, and
explains `(9.6)` as a recursively consistent triangle/square potential.

Now switch a random six-tournament so its last vertex dominates all others.
Its remaining five-tournament `U` is uniform, and exact boundary comparison
using `(9.9)` gives

```text
T_6 is f_5-directed  iff  U is f_4-directed.            (9.10)
```

The right side holds for `144` of the `1024` choices, giving `9/64` without
recounting all `32768` lifts.  The recursion stops here: every monomial in
`f_5` is Eulerian, so `f_5` is switching-invariant; a further apex potential
would be consistent only on the already directed six-sets, not on every
six-tournament.

The favorable locus has an exact external identification.  The first
successful six-vertex edge mask is `10`.  Its closure under all `S_6`
relabelings and all vertex switches has size

```text
(6!/5)2^5=4608,                                        (9.11)
```

and every member succeeds.  Its triangle-flux oriented two-graph has
automorphism group of order `5`.  Its switching class contains eight ordinary
tournament types, with representatives and automorphism orders

```text
representatives: 10,21,81,140,147,148,313,314
Aut orders:       1, 1, 1,  5,  1,  1,  1,  5.         (9.12)
```

These data exactly match Gunderson--Semeraro's special oriented two-graph
`g`, whose switching class gives their cited bound
`pi(H(6))>=9/64` for the six-uniform three-edge Turan problem.  The new bridge
is that the same class is exactly the directed-six-set locus of the local
five-face potential `(9.4)`.

The invariant formula `(9.6)` replaces the opaque table; the next proof target
is a classification-free derivation that its directed locus is `g`.  For
`d=6`, the same
optimization has `56` Boolean type variables, so a raw `2^56` search is no
longer sensible.  Each seven-vertex edge tournament imposes a small affine
constraint on its six-card type variables, turning the problem into a
symmetry-compressed weighted Max-CSP/Ising optimization over the `456`
seven-vertex types.

## 10. Oriented matroids: four-sets carry the first determinant law

Independent higher tournaments are far more general than determinant sign
data.  A uniform rank-`r` oriented matroid has a chirotope `chi` satisfying
Grassmann--Pluecker.  In rank two, for four nonparallel vectors,

```text
Delta_ab Delta_cd-Delta_ac Delta_bd+Delta_ad Delta_bc=0. (10.1)
```

Therefore the three signs

```text
chi(ab)chi(cd),
-chi(ac)chi(bd),
chi(ad)chi(bc)                                           (10.2)
```

cannot all agree.  This is a genuinely four-vertex constraint: an edge-only
tournament does not expose the three perfect-matching channels in `(10.1)`.

In rank three, a flat oriented-two-graph triangle field need not be a
chirotope.  A minimal five-vertex hostile assigns `+` to every increasing
triple except `chi(135)=-`; for `lambda=1` and the remaining indices
`2,3,4,5`, its three Grassmann--Pluecker terms all have the same sign.

Finite exact counts produced a clean conjecture.  The numbers of flat
rank-three chirotopes on labeled sets of sizes `5,6,7` are

```text
24,120,720=(n-1)!.                                     (10.3)
```

These are exactly the counts of oriented circular orders.

> **OPEN conjecture.**  For `n>=5`, an oriented-two-graph triangle field that
> is also a uniform rank-three chirotope is exactly a circular order.

Cameron's characterization of circular two-graphs suggests a short route:
show that each noncircular flat `K_4` pattern violates circuit elimination
when extended by one fifth point.

### Planar-Jacobian consequence contract

For exponent vectors in a Jacobian coefficient fibre,

```text
[x^w] Jac(P,Q)
 =sum_(u+v=w+(1,1)) det(u,v)p_u q_v.                    (10.4)
```

The existing gain graph treats contributors pairwise.  Equation `(10.1)`
shows that the first intrinsic determinant syzygy is instead supported on a
four-set and couples its three perfect matchings.

| Contract field | Exact content |
|---|---|
| Source | four exponent vectors in one or adjacent Jacobian coefficient fibres |
| Target | the three-term Pluecker circuit `(10.1)` |
| Preserved | alternating determinant signs, zeros, and circuit compatibility |
| Lost by sign-only higher tournament | determinant magnitudes, coefficient phases, and the global `p_u q_v` Segre constraints |
| Required sidecar | actual minors/magnitudes, zero pattern, target fibre, and coefficient factorization |
| Cheapest test | add all four-set Pluecker faces to the first live leaf-peeled response passport and compare its surviving dimension with the edge-only gain graph |

This does not prove or disprove the planar Jacobian conjecture.  It identifies
a precise reason that a square/face complex can prune candidates that an
edge-only cancellation graph cannot see.

## 11. Cross-domain connections with explicit loss ledgers

### 11.1 Johnson/Kneser walks and dephasing networks

For `n=2r+1`, the middle incidence operator is the adjacency of the odd graph
`O_(r+1)`.  A coherent walk on that graph and its strongly dephased resistor
walk therefore propagate data between the two middle subset layers.  The
`S_(2r+1)` harmonic sectors have eigenvalues `(4.4)`, and the exact inverse is
a signed spectral filter.

- **Preserved:** subset intersection harmonics and all additive face data.
- **Destroyed by dephasing:** the signed interference needed to implement
  negative coefficients in `(4.3)` directly.
- **Sidecar:** coherent phase control or a signed post-processing filter.
- **Test:** for `n=7`, implement the cubic filter `(4.6)` and recover a random
  triangle field from its four-set aggregates.

### 11.2 Flag algebras and formal zeta dimensions

Unlabelled tournament profiles form the degree pieces of the tournament flag
algebra; overlap-labeled cards supply its multiplication/gluing data.  This is
a controlled analogue of the shuffle/stuffle distinction for formal multiple
zeta values: graded marginal counts alone do not settle the kernel of the
evaluation map.

There is no map from tournament densities to zeta periods here.  What is
preserved is only the proof architecture—graded ranks, products, and hidden
kernel audits.  The missing MZV sidecar is period injectivity, exactly the
coordinate that prevents a formal upper-bound theorem from proving Zagier's
actual-number dimension conjecture.

### 11.3 Exceptional-prime audits

The odd-middle transform has completely explicit bad primes through `(4.4)`.
This suggests a generic validity gate for any reported “all primes except one
huge prime” theorem: identify the integral relation matrix, exhibit the exact
minor or Smith factor containing that prime, and show why all other primes
are units.  Without such a certificate, a lone decimal or social-media report
has no analogue of the Kneser determinant mechanism.

### 11.4 Designs, trades, and Hodge pieces

Cochain flatness `delta z=0` tests lower-face origin.  Chain balance
`partial z=0` instead describes signed designs, trades, or cycles.  They must
not be conflated.  On an incomplete complex the Hodge decomposition adds a
genuine harmonic sidecar; on the full simplex positive-degree real homology
vanishes, but Boolean `+/-1` feasibility remains nonlinear.  Boundaries of a
tetrahedron or a four-simplex are the first local cycle carriers.

## 12. Best next frontier

1. **Direct `g`-locus proof (`PROOF TARGET`).**  Starting from the bowtie
   Pfaffian formula `(9.6)` and recursion `(9.9)`, prove without catalogue
   classification that the directed locus is the special oriented two-graph
   `g` and recover its order-five automorphism.
2. **Flat chirotope iff circular (`OPEN`).**  Reduce to the noncircular
   four-vertex restrictions and one fifth point; certify with one
   Grassmann--Pluecker/circuit-elimination argument.
3. **Pair-flag reconstruction depth (`OPEN + COMPUTATIONAL`).**  For the
   Stockmeyer--Kocay families, add joint types of two overlapping deleted
   cards and find the least flag order that repairs reconstruction.
4. **Connected motif cumulants (`PROVED FRAMEWORK, OPEN CLASSIFICATION`).**  Use
   `(5.1)--(5.2)` for three or more cyclic triangles and split by intersection
   hypergraph, not only union size.
5. **Square consistency presentation (`OPEN`).**  Test whether the local
   `K_4` relations plus `K_(2,3)` three-square relations generate the full
   even-cycle relation space through `n=12`.
6. **Higher local-rule optimization (`OPEN`).**  Orbit-compress the `d=6`
   problem to a weighted Max-CSP on the `56` six-tournament types; use exact
   branch-and-bound or an SDP/Fourier upper certificate.
7. **Jacobian Pluecker face complex (`OPEN`).**  Add four-set determinant
   circuits to the current response gain graph before phase solving, retaining
   all coefficient magnitudes and the target-range sidecar.

## 13. Sources and scope

Primary external sources used in this session:

- [Leader--Tan, *Directed Simplices in Higher Order Tournaments*](https://arxiv.org/abs/0906.4027).
- [Babai--Cameron, *Automorphisms and Enumeration of Switching Classes of Tournaments*](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v7i1r38).
- [Gunderson--Semeraro, *Turan Numbers and Switching*](https://arxiv.org/abs/2204.10775).
- [McKay, *Reconstruction of Small Graphs and Digraphs*](https://ajc.maths.uq.edu.au/pdf/83/ajc_v83_p448.pdf),
  with the [seven-tournament catalogue](https://users.cecs.anu.edu.au/~bdm/data/tourn7.txt).
- [Stockmeyer, *The Falsity of the Reconstruction Conjecture for Tournaments*](https://doi.org/10.1002/jgt.3190010108).
- [Kocay, corrected constructive proof](https://doi.org/10.1002/jgt.3190090406).
- [Linial--Morgenstern, *On High-Dimensional Acyclic Tournaments*](https://arxiv.org/abs/1302.1684).

Status boundaries:

- `PROVED`: identities `(3.3)--(7.7)`, the odd-middle inverse, the four-card
  packet, cochain dimensions, and the implication from an exhibited local
  rule to a `c_5` lower bound.
- `FINITE-EXACT`: the seven-vertex catalogue census, affine hostile,
  local-rule optimizations, switching-orbit identification, square controls,
  and low-order flat-chirotope counts.
- `CITED`: reconstruction history, Leader--Tan bounds, oriented-two-graph
  background, and Gunderson--Semeraro's `pi(H(6))>=9/64` bound.
- `OPEN`: a classification-free proof of the `g` locus, the circular-order
  conjecture, `d>=6` local-rule optima, and every proposed planar-Jacobian
  pruning effect.
- Not claimed: a full tournament reconstruction theorem, a new value of
  `c_5`, a formal literature-priority claim, a Jacobian counterexample, or a
  transfer theorem between tournament flags and multiple zeta values.
