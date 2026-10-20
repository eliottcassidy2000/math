# First-exit charge, complementary path games, and hinged colour boundaries

2026-10-20. **PROVED scoped deductions; FINITE-EXACT symbolic and geometric
certificates; same-agent review, not an independent-agent audit.** Collatz
convergence, the unrestricted first-exit positive-charge inequality, and a
six-chromatic unit-distance construction remain **OPEN**. No priority claim
or formalization is made. [Inheritance and connections](INHERITANCE.md),
[validation](AUDIT.md).

## 1. A sharp all-length bound for the positive-charge threshold

Use the odd map `U(n)=oddpart(3n+1)`, positive valuation word
`w=(a_1,...,a_m)`, `K=sum a_i`, `d=2^K`, and `A=3^m`. Set

    B = sum_i 3^(m-i) 2^(a_1+...+a_(i-1)),
    C = sum_(a_i=1) 3^(m-i) 2^(a_1+...+a_i).             (1)

Indices in formulas run from1 to m. The actual affine endpoint is
`(An+B)/d`. In `z=n-1` discard only the negative constants at steps
`a_i>=3`. The resulting positive-charge upper bound is

    Uplus_w(n) = [A(n-1)+C]/d.                           (2)

This is the earlier packet's `A_slope(n-1)+P(w)`, with `A_slope=A/d`
and `P=C/d`. The actual endpoint minus1 is `Uplus-D` where

    D = sum_(a_i>=3) (1-2^(2-a_i)) 3^(m-i)/2^(K-K_i),
    K_i=a_1+...+a_i.                                    (3)

For `d>A`, payment is exactly `n>R(w)`, where

    R(w)=1+C/(d-A).

**Proposition 1 (PROVED, inheriting the stated logarithmic input).**
For every contracting positive valuation word,

    R(w)/d <= 161/256.                                  (4)

Equality holds only at `w=(4,1,1,1,1)`, for which
`A=243,d=256,C=2080,R=161`. This bound concerns every contracting word,
not just first exits. The maximizing word itself is not a first exit at
length5 for any positive odd source greater than1.

**Proof.** Move every exponent1 to the right of the other exponents.
Swapping adjacent `(1,b)` to `(b,1)`, `b>=2`, increases C: before
the common suffix the difference is `2^H(2^(b+1)-6)>0`, where H is
the preceding clock. This is an algebraic maximization, not a legal change
to the source's trajectory. If r exponents equal1, the resulting numerator is

    C = 2^(K-r+1)(3^r-2^r).

At fixed K this increases strictly with r. A contracting word cannot
consist entirely of ones, so `r<=m-1`. Consequently

    C <= Cmax = 2^(K-m+2)(3^(m-1)-2^(m-1)),              (5)

with equality at `(K-m+1,1,...,1)` (the trivial m1 case has C=0).
At fixed m, `1/d+Cmax/[d(d-A)]` decreases strictly with d. It is enough
to take `K=bit_length(3^m)`. Exact rational arithmetic for `1<=m<=64`
gives the unique maximum `161/256` at `(m,K)=(5,8)`.

For the infinite tail we reuse, rather than re-prove, the actual carry
bound from [THM-4512, coefficient-descent classes](../../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md)
and [BC1, steps2--3](../collatz_boundary_compiler_20261004.md#2-the-sharp-all-length-carry-bound):
for every `m>=65`, `B/[d(d-A)]<1/1024`. That proof includes its rational-log
bridge through `2^42` and cited explicit Matveev input beyond it. Formula(1)
gives `C<=2B` term by term. Since `d>=2^m`,

    R/d < 2^-65+1/512 < 161/256.

This covers all lengths and proves uniqueness. The new finite sharp
calculation is in `coarse_charge_audit.py`; the infinite-tail input is an
explicit inherited dependency, not an unexecuted new finite search. QED.

**Corollary 1 (one remaining source).** Let
`r=-B A^(-1) mod d`, `0<r<d`, be the least positive member of the coarse
cylinder. Every source `n=r+td`, `t>=1`, pays (2). Only r can fail, and
failure requires `r<=161d/256`.

The exact word has modulus `2d`, with residue `(d-B)A^(-1) mod2d`.
The coarse class fixes the earlier valuations and only a lower bound on
the last. Since contraction forces the last exponent to be at least2
when this is the first coefficient crossing, extra final divisions
decrease (2). For a general contracting word whose last exponent is1,
an extra final division also removes its positive unit emission and
decreases the bound. Thus the sufficient payment conclusion is valid on
the coarse class; no exact-endpoint equality on that larger class is claimed.

This sharpens the **charge** question using the earlier one-exception
compiler. It does not prove that the least representative pays, or that
every source ever encounters a contracting word.

## 2. First-exit certificates, and a genuine history obstruction

**Proposition 2 (PROVED).** Suppose the actual word first descends from odd
`n>1` to `h<n`, and contains r exponents at least3. Then `D<r` if `r>0`
and `D=0` if `r=0`. In particular, payment holds whenever `r<=n-h`.
Every first-exit word with at most two such resets therefore pays.

**Proof.** A reset at the final step contributes less than1. A reset at
an earlier index i has actual source for its suffix `x_i>=n`. The
unshifted suffix has positive affine constant, so its slope is
`(h-positive_constant)/x_i < h/x_i <1`. Each reset contribution in(3)
is consequently less than1. Sum, then use
`Uplus_w(n)=h-1+D<n-1`. Odd endpoints make `n-h>=2`. QED.

The more informative, still source-dependent version is

    D <= eta_m + h sum_(i<m,a_i>=3) eta_i/x_i,
    eta_i=1-2^(2-a_i), eta_m=0 if a_m<3,                 (6)

with strict inequality whenever there is an internal reset. This replaces
the exact suffix discounts by a reciprocal-value budget. A uniform bound
on that budget is a possible stronger target, not a proved consequence.

**Proposition 3 (FINITE-EXACT symbolic certificate).** For every positive
odd source greater than1 whose first coefficient crossing occurs within
18 odd steps, its positive-charge bound is paid at that crossing. Hence
its actual first descent occurs there. In particular the original
positive-charge hypothesis holds for **every source with actual first
descent of at most18 odd steps**, without a bound on the source magnitude.

Two different enumerations establish this finite language statement:

* `symbolic_charge.py` retains `(A,B,C,d,residue,source interval)` and
  splits exact cylinders according to actual first descent. It visits
  1,166,058 non-exited prefixes. A proved monotone bound handles each
  infinite tail of possible next valuations. There are no failures;
  1,900,470 non-exited depth18 cylinders remain, with least source27.
* `coarse_charge_audit.py` instead enumerates first **coefficient** crossings,
  joining every larger last exponent into a single coarse cylinder.
  It reconstructs C by suffix sums, rather than the producer's recurrence.
  All 1,166,058 exit cylinders pay at their least eligible source. Since
  payment increases with the source, this proves it for every source in
  each cylinder. The smallest integer payment margin is4.

The count by exit length is
`1,1,2,3,7,12,30,85,173,476,961,2652,8045,17637,51033,108950,312455,663535`.
This extends the inherited length16 **actual-carry** census to a stronger
positive-charge test at length18; it is not a claim to exceed published
verification ranges for the different coefficient-stopping conjecture.

The exact producer's infinite-tail rule is particularly reusable. If
the least source of a retained cylinder is L, then

    2^a > 3[A(L-1)+C]/[d(L-1)], a>=2                  (7)

pays every source in that cylinder at the next step of exponent a.
The right side decreases with the source. All lower exponents are
explicitly split; there is no guessed maximum valuation.

**Hostile boundaries.** A floating scout through `2^28-1` found no failure,
but that run is only **NUMERICAL screening**, not an exact certificate.
Its candidates are replayed in exact fractions:

* Source63 remains the largest screened dividend/descent ratio:
  `215402584767891/281474976710656`, at its first exit to61.
* Source4575 first exits at4465 with `D>2`; source123027327 first exits
  at120238921 with `D>4`. Thus dropping the endpoint gap and asserting
  a uniform `D<2` is false.
* Source55 has a later length39 word ending at53 which is **not paid**,
  although its endpoint is smaller. Its first strict descent was at step3.
  The later word's threshold is
  `33176444486369026739/559130865408411637 >55`.

Thus even actual descent of an authenticated word does not allow us to
discard reset dividends. The first-exit history is substantive information.
The old question is still **OPEN** at unbounded first-exit depth.

### A uniform logarithmic-depth region

**Proposition 3a (PROVED, with CITED Ellison input).** Let `n>=65536` be
odd. If its first coefficient crossing occurs at odd-step depth
`m<=4 log(n)` (natural logarithm), the positive-charge bound pays there.
Consequently actual first descent occurs at that crossing. The
integer-only sufficient condition `m<=2 floor(log_2 n)` also works.
This allows unbounded source magnitude **and** unbounded depth.
It also covers every actual first exit within that depth bound, since
the first coefficient crossing cannot occur later than actual descent.

For m<=18 use Proposition3. Otherwise every proper coefficient prefix
has `2^K_i<3^i`, so (1) gives `C<=2B<=2m*3^(m-1)`. Let
`k=bit_length(3^m)`; replacing the final clock K by k only increases the
threshold estimate. Since m>=19, k>=31. Ellison's
[Theorem3, printed p12-05](https://www.numdam.org/article/STNB_1970-1971____A10_0.pdf)
gives `2^k-3^m>2^k exp(-k/10)`; its exceptional exponents are all<=27.
Thus, using `log_2 3<8/5` and `exp(1/10)<10/9`,

    R < 1+(2m/3) exp(k/10)
      < 1+(20m/27) exp(4m/25)
      <= 1+(80/27) log(n) n^(16/25).                   (10)

The function `n^(9/25)/log(n)` increases for `log(n)>25/9`.
At `n=2^16`, it exceeds `56/15`: use
`2^(144/25)>2^(11/2)>224/5` and `16 log2<12`.
Consequently (10) is below `1+(50/63)n<n` throughout the claimed range.
Finally `log2>1/2` implies `2 floor(log_2 n)<4 log(n)`. QED.

The proof does **not** establish this depth bound on every orbit. It
identifies a concrete remaining direction: a counterexample to positive
payment at the first coefficient crossing, with n>=65536, must survive
past `4 log(n)` odd steps before that crossing. This stronger
coefficient-crossing question is distinct from the original hypothesis,
which assumes an actual first exit. No clock equality is silently assumed.

## 3. The four complementary graph games

Interpret the game as an undirected spanning graph G of Player1's edges
and its complement H of available edges. Player2 may choose the starting
vertex. A path visits every vertex once; a spanning walk may repeat edges
and vertices and need not use every available edge. A spanning walk exists
exactly when the graph is connected. If every remaining edge must be used
exactly once, that is a separate Euler-trail game.

**Proposition 4 (PROVED).** The complete existence table is:

| Forbidden in Player1's G | Forbidden in Player2's H | Winning choice exists |
|---|---|---|
| Hamilton path | Hamilton path | exactly when n>=4 |
| Hamilton path | spanning walk | exactly when n>=4 |
| spanning walk | Hamilton path | exactly when n>=4 |
| spanning walk | spanning walk | never |

For n>=4 take the complementary pair `K_(1,n-1)` and `K_(n-1) disjoint K1`,
in the needed order. Neither has a Hamilton path. The former is connected
and the latter disconnected. If G is disconnected, every two vertices
in different components are adjacent in H; vertices in one component
connect through any other component. Hence H is connected. Orders1--3
are immediate, also covered by the complete exact census.

**Sharp edge cost.** For the first two rows Player1 needs at least `n-1`
edges, and the spanning star attains it. For n>=5 it is the unique
minimum-size winning G up to isomorphism. At n4 the Hamilton/Hamilton
game also admits the triangle plus isolate at that cost.

Here is a self-contained degree-sum proof. Adding a nonedge uv to an
n-vertex graph preserves existence of a Hamilton path if
`deg(u)+deg(v)>=n-1`. Adjoin a universal vertex and use the following
cycle argument: if a newly added edge closes a spanning path
`u=v_1,...,v_N=v`, the two endpoint-neighbour index sets intersect when
their degree sum is at least N; the intersecting index reroutes the path
into a spanning cycle without the new edge. Delete the universal vertex.
Repeatedly adding nonedges with this degree condition gives a complete
graph whenever every nonedge qualifies.

If G has at most n-2 edges, each of its edges uv satisfies
`deg_G(u)+deg_G(v)<=|E(G)|+1<=n-1`. Thus every nonedge of H qualifies,
and H has a Hamilton path. If G has exactly n-1 edges and H is
nontraceable, every edge uv of G must instead have degree sum n:
otherwise remove uv and use the same closure lemma plus the preceding
bound. Degree sum n says that every edge of G meets uv. All edges of G
are therefore pairwise intersecting. Such a family is a star or a
triangle: two edges ab,ac force every edge avoiding a to be bc, which
then prevents any fourth edge. This proves the equality classification.

### Both colours can be connected and still both fail

**Proposition 5 (PROVED family; FINITE-EXACT minimum order).** The first
such example has n7, and examples exist for every n>=7.

Let k=n-4>=3. Take a triangle u,v,w, a vertex z joined only to w, and
k independent vertices each joined to u and v. This is connected G.
Deleting u,v leaves k isolated vertices and the component wz: `k+1>3`
components. A spanning path, after deletion of two vertices, has at most
three components, so G is not traceable. In H, u and v are two leaves
with the same neighbour z. Both would have to be path endpoints and
would exhaust z's two path incidences, omitting all other vertices.
H is connected because z and the k twins form a clique, w meets every
twin, and u,v meet z. Hence H is not traceable either.

The exhaustive census of every labelled graph through n6 gives no
both-connected winning pair, proving the stated minimum. For n4,5,6,
the Hamilton/Hamilton winning counts are8,70,1092; each consists of the
two complementary mixed walk/path families.

When k is even, G has exactly two odd-degree vertices, w and z, so it
has an Euler trail of length `2k+4=2n-4`, but still zero Hamilton paths.
For odd k it has four odd-degree vertices and no Euler trail. This gives
an infinite, concrete separation of the user's two measurements.

**Incoming stronger classification, integrated after the first checkpoint.**
The independently reviewed [complement-path packet](../complement_path_games_20261020/README.md)
and [sharp-budget packet](../kernels_games_detours_20261020f/game.md)
close the edge-cost question for the third row too: its minimum is
`2n-5`, except6 at n6 and10 at n8. Their proofs retain a Pósa-type degree
obstruction or the classical Erdős extremal bound, not only our small census.
The incoming connected-both family is isomorphic to Proposition5.

The former packet also resolves a prescribed-start reading: if Player2
must succeed from **every** start, defeating that guarantee in the
Hamilton column is possible from n3 and costs n-2 red edges. This does
not change the free-start table above. Its exact join criterion is useful
for future composition: for `G=join G_i`, a Hamilton path exists iff
`pc(G_i)<=n-|G_i|+1` for every i, with a rooted path-cover refinement for
prescribed starts. These questions are now closed in the cited scopes,
and are removed from this packet's open list.

## 4. A non-lattice nineteen-point gadget with an exact six-colour relation

Put `s=sqrt(-3)`, `t=sqrt(-11)`, `zeta=(1+s)/2`,
`r_j=zeta^j`, and `t_i=r_(2i)+r_(2i+1)`, i0,1,2.
The previous ten-point gadget is `V={0,r_0,...,r_5,t_0,t_1,t_2}`.
Translate it by `-t_0`, so one tip is the origin. Adjoin the rotated
copy `rho(V-t_0)`, where

    rho=(5+t)/6, |rho|=1.

The two copies share exactly their origin tip. Their hubs are `-t_0`
and `-rho*t_0`, and their distance is1 because `|t_0|^2=3` and
`|1-rho|^2=1/3`. In the basis `1,s,t,st`, exact rational arithmetic
checks all171 pairs. The induced graph has19 vertices and37 edges.
Its **only additional edge** beyond the two copies is the hub--hub edge;
the other four cross-copy unit pairs are already tip--rim edges incident
with the shared point. No assumed drawing or tolerance is used.
Adding t_0 back puts the first copy at its original coordinates: the
second is rotated about its named tip, with every original port fixed.

Keep the17 non-hub vertices as named ports. Let M,N be the subsets of
the six colours absent from the two six-vertex rims.

**Proposition 6 (PROVED exact extension relation).** For every proper
assignment to all17 ports, the number of extensions to the two hubs is

    K(M,N)=|M||N|-|M intersection N|.                     (8)

The proof is literal: choose one missing colour at each hub, excluding
equal choices across their edge. Thus both individual wheels can extend
while their union cannot: take `M=N={c}`. Missing palettes of the same
size with `M={c},N={d}`, c!=d, give one extension instead.

With the common tip fixed to colour0, an example rim missing0 is
`(1,2,1,3,4,5)`, and one missing1 is `(2,3,0,2,4,5)`.
All unshared tips have available colours. These give a concrete
blocked/accepted pair with identical local extension counts.

The full count is computed by factoring through the missing-colour masks:
there are166720 proper port assignments per copy with the shared tip
fixed, hence27795558400 joint assignments. Exactly836075520 are blocked
despite **both** copies admitting individual extensions. Another2501836800
are already blocked in a single copy. All mask frequencies and the
extension-count histogram are retained in `geometry-games.json`.

If v_M is the six-bit incidence vector of M, then

    K(M,N)=v_M^T (J_6-I_6) v_N.                         (9)

The kernel has rank6 on the realizable palette masks, because every
singleton occurs and their submatrix is `J_6-I_6`. It is indefinite
(eigenvalues5,-1,-1,-1,-1,-1 on that submatrix), despite nonnegative
entries. A count-only or positive-Gram replacement loses the common
colour label that makes (8) vanish.

**Boundary.** The graph contains the classical seven-vertex Moser spindle
on the common tip, the two hubs and their four shared-tip rim neighbours.
In a three-colouring, each adjacent pair of rim neighbours forces its
hub and the common tip to share the third colour; the hub edge then
contradicts that equality.
The exact checker rejects three colours and supplies a four-colouring.
So this is a genuine escape from the single triangular lattice, with a
new six-colour **boundary relation**, not a six-chromatic graph.
Moreover [THM-4558, field-plane residue colouring](../../../01-canon/theorems/THM-4558-residue-colourings-of-field-planes-the-heegner-tower-is-at-most-four-chromatic.md)
already bounds the entire field `Q(sqrt(-3),sqrt(-11))` by four colours.
Repeated constructions remaining in that field cannot reach six colours.
The next target needs both a further carrier escape and a context forcing
the same missing singleton at adjacent hubs. Keeping only palette sizes
will not formulate that target correctly.

The new incoming sixteen-point gate yields a further exact composition:
[two copies sharing their output](PALETTE-GLUE.md) have31 vertices and62
edges, with no unintended cross contacts. Two prescribed three-colour
input palettes A,B give exactly `576*|A intersection B|` extensions;
`A={0,1,2},B={0,3,4}` forces output0. This remains a conditional boundary
gate. The incoming18-point construction is smaller when only a single
externally prescribed output is required. The31-point result supplies
an exact two-input composition and a general finite rotation-avoidance
lemma, with all first-copy ports fixed.

## 5. What relates Hamilton paths, Euler trails and torus parameters

There is an exact graph map: a directed edge-word traversing every edge
of D once is a Hamilton vertex-word of its directed line graph L(D).
Vertices of L(D) are the **labelled edges** of D and adjacency means
head-to-tail compatibility. This preserves the ordered witness and
length, not the original graph's vertex set. For undirected line graphs,
pairwise incidence alone is insufficient: the claw's line graph is K3,
but its three pairwise incident edges cannot be one Euler trail. Endpoint
orientation is the missing sidecar, consistent with the earlier K4 atlas.

For a genuine torus construction, retain p ordered strands and q full
applications of `sigma_1 ... sigma_(p-1)`. Its strand permutation is
`i -> i+q mod p`, with `gcd(p,q)` cycles of length `p/gcd(p,q)`.
The closure is a knot precisely when that gcd is1. The positive braid
word has `q(p-1)` crossings; this is a **presentation length**, not a
universal minimal-crossing formula. This is the standard torus-braid
model, also retained in [THM-2646, central braid height](../../../01-canon/theorems/THM-2646-braid-three-modular-central-pullback-and-full-twist-knot-fibre.md).

The permutation graph only remembers q modulo p. In particular
`T(2,1)` and `T(2,3)` have identical two-strand permutation graphs,
Hamilton-path counts and Euler-cycle lengths, while one is an unknot
and the other a trefoil. Thus no recovery of the pair(p,q) from those
two graph statistics follows. A cyclic strand order and **integer twist
height** restore a meaningful map. This is the same representational
discipline as retaining the Collatz integer clock and the actual missing
colour names; it does not identify their dynamics.

The incoming [integer-winding construction](../kernels_games_detours_20261020f/walks_knots.md)
makes the sidecar local: attach `c_i=floor((i+q)/p)` to the arc from i to
`(i+q) mod p`. Telescoping gives total carry q and carry q/gcd(p,q)
on each component. A finite graph can therefore preserve an unbounded
integer exactly. This restores representation, without supplying a
decreasing Collatz potential.
