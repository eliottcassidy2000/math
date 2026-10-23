# Endpoint-marked composition gives guarded cross-orbit cubes

23 October 2026. **PROVED:** the composition, cube, routing and conditional
weight bounds below. **FINITE-EXACT:** the eight-owner seed and the complete
409-path pointed atlas. **OPEN:** coverage or mixing of unrestricted guarded
fibres. The local native event here is an authenticated two-owner path
replacement, not necessarily an adjacent swap or an old metagraph tile edge.

## 1. Inherited mechanism and boundary data

The source is the actual rigid order-eight tournament of merged class1238,
ordinary canonical code29558785, as saved in `native_cube.py`. Its two paths are

    P0=(3,0,6,7,1,5,4,2), P1=(3,1,6,7,0,5,4,2).

The literal map `g=(4,0,2,3,1,5,6,7)` is an isomorphism from K_P0(T) to
K_P1(T), fixing owners3 and2. Here K_P keeps the consecutive path arcs and
reverses all other arcs. The seed is inherited from
[Session44 stable padding](../session44_guarded_native_moves_20261023/stable-padding/THEOREM.md).
The more general warning that complement/substitution needs actual exit and
entry arcs is already in
[Session25 spectrum wiring](../session25_spectrum_wiring_20261021/session25-wiring/RESULT.md).
This work extends that marked mechanism; unpointed partner equality alone
does not justify the construction.

**Composition lemma.** For tournaments T_i and Hamiltonian paths P_i,Q_i,
suppose g_i:K_Pi(T_i) -> K_Qi(T_i) sends each actual first endpoint to the
corresponding first endpoint, and likewise for the last. The blockwise map
g is an isomorphism between the partners of the concatenated paths in the
ordinal source T_1 -> ... -> T_r.

**Proof.** Internal arcs obey the given isomorphisms. Between distinct blocks
every arc reverses, except the one actual last-to-first junction between
successive blocks. The endpoint conditions preserve exactly these exceptions.
Every Hamiltonian path of the ordinal source visits the blocks consecutively:
after leaving one it cannot return. This also proves the exact path-count
product. Empty external padding simply removes a junction condition. QED.

## 2. A complete pinned fibre crossing exponentially many source orbits

For every r>=1 take T_r to be the ordinal sum of r copies of T. Pin clocks
0,2,3,5,6,7 in each eight-owner block to the corresponding owners of P0.
In block i the two remaining clocks contain precisely owners0 and1 in either
order. Both orders are legal. Thus the **complete** path family satisfying
these pins is

    X_r={P_(b1)...P_(br): b in {0,1}^r},       |X_r|=2^r.

No extra target filter is needed to obtain completeness. The composition
lemma proves all partners isomorphic. Therefore every target-class-dependent
guard, at any modulus, is unchanged on X_r with the source fixed. This is
not a numeric assertion about a newly computed higher-order electrical atlas.

The seed is strongly connected and has trivial automorphism group (all8!
permutations checked). Hence the strong components of T_r are its ordered
eight-owner blocks, and Aut(T_r)=1. The 2^r paths occupy **2^r different
source-automorphism orbits**. A cube edge changes two actual owner clocks,
8i+1 and8i+4, and has an explicit partner map g or g^-1 in that block.
The common indexed arcs(6,7) and(5,4), and the actual outer endpoints,
remain fixed. Converse reverses the source and path, transporting clock j
to 8r-1-j; it does not silently keep the old printed clock.

These are literal **source-path** pins. The partner isomorphism need not fix
every interior target owner (for example g(4)=1). Arbitrary extra target-side
attachments at those owners are not preserved by this lemma; they require
additional marked transports. Inverse target ports are retained as complete
actual records, not identified merely because their target classes agree.

Literal normalized source words are distinct, since a coincidence would
give a source automorphism. Complementing their non-path bits gives distinct
inverse target ports. Their converse copies are distinct as well. In fact
the full four-way reservoir has4*2^r ports. For the seed, the source and
partner degree multisets are respectively
`(3,3,3,3,3,3,4,6)` and `(1,3,3,4,4,4,4,5)`; the source is not self-converse,
and the partner is strong. For r>1, the source is decomposable whereas its
partner is strong: each strong partner block has arcs in both directions
to its next block, including its retained junction. Thus source and inverse
types remain disjoint at all r. Complement commutes with converse, so the
remaining disjointness follows too. Any unchanged background must avoid
the actual reservoir consumed by a proposed route.
In particular these are2^r different normalized K/converse tile resources,
not merely different raw lifts of fewer tiles: their four-word sets are
pairwise disjoint. This stronger fact authenticates the unit-resource menus
used in PORTS.md under either carrier.

## 3. Exact source-weighted cost and gap

Let P choose a coordinate uniformly and resample it fairly. Thus a specific
cube-edge transition has probability1/(2r) and the hold probability is1/2.
Walsh characters have eigenvalues `1-|S|/r`, so

    gap(P)=1/r.

This is a statement about the conditional uniform measure on the literal
pinned family. It never assigns one unit of mass to each unweighted quotient
class. For arbitrary positive literal weights w, define the Metropolis chain
P_w with stationary distribution pi=w/Z. If max(w)/min(w)<=kappa, then

    gap(P_w)>=1/(kappa*r).

Indeed its edge conductances are min(w(x),w(y))/(2rZ). Relative to uniform
Dirichlet energy they are at least N*min(w)/Z; its variance is at most
N*max(w)/Z times uniform variance. Taking the ratio proves the claim.
The source weights have not been changed; the ratio is an explicit hypothesis.
Unbounded or zero weights are outside this estimate.

Flipping differing coordinates in increasing order gives a route of at most
r events, with exactly two changed owner clocks per event. For any fixed
oriented edge, exactly N/2 ordered endpoint pairs use it. Their mean route
length is(r+1)/2. Thus the exact uniform all-pair length-weighted congestion,
relative to conductance1/(2rN), is

    A=r(r+1)/2.

More generally let R be any return kernel on X_r stationary for pi (it need
not be reversible). Route x to a uniform intermediate z and then to y,
where the marked demand is pi(x)R(x,y). Every return uses at most2r native
events. Each segment's endpoint-pair flow is at most max(w)/(ZN); each
oriented edge has N/2 endpoint pairs on each segment. With full route length
at most2r, congestion relative to P_w is at most

    4*kappa*r^2.

This pays for arbitrary changing terminal controllers **on this same complete
pinned family**, retaining their stationary source demand. It does not route
a controller whose endpoint lies outside the family or whose source marginal
is silently changed. An original integer time tau remains in the record;
the compiled return adds an integer event count ell<=2r and owner-clock work
2ell. It does not equate this count with an unrelated original dynamical time.

## 4. Universal coverage by weaker pointed charts

For any finite source block, partition its actual Hamiltonian paths by

    (actual first owner, actual last owner,
     isomorphism type of partner with those two owners individually marked).

Every path has exactly one such type. On an ordinal sum the product charts
cover **all** paths. The composition lemma authenticates every coordinate
replacement within a chart, including the literal endpoint transports.
Choosing a coordinate and uniformly resampling its whole local fibre has
gap1/r on any non-singleton product chart. The same kappa comparison applies
to positive actual weights. A chart with no varying coordinate is a singleton.
Such a full-core replacement can change up to eight owner clocks; it is a
different event alphabet from the two-clock cube edges of section2.

The complete original-source census, computed by enumerating all8! paths
and all6! interior relabelings of each partner, is

| Local chart size | 1 | 2 | 3 | 4 |
| --- | ---: | ---: | ---: | ---: |
| Number of charts | 320 | 38 | 3 | 1 |
| Actual paths | 320 | 76 | 9 | 4 |

Thus there are362 charts and409 actual paths, of which89 lie in nontrivial
charts. In T_r there are362^r complete product charts covering409^r paths.
Under the actual uniform path measure, the number of varying factors is
Binomial(r,89/409), and the probability of a singleton chart is(320/409)^r.
These are exact counts for this ordinal family, not a random-tournament model.

The selected two-state cube covers only(2/409)^r of all source paths.
Even the full atlas preserves its chart label and has no cross-chart mixing.
The universal weaker coverage is therefore real, but neither universal
connectivity nor a well-founded descent follows. This pinpoints the next
obligation: authenticate interfaces between pointed charts, with their
actual source masses, instead of trying to enlarge a single convenient chart.

## 5. Reproduction and scope

`native_cube.py` checks all cube states and edges for1<=r<=6, exact Walsh
eigenfunctions and all ordered canonical routes, plus converse and four-way
resource injection. `pointed_atlas.py` saves every member of every pointed
chart and checks its canonical witness as a literal partner map. Normal and
optimized runs must agree. All-size claims follow from the proofs, not from
the finite range. Single-cell source automorphism kernels are trivial in
this construction; they are not quotiented away in the carry construction
to which this packet is compared.
