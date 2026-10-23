# The original electrical guard: zero kernel footprint, large path footprint

23 October 2026. **PROVED:** the transport and factorization statements
below and the exact conflict-set criterion. **FINITE-EXACT:** complete
native path and harmonic-character controls through order five; the stated
six-state order-six fixture; and four complete order-eight source fibres
with the original frozen character. **OPEN:** deriving those characters
without knowing the native network, and all-source guarded wiring.

## Inheritance and the actual carrier

The closest proved mechanism is the inverse-free marked charge rule in
[`metagraph_hall_bott.../README.md`, sections 1--3](../../metagraph_hall_bott_20261023/README.md).
The canonical hostile is its original order-eight switch at merged nodes
`1238,2242,2246,2315`. The corrected near miss is replacing this guard by
the order of a ternary kernel, a Hamiltonian count residue, or a clock
residue. The least-used sidecar is the actual target owner of each
Hamiltonian presentation. The concept board is: frozen character,
transported character, actual target owner, indexed path footprint,
shared owner pins, and aggregate versus per-port conservation.

The electrical vertices are **converse-merged tournament classes**, not
the vertices inside one tournament. The full literal complement-line
network is `G=B+2Q`; loops retain two endpoint stubs but contribute zero
to its Laplacian `L_G`. This is the carrier of
[`metagraph_network_invariants.../PROOF.md`](../../metagraph_network_invariants_20261021/PROOF.md),
not an auxiliary edge-cube network. The old Session38 dyadic tower acts on
this integer Laplacian; it does not provide a theorem that an arbitrary
small port record determines its harmonic characters.

Fix a root and a known native potential

    x in ker(L_G modulo 8),    x_root=0.

For an actual directed Hamiltonian path `P` in `T`, let `K_P(T)` keep
its path arcs and reverse every other arc. Write `[T]` for the merged
owner and define the actual directed port phase

    phi_x(T,P) = x_[T] - x_[K_P(T)] modulo 8.       (1)

The original necessary guard on candidate wiring is, at every merged node,

    sum_(literal outgoing stubs p) phi_x(p)=0 mod8. (2)

Equivalently `L_candidate x=0 mod8`. Actual edge multiplicities, source
ownership, target ownership and loop stubs remain part of the input.
Equation (2) is not a test of one arbitrarily colored vertex inside `T`.
The phase can be evaluated only after authenticating the target owner or
an equivalent sufficient target observer.

## 1. The genuine kernel-lift guard has footprint zero

**Theorem.** For every actual source `T`, Hamiltonian path `P`, and
`g in Aut(T)`,

    K_(gP)(T)=g K_P(T),
    phi_x(T,gP)=phi_x(T,P).                        (3)

The first identity follows arc by arc: `g` carries precisely the actual
consecutive arcs of `P` to those of `gP`, and carries their complements
to complements. Thus both merged owners in (1) remain the same. The
electrical potential is indexed by those owners, so its phase is unchanged.
No invariance of the chosen coordinate values under unrelated permutations
of merged-node labels is needed.

Consequently the real phase (1), and any guard depending only on these
fixed endpoint phases and their owned sums, have **zero exposed free
coordinates** on the native automorphism-orbit compiler of
[`hamiltonian_observer_restrictions.../compiler`](../../hamiltonian_observer_restrictions_20261023/compiler/README.md).
That compiler's Hamming lift moves preserve this actual guard automatically.
This discharges its formerly unspecified electrical-guard obligation on
that precise orbit. It does not preserve additional conditions on literal
vertex labels unless those conditions are also pinned or transported.

Converse plus reversal gives
`K_(reverse P)(T^op)=K_P(T)^op`, so (1) is unchanged in the merged carrier.
Exchanging the two ends of the complement line, by applying `K_P` again,
negates the phase. These are different operations.

This proof does not identify all Hamiltonian paths of `T` with one
automorphism orbit. Across different orbits, the target owner and the
phase may change. The finite examples below show that change sharply.

## 2. Transporting a character is different from freezing its coordinates

For a permutation matrix `R` on the **merged network**,

    L'=R L R^t,  x'=R x  =>  L'x'=R Lx.

Hence transporting both the network and its marked character preserves
harmonicity and transported port phases. Freezing `x` while relabeling
only the network instead tests `R L R^t x`, which need not vanish.

The checker exhibits this distinction on the original order-eight network.
Swapping merged labels `1238` and `2242`, while freezing the original
character, gives **217 nonzero node defects modulo eight**. Transporting
the character with the same swap restores all equations. The graph is
isomorphic throughout; the frozen sidecar is a genuinely marked datum.

The original primitive lift is modulo 256 on 3,528 merged nodes, with
root 3,527. Its four selected entries are `(97,18,65,212)`, hence the
mod-eight residues are `(1,2,1,4)`. The checker replays every edge of the
original 438,216-pair export and verifies `L_G x=0 mod256`. For the
published switch `L_t=L+2t(uv^t+vu^t)`, its exact mod-eight defect at
those four nodes is

    t*(2,6,6,2).

Thus among the inherited feasible amplitudes `t=-1,0,1,2`, only zero
passes. The new work consumes this exact original character and gauge;
it does not replace them with an abstract harmonic example.

## 3. Exact fixed-clock footprints

Fix an actual labelled source tournament and its finite set of all
Hamiltonian paths. For a function `F(P)` and each pair with different
values, retain the disagreement set

    D(P,Q)={j : P[j] != Q[j]}.

**Theorem.** A set `S` of indexed owner clocks determines `F` iff it
intersects every `D(P,Q)` with `F(P)!=F(Q)`.

Indeed a missed disagreement set gives two agreeing observation records
with different values. Conversely, meeting all such sets ensures every
equal observation record has one value. This is a complete finite
criterion, not merely a sufficient feature score. An indexed arc flag
at clock `j` exposes the two real owners at `j,j+1`; use that union of
clocks in the same criterion.

The receipt retains every minimum clock set, every minimum arc-flag set,
and an actual conflicting path pair for **every** clock set one smaller
than the optimum. Upper bounds are checked independently by grouping the
actual paths by their observed owner tuples.

There is a useful all-size upper mechanism. If exactly two positions are
unmarked and their distance is at most two, the full Hamiltonian path is
determined. At distance one, interchanging the two remaining vertices
would require both arc orientations. At distance two, their common marked
middle vertex would require both orientations of an incident arc. Thus
the `2n-3` pairs of holes at distance one or two give sufficient banks of
`n-2` marked owner clocks. This is not asserted minimal for every source
or every observable.

It also gives an all-size upper bound of `floor((n-1)/2)` actual arc flags
for identifying the whole path when `n>=2`: for even `n`, mark the first
`n-2` positions by disjoint consecutive flags, leaving two adjacent holes;
for odd `n`, mark the first `n-1`, leaving one owner determined by the
remaining vertex. This is a sufficient bank, not a universal minimality
claim.

## 4. Original order-eight phase really needs six clocks in four source fibres

The four sources are the actual nodes of the old switch. Their full
Hamiltonian fibres are enumerated from their pinned canonical matrices,
and each normalized mask is checked against the original owner map.
All eight residues occur in every row:

| Merged source | Ordinary class | Actual Hamiltonian paths | Distinct target owners | Minimum owner clocks for original phase | Minimum actual arc flags |
|---:|---:|---:|---:|---:|---:|
|1238|1763|409|225|6|3|
|2242|3481|407|204|6|3|
|2246|3490|523|313|6|3|
|2315|3676|527|323|6|3|

For each source, exactly thirteen six-clock sets work. They are precisely
the complements of two holes at distance one or two, matching the preceding
structural upper bound. Each of the 56 five-clock sets has a saved literal
counterexample. The four minimum arc-flag banks are

    {0,2,4}, {0,2,6}, {0,4,6}, {2,4,6}.

The full merged target-owner observer has the same six-clock minimum on
these sources. This is not a claim that the three-bit phase identifies
that target owner; each phase value merges many owners. Rather, compressing
an output to three bits need not make it depend on fewer input positions.
These are nonadaptive minima in the declared indexed-clock menu.

For an explicit endpoint hostile at source 1238, the actual paths

    (0,4,2,1,6,7,5,3),   (0,4,1,5,2,6,7,3)

have the same first and last owners but phases two and three. Endpoints
alone do not determine even this one electrical character.

At smaller native orders, the checker reconstructs a complete grounded
mod-eight harmonic bank independently by binary lifting, and also generates
it from exact inverse-column current certificates. The minimum fixed-clock
cost distribution over all ordinary source classes is:

| Native order | Cost zero | Cost one | Cost two |
|---:|---:|---:|---:|
|3|2|0|0|
|4|3|1|0|
|5|4|2|6|

The order-five bank records every mod-eight character, not just a selected
one. The first same-endpoint phase-bank hostiles occur at ordinary classes
2 and 6, a converse pair. All these are full actual path enumerations;
the frozen owner maps remain explicit inherited inputs.

## 5. A small actual electrical split with a retained shared phase

Take two directed triangles on labels `3i+u`, with internal arc
`u->u+1 mod3`. Direct the first cell to the second iff `v-u` is one or
two modulo three. The six-vertex source is **not regular**. For phases
`a,b`, concatenate the two cyclic orders

    P(a,b)=(a,a+1,a+2, 3+b,3+b+1,3+b+2),

with each cell coordinate read modulo three. This is Hamiltonian exactly
when `delta=b-a` is zero or one. The actual source is old atlas ordinary
class 33, merged node 25. The native character is

    (0,4,0,0,0,0,4,4,4,4,4,4,4,0,0,0,0,
     0,4,4,4,0,4,4,0,0,0,4,4,4,4,0,0,4).

It is on the full 34-node native order-six network, rooted at zero.
The receipt supplies both a short exact certificate `Lx=8j` and the
integer inverse-column lift from which it was obtained. Independent
literal target isomorphisms are retained, and all 720 permutations verify
the two targets are neither isomorphic nor converse.

| Delta | Actual path valid? | Old backward-tile mask | Merged target | Original phase |
|---:|---:|---:|---:|---:|
|0|yes|684|32|0|
|1|yes|841|11|4|
|2|no|undefined|undefined|undefined|

The full nine-state table retains validity separately from `valid*phase`.
An invalid record with a zero placeholder is not a phase-zero port.
The source's complete automorphism group has order three: common rotation
of the two cells. Thus each phase sector is one three-element common-phase
orbit. The electrical split is between two different target-owner sectors,
not between hidden lifts of the same sector.

### The original order-eight character also sees the pinned extension

Adjoin a source vertex zero and sink vertex seven around that six-vertex
core, whose vertices now have labels one through six. The paths start at
zero and end at seven, so these two external owner-clock pins are real.
Using the **original** order-eight character of section 2 gives:

| Delta | Native merged source | Native merged target | Original forward-bit mask | Phase modulo eight |
|---:|---:|---:|---:|---:|
|0|3517|532|1765183|1|
|1|3517|1290|1799615|5|
|2|undefined|undefined|undefined|undefined|

This is a literal padded tournament and its actual Hamiltonian paths,
not a claim that harmonic potentials extend naturally under padding.
The new phases were read from the old full order-eight character. The
complete adjacency matrix, path words and actual pins are saved.

## 6. Two observer slots: aggregate conservation permits a sharp paired move

On the same six-vertex source, take two labelled observer slots with paths
`P(a,b)` and `P(a,c)`. Fix `a=0`; their shared first-cell owners at clocks
zero, one and two are therefore exactly `(0,1,2)`. Require the target
multiset `{32,11}`. There are exactly two assignments:

    (b,c)=(0,1), (1,0).

Both retain the actual shared owners, target incidence multiset, phase
multiset `{0,4}`, and hence the source's aggregate electrical contribution.
The simultaneous swap changes two phase coordinates. Every one-coordinate
intermediate either is non-Hamiltonian or has the wrong target multiset.
Budget two is therefore sharp **for this labelled assignment menu**.

The pinned order-eight extension has the same control, with shared clocks
`0,1,2,3,7`, targets `{532,1290}`, and original phases `{1,5}`. Every clock
and path endpoint is explicit in the receipt.

If per-slot phases, per-slot targets, or the second-cell owners are
individually pinned, this swap is forbidden. Moreover its unlabelled
multiset of actual paths and complement lines is unchanged. It is a
coupled repair of labelled observer-slot assignments, not a switch that
changes the native merged wiring. The distinction prevents aggregate
conservation from silently replacing individual port authentication.

## Reproduction, gauges, and the remaining obligation

    python -B 05-knowledge/results/native_guard_moves_farey_20261023/guards/check.py
    python -O -B 05-knowledge/results/native_guard_moves_farey_20261023/guards/check.py

The small atlas uses its original triangular tile order with a one bit
meaning a backward arc. The order-eight producer uses lexicographic
nonconsecutive chords with a one bit meaning a forward arc. Its canonical
full tournament words use a third, most-significant-bit-first pair order.
The checker implements each gauge separately and pins input hashes. Their
class IDs and masks are never interchanged.

The decisive unresolved obligation is now precise: derive a useful native
character or its needed target values from incomplete port data **without
already knowing the target owner/network**, and prove that a proposed
move preserves the actual retained per-port or aggregate guard. The
automorphism-orbit case is closed by (3); the original order-eight fibres
show why that case does not settle arbitrary native path changes.
