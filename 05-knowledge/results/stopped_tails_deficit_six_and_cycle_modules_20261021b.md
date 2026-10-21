# Stopped Collatz tails, deficit-six centers, and cycle-count modules

Research checkpoint, 2026-10-21. **PROVED:** a monotone all-future energy
envelope, an exact four-point flow model, a quadratic center cutoff, and an
all-order change of enumeration state. **FINITE-EXACT:** the declared native
scales, order-ten frames, and coefficient comparisons. Uniform all-N 153N,
arbitrary-source Collatz entry, and full tournament-metagraph classification
remain **OPEN**.

## What changed

The most useful change is that a long-horizon certificate need not finish
every source path. Advance the sources to a finite native depth, retain each
unresolved current node and its separate old/new charge vectors, and bound
every possible later placement of those blocks. The resulting upper bound
can only decrease as the actual prefixes are extended. It succeeds at three
declared scales, including the new scale N=32768, while 1401 source paths
are still unresolved. This is a sufficient stopping certificate; its eventual
success at every scale is not proved.

The tournament work reduces the maximum possible order of a strong
self-converse deficit-six center from4096 to248. A complete necessary
128-frame order-ten universe has no such center. Consequently the unknown
center polynomial has known beginning `z^5+z^6+z^7+O(z^12)` and91 possible
remaining orders. Separately, a cycle-count polynomial changes the exact
Burnside state space and substantially accelerates the tested signed-core
enumeration. For ordinary A000568 the incoming scalar algorithm beats the
full polynomial method and is essentially tied with the new hybrid tail.

## Inheritance, portfolio, and connection board

Anchor: long-horizon source-owned coverage. Niche: deficit-six centers.
Wildcard: exact normalized enumeration. The repository's LRC(14) priority
is retained; no unsupported transfer to it is claimed.

The session began from the [previous scale checkpoint](scale_credit_and_tournament_compression_20261021a.md),
integrated the [recovered cross-scale credit proof and receipt](cross_scale_frontier_credits_20261021/README.md),
and read the incoming [degree-local enumeration](degree_local_tournament_enumeration_20261021/README.md).
The last gives a stronger performance comparator, not a reason to claim the
new ordinary-count representation is faster. Explicitly reconstructed
Session21 documents supplied questions and formulas; their missing old
receipts are not treated as fresh evidence.

Before close-out, incoming commit `777c7240f7` supplied the independently
reviewed [fixed-degree rationality theorem](fixed_degree_metagraph_rationality_20261021/README.md).
Its exact connection to the new cutoff is recorded below. Both incoming
commits were integrated from a clean tree.

Closest proved mechanisms: synchronous one-edge cancellation, source-owned
interval banks, regular rotation-orbit connectivity, half-supported path
cliques, and future-gcd-response enumeration. Hostiles: the native N7 cut8
failure; opposite charges at different current nodes; the pure-blue C5;
and histories with equal residual responses but different actual cycle count.
Corrected near misses: optimizing a zero packet does not reduce physical
energy, a smaller state-key count need not reduce runtime, and a cycle
without a legal common suffix does not certify black path orbits.
Least-used sidecars: unresolved current nodes, original clock words,
the literal common suffix, and the cycle-count polynomial.

The live board is: source ownership / unresolved continuation / four-point
flow defect / reflected prefix support / normalized cycle count / original
integer clocks. The productive common operation is to remove a contribution
whose future action is known, while retaining a coordinate that reconstructs
the observable. This is a mechanism comparison, not an isomorphism between
Collatz dynamics and tournament enumeration.

## 1. A monotone envelope that does not assume future arrival

The full proof, implementation, receipt, hostile controls and audit are in
[future_tail_blocks_20261021b](future_tail_blocks_20261021b/README.md).
At a finite cut D, let `c_t=(W_t,H_t)` be the arrived old/new vectors at
terminal t and `b=(w,h)` a pending common-current-node block. Put

\[
\delta(w,h)=Q(h)+\beta(w,h),\qquad
\kappa(b,b')=\beta(h,h')+\beta(w,h')+\beta(w',h),
\]
\[
L_b=\max\{0,\delta(b)+\max_t|\kappa(c_t,b)|\},\qquad
\mathcal U_D=\Delta_D+\sum_bL_b+\sum_{b<b'}|\kappa(b,b')|.
\]

The terminal maximum includes zero. Every pending block may later be absent
or arrive with either sign at any terminal. Bounding this enlarged choice
space proves `Delta_h<=U_D` for every h>=D. Signed transport, merging,
arrival and both-channel-zero erasure each preserve or decrease U. Every
cross-block term appears once. This proof uses no hypothetical terminal or
ROOT suffix for an unresolved source.

| N | stopped depth | pending nodes | unresolved nonzero-charge sources | U_D | exact153N budget |
|---:|---:|---:|---:|---:|---:|
|4096|105|195|373|557164|626688|
|8192|110|348|689|1138407|1253376|
|32768|130|760|1401|4241482|5013504|

All earlier cuts also pass. Thus these are **all finite horizon** certificates
for these three scales, not merely checks through the displayed depths.
The first two scales had already been proved by completing every trajectory;
the new certificate avoids that requirement. N32768 extends the supplied
scales beyond the previous contiguous4096..16384 census.

For any authenticated ownership partition,

\[
\Delta_h=A_h+I_h,\quad A_h+M_h=153N,\quad I_h=R'_h+J'_h
\]

gives `R'_h<=U_D-A_h-J'_h`. Once U_D<=153N this is precisely
`R'_h<=M_h-J'_h`. It does not assert the stronger `R'<=M`, or deposit an
old cancellation twice. All source addresses, mod15/mod35 and performed
integer clocks remain decodable. The observable is native-depth Nyquist
energy; no all-frequency or all-integer-clock-cut theorem is inferred.

**Remaining obligation:** prove a successful cutoff D(N), with the earlier
cuts controlled, for every N>=4096. A scale lift of a selected zero block
does not cover the other sources in the enlarged annulus. The prior theorem
`N>=2^(8d+12)` controls fixed early horizons; the new stopping predicate
controls a whole tail after success. Their quantifiers do not yet join into
a uniform all-N proof.

## 2. Larger zero blocks via a four-point flow, with no recycled gain

[long_horizon_pairs_20261021b](long_horizon_pairs_20261021b/README.md)
identifies the twelve native charge types with the twelve directed differences
between `(0,0),(0,1),(-1,-1),(1,0)`. A selected integer multiset sums to zero
exactly when its four-vertex flow divergence is `t*(-3,1,1,1)` for some
integer t. This includes a defect coordinate absent from ordinary circulation.
The carrier is a bidirected graph, not a tournament with artificially broken
ties. Capacities and source identities survive the map.

Network cuts characterize feasible t. Integral residual potentials certify
the maximum selected norm at fixed t, and concavity plus neighboring optima
certifies the global integer optimum. The irreducible quartet
`{(-1,-1),(1,-1),(1,2),(-1,0)}` demonstrates a zero packet which no pair or
triple decomposition can extract. At N4096, terminal711, cut47, an exact
3747-source packet cancels norm6317 of6378; antipodal pairing cancels5930.
All selected sources have actually arrived at that cut.

This is stronger information compression, not an additional energy saving.
The residual vector is unchanged. Writing `I=Gamma-Z`, increasing the removed
zero norm Z increases Gamma by the same amount. Applying the optimizer within
a pending node likewise leaves its aggregate `(w,h)` and U unchanged. Also,
optimal packets need not grow by inclusion as new sources arrive: the receipt
has a minimal temporal example. Ownership must remain cut-specific, or a
newly chosen bank must replace the old one rather than augment it.

## 3. Deficit-six centers and a finite recurrence

[deficit_six_centers_20261021b](deficit_six_centers_20261021b/README.md)
proves the following all-order bound. For a strong self-converse tournament
with black path-orbit deficit Delta<=d, d even, take the least power of two
r>=d+4. Then

\[
n\le\max\{3r-1,(d+2)(2r-1)\}<4(d+2)(d+4).
\]

The improvement averages exits of an r-cycle inside one regular rotation
orbit, where the average is exactly r/2. Inherited connectivity provides a
common actual Hamilton suffix. Relative even-order permutations and disjoint
reflected supports force each selected black orbit to have a distinct mate
outside the selected list. This is the factor-two deficit bound; merely
counting the selected words loses it.

At d=6 the cutoff is248. The complete normalized two-shore order-ten census
has128 frames,64 self-converse, and none of deficit six. A nonautomorphic
Hamilton-cycle successor also excludes order11. The three known centers
retain every literal black path port and target; their orders are5,6,7.
Ninety-one candidate orders from12 through248 remain unclassified.

With `p=1-z-z^3`, the trivalent vertex series has exactly reduced denominator
`p(z)^2 p(z^2)^2`, hence minimal eventual recurrence order18. The center cutoff
gives a safe recurrence onset n>=261. The reduction is independent of the
unknown center coefficients; it does not supply those coefficients.

The incoming rationality theorem extends the counting framework to every
fixed fundamental degree k. Our new cutoff applies to all its self-converse
center polynomials: `x(x*a-b)=2k` implies `0<=a-b<=2k`. Consequently their
degrees are at most `max(3r-1,(2k+2)(2r-1))`, with r the least power of two
at least2k+4. This is a quadratic finite-center bound; its denominator law
is unchanged. Ordinary non-self-converse atoms retain their separate inherited
bound. Neither result recovers adjacency from the counting series.

## 4. Faster signed-core enumeration by retaining a polynomial sidecar

[tournament_tail_kernels_20261021b](tournament_tail_kernels_20261021b/README.md)
removes `binom(t,2)` from Burnside energy, uses future responses `S_j=R_j-t`
as the structural key, and retains actual cycle count t in a polynomial
payload. Inserting k cycles of odd length ell multiplies the polynomial by

\[
\frac{2^{(\ell-1)k^2/2+kS_\ell}}{\ell^k k!}\,y^k.
\]

The final linear functional `y^t -> 2^binom(t,2)` restores the removed energy.
It is not a ring homomorphism. In the signed-core calculation it can be
applied before distributing terminal excess, so the same calculation is
not repeated for every excess degree. Full factorial normalization supplies
integer implementations; the earlier odd denominator is not silently reused
after the normalization changes.

Ordinary bound200 structural transitions fall from2994814 to1515929, but the
incoming degree-local scalar routine is faster in the matched timings.
Signed-core bound100 transitions fall from1697975 to289706, with matched
median runtime falling from0.8517s to0.2486s. Ordinary incoming2.1967s versus
hybrid2.1777s is treated as tied, without an ordinary speedup claim.
The packet retains paired timing samples,
coefficient-operation counts, memory controls and exact output hashes. Counts
do not reconstruct tournament sources, automorphisms or metagraph edges.

## Next decisive targets

1. Prove a scale-uniform bound on the stopped envelope's arrived, local and
   pending-pair terms, using actual synchronous fibres. The missing statement
   concerns every annulus source, not just an infinite lifted subset.
2. Classify or exclude the remaining cyclic deficit-six frames, beginning at
   order12, by legal prefix certificates or full orbit deficits. Keep literal
   complement targets when asking about black-core connectivity.
3. Reduce cycle-polynomial coefficient work while retaining the final linear
   functional. Continue comparing against the incoming scalar algorithm;
   use the polynomial tail where its signed-excess savings actually pay.

The packets provide their proof boundaries, reproducible commands, immutable
receipts and independent proof/API audits. Successful finite stopping,
lossless observable compression, and complete classification are different
claims throughout.
