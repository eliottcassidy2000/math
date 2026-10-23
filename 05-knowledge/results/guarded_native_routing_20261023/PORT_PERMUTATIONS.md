# Distinct observers: permutation schedules on one complete accepted fibre

23 October 2026. **PROVED scoped.** Let X be a finite nonempty accepted
observer family of ONE fixed actual source. All its states share every
required indexed pin and guard. Suppose a connected undirected graph of
authenticated native moves on X has diameter at most D. These hypotheses
hold for the count-residue carrier in PROOF.md with D as stated there.
They are not inferred for an arbitrary electrical filter.

Put one distinct token, labelled by its original actual observer, at each
state of X. A time slot exchanges the tokens at the endpoints of a native
edge, leaving all other tokens in place. Both participating observers then
make a legal native move in opposite directions. Every slot is a bijection
of the occupied actual states. The two endpoint records and all fixed pins
are retained, rather than replacing two tokens by two copies of one path.

## 1. Any uniform-stationary demand has a distinct-terminal coupling

Let R be a Markov kernel on X with uniform stationary law, N=|X|. Its row
and column sums are one. Then R is a convex combination of permutation
matrices. Here is the finite argument, including the relevant capacity test.

For a nonnegative residual matrix with common positive row and column sum
t, let S be any collection of rows and Gamma(S) the columns meeting its
positive entries. Its total mass |S|t is supported in columns of total
capacity |Gamma(S)|t. Thus |Gamma(S)|>=|S|. Hall's matching argument gives
a permutation of positive entries. Subtract their minimum epsilon from
the corresponding permutation matrix. Every row and column sum becomes
t-epsilon, and at least one positive entry disappears. Iterate. There
are at most N^2 iterations; the coefficients sum to one. Rational input
gives exact rational coefficients without a limit argument.

Sampling a whole permutation from this decomposition gives each initial
observer x exactly its demanded terminal distribution R(x,.), while every
actual terminal is used exactly once on every trial. Independent terminal
sampling would not give this stronger ownership property. This is the
classical doubly-stochastic decomposition, with a proof included here;
no novelty claim is made for it.

## 2. Implement the sampled permutation by legal native exchanges

Any permutation of N states is a product of at most N-1 arbitrary pair
transpositions. For endpoints u,v choose a simple native path

    u=v_0, v_1, ..., v_d=v,       d<=D.

Exchange consecutive positions forward along all d edges, then backward
along the first d-1 edges. This uses 2d-1 slots. It exchanges the two
endpoint tokens and restores every intermediate token. Consequently every
sampled permutation has a collision-free sequential schedule of at most

    (N-1)(2D-1) slots,
    2(N-1)(2D-1) individual non-idle native events.       (1)

Every intervening state is accepted. Keep the complete integer slot ledger
and, for each token, its individual move/idle ledger. A slot can be padded
by idling all tokens, but it cannot be silently identified with one old
ambient excursion. For the ternary carrier, D<=4r, so the bounds are
(N-1)(8r-1) slots and twice that many individual events. Each event still
changes at most two cell phases. The dependence on N is exponential in r
in the large families; (1) is not a replacement for the polynomial
single-demand congestion theorem.

The construction is an existence and exact finite compilation result. It
does not assume that a dense R or its permutation decomposition can be
processed in time polynomial in the short word length r. The transfer
sampler in PROOF.md avoids listing X for ONE demand; it does not eliminate
this global schedule's output size.

## 3. What this does and does not repair

This resolves distinctness for unit-capacity permutations of the COMPLETE
accepted family, subject to authenticated moves. It leaves the multiset
of occupied actual states unchanged. Hence it preserves any aggregate
observable depending only on that multiset; it does not create a new
source/target incidence multiset or a new metagraph edge.

For a selected subset of ports, different terminal capacities, or prescribed
cross-source rectangles, the equal-margin hypothesis must be checked on
the actual owned states first. The incoming amplitude-two rectangle fails
before scheduling: one needed target has only one available actual path.
A decomposition cannot manufacture another. Likewise, a full family
permutation cannot repair its nonzero higher-dyadic target current.

The subsequently received [native coupled-packet theorem](../session44_guarded_native_moves_20261023/coupled-moves/THEOREM.md)
gives an especially useful boundary: three slots with distinct pinned
eligibilities {p0,p1}, {p1,p2}, {p0,p2} have exactly two injective assignments,
which differ in all three slots. No two-slot exchange preserves eligibility.
Our pair-exchange compiler requires a COMMON accepted family for all tokens.
It does not contradict or bypass that authenticated three-slot obstruction.
The incoming alternating-cycle compiler is the applicable general mechanism
when slot eligibility differs; our new contribution there is no all-size
small-cycle bound.

The weighted routing statement in PROOF.md keeps arbitrary positive
weights. It does NOT automatically give one-token-per-state permutation
couplings: a nonuniform stationary kernel need not be doubly stochastic.
Rational weights may suggest a capacity lift, but the extra copies need
authenticated actual owners and native lifted moves. They cannot be
invented as duplicate witnesses. No such general weighted lift is claimed.

The next useful strengthening is a polynomial-depth schedule of disjoint
native exchanges, with an exact capacity/lift description. Short routes,
finite-state guard encodings and the present fractional congestion estimate
alone do not supply that schedule.

## 4. Exact controls

`ports.py` constructs the native bit-move graph independently of the route
producer. It tests all 120 permutations of the five-state r=4 ternary
endpoint fibre and selected permutations of the 22-state r=6 fibre. Each
intermediate occupancy is a permutation, each actual exchange is a native
edge, and the entire terminal permutation and integer slot bound are
checked. A nonreversible doubly-stochastic kernel is decomposed by exact
augmenting matchings; reconstruction of every matrix entry is checked.
A two-token/one-target demand is the negative capacity control. These are
finite checks of the proof, not an enumeration of all large permutations.
