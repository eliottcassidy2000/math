# Shared order-seven construction and incoming credit

The order-seven pair and four-arc route are shared with parallel repository work by **opus-2026-11-04-S76**. They are not presented here as a Session190 priority claim. The incoming commit [9036e4dcb9e1016ec42c11802ca193e6a09e6f34](https://github.com/eliottcassidy2000/math/commit/9036e4dcb9e1016ec42c11802ca193e6a09e6f34) was already on main before this checkpoint and is preserved in its ancestry.

Its [HYP-9465 preregistration](https://github.com/eliottcassidy2000/math/blob/9036e4dcb9e1016ec42c11802ca193e6a09e6f34/05-knowledge/hypotheses/HYP-9465-unit-start-plus-minus-two-port-pairs-and-the-symmetry-split-of-three-cycle-tournaments.md) distinguishes completed lower-order evidence from then-unrun order-13/14 predictions. The [exact route output](https://github.com/eliottcassidy2000/math/blob/9036e4dcb9e1016ec42c11802ca193e6a09e6f34/04-computation/experiments/plus_minus_two_port_routes_20261104.out) already records the (29,31; F=157) pair, its shortest four-arc routes, 38 all-strong shortest paths, and two below-start routes. These are incoming shared results, not claims of new priority in the present publication.

There is an exact comparison, not merely agreement of scalar counts. The single rooted relabeling

    p = (2,1,0,3,6,4,5)

maps every one of the five physical tournaments in the retained Session190 route to the incoming output's second below-start route. Here p[i] is the incoming vertex name corresponding to Session190 name i, so the Session190 root 2 maps to incoming root 0. The incoming flips are {0,1}, {1,2}, {0,3}, {3,4}. `check_incoming_credit.py` exhausts all root-compatible permutations and `INCOMING_ROUTE_MATCH.json` records the unique match. `INCOMING_SPECTRUM_SOURCES.json` retains the two exact permitted source-file reads, their Git blob IDs, and their immutable commit.

The frozen source proof and independent audit are preserved unchanged. Their scoped all-order wrapper proof, exact receipts, shortest-path analysis and fixed-parent exclusions remain the supported content of this checkpoint. This comparison neither independently reproduces the incoming higher-order census nor treats its order-13/14 predictions as completed computations. It makes no assertion about the eventual status of that separate programme.
