# Independent audit: one action, common versus individual fixed spaces

2026-10-25. **PASS / PROVED.** Reviewed `COMMON_ACTION_BRIDGE.md` and
`check_common_action.py`. The bridge transports the same permutation
matrices into two constructions; it does not identify their target
obstructions or make the added graph intrinsic to the tournament.

For any field `K` and finite permutation action on `n` points, a vector is
fixed by the whole group iff it is constant on each whole-group orbit.
Thus `dim (K^G)^A=r`. For one element, the analogous dimension is its
number of cycles, in every characteristic. In particular the clone
quantity `d=min_g dim V^g` and the graph quantity `r=dim V^A` are different
fixed-space questions; semisimplicity is unnecessary for this statement.

For a connected graph with cycle rank beta and this specified holonomy,
the space of global flat sections has dimension r. The graph cochain
complex has dimensions nv and ne. Therefore

```
dim H^1(Gamma;V) = n(beta-1)+r.
```

Over `F2`, when n is odd, augmentation splits equivariantly by the constant
all-ones vector. Hence, more generally than the regular example,

```
dim E^A = r-1,
dim H^1(Gamma;E) = (n-1)(beta-1)+r-1.
```

The ordinary subset cancellation theorem still requires the odd-action
hypothesis. No claim about cancellation over coefficient fields of positive
characteristic is implicit in the characteristic-free fixed-space counts.

For the regular action `C3^k`, a bouquet of k standard generator loops has
`beta=k`, and subdivision into triangles gives the advertised simple graph
without changing holonomy or cycle rank. Transitivity gives `r=1`; each
nonidentity translation has `3^(k-1)` cycles. Thus

```
d=3^(k-1),
dim H^1(V)=3^k(k-1)+1,
dim H^1(E)=(3^k-1)(k-1).
```

The independent audit explicitly enumerates all **512** nine-bit vertex
assignments at k=2. Their two-generator coboundaries have exactly 256
values, each with two preimages; restricting to the 256 even assignments
is injective with the same image. This yields `H^1(V)=18-8=10` and
`H^1(E)=16-8=8` without repeating the candidate's matrix-rank code. The
same nine-point action independently has `r=1,d=3`.

The candidate's direct binary coboundary calculations at k=2,3,4 were also
replayed in ordinary and optimized mode, with identical outputs. The
independent suite's corresponding evidence is under `common_action_bridge`.
The finite checks support the rank-nullity proof, not an all-k extrapolation.

One qualification is essential: the three-pair holonomy table also includes
`C2` and `S3`, which are not odd groups. Its local-system formula remains
valid there; it does **not** inherit the odd-action subset-root formula.
The bridge text correctly keeps these roles separate. The full three-pair
table follows from `f=r-1`, giving `2beta`, `2beta-1`, and `2beta-2`.

The added graph, coefficient field, and cycle generators remain mandatory
sidecars. No plane-coloring, tournament-move or counting-complexity conclusion
is certified.
