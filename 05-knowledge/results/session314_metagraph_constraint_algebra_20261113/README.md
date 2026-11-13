# Session 314: Exact recurrences and obstructions for the merged tournament metagraph

The central distinction is between the **actual merged complement-line multigraph**, its **blue score quotient**, the **labelled cycle-space shadow**, and the auxiliary **single-arc-flip graph**. Results in one are not assumed to carry to the others.

## The actual graph

- [Alignment recurrence](ALIGNMENT_RECURRENCE.md): given the two class representatives, a line is an alignment whose disagreement arcs form one spanning directed path. The degree guard and the determinant cycle guard give a clean exact recurrence for every colored incidence and loop. Keep the witnessing alignment orbits and the construction retains the actual half-edges, not just their totals. The two independent enumerators agree through order five.
- [Cycle carriers](CYCLE_CARRIERS.md): closed ends are grouped by directed Hamilton-cycle orbits. The exact divisor polynomial tracks seam period and reflection type, with one blue seam at odd order and two at even order for a reflective carrier. Successive seam choices give open partners differing by just the two seam arcs.
- [Line flux](LINE_FLUX.md): product formulas give the joint triangle-level generating function, the exact blue and black flux distributions, and the open/closed covariance. These are constraints on the full graph, not replacements for its edges.
- [Marked stabilization](HISTORICAL_BRIDGES.md): symmetric source/sink padding stabilizes the actual bipartite incidence kernel only after keeping both endpoint names. A blue self-line already splits at order four if those marks are recovered. The [radial stabilization](RADIAL_STABILIZATION.md) theorem concerns the auxiliary single-flip graph instead, where the padded copy is geodesically convex.

## The blue score quotient is classified

[Complete profile support](BLUE_PROFILE_SUPPORT.md) proves, at every order, that distinct feasible score profiles are adjacent exactly when they differ by one Young-lattice box and at least one is strong. At even order same-profile lines occur exactly at strong profiles containing the two central scores. [The cut and rounding theorem](BLUE_CUTS.md) supplies the physical lift from a sorted-interior score word, and the corner recurrence computes the exact unweighted support generating functions. [The Laurent product](OBSTRUCTIONS.md) supplies the actual quotient multiplicities.

The same Laurent product is a type-$B$ positive-root product with the simple roots removed. Inserting the center at $4\to5$ or $8\to9$ adds an exact centered spin-weight factor, not a claim of Bott periodicity. The finite order-19 control had 25,826 odd profiles and 84,222 vertical support edges; the proof does not depend on that ceiling.

## Three roles for the small graphs and the number 19

[Shape obstructions and witnesses](OBSTRUCTIONS.md) distinguish the actual support from its shadows. The first combined-support $K_5$ and $K_{3,3}$ occur at order five; the first Petersen subgraph occurs at order six, and the first induced Petersen at seven. The ten-node order-five support has no nontrivial $A_5$ action, although the ten labelled triangle vectors carry an actual $A_5$-equivariant Petersen relation. The [historical bridges](HISTORICAL_BRIDGES.md) also give the six-point Sylow model: ten $K_{3,3}$ frames have identical pairwise overlap, but their order-three generators recover Petersen adjacency. The [orientation lift](ORIENTATION_LIFT.md) exhibits a literal $A_5$ label action and a period-four effect on the lifted walk determinants. None is an eightfold periodicity of the full metagraph.

The historical note traces the $4,5,8,9$ square-sum forced corridor, the short $Q_{19}$ obstruction, Paley order 19, the old Collatz 19-adic observers, and the Session 313 physical 19-block. Neutral padding there gives one fixed-cut family reaching every unit modulo every power of 19 by length alone. This remains modular, not absolute Hamiltonian-count coverage.

**Boundary.** The class-level incidence is exactly computable from representatives, and the blue score quotient has a complete support theorem and generating functions. A compact recurrence that generates the full class-level wiring without constructing representatives or retaining path information is not established. No Chan–Pak proof or metagraph periodicity is claimed.
