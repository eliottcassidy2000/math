# Hamilton-cycle packing closes beta(P11)=10

Status: PROVED. This is a separate correction/supplement to the earlier P11 stress packet: its previously unresolved global blocking number is now closed. The full number and classification of minimum blockers remain unresolved.

## The elementary packing bound

Let D be a directed graph with r pairwise arc-disjoint directed Hamiltonian cycles. Any arc set B whose deletion destroys every Hamiltonian path must meet each cycle in at least two arcs. Indeed, a Hamiltonian cycle with zero deleted arcs contains a spanning directed path; with exactly one deleted arc, removing that arc from the cycle itself leaves a spanning directed path. Therefore

    beta_HP(D) >= 2r.

Only arc-disjointness is used, not an assumption that the cycles are all Hamiltonian cycles or that the whole digraph is regular.

If a tournament T on n>=3 vertices has a decomposition into directed Hamiltonian cycles, the number of cycles is r=(n-1)/2 by counting arcs. Isolating any vertex deletes exactly n-1 arcs and destroys every Hamiltonian path. Consequently

    beta_HP(T)=n-1.

Moreover every minimum blocker contains EXACTLY TWO arcs from each cycle of the specified decomposition. This is necessary, not sufficient.

## Prime-order circulants

Let q be an odd prime, and let Q contain exactly one member of each pair {d,-d} of nonzero residues modulo q. Orient x->x+d for d in Q. For each d, the sequence

    0,d,2d,...,(q-1)d,0

is a directed Hamiltonian cycle because d is invertible modulo q. Different d give disjoint arc sets, and these cycles partition all arcs. Thus every prime-order circulant tournament satisfies

    beta_HP=q-1.

In particular, Paley P7 has beta=6 and Paley P11 has beta=10. This is a proof with explicit short certificates, not a large deletion-set enumeration. It also upgrades the existing finite verification of the32 circulants on Z11 and64 circulants on Z13 to instances of an all-prime theorem.

This does not invoke or prove Hamilton decomposability of every regular tournament. For composite cyclic orders a nonunit difference need not form a Hamiltonian cycle: on Z9, difference3 makes three cycles of length3. The checker rejects that false extension. Other Hamilton decompositions of composite-order regular tournaments are not ruled out.

## Consequences for the P11 packet

The121 known ten-arc blockers are now genuinely GLOBAL MINIMA, not merely inclusion-minimal examples. The110 Hall minima remain completely classified by the earlier Hall theorem; the eleven stars are certified exotic minima. The exact total of all minima is still open here, and the bare possibility1023 is not refuted by this result.

Every P11 minimum must have two deletions in each of the five difference cycles, reducing the candidate universe from binomial(55,10) to binomial(11,2)^5=503284375. This remains a large universe; it is not enumerated here. The P11 fixed-point obstruction now says directly that the minimum-blocker G-set has exactly one point fixed by the order-five multiplier, while the proposed F1024× action has three. The equivariant field extrapolation remains refuted for the sharper, fully proved minimum family.

At P7 the two-per-cycle candidate universe has21^3=9261 elements; only63 are blockers. The checker gives explicit P7 and P11 deletion masks satisfying the necessary two-per-cycle profile but still retaining Hamiltonian paths. Thus this profile is not mistaken for a full classification.

## Recovery and novelty scope

The checked repository source `HYP-9168-hp-blocking-number-equals-hall-bound.md`, copied under inputs, states the general beta=hall conjecture as OPEN, the regular beta=n-1 consequence as verified in finite ranges, and the Z11/Z13 circulants as FINITE-EXACT. THM-4553's P7 census was already known. Targeted searches of canon/results/hypotheses for Hamiltonian-cycle decomposition, cycle packing, and circulant blocking did not locate this simple structural lower-bound argument. This is a repository-level recovery and strengthening, not a literature-priority claim.

## Verification and reproduction

verify.py constructs and verifies every difference-cycle decomposition for all circulant orientations at primes3,5,7,11,13 (2,4,8,32,64 cases). For every cycle and every possible single missing cycle arc it checks the explicit surviving spanning path. This verifies the entire local certificate needed by the proof; the universal lower bound is then the pigeonhole argument above, not an unperformed brute-force check.

It also checks that all63 inherited P7 blockers and all121 known P11 blockers have exactly two arcs per cycle; supplies surviving paths for profile-only hostiles; and rejects a composite nonunit-cycle misuse. The inherited blocker lists are included as pinned data, not recomputed. Normal and optimized results are identical.

Run replay.py --self-test from any directory for exact file-set/digest verification, normal/-O replay, and file-add/content-mutation hostiles. Prior frozen packets are unchanged; this supplement explicitly supersedes their former beta(P11) uncertainty.
