# Physically implementing optimal port marks

Use Session 311's cyclic ternary tower of depth d≥1 with N=3^d leaves, each a copy of the same rigid eight-vertex tournament B. Marking leaves means pointwise fixation, not merely setwise preservation. In each copy use the same designated vertex via the unique identification with B. Let Γ_h be the full recursive rotation group of a depth-h uncolored tree and E_d=(3^d−1)/2.

## Two probes lift every mark set

Add two vertices x→y. Normally x beats all tower vertices and all tower vertices beat y. For a selected leaf choose either red (reverse the x edge to its designated vertex) or blue (reverse the y edge to its designated vertex). All other edges remain unchanged. For any prescribed set S of selected leaves there exists a choice of red/blue with automorphism group isomorphic to the pointwise stabilizer of S in Γ_d. In fact its restriction to the old tower is exactly that subgroup.

The probes are intrinsically distinguished even without names. Every core vertex has outdegree between (8N−8)/2 and (8N−8)/2+9 after attachment: it beats exactly half the vertices outside its own leaf, and at most 7 inside its leaf plus at most 2 probes. Meanwhile x has outdegree at least 7N+1, and y has outdegree at most N. Since N≥3, both probes are fixed by all automorphisms. The full automorphism group is therefore precisely the color-preserving recursive rotation group on the leaves; their interiors are rigid. Reversing the incidences at x or y gives distinct colors.

Now color bottom-up. Each occupied subtree admits a pair of red/blue colorings exchanged by global red–blue complement, nonisomorphic, and having no rotations at any occupied node. For an occupied leaf choose the red and blue alternatives. At an internal node with one or two occupied children, the child tuple is automatically nonconstant since an unoccupied child is uncolored. If all three are occupied, choose an alternative in one child to avoid a constant triple. The resulting parent has no top rotation, and its global complement is distinct as a class: any rotation that matched the child tuple to its complement would, on iterating three times, force every child to be complement-fixed. At least one child is occupied, a contradiction by induction. Thus precisely the rotations at occupied ancestors are lost, just as for pointwise marks. All unoccupied subtrees keep the full group. Conversely any physical coloring with this touched support retains the groups on unoccupied subtrees, whether or not extra symmetries survive.

This physically realizes the entire sharp stabilizer exponent interval from the marking budget notes by adding only two vertices and exactly |S| selective arcs. If a desired count c has v₃(c)=a, every construction in this exterior-probe family must satisfy

    Σ_(j=0)^(d−1) max(3^j−|S|,0) ≤ a.

Indeed the unoccupied subtrees give that many factors of 3 in the automorphism order, and automorphisms act freely on directed Hamilton paths (fixing a directed path fixes every position). Conversely greedy spread and the red/blue coloring attain the smallest group exponent allowed by |S|, so the obstruction from this automorphism group is sharp. This does NOT assert that the produced Hamilton count has this valuation or equals c.

## One probe costs more

In the same protocol with just x, each leaf is either untouched or has its x incidence reversed at the designated vertex. Up through E_d reversals the probe is uniquely recognized by its high outdegree. Its extension group is precisely the binary color stabilizer. Rigidity first occurs at E_d touched leaves, by Session 311's binary neutrality theorem. With two probes it first occurs at N/3, by the first theorem and the fact that each bottom triple needs at least one touch. The two budgets coincide only at d=1. They are minimal for this single-vertex-per-leaf protocol, not for arbitrary extensions.

## Mixed neutral replacements and point marks

Before adding exterior arcs, suppose we may use both named pointwise marks and neutral B→D leaf replacements. Keep every marked leaf unchanged; D is the rigid equal-polynomial partner from Sessions 310–311. If O internal nodes have a marked descendant, the exact minimum number of replacements to get a trivial pointwise stabilizer is E_d−O, for this specific mark set. This follows by splitting into maximal entirely unmarked subtrees. A subtree of height h must be made rigid independently, costing (3^h−1)/2 replacements, and the sum is exactly the number of wholly unmarked internal nodes E_d−O. Occupied ancestors need no neutral replacement, because their rotations are already blocked by the named marks. Thus the Pareto curve with k marks is

    minimal additional replacements = Σ_(j=0)^(d−1) max(3^j−k,0),

and the greedy spreading theorem realizes the whole curve incrementally. The same placements minimize the required replacement cost even though a single binary recolor, unlike a named mark, may create new upper symmetries. For a nonempty fixed mark set, Session 312's maximum-spanning-tree formula computes O without choosing a lexicographic gauge.

The curve is a theorem about the unchanged path-slot polynomial of the core. Adding physical probes need not preserve that polynomial or a specified total. The next step for an exact target selector is to find an exterior with an exact count response to these red/blue channels, or a second legal route that cancels that response while retaining the selected state.

## An actual 3-coprime output and a control

At depth one the all-$B$ core has 24 vertices and $H=127749360032069607$. If a single red incidence is flipped at vertex $v$ and the other probe is left as a sink, the exact new count is
$$
H'=H-s_v+e_v+D_v ,
$$
where $s_v,e_v$ count paths starting and ending at $v$, and $D_v$ counts ordered spanning covers by two nonempty paths whose first ends at $v$. The lost paths started at the added source and then at $v$; the new ones end at the added source or insert it after $v$ between the two pieces. For blue, reverse all paths to get the corresponding start-rooted formula. Each number can be computed by summing the quotient walks of the three leaves, multiplying the ordered-cover counts for the indicated runs in each copy of $B$.

At designation 1 (zero-based $B$ order), the red and blue counts are respectively $378758117678610131$ and $326103131557254703$, congruent to 26 and 22 modulo 27. Designation 0 has both counts congruent to 24 modulo 27, despite the same trivial group. This is an exact finite example of eliminating divisibility and an exact counterexample to automatic elimination by rigidity. It does not certify an all-depth or absolute-target selector. See [replay](probe_replay.py).

## Path packets retain symmetries that the extended graph can lose

There is a stronger obstruction that does not depend on the red/blue protocol. Add $r$ arbitrary exterior vertices to the monochromatic tower. Keep the original tower arcs, but allow arbitrary directions both among exterior vertices and between them and the tower. Then the full Hamiltonian path count is divisible by
$$
3^{\,\sum_{j=0}^{d-1}\max(3^j-2r,0)} .
$$

Partition the paths by the order of the exterior vertices, by which of the $r+1$ intervening tower segments (including the initial and final ones) are empty, and by the **actual** tower vertex on each side of each exterior vertex whenever that side is in the tower. This freezes all used arcs involving the exterior and names at most $2r$ tower leaves in a packet. The pointwise stabilizer of those named leaves inside the original tower acts on the packet: it preserves every within-tower arc of the path and fixes all its used exterior arcs. It need not preserve unused arcs of the extended graph. The action is free, since fixing a path as a directed sequence fixes every vertex. The stabilizer order is divisible by the displayed power, since at depth $j$ at most $2r$ nodes can contain a named leaf. Sum the packet sizes. The $r=0$ case is the ordinary free action on all paths.

For instance a depth-three tower extended by two probes can be rigid, yet **every** such extension, even with arbitrary incidences, has $3^5\mid H$. Targets coprime to three require at least $\lceil 3^{d-1}/2\rceil$ exterior vertices if the core is kept unchanged. This is a hierarchical version of the triangle-contact argument: name used endpoints, not every arc in the full graph. The result is a lower bound, not a sufficiency claim at that threshold. Any proposed fixed-probe exact selector must replace or recolor part of the growing core, or abandon this tower.
