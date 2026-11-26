# Physical synchronization eliminates the last tight 128-matching case

2026-11-26. Continuation of [activity transfer](TRANSFERRED-BOUNDARIES.md). All statements about full degree-57 graphs are conditional on two distinct completions sharing the three actual frozen stars. They neither exhibit nor disprove a Moore graph.

## 1. All four zero-slack cases are impossible

The last zero-slack row was $(h,s,l)=(40,4,0)$, $p\in\{12,14\}$. The four closed branches are all high and have changed degree five. All their 20 changed incidences go to outside lows. Let $J$ be the number of such lows. Their active sizes are exactly $p$ by activity transfer. Of the 14 outside lows, at least seven have cap 14 by the hub bound. The integer second-moment test uses $J$ caps $p$, another $\max(0,7-J)$ caps 14, remaining low caps 28, and 39 caps 56. It gives $J\le7$ at $p=12$ and $J\le9$ at $p=14$: the gaps at $J=8$ and $J=10$ are +5,376 and +3,136 respectively, and are nondecreasing thereafter.

Regenerate the 1,057 four-contact trades in the 152 branch quartets. With the link to branch 56 gauged to identity, a target row includes the actual target sheet. The replay finds at most one trade at each triple $(Q,a,x)$ of quartet, target branch in the quartet and antipodal pair label $x$ at that branch. Orient it by putting the actual old sheet $x\in[0,27]$ in the old side; the other orientation shifts all four old contacts by 28.

If two distinct outside lows hit the same active sheet in a closed branch, their two oriented trades cannot have a common old contact away from that branch: this would give a 4-cycle through the two lows. They cannot use the same quartet either, because the oriented trade at the shared target sheet would be identical. Therefore each used quartet describes one particular outside low and the number of distinct utilized quartets is exactly $J$.

For every potential closed branch $a$, the replay enumerates five incident quartets and intersects their actual pair masks, and the ten pairwise masks requiring disjointness of their other old contacts. In every selected quartet the mask has to support the *entire* active set at $a$, not just part of it, because the relevant low-to-closed matching is a bijection on active sheets. Exactly 24 labels are eligible at the six-pair threshold. The same 24 are eligible at seven pairs, with 573 local five-quartet assignments in total and exactly seven usable pair labels on each.

* **$p=12$.** On all 10,626 four-sets $S$ from eligible labels, give a quartet weight $|Q\cap S|$ if it supplies six pairs. The seven largest weights sum to at most **18**, short of the required 20.
* **$p=14$.** A nine-quartet incidence screen at the seven-pair threshold leaves **115** four-sets. For any two closed branches the local lists must agree on exactly which selected quartets meet both. A common quartet also transports the seven active pair labels at one branch to the seven at the other along the actual catalogue rows. Test these equalities and require the union of quartet lists to have size at most nine. **None** of the 115 sets passes.

This closes the last zero-slack row. The conditional lower bound remains **128 changed whole branch matchings**: the six positive-slack $(h,s)$ rows have not been excluded. Those cases can have low closed branches, so the cutoff is not a global no-go. Replay: `node sync_tight.js ../moore57_quadratic_tensors_20261125/period_fourteen_site.txt`.

## 2. Another edit layer excluded on the marked partial board

The 53-row marked contact problem is separate from the full-chart comparison. Starting from the 50 fixed reference contacts and the two holes, the original centerwise exact search excludes *all six* allocations of at most two total collisions within **11 reference edits**, improving the earlier 10-edit exclusion. There was no node-limit truncation. The largest search was the `110` allocation, at 343,746,028 nodes. Since any total at most two lies under one of those six allocations, every solution with at most two squares would require at least 12 reference edits in this fixture. This does not say whether one exists. See [the receipt](anchored_cap11.out).

Rebuild the input as in the [previous replay](README.md), compile `anchored_centers.cpp` with `g++ -O3 -std=c++17`, and run

    for c in 002 011 020 101 110 200; do
      ./anchored ${c:0:1} ${c:1:1} ${c:2:1} 11 11 2000000000 < board.txt
    done

## 3. Relation to regular-group obstructions

[Derek H. Smith and Roberto Montemanni, Theorem 3](https://www.mdpi.com/2227-7390/14/23/4318), independently rule out construction by a group of derangements. In the one-center identity gauge, their group is regular on 56 sheets. Their proof sends Cayley transversal equations to the abelianization. The [spectral transport argument](../moore57_quadratic_tensors_20261125/PRIMITIVE-HOLONOMY-AND-EDGE-SEPARATION.md) uses a different invariant: the dimension of every real centered invariant subspace is divisible by five. It excludes every imprimitive sheet group, including all regular groups of this degree, since a proper subgroup gives a nontrivial block system in the regular action. This extends the exclusion to nonregular imprimitive transport, without deciding existence.

## 4. Small extension of the tournament rigidity result

The [all-order uniform spine theorem](ISOLATED-SPINE-REVERSALS.md) in fact rules out every other edge in the fixed binary left arm, including a root spoke. The bounded-shift proof forces $\Delta R=0$ independently of the edit type. A fixed binary root spoke changes the address, and hence $R$, by a nonzero signed power of two. Thus it cannot give a bounded shift across all odd binary right addresses. In both problems this is why retaining the actual contact profile matters: a conserved boundary channel forces the entire surviving pattern, not just its total occupancy.
