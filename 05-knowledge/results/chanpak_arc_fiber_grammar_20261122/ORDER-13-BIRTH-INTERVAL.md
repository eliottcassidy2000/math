# A million consecutive odd heights first occur on thirteen vertices

**FINITE-EXACT.** Every admissible odd height through **3,161,043** has a tournament on exactly 13 vertices. Here admissible excludes 7 and 21. Every one of the **1,060,937** odd heights from **1,039,171 through 3,161,043** has minimum order exactly 13. All minimum representatives in that interval are strongly connected.

Let M(n) denote the maximum Hamiltonian path count on n vertices. Split each path on 12 vertices after its first two vertices. Their unordered pair has a unique internal order, and its ten-vertex complement supplies at most M(10) possible orders. Thus

\[
                  M(12)\le {12\choose2}M(10)
                         =66\cdot15,745=1,039,170.
\]

Oddness tightens the numerical upper bound to 1,039,169. A smaller-order witness pads to order 12. This bound improves the first-vertex bound of 1,141,140 used in the component-pair search. The exact maximum data required here stop at order ten. As of this certificate, the public extremal sequence tabulates exact maxima through order eleven. [Underlying census](../path_spectrum_two_angles_metagraph_20261106/data) · [Gordon Royle’s order-11 computation](https://mathoverflow.net/questions/415359/maximum-number-of-hamilton-paths-in-a-tournament-on-n-vertices) · [OEIS A038375](https://oeis.org/A038375)

The order-13 upper witnesses and the bound at order 12 give the asserted exact births. A smallest witness of a height greater than one has no singleton strong component. By the [component-size obstruction](BIRTH-CUT-COVER.md), two nontrivial components would need combined order at least sixteen. Hence all such minimum witnesses are strong.

## What the parent sample covers

Use the same universe as the [fourteen-vertex interval](ORDER-14-INTERVAL.md), every admissible odd height through 3,423,417, and again remove the 43,911 heights already represented by order 11. Generate 55,000 labelled twelve-vertex parents and evaluate all 4,096 apex contact masks for each by the [two-piece path formula](BIRTH-CUT-COVER.md#one-parent-supplies-a-quadratic-table-of-counts). Each equality is a literal thirteen-vertex tournament. Smaller known witnesses can be padded to thirteen by ordinal singletons.

One third of the parents start from the cyclic regular tournament on vertices modulo 13 with positive steps 1 through 6. Make 170 deterministic proposals to reverse a directed triangle, accepting precisely those triples cyclically oriented in the proposed order. Reversal preserves every outdegree. Delete vertex 12 to obtain the twelve-vertex parent. For the other two thirds, independently orient pairs using a varying bias toward increasing order. The source fixes the integer generator, seed, all edge-loop orders, and the schedule. No independence assumption enters the proof.

Of the 1,667,796 targets remaining after the small census, this family misses exactly 22,362. Its first residual is 3,161,045 and its last residual is 3,354,713. These sample omissions have known fourteen-vertex upper witnesses; the sample cannot exclude different thirteen-vertex witnesses. The same computation therefore also gives exact birth 13 to every odd height from 3,354,715 through 3,423,417. Over the full sampled window above 1,039,170, it certifies **1,169,762 distinct heights** of exact birth 13, including the consecutive interval in the title.

For the height 3,161,043 at the top of the first interval, one strong witness has outgoing masks **2712, 4789, 7569, 358, 1896, 5061, 1543, 6488, 3651, 7308, 171, 5242, 1369**.

Run `python3 -B birth_cut_cover.py interval13` from this directory. Success requires the stated residual count and both boundary values over the entire window. Add `--witness 3161043` to emit the displayed graph and independently recount it by subset recurrence. The highest observed contiguous prefix is a lower bound on the full order-13 prefix. [Driver](birth_cut_cover.py) · [C++17 oracle and schedule](birth_cut_cover.cpp)
