# Two-edit prime ports: finite witnesses, joint cuts and circular avoidance

**PROVED** for the reduction lemmas and wheel; **FINITE-EXACT** for the indicated censuses. The prime/\(7p\) [reduction of the spectrum conjecture](../chanpak_two_holes_atoms_20261123/README.md) is unchanged. Nothing here asserts a uniform two-edit rule for every prime or resolves Lonely Runner or planar unit-distance coloring.

## Two-edit follow-up

- [Collatz-27 itinerary and dyadic/ternary separation](FILLER-COLLATZ-27.md): full parity lifts, a compact literal sevenfold witness, and a dimension drop from 27 to 19 for prime candidates.

- [Joint filler selector and sharp capacity](FILLER-JOINT.md): physical same-phase straddling under 0.164 ppm, plus sharp dyadic sparsity for the full fixed bank.

- [Unrestricted flip squares and planar Jacobians](UNRESTRICTED-JACOBIAN.md): incident-port nondegeneracy, real-square fold classification, and certified integral and disjoint-port examples.
- [Scalable sevenfold module families](TWO-APEX-SCALABLE.md): two literal exact-ratio squares for arbitrary tournament module, with joint residue selection of the cofactor and its curvature.

- [Complete one-apex two-contact inverse](ONE-APEX-INVERSE.md): exact admission in all mark patterns, only \(O(\log p)\) candidate words, a global sparsity bound for the chart, and three all-order nonmonotonic targets within it.

- [Two-adic support compression](TWO-ADIC-SUPPORT.md): use the short-gain arm to cut the leading finite-search order bound in half; explicit nonnested square at 43.
- [Prescribed cofactor phase and common content](PRESCRIBED-PHASE.md): straddle a chosen prime target in the same phase with a simultaneous physical wheel, and isolate the exact shared gcd across all readings.
- [Transitive-spine one-apex seam](ONE-APEX-SEVENFOLD.md): all-order single-contact obstruction, a sparse exact inverse in the same-side chart, and a finite first two-contact realization.

## What actually transfers

The repository's [Lonely Runner anchor-floor work](../../../07-reflections/the-anchor-floor-turns-the-maxgap-crux-into-finite-order-statistics-klein-S153.md) replaces a gap question by exact finite anchored reads, with a quantified sampling loss. The corrected note expressly does not claim AP configurations minimize the mean gap. We use the discrete version of that principle: move from a continuous picture to an exact finite certificate, keeping the same choice of clock for every constraint. The [unit-distance two-port construction](../two_port_unit_distance_20261103/geometry/PROOF.md) is a second model: a private ear and its boundary extension tensor are controlled jointly, not by independent port tallies. Here the bridge is an explicit Hamming square of flips and an exact joint-cut correction, not an identification of graph coloring with Hamilton paths.

## 1. Localize a short walk by its two endpoints

Write \(d(a,b)\) for the minimum number of physical arc reversals from some tournament of height \(a\) to one of height \(b\), allowing arbitrary finite order. For an attainable atom \(a\) in the count monoid, every realization has just one nontrivial strongly connected component (SCC). Let \(a,b>1\) be attained atoms. For any two endpoint orientations on the same vertex set, form their directed union, making each changed pair bidirected. Its SCCs are a common ordinal partition for both endpoints and every intermediate state obtained by flipping a subset of the changed arcs: between union SCCs the orientation is fixed and uniform. No changed pair straddles those SCCs.

Let \(d_{\mathrm{same}}(a,b)\) be the minimum with the nontrivial cores of both endpoints in a common SCC of that union. Then
\[
 d(a,b)=\min\big\{d_{\mathrm{same}}(a,b),\ d(a,1)+d(1,b)\big\}.
\]
In the common-carrier case all other union blocks have height one at both ends and can be discarded. In the different-carrier case the block supporting \(a\) must make a transition \(a\to1\), and the block supporting \(b\) must make \(1\to b\); these consume disjoint marked pairs. Conversely ordinally join two minimizing witnesses for these transitions. The decomposition depends on actual incidence in the directed union, not on two separate component-count marginals.

For \(a=p,b=7p\) with \(p\ge5\) an odd prime other than 7, the different-carrier term is at least three. A single flip from transitive has height \(1\) or \(1+2^j\) (reverse an interval and choose a subset of its interior to go before the reversed arc); neither can be \(7p\), since \(-1\) is not a power of two mod seven. Thus every genuine two-edit witness can be confined to one union SCC, and when there is no one-edit witness both marked arcs must be in that carrier.

A second compression applies to the **ordinary** endpoint SCCs, with no extra vertices needed for intermediate steps. Their nontrivial supports intersect in at least three vertices in a two-edit \(p\to7p\) witness. If neither arc is internal to the lower strong core \(C\), the induced copy of \(C\) persists. If both are internal, the upper nontrivial core stays within \(C\). If exactly one is internal, commute it first. The induced \(C\) after that edit cannot be transitive, since the whole host would then be transitive and its remaining one edit could not reach \(7p\). A directed triangle survives in \(C\), and in the upper endpoint it lies inside the unique nontrivial SCC. This proves the intersection bound in all cases.

Delete every vertex singleton in both endpoints. The endpoints' heights stay exact, and the induced arc distance can only decrease. If one edit is globally impossible, neither of the two marked arcs can disappear; exact distance two is still witnessed. Therefore, with \(B(h)=\lfloor 1+3\log_5 h\rfloor\) from [Busch's strong-order bound](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf), the existence of a two-edit witness in that case can be decided at a single finite order
\[
B(p)+B(7p)-3,
\]
padding smaller witnesses ordinally. This is a complete per-target test, not a bound independent of \(p\) or a promise that the test succeeds.

## 2. Exact two-edit distances at the first seventeen eligible primes

The complete upper-endpoint census anchors an arbitrary Hamilton path, tries all remaining orientations, keeps \(H=7p\) of eligible prime cofactor, and tests every unordered pair of distinct arc flips, including frame arcs. Trying the arcs unordered is complete for the endpoints since their toggles commute. For a negative, every eligible upper frame for its \(p\) is tested; after finding a positive, further frames at the same \(p\) need not be tested. Relabelling and singleton padding preserve the applicability of this finite order test.

| Host order | All eligible \(p\) reached from their own \(7p\) in two flips |
|---:|---|
| at most 6 | none |
| 7 | \(5,11,13,19\) |
| 8 | \(5,11,13,17,19,23,29,31,37,41,43,47,53,59,67\) |

Order 8 has 156422 eligible upper frames and the replay makes 21907113 pair tests before all decisions for that order. Explicit order-nine examples add \(61,71,73,79,83,89\) with no claim that this is the full order-nine list. The [witness replay](two-port-witnesses.js) includes all 21 rows and **each has a strictly increasing choice of intermediate height**. For example the order-nine witness for \(p=61\) starts at \((414,284,280,272,288,335,415,318,0)\), and flips \((0,1),(7,8)\), giving \(61\to93\to427\).

**Every eligible prime \(5\le p\le71\) has exact edit distance two to its own \(7p\), over all orders.** The lower bound uses the [previous one-edit obstruction](../chanpak_planar_prime_ports_20261123/README.md) and its full core-order-eight rooted census: the strong lower core is at most order eight by the strong order-nine minimum 75, and none of these primes occurs in the rooted hit lists. An internal edit is ruled out by the earlier unrestricted one-edit census through order eight. Thus the previously impossible one-flip range is now filled by two flips, not by a universal one-flip repair. For \(p=73\), the earlier one-edit witness uses ten vertices, whereas two edits fit on nine: minimum edit length and minimum host order need not select the same witness.

## 3. Literal unit squares and the missing joint read

Changing two distinct labelled arcs yields a physical square with vertices \(00,10,01,11\). The map to \((0,0),(1,0),(0,1),(1,1)\) is a faithful Euclidean unit-distance representation of this **labelled** flip subgraph, with the two diagonals not unit. It preserves physical arc adjacency, not Hamilton-height differences. In fact the two ends of a \(p\to7p\) path have opposite split/inert signs, so any two-edit route has exactly one sign-changing edge; it cannot consist of two edges of the old sign-changing unit-height model.

At \(p=89\) the order-nine square from rows \((510,508,200,48,36,68,152,312,124)\) and flips \((0,2),(0,6)\) has heights
\[
 H_{00}=89,\qquad H_{10}=H_{01}=373,\qquad H_{11}=623.
\]
The two intermediate tournaments are **not isomorphic** even though they have equal height and even the same sorted score sequence. Each has one vertex of outdegree three. Its outneighbors have outdegrees \(\{2,2,4\}\) in the first intermediate and \(\{2,4,4\}\) in the other. The equal count therefore forgets which physical port was used, while the marked square retains that information.

There is an exact two-port analogue of a boundary extension tensor, using the inherited [apex cut gauge](../chanpak_mirimanoff_unit_lifts_20261122/CUT-GAUGE.md). Append one singleton sink to a core \(C\). Reverse its contacts at \(i,j\), alone or together. Let \(D_i,D_j\) be the first individual gains. The polarization is
\[
H_{11}=H_{00}+D_i+D_j-K_{ij},
\]
where \(K_{ij}\) counts ordered two-path covers of \(C\) with one marked vertex as the first path's end and the other as the second path's start, summed in both orders. This follows by expanding the inherited weighted directed cut at the two contacts; joining both contacts removes the internal cut pair. In particular \(K_{ij}\) is positive even. For sevenfold output the exact condition is \(D_i+D_j-K_{ij}=6H(C)\); summing the gains would accept the wrong endpoint.

The failure of the two marginals is visible **on the same five-vertex host**. Start at \((18,20,25,19,0)\), height 5, with sink 4. The contact pairs \(\{0,2\}\) and \(\{1,2\}\) both give single-contact heights \(9,11\); the joint heights are 13 and 11 respectively (\(K=2,4\)). In the 89 square both edits also use a common apex, initially a source; taking the converse puts it in the sink gauge. Both individual gains are 284 and \(K=34\). We now have a literal geometric square, explicit irrecoverable marginal loss, and an exact integer joint read, without implying either color forcing or geometric height.


There is a useful inverse interpretation. For the old core arc \(i\to j\), let \(M_{ij}\) and \(M_{ji}\) be the two ordered-cover counts in \(K_{ij}\), with directions as subscripts. Concatenation identifies \(M_{ij}\) with the number of internal Hamilton paths lost by reversing the arc and \(M_{ji}\) with the number gained. Hence a core flip gives signed difference \(\Delta=M_{ji}-M_{ij}\), while the appended-apex square gives their sum \(K_{ij}\). Taken together the two experimentally accessible reads recover both loss and gain as \((K_{ij}\mp\Delta)/2\). The joint read therefore restores information that a net height difference cancels.

## 4. The other exact transfer: find one point outside all forbidden phases

The physical factor-seven fork has parent \(7(40Q+D)\) and marked exits \(840Q+19D,\ 840Q+23D\), where \(D=H(M)\) and \(Q=Q_2(M)\) are the readings from **one** filler \(M\). For \(\ell\ge11\), fix \(D=1\bmod\ell\). The three forbidden residues of \(Q\) are distinct: the pairwise determinant magnitudes of the affine forms are \(80,80,3360\). Thus exactly \(\ell-3\) of \(\ell\) phases avoid division of all three readings by \(\ell\).

Place the residues at the vertices of a regular \(\ell\)-gon of side one. These are genuine adjacent unit-distance contacts on this auxiliary circle, but not assumed to be physical flips. If the largest cyclic gap between forbidden vertices is \(G\) grid steps, an integer midpoint in that gap has maximal minimum circular clearance \(\lfloor G/2\rfloor/\ell\). For three forbidden residues \(G\ge\lceil\ell/3\rceil\). This is the exact discrete anchor selection; e.g. at 13 the bad residues are \(2,9,12\), and 5 and 6 attain the optimal clearance \(3/13\). Individual nonzero choices of \(Q\) for different forms do not suffice.

The [full wheel proof](WHEEL.md) handles the small primes and the normalized readings of the entire actual five-edit route. It uses the existing joint modular controller, so the CRT prescription is realized by a single physical filler rather than by independent residue witnesses. It also makes the wheel depth grow like \(\log\log X\) in a sublinear height window. That is a simultaneous roughness result, **not primality**, and not exact matching to a prescribed \(7p\).

## Replay

- `c++ -O3 -std=c++17 two-edit-census.cpp -o /tmp/two-edit-census && /tmp/two-edit-census`
- `node two-port-witnesses.js`
- `python3 wheel-check.py`
