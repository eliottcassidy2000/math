# Partial matching insertions and an unmarked satellite chart

All degree-57 bounds here are conditional on existence and hold for both members of the marked pair. This improves the [edge-rooted grid and native windows](EDGE-ROOTED-GRID-AND-NATIVE-WINDOW.md).

| Encoding | New lower bound | Upper bound |
|---|---:|---:|
| Original, every complement spine | \(2^{19503}\) | \(\tfrac9{10}2^{26853}\) |
| Original, selected full spine | \(2^{20445}\) | \(2^{26725}\) |
| Monotone grid | \(2^{20445}\) | \(2^{25627}\) |
| Grid with two cyclic halves | \(2^{21386}\) | \(2^{26795}\) |
| New: monotone between rows, near-regular inside each row | \(2^{20443}\) | \(\tfrac9{10}2^{25628}\) |

All inequalities are strict. The new uniform native floor excludes \(E_p\pm1\) for \(2\le p\le9752\); it does not assert existence.

## Retain the matching in the completed path

On the grid's balanced cut take a full crossing matching of size \(N=1568\); the Schrijver lower bound is \((27^{27}/28^{26})^N\). Delete \(k\) matched edges and use the remaining \(r=N-k\) backward \(R\to L\) arcs as consecutive blocks. The cyclic switching lemma, with an unconstrained cut anchor, gives
\[
 r!\left(\frac{r-110}{r-109}\right)^{27r+1}
\]
block orders as in the previous grid proof. The degree 27 excludes each matched edge inside its block; the isolated extra mate prohibition can be processed last at the same branching budget. Insert the \(2k\) unmatched sites, one per internal gap. Each uses one within-shore step and one forward cross-shore step. Therefore all \(R\to L\) steps of the finished path recover precisely the retained matched edges and their order. Each singleton's slot is also recoverable.

Rather than divide by the worst possible \(28\)-regular number of completions, choose sparse deletions. For a given full matching use a deterministic pool of \(L\) admissible edges. A uniformly chosen \(k\)-subset has deleted-endpoint crossing-edge count \(X\) of expectation at most
\[
 \mu=k+27k(k-1)/(L-1).
\]
The diagonal contributes \(k\); at most \(27L\) off-diagonal arcs have both endpoints in the pool, and each survives with probability \(k(k-1)/(L(L-1))\). Thus a fraction at least \(1-\mu/(ck+1)\) of the choices have \(X\le ck\). For each resulting output, completions are perfect matchings on those deleted endpoint sets. Brégman and the [concavity of \(\log(d!)/d\)](DEGREE57-HEIGHT.md) bound their number by \((c!)^{k/c}\) because the average row degree is at most \(c\). Any additional pool restrictions only reduce multiplicity. With a fixed order in which to insert the singletons, any initial \(B\) eligible slots for each gives at least \(B!/(B-2k)!\) choices, with no ordering multiplier to divide out. The replay combines these factors with integer powers.

- **Two cyclic halves:** each shore has minimum internal indegree and outdegree 729. Lose at most one to a mate spine change, up to \(k\) unmatched same-shore sites, one block endpoint, and at most 29 forbidden cross joins. This leaves \(698-k\) slots. Use \((k,L,B,c)=(211,1568,487,5)\).
- **Monotone:** require deleted edges to join an early right site within the first \(N-q\) on its shore and a late left site within the last \(N-q\). There are at least \(N-2q\) such edges in each matching. Each deleted site has \(q\) favorable same-shore positions before exclusions. Lose at most 27 graph edges, one spine change, \(k\) unmatched sites, one endpoint and 29 cross failures, leaving \(q-58-k\) slots. With \(q=567\), use \((128,434,381,9)\). This also lifts to the selected full spine: all insertions are internal, so the first site stays on the right and the fixed native prefix can still attach.
- **Near-regular row interiors:** the same monotone argument uses \(q=560\), taking entire rows. This ensures the favorable same-shore positions are outside the singleton's row, so changing its row interior is harmless. Use \((131,448,371,9)\).

## Transfer to an arbitrary native complement spine

Start from the [uniform lower certificate](EDGE-ROOTED-GRID-AND-NATIVE-WINDOW.md#2-uniform-native-lower-bound-spread-the-auxiliary-deficit), but take only its factor/entropy lower bound for the distinct 1547-edge matchings on the balanced native cut. It applies before the factorial and switching factors are multiplied in. There are 78 unmatched vertices per side. Keep their fixed paths as the original prefix and suffix, and select another \(k\) crossing matching edges to delete and insert into the internal block gaps. A completed path can be trimmed by exactly 78 sites at each end; this recovers the original unmatched sets separately from the newly deleted singleton sets. This prevents counting the old unmatched sites as deletion choices.

Take \(q=583\) favorable positions per shore, \(L=1547-2q=381\) candidate edges, and \(k=88\). From the \(q\) favorable positions subtract 57 graph neighbors, one same-shore mate spine change, \(78+k\) unmatched sites, one endpoint, and at most 58 cross-join failures. This leaves \(B=q-195-k=300\). The same deletion expectation uses 56 instead of 27, so take \(c=14\). For \(r=1547-k\) blocks the prior switching bound with constrained prefix/suffix anchor is \(r![(r-238)/(r-237)]^{52929}\): the full cut still bounds forbidden pairs and the maximum indegree/outdegree is still 59. The insertions leave both joins to the anchor intact. The integer product with the inherited entropy bound is above \(2^{19503}\), an improvement of **478 bits** for every original spine.

## Near-regular interiors are decodable without a marked spine

Keep monotone arrows between rows, reversing the satellite graph edges as before. Inside each row put a near-regular 56-site cyclic tournament: send each site to its next 27 cyclic positions and orient antipodal pairs from the lower to the higher index. Use the successive-position cycle, rotating its starting site as needed to avoid the one forbidden graph neighbor at the boundary from the previous row. With no graph edges inside a row, those cyclic paths are legal.

Row \(i\) has a disjoint score band \(3025-54i+\{27,28\}\). In the converse of the mate only the two full-spine endpoints change scores, by one each, still in the disjoint bands enlarged by one. Thus the rows and their order are intrinsic from scores. In the source, the arrows opposed to the monotone template are the matching between each row pair. In the mate's converse, consecutive rows additionally have one spine edge, distinct from the matching. The extra edge is uniquely the edge joining the sole degree-two site on each shore of that matching-plus-edge. Discard it, recover columns by the rook zeros and restore the roots and heads using the exact grid equation. The arrows themselves distinguish the two interpretations by the 56 versus 57 against-template edges at consecutive rows (try the converse for the marked mate). This is not recovery from the scalar Hamilton count. The source maximum score 3053 also differentiates it from the original two grid templates.

For the upper bracket retain the rational row weights of the monotone grid. Instead of the internal outdegree profile \(55,54,\ldots,0\), use 28 copies each of 27 and 28. In the mate the same two boundary corrections go on the first and last sites. The replay maximizes each row product over the four possible assignments of those two local scores and rounds upward exactly before multiplying. Thus the upper bound holds for every legal rotation.

Replay: `python3 grid-insertion-replay.py` (uses the neighboring native deficit certificate).
