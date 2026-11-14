# A three-square reflection hexagon and the parity memory of 25-block returns

Let \(Q_n\) join distinct positive integers at most \(n\) when their sum is a square. The template and first base path below are inherited, from Gerbicz through the repository's [25-fold blow-up](../collatz_mod6_20260922_w6_square_sum_hamiltonicity.md#6-gerbiczs-25-fold-blow-up-s6-proved-by-a-finite-junction-check). The reflection-height certificate, second base path, independent-switch theorem and template fibre enumeration are the present claims.

## Three cutoffs: motif, readiness, repetition

Write \(R_s(x)=s-x\). For three distinct square labels \(a<b<c\), the composition \(R_cR_bR_a\) is itself a reflection. Traversing the labels \((a,b,c,a,b,c)\) from \(x\) produces the cycle candidates
\[
x,\ a-x,\ b-a+x,\ c-b+a-x,\ b-c+x,\ c-x.
\]
Their positivity is exactly \(c-b+1\le x\le a-1\); their height is \(\max(x+b-a,c-x)\). The six candidates must also be distinct. Arbitrary permutations of the three labels or orientations of the cycle reduce to the displayed sorted order by rotating/reversing. Hence the least possible height of a simple three-label repeated-word hexagon is **23**: below or at 23, the available squares are 4,9,16,25,36, and only \((16,25,36)\) satisfies \(a+b-c\ge2\). Its two balancing choices \(x=13,14\) give the same unoriented cycle, \(2,14,22,3,13,23\) in one orientation. Each alternating matching uses exactly the same three square labels.

This motif is present at 23, but is only later usable as an exchange in a spanning structure with specified endpoints. At **35** the inherited path
```
1,8,28,21,4,32,17,19,6,30,34,15,10,26,23,13,12,24,25,11,5,20,29,35,14,2,7,18,31,33,16,9,27,22,3
```
has a new partner
```
1,8,28,21,4,32,17,19,6,30,34,15,10,26,23,2,7,18,31,33,16,9,27,22,14,35,29,20,5,11,25,24,12,13,3
```
with each label at the same position parity. The switch replaces \(\{2\!-!14,13\!-!23,3\!-!22\}\) by \(\{3\!-!13,2\!-!23,14\!-!22\}\). Both paths close via 1–3. This is a real three-edge exchange, not three arbitrary edges preserving degrees: both literal outputs are connected spanning paths. The repository's [four-cycle cutoff](../crossroads_family_20260926_squares.md#32-a-sharp-hamiltonian-switch-threshold) is 46, so an endpoint-preserving three-edge switch really can activate before any two-edge switch. The standalone fixed-port exact search gives no 1-to-3 spanning path at earlier odd sizes, and five parity signatures at 35 of multiplicities 1,1,1,2,8.

## Full parity fibre of one subsequent return

For an odd length \(n\), choose a separate 1-to-3 Hamilton path \(P_c\) for each \(c\in[-12,12]\). Its block replaces each \(x\) by \(25x+c\) in odd positions of \(P_c\) and by \(25x-c\) in even positions (before any block reversal). The fixed glue order and valid square junctions reside in the verifier; they depend only on the endpoint labels, not on the chosen paths. Let \(\sigma_c(x)\in\{0,1\}\) record the position parity in \(P_c\), with 0 the initial parity. Since \(n\) is odd, reversing a block does not change position parity. Let \(t_c\) be the starting parity of the block in the entire output. It is the index parity of its token in the fixed glue word, as all blocks and singletons have odd length.

The template fills the entire interval if and only if \(\sigma_c=\sigma_{-c}\) for every \(c>0\). Indeed each label \(25x+c\) and \(25x-c\) must occur exactly once and these two offsets are available only from the two signed blocks for \(|c|\). Unequal input parities duplicate one label and omit the other. All junctions are automatically valid when this partition holds. Thus for an input family with parity multiplicities \(m_\sigma\), the template produces exactly
\[
\left(\sum_\sigma m_\sigma\right)\left(\sum_\sigma m_\sigma^2\right)^{12}
\]
distinct fixed-port outputs; injectivity follows by reading the blocks back at their known positions. Fixed endpoints 1 and 3 also ensure distinct closed cycles.

For a matched pair with common input bit \(b=\sigma_c(x)\), the two output parities at \(25x\pm c\) are respectively \((t_c,t_{-c})\) if \(b=0\) and \((t_{-c}\oplus1,t_c\oplus1)\) if \(b=1\). Opposite starts erase that entire input signature, while same-parity starts expose it. The visible pairs are exactly \(|c|\in\{3,5,8\}\), the other nine are hidden. The central block exposes its input signature even after reversal. Hence a fixed output signature with central signature \(\sigma_0\) and visible signatures \(\sigma_3,\sigma_5,\sigma_8\) has precisely
\[
m_{\sigma_0}m_{\sigma_3}^2m_{\sigma_5}^2m_{\sigma_8}^2
\left(\sum_\sigma m_\sigma^2\right)^9
\]
preimages and output paths from this template. Applied at 35, it supplies exactly \(13\cdot71^{12}\) output cycles at **887** in \(5^4\) output signatures; the largest such fibre has \(8^7\cdot71^9\) cycles. These are counts for the template, not exhaustive counts of all 887-vertex cycles.

## Unbounded identical-interface branching, not cofinite size coverage

If all inputs start in one parity class of size \(m\), every block may choose freely. Every output has the same full parity signature and still starts at 1 and ends at 3. Iteration therefore gives at least \(8^{25^d}\) same-signature Hamilton cycles at \(N_d=(71\cdot25^d-1)/2\). Using just the explicit pair gives a cube of \(25^d\) independent three-edge switches: each leaf block has disjoint vertex support, and changing a leaf does not change any junction edge or position parity at any higher level. The embedded cube has edge symmetric-difference size exactly six times Hamming distance. The sizes themselves are sparse. A genuine cofinite-size argument needs an analogous family of returns covering all possible remainder digits, as in the cited but distinct square-sum 49-residue construction. This difference is the analogy to tournament cutoff proofs: a motif can be affordable, reusable and even carry exponentially many hidden alternatives without being able to choose the absolute output size.

Reproduce: `python replay_square_switch_cube.py`. The finite census enforces fixed endpoints and records position parity. The proof above gives the unbounded statement; the replay checks literal switched paths through the first two levels.

## Why 25 is the first nontrivial scale for this exact interface

Consider the same centered odd-square template at (q^2): blocks with offsets from (-(q^2-1)/2) to (+(q^2-1)/2), each input odd with ports 1 and 3; the remaining labels are singleton vertices (1,dots,(q^2-1)/2), and the output is again to have endpoints 1 and 3. For (q=3), the block with offset zero has ports 9 and 27. The available vertices at junctions are exactly ({1,ldots,13}cup{23,ldots,31}). The only available square-sum partner of 27 is 9 (the next candidate 22 is not a port). But 9 is the other end of the same block, and 27 is not a global endpoint. Thus no ordering and orientation of the blocks can work. Since the supplied (q=5) template does work, 25 is the first nontrivial scale in this exact named grammar, not a lower bound for other 9-based or changing-port constructions.
