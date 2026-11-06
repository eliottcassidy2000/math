# Independent audit of freely prescribed positive even count differences

6 November 2026. **PASS for the positive result in Sections 1–2 and its controls.** The exploratory nine-core negative companion in Section 3 is outside this audit. The result constructs a pair of counts with a prescribed positive even difference; it does not hold the starting count fixed or prove eventual coverage.

## Reviewed inputs and scope corrections

The reviewed ADDITIVE_NOTES.md has SHA-256 382319cb2e26c2d9a617f120fe09218baec638a275fe53977435f14a6a92a0c1. The original da7172d1e314b439154ea54be9fb6a394655e697d86de1d56583a281fce79900 version was preserved by the author. Its only changes narrow “any even integer” to “any positive even integer” and replace the ambiguous length-independence phrase by “freely prescribed with logarithmic order.” No formula or code changed.

The reviewed free_additive.cpp has SHA-256 c808c5e4845658e399f0ca3204ac4d70f683a9d0307eba57fde94c2f337609d5. Its data types and shifts are safe in the tested range, whose largest output order is 14. The graph builder preserves the input root incidences and all edges among the input's other vertices; every added nonroot core vertex dominates those vertices. The independent checker additionally validates every pair of every generated tournament.

## Physical additive identity

Let h=H(U), s=start_v(U), and L be the ordered two-path-cover response with the designated nonempty path ending at v. The rooted edge v->a with a dominating U-v gives H(E(U))=h+L-s by the inherited ordered-root-amalgamation law, including singleton U.

The five-vertex parent core (2,0,11,19,7) is exactly a directed triangle dominating the rooted edge. Hence the complete parent is the ordered join C3->E(U), and its path count is 3(h+L-s). This parent is weak, not strong.

The four-vertex child core (2,8,3,5), rooted at vertex 0, has h_G=5, e_G=2, R_G=6 and d_G=3. Expanding the inherited one-sided amalgamation law gives

    H(G*U) = d_G*h + (R_G-d_G)*L + (e_G-d_G)*s
             + (h_G-e_G-R_G+d_G)*H(U-v).

The final coefficient vanishes for this core, leaving 3h+3L-s. Therefore H(child)=H(parent)+2s for every rooted nonempty U. The core is strong; when U is strong, the two strong induced subgraphs meet at the retained root, so their union, and hence the child, is strong.

The output orders are |U|+4 and |U|+3. The smaller child is not a birth extension of the whole parent. The construction supplies both graphs from the same rooted input, not an operation that accepts an arbitrary tournament with a specified scalar count.

## Arbitrary positive endpoint source

For m>=1, W_m has a transitive private chain, followed cyclically by x,y: chain->x->y->chain. A path starting at p_i must use x immediately before y. The earlier private vertices must lie after y, while each later private vertex independently chooses its increasing block before x or after y. This gives exactly 2^(m-i-1) paths. Paths starting at x or y each have one forced order. Thus H(W_m)=2^m+1 and its start-count vector is the claimed binary weight vector with two additional ones.

Given s>=1, choose m=bit_length(s). Add root v, with out-neighbors precisely the private vertices selected by the m binary digits of s, and with x,y beating v. A Hamilton path starting at v is exactly v followed by a Hamilton path of W_m whose initial vertex is an out-neighbor of v. The sum of its binary start weights is therefore s. There is no extra factor or dependence on subsequent root incidences.

W_m is strong, and the new root has an out-neighbor (the leading binary bit is one) and in-neighbors x,y. Thus U_s is strong and has order m+3. Substituting it in the preceding literal pair gives

    H(C_s)-H(P_s)=2s,
    |P_s|=bit_length(s)+7,
    |C_s|=bit_length(s)+6.

Both counts are odd, as they are tournament path counts; the child is strong. This proves every positive even difference. The audited statement does not assert this construction for negative s or an independently chosen output order.

## Coupling and iteration limits

The parent count remains 3(h+L-s), determined by the actual U_s. In particular every parent count is a multiple of 3. There is no freely chosen base count and no fixed-parent interval of children.

The two cores agree in (d,R,h-e)=(3,6,3), but their ending counts are 0 and 2, and their root-deletion two-bin counts are 24 and 12. The scalar difference does not identify their complete rooted profiles. An additional literal control makes the changed endpoint response concrete: at the retained root, P_s has start count 0 and C_s has start count 2^m+1. In the one-sided construction this follows from start(G*U)=start(G)*H(U-v), with the two core start counts 0 and 1. These output responses are not the freely selected input s.

This does not prove that some future representation change cannot normalize the base or preserve a useful boundary. It establishes that the present pair does not itself supply that extra operation. No base-normalization, unrestricted addition, full-profile rebase, iterability, or threshold conclusion is approved by this audit.

## Independent controls

The standalone check.py imports no producer code. It uses nonmemoized literal directed-path enumeration, rather than the producer's subset dynamic programming. For needed endpoint and cover responses it counts all literal paths on every subset, then pairs their actual path choices across complementary vertex sets.

- Every labelled root-zero input of orders 1 through 4, 75 inputs, passes the edge formula, triangle-product parent, child formula, physical orders, and the strong-input child control.
- Every prescribed s from 1 through 127 passes the exact endpoint count, strong source, positive even difference, oddness, output orders, and changed-root-response controls.
- Independent binary start weights agree for every chain length 1 through 9.
- Both core profiles are reproduced exactly, including the higher-coordinate mismatch.

There are 683 source/output/core graph calculations, plus the nine all-start weight controls. The main graph calculations visit 84,556,997 path prefixes. The largest child count is 21,773, matching the producer's receipt. Normal and optimized Python receipts are byte-identical and use explicit exceptions rather than removable assertions.

The inherited rooted-amalgamation source path session197/path-polynomials/ROOTED_ORDERED_AMALGAMATION.md exists in the repository; its general formula, not finite tests, establishes the all-input identity. The binary-source combinatorial proof establishes the arbitrary-s quantifier.

## Reproduction and publication boundary

    python3 check.py > NORMAL.json
    python3 -O check.py > OPTIMIZED.json
    cmp NORMAL.json OPTIMIZED.json

REVIEWED_ADDITIVE_NOTES.md preserves the reviewed note, including its explicitly exploratory Section 3. REVIEWED_VERIFIER.cpp preserves the producer code. A status promotion should identify Sections 1–2 as proved and leave the exploratory companion outside the independently audited result. No new claim of literature priority is made for the inherited one-sided specialization.
