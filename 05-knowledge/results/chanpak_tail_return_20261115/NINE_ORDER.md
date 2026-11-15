# Nine consecutive odds versus exactly nine vertices

The [mixed-contact slope-255 lane](../chanpak_carry_admission_20261115/MIXED_CONTACT_LOCKS.md) gives
\[
255\,2^j+\{235,237,\ldots,251\},\qquad j\ge2.
\]
Its nine arms use seven-vertex spines; the final construction has order j+8 or j+9, depending on whether the pre-suffix coefficient S+h is 510 or 255. In particular the original recipes have order ten or eleven at the lowest height j=2, not nine. This note checks separately how many vertices the actual counts require and how much of the nine-run can share a literal continuation state.

## 1. Both first heights really are nine-vertex phenomena

An exact isomorphism census through order nine gives maxima H=661 at eight vertices and H=3357 at nine. Every member of both first two heights occurs at order nine:

| Height | Count interval, odd increments | Number of order-nine classes for each successive count |
|---|---|---|
| 2 | 1255 through 1271 | 212,198,193,206,210,221,211,198,178 |
| 3 | 2275 through 2291 | 129,135,149,131,124,145,115,119,147 |

Hence **all eighteen counts first occur on precisely nine vertices**. Moreover every order-nine witness for these values is strong: a nonstrong tournament splits into ordered components, and the path count is a product over the split. The computed smaller maxima bound the split product by max(1·661,1·189,3·45,5·15)=661. Every higher rung needs more: already at height 4 the lowest target is 4315, greater than 3357. Each table entry's representative is in [the exact census receipt](replay_nine_census.txt). A single tournament has only one H value, so this is about a family at fixed order.

The census starts with the singleton and adjoins one vertex with every possible out-neighbor mask to each retained representative. This is complete by deletion of any vertex. Exact canonical form refines an ordered partition of the vertices by the vector of out-neighbor counts into its cells; when refinement stalls it individualizes all possible vertices in the first nonsingleton cell and minimizes leaf encodings. Isomorphism preserves partition refinement and the entire leaf set. Two equal leaf encodings actually relabel to the same tournament. The resulting class counts are 1,1,2,4,12,56,456,6880,191536 for orders 1–9, and the endpoint-subset DP yields H.

## 2. The nine first-rung counts are not already one matched-window family

Fix a distinguished apex on each order-nine tournament. Enumerate all 6880 nonisomorphic eight-vertex spines, and for each all 256 apex masks. This covers every marked isomorphism type, with some repeat presentations from spine automorphisms. Let E_x,O_x sum the Hamilton counts of the complements of all nonempty even/odd-order spine paths beginning at a vertex beaten by the apex. Then h is the spine path count, S=E_x+O_x and d=O_x−E_x. These are the inherited exact right interface; an ordinal singleton suffix b updates the marked count by (1−b)(S+d) and sends (S,d) to (2S+bh,bh).

The masked H computation in the second replay is separately derived and checked against literal whole-tournament DP: a cycle containing the apex consists of an even-order spine path from a marked start to an unmarked end. Grouping by the cycle containing the apex in the [odd-cycle packing formula](https://arxiv.org/abs/2412.10572) gives twice the Hamilton count of the complementary spine per such path; packs not using the apex give h. Thus no unresolved search over right tails is used for the finite assertions.

- For the first rung (1255 through 1271), there are 16,279 distinct literal state triples (h,S,d). The largest set of first-rung counts in any one triple has size three, but the longest adjacent run has only two entries.
- This maximum two is attained at 1267 and 1269 with common (h,S,d)=(81,1717,23), using rows (2,12,49,21,67,27,47,127), mask 232, and (2,36,25,67,11,29,55,127), mask 228. A common 0 gives 3007 and 3009, and leaves shared (81,3434,0). Thus the entire pair persists for all canonical common transitive tails as 3007+3434t+81p(t)+{0,2}.
- A formal reason to consider a different grouping is the neutral-one law: a run of ones leaves N unchanged, doubles K=S+h once per bit, and resets d=h. Hence exact alignment by neutral runs is possible **if and only if** two seeds have equal h and equal oddpart(K). This criterion is from [matched windows](../chanpak_carry_admission_20261115/MATCHED_WINDOWS.md); the new census applies it without an a priori tail-width bound.
- Under that enlarged criterion the longest adjacent run within the first rung is three: 1261,1263,1265 with h=239 and K respectively 1395,2790,1395. Witnesses and masks appear in the receipt. Append g,g−1,g ones with g≥2 and a common 0. The three now share h=239, S=2(1395·2^g−239), d=0, and ordinary common growth gives the exact two-parameter window
\[
1261+1395\,2^g(2t+1)-478t+239p(t)+\{0,2,4\},
\qquad g\ge2, t\ge0.
\]
The post-alignment h and S are relatively prime. This does not give all absolute targets; it gives a wider matched branch than any single rigid apex state on that rung.

## 3. The best nine-vertex matched widths anywhere, not just this rung

Across **all** order-nine marked choices, exact equality of (h,S,d) gives at most four adjacent odd counts, attained by 89,91,93,95 with state (9,245,9). If neutral-one alignment is allowed, the exact maximum is **six**, at 129,131,…,139. The six all have (h,S)=(45,129), but their biases need not match. The same 1 resets all biases to 45 and takes S to 303; the same 0 then adds 348 and takes S to 606, resetting bias to zero. Thus one obtains six literal counts at 477,479,…,487 that survive arbitrary common transitive tails:
\[
477+606t+45p(t)+\{0,2,4,6,8,10\},\qquad t\ge0.
\]
No strong-tournament claim is made for all six; a dominating spine component not reached by the apex can persist through ordinal growth. Here p is the [exact binary suffix correction](../chanpak_collatz_exact_addresses_20261114/BINARY_SUFFIX_DECODER.md): p(0)=0, p(2u+1)=p(u), p(2u)=p(u)+u+(u mod 2). These maximality results are for marked **order nine** and the named form of alignment, not a global theorem that matched windows stop at six.

In fact the maximum stays six even if all seeds must have a **strong** completion under a neutral-one tail ending in 0. This condition is equivalent to the original apex reaching all eight spine vertices, which the marked census checks by transitive closure. A final 0 then gives a return route to the apex from every vertex. For the same nine-vertex counts 129,131,…,139, the strong-completable witnesses have K=696,696,174,696,696,696 and h=45. Append g−3 ones when K=696 and g−1 ones when K=174, then a common 0, for any g≥4. This gives the strong two-parameter family
\[
129+87\,2^g(2t+1)-90t+45p(t)+\{0,2,4,6,8,10\},
\qquad g\ge4,\ t\ge0.
\]
In particular g=4,t=0 gives 1521 through 1531 by twos, with common h=45,S=2694,d=0. This grows the four-arm interface to six, but does not yet supply a nine-arm return.

The same order-nine maxima hold if one allows **any count-neutral transitive preparation**, not just an explicitly all-one one. A nonempty mark set has positive total odd-start weight: the singleton path at any marked vertex contributes a positive complement count. A 0 would therefore increase N by twice that weight. A 1 keeps the mark set nonempty, so no later 0 is neutral either. With an initially empty mark set, neutral leading zeros leave the whole state unchanged and can be stripped before the first 1. This explains the scope of the finite normalized maximum without imposing a bounded number of neutral operations.

## Replay

Compile and run `replay_nine_census.cpp` and `replay_nine_marks.cpp` using C++17. The second includes the first's canonical code. The first returns the multiplicities, and the second covers exactly 1,761,280 spine/mask presentations, independently counts the marked-interface widths, and outputs physical row data. The named continuations also have the separately coded literal Python replay `replay_nine_tails.py`. The infinite suffix statements follow by the exact transfer identities, not by extending the finite replays.
