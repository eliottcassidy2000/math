# Retrofitting the first edge of the inherited fold

This is a stricter test than searching for an unrelated matched window. Retain the four prepared spines and their existing apex marks from the [matched windows](../chanpak_carry_admission_20261115/MATCHED_WINDOWS.md) at preparation height (j). The original neighboring-address fold changes the first tail mark between 0 and 1; its second chronological mark is 1 in the three first-height examples (j=3,4,5). These examples use fold step (2,4,8), respectively. Allow the first new tail vertex (x) to have **arbitrary arcs back to the entire old prepared spine** and retain its two possible apex marks. Every later tail vertex still follows ordinally.

## Exact search with entire first-vertex contact language

The comparison key after the new vertex is ((h,S)), *not* ((h,S,d)). The common next painted 1 resets the signed bias without changing the count, so discarding cases with differing old (d) would be an overly strong test. Enumerate every subset of old vertices beaten by (x) in each old arm, and both apex marks of (x). Direct subset DP gives both (h,N) and the weighted start (S); no interval approximations or finite-modulus tests are used.

| Preparation height | Marked candidates, including redundant presentations | Distinct ((h,S)) keys | Longest consecutive odd run | Longest such run using both first marks |
|---|---:|---:|---:|---:|
| 3 | 3072 | 2959 | 4 | 2 |
| 4 | 6144 | 5986 | 4 | 2 |
| 5 | 12288 | 12039 | 4 | 2 |

Thus no width-eight closure is possible by editing only contacts of this one new vertex, at the three exact configurations replayed, with all older masks frozen and a common prescribed continuation. The full search has no bound on how many of that vertex's contacts can reverse. The restriction is its number of returning vertices and the three fixed old seed families, not a statement about other seams or greater heights.

## A positive compensated pair involving two actual original arms

Let (T=(18,20,17,7,8)) be the third old unprepared spine, and (R=(48,29,25,1,8,30)) the fourth, with old apex mask 1 each. The third original preparation inserts (j-1) ones then the common zero, and the fourth does the same. In the third arm, mark (x=1). Let (x) beat vertices 0 and 4 of (T) and all old preparation vertices **except the first one**. In the fourth arm, mark (x=0). Let (x) beat vertices 0, 1, 4 of (R) and all old preparation vertices. All unlisted contacts retain the ordinal orientation into (x).

The prepared spines have a simpler exact description after relabeling. In the first one, place (x) directly after the first preparation vertex; in the second put it directly before all prep vertices. The remaining prep words follow ordinally, so the two post-move bases are the following, with apex masks in the new ordering:
```
U=(50,116,113,103,40,64,17),  mask 97: (h,N,S,d)=(77,173,285,51);
V=(48,29,89,65,8,94,19),      mask  1: (h,N,S,d)=(77,175,104,6).
```
Notice (S_U+h=362=2(S_V+h)). After the respective remaining neutral-one widths (j-2,j-1) followed by the guard 0, both have (S'=181\cdot2^j-154) and (d'=0), and their counts are
\[
173+181\cdot2^{j-1},\qquad175+181\cdot2^{j-1},\qquad j\ge3 .
\]
This is an all-height identity: each appended 1 doubles (S+h) with (N) unchanged, then the guard adds the resulting (S+h) to (N). Any common remaining word preserves the +2. In particular it is a literal diagonal cross-contact closure for two of the actual unfolded/corrected arms, at a new (h=77), rather than an extrapolation of the table. It does not recover the complete old width-eight window or keep its original absolute position.

The abstract move resembles the randomly sampled [retained certificate routing](../multiplication_certificate_routing_20261011.md): retiming a constituent to an already shared interface only helps after checking what its consumers actually read. Here part of the apparent long-distance backward motion is removed by relocating (x), but the contacts into the active old body survive and create the new common count.

Replay: `g++ -std=c++17 -O3 replay_fold_first_seam.cpp -o replay && ./replay` gives the exhaustive table with witnesses; `python replay_relocated_fold.py` directly checks the original orientation construction against the relocated bases.
