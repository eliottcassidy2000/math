# Exact equal-height alternatives in one Hamiltonian frame

This complements the [physical inverse atlas](INVERSE_METAGRAPH.md): replacing the source need not mean changing its count or order. An infinite ordinary family admits a chain of genuinely different representatives in each count fiber, with literal two-arc squares joining consecutive representatives. These are common-frame arc edits, not a claim of two complement lines.

## Construction and height

Take an increasing transitive spine S of length m≥3, and after it a core consisting of a transitive clone list C of length k+t+1, a transitive two-vertex block V=(u,v), and one vertex z. Orient the three core classes cyclically, C→V→z→C. All spine vertices beat all core vertices except that the (k+1)-st clone r beats the first two spine vertices. Assume k≥2, t≥1. The k and t are the numbers of clones before and after r.

In the deletion of r inside the core the ordinary count is h=3·2^(k+t)−1, by the balanced C3 formula with a singleton; the whole core count is κ=3·2^(k+t+1)−1. Let W count ordered two-slot covers of the core with a distinguished path ending at r. Then
```
W=6·3^k·2^t.
```
To see it, separate which path contains z. If it is the marked path, either neither V vertex appears before z (giving 2^k choices), or a nonempty subset of V does (three choices), and each predecessor clone can occupy the marked prefix, marked suffix or ordinary path; each successor clone can occupy the marked prefix or ordinary path. This gives 2^k+3·3^k2^t. If z is in the ordinary slot, the marked slot contains only predecessor clones and r. When both V vertices follow z the ordinary path has a nonempty following clone run since t≥1: 2^k choices. When both V vertices precede z there are 3^k2^t choices. With a split of V there are two choices of its order across the runs, and each requires some clone after z, giving 2(3^k2^t−2^k). The sums give W.

Now use the transitive-spine apex count: if the root beats just the first two spine vertices, H=3·2^(m−2)W+κ−2h. This follows from odd cycles over selected even middle supports; for a mixed cycle either spine start offers a free parity bit, and the separate spine cycles give h(3·2^(m−2)−2). Thus
```
H=3^j 2^w+1,     j=k+2, w=t+m−1,
order=j+w+3.
```
For fixed j≥4,w≥3 this gives w−2 different classes (t=1,…,w−2). The spine degrees occupy every value from N−2 down to N−m, and N−m−1 is absent because all core degrees are at most N−m−2. So the first gap from the top recovers m. Converse gives exactly the reflection t↔m−2. In the converse graph put the two old V vertices in reverse order at the beginning of the new spine, followed by the old successor clones in reverse order. Its clone list is the old predecessor clones reversed, the old z as its new marked vertex, then the old spine with its first two vertices removed in reverse order. Its new V consists of the first two old spine vertices reversed, and its new z is the old marked vertex. Checking the three core directions and two exceptional spine arcs gives G(k,t,m)^op≅G(k,m−2,t+2) literally. Thus the merged family contains exactly ceil((w−2)/2) classes, with precisely one self-converse class in the family when w is odd and none when w is even.

## A literal square and exact edit distance

When t≥2 let x be the first successor of r and y the first predecessor. Independently reverse the disjoint arcs r→x and z→y. All four tournaments share the directed Hamiltonian frame consisting of r, the whole spine, the predecessor clones, the successor clones, V, and z. Neither toggled arc is a frame edge. The corners are:
- no flips: (k,t,m), height H;
- just r→x reversed: (k+1,t−1,m), height 1+(3/2)(H−1);
- just z→y reversed: (k−1,t,m+1), height 1+(2/3)(H−1);
- both: (k,t−1,m+1), height H.

The first move reassigns x to the end of the predecessor-clone list, and the second independently reassigns y to the end of the spine. Outside those roles no arcs change. So one route overshoots, but the other stays at or below the desired starting and ending height. It is a true square in the anchored Hamming / arc-edit graph. It does **not** give a constant-count path between the two neutral ends.

The neutral corners are exactly two single-arc moves apart up to orientation-preserving isomorphism. Passing from (k,t,m) to (k,t−1,m+1) changes the multiset of degrees by deleting [k+t+2,t+1,t+4,k+t+1] and inserting [t+2,t+3,k+t,k+t+3]. After at most one cancellation (check k=2,3,4 and then k≥5), at least three degree entries change. A single arc reversal can change only two, even after relabeling. The square supplies distance two; no corresponding lower bound after merging with the converse is asserted.

This gives a family of actual equal-height alternatives that can be reached with controlled intermediate counts. A future inverse still has to specify which marked or framed coordinate the next operation can read. The square alone does not force any one of these representatives to have a prescribed decreasing line.

Replay: `python3 neutral_family.py`.
