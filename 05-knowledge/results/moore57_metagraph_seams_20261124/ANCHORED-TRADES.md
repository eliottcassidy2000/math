# The anchor forces large, interacting trades

**PROVED conditional.** Suppose two distinct degree-57 Moore completions are expressed in one common 57-branch transport chart with the same labeled anchor matchings (identity). Let H be the graph of changed matching seams on its active branches, r its order, t=57−r, and δ=t+1. The anchor is untouched. A leaf is moved if any incident matching changes; its deleted and added edge counts are equal.

## Exterior witnesses force δ edits at a moved leaf

Take a gained edge x–y between branches i,j. For each exterior (untouched) branch k, let z_k be y's unchanged neighbor in k. In the new graph x–y–z_k is the unique two-step route. The pair x,z_k is nonadjacent in both graphs. Its old unique route x–b_k–z_k must have lost the first edge, since incidence at z_k is unchanged and x–y was absent. The old intermediate b_k is not in branch j, because z_k already has its unique j-neighbor y. Also b_k cannot equal b_l for different exterior branches k,l, since b_k,y would share two old neighbors. Together with the deleted matching edge of x toward j, there are at least t+1 deleted edges at x. Every moved leaf has a gained edge. In particular every active branch has degree at least δ in H, giving r≥30.

## Two minimum-degree branches cannot share an edited seam

If the moved leaf x has deletion-degree exactly δ, the t witnesses exhaust all its deletions other than the one toward j. Let w be the unique old common neighbor of x,y. It cannot be in branch j or equal any b_k (otherwise w,y,z_k would form an old triangle). Thus x–w is unchanged if x has deletion-degree δ, and similarly for y–w if y has deletion-degree δ. A gained edge between two such leaves would create a new triangle.

Thus the degree-δ branches of H form an independent set; a changed seam between them would supply a gained edge violating the previous argument. If a is their number and m is the number of changed seams,

    m ≥ δa,    2m ≥ (δ+1)r−a,
    m ≥ ceil[r δ(δ+1)/(2δ+1)].

For r≤55, even the weaker m≥r(58−r)/2 gives at least 83 over 30≤r≤55. For r=56 the independence bound gives m≥68. Hence exact distinct anchored completions are separated by at least 68 matching edits, even if the unlabelled graphs happen to be isomorphic.

## At the seam bound, most labels would have to move too

If m=68 then r=56 and 32≤a≤34. The total excess beyond degree 3 at higher-degree branches is a−32≤2, so max degree of H is at most five.

More generally if the maximum number of old edges deleted at a touched vertex is D, the old graph induced by the touched vertices T has min degree at least 57−D. For any touched v choose a deleted edge uv. Each untouched old neighbor z of v must now reach u via an added edge u–w followed by unchanged w–z: if u–w were old the old unique route u–v–z had a competitor. The w must be distinct as z varies because any two fixed vertices v,w have at most one old common neighbor. Hence v has at most D untouched neighbors.

At m=68 all touched vertices are leaves of the 56 active branches. The graph induced by the 3136 active-branch leaves is 55-regular and its second eigenvalue is at most 7 by interlacing from the full leaf graph (fiber-constant eigenvalues 56 and −1, remaining eigenvalues 7 and −8). The average old internal degree of T is at most 7+48|T|/3136 by splitting its indicator into constant and orthogonal components. As D≤5, |T|≥2940. Each leaf in T deletes at least two edges, so at least 2940 actual old matching edges would be removed and the same number added in such a minimum-seam trade. No such trade is known.

Counting active seams alone is inadequate: each old path must be repaired at its own label and its own exterior branch.
