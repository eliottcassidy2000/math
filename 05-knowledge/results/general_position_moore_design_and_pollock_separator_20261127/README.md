# Transfer from oriented trees: geodesic constraints and shared separator data

The linked [oriented-tree proof](https://github.com/1337isnot1337/solved-conjectures/blob/main/results/general-position-oriented-trees/PROOF.md) treats the tree clause of [Chandran et al., Problem 5.2](https://arxiv.org/html/2604.15909#S5). Its two useful moves are distinct. A directed path in an oriented tree is already a geodesic, so excluding three selected points on a geodesic becomes a reachability-chain constraint. And a finite label at a parent is enough to make the child branches independent. Below, the first move yields an extremal structure for Moore graphs and an exact triangle criterion for regular tournaments. The second yields a complete regular four-moment selector for a fixed Pollock fiber. None of the transfers assumes reachability in a dense graph is a tree order.

## 1. Moore graphs: general position is an induced matching with optional singletons

For an undirected Moore graph of diameter two the only geodesics with three vertices are induced paths of length two. It has no triangles. Thus a general-position set \(S\) induces maximum degree at most one: a matching plus isolated vertices. In particular it forbids a local two-contact collision, not an unrestricted long path. This differs from the tree reachability test precisely because distances greater than two are absent.

If \(G\) is \(d\)-regular on \(N\) vertices and has smallest adjacency eigenvalue \(-h\), put \(s=|S|\). Centering the indicator of \(S\) against the all-ones vector gives
\[
2e(S)\ \ge\ (d+h)s^2/N-hs.
\]
Since \(2e(S)\le s\), it follows that \(s\le N(h+1)/(d+h)\). Equality needs **both** a perfect induced matching and equality in the eigenvalue bound. The exact replay exhibits sharp sets of sizes **6** and **20** in Petersen and Hoffman–Singleton, respectively, using their ordinary constructions. Their exterior vertices see respectively three and four chosen vertices. In the Petersen case the four exterior neighborhoods select one end of each of the three matched pairs, and realize all four settings of a binary parity constraint (every prescribed pair of endpoints from different matching edges occurs once).

### Conditional degree 57: a 450 bound and a rigid equality template

Here \(N=3250,\ d=57,\ h=8\), from the [Moore adjacency equation](../odd_pair_matrices_moore_gram_and_pollock_fibers_20261126/README.md). The bound is **450**. If attained, the equality vector satisfies
\[
(A+8I){\bf1}_S=9{\bf1}.
\]
Thus the 450 selected vertices form 225 induced matching edges, and each of 2,800 exterior vertices sees exactly 9 selected endpoints, from 9 distinct matching edges. Each endpoint occurs in 56 exterior neighborhoods. Two endpoints from different matched edges are nonadjacent and have a unique common neighbor, necessarily exterior. Hence the exterior neighborhoods form a **group-divisible design of type \(2^{225}\) with block size 9 and pair multiplicity 1 across groups**. Unlike the earlier coclique designs, the matched pairs are groups, and the design excludes both points of a group in one block. This is conditional on saturation, not a requirement for arbitrary candidate Moore graphs.

One can retain the entire matrix readout. In the block decomposition
\[
A=\begin{pmatrix}M&K\\K^\top&B\end{pmatrix}
\]
with \(M\) the matching, \(K\) the selected-to-exterior incidence, and \(B\) the exterior adjacency, the equation implies \(KK^\top=55I+J-M\) and \(KB=(J-(M+I)K)\). Thus \(KK^\top\) has eigenvalues \(504\) once, \(54\) on the 224 zero-sum directions constant on each matching edge, and \(56\) on the 225 directions antisymmetric on their edges. On their images under \(K^\top\), the eigenvalues of \(B\) are respectively \(48,-2,0\). On the remaining \(2350\) dimensions, the exterior block equation forces \(7\) or \(-8\); the zero trace gives multiplicities \(1280\) and \(1070\). The whole exterior spectrum is therefore
\[
48^1,\ (-2)^{224},\ 0^{225},\ 7^{1280},\ (-8)^{1070}.
\]
These are constraints on the required exterior completion of the design, not a sufficiency statement.

The same local-budget viewpoint works below saturation. For a general-position set of size \(s\) let \(s_0\) be the number of singleton vertices in it and let \(t_y=|N(y)\cap S|\) outside it. Summing the Moore equation on the indicator of \(S\), then subtracting the inside contributions, gives the exact deficit identity
\[
\sum_{y\notin S}(t_y-9)^2=(450-s)(585-s)-16s_0.
\]
For example \(s=449\) allows at most eight singletons and bounds the total squared deviation of the outside local contacts. This turns the scalar extremal bound into a test on actual routes.

## 2. Tournaments: the forbidden local triple is a directed triangle

Every regular tournament has directed diameter at most two. If \(v\to u\) and there were no \(u\to w\to v\), every outneighbor of \(u\) would be an outneighbor of \(v\), which also beats \(u\); that contradicts the equal outdegrees. A directed two-step path is geodesic exactly when the back edge closes a directed triangle. Consequently for **any diameter-two tournament** general position is exactly the absence of a selected cyclic triple, equivalently the induced tournament is transitive. This transfers the local geodesic obstruction, but not the oriented-tree chain-partition formula. It also reuses the [marked cyclic-triple readout](../chanpak_uniform_and_rooted_readout_equivalences_20261126/README.md): that hypergraph is the complete local constraint system for this observable on the diameter-two class.

For the [Fermat circulants at 17](../fermat_dihedral_tournaments_and_lonely_chord_grids_20261127/README.md), the carousel with differences \(1,\ldots,8\) has general-position number 9, while replacing just the chosen difference 1 by 16 gives 6. This replacement reverses the entire orbit of arcs at distance 1, not one physical arc. The replay verifies both maxima by testing for a source and recursing on each induced subtournament. It is an example of a secondary observable varying within the same dihedral symmetry type.

## 3. Pollock: one shared level really does decouple the branches

The second move from the tree proof, conditioning on the label that makes subproblems independent, applies to the [four-moment Gram collision](../odd_pair_matrices_moore_gram_and_pollock_fibers_20261126/README.md) in a different way. Fix \(c,Q,B,p_4\), where the centered offsets have sum zero and power sums \(2Q,6B,p_4\) in degrees two through four. Put
\[
f(x)=x^5-Qx^3-2Bx^2+(Q^2/2-p_4/4)x .
\]
By Newton identities every admissible multiset consists of the five roots of \(f(x)-\lambda\), with only the common constant \(\lambda\) free. For the **regular** stratum, meaning five distinct roots, the following selector is complete:
- If \(f'\) does not have four distinct real roots, reject by Rolle's theorem.
- Otherwise split the real line into its five open monotonicity intervals. In each interval collect the values \(f(y)\) for integers \(y\ge1-c\) with \(5y^2\le8Q\). This bound follows by using the sum-zero condition and Cauchy–Schwarz on the other four offsets.
- Intersect the five value sets. Each common level gives exactly one increasing tuple, since every interval is strictly monotone. Conversely any regular tuple must meet every interval exactly once by interlacing. Newton identities verify all the prescribed moments for the reconstructed roots.

This is a separator computation with one central label \(\lambda\) and five independent branches conditional on it. It is not an independence approximation, nor does it resolve the repeated-root stratum or arbitrary choices of moments. It scans only \(O(1+\sqrt Q)\) bounded integer offsets, rather than tuples of them. The exact replay uses polynomial root counts to place each integer in its interval, not floating critical points.

For the earlier triple \((Q,B,p_4)=(335,-105,143554)\), this proves there are **exactly two** regular centered integer multisets even with all allowed signed offsets: the common levels are \(-57120\) and \(53760\). At center \(c\ge18\) both are positive and represent \(5F(c)+335c-105\), \(F(x)=(x^3-x)/6\); at \(c=17\) only the second remains and at \(c\le16\) neither remains. This pins down what the shared Gram omits in a specific family without asserting that Pollock always has a representation.

Run `python3 replay.py` (SymPy for the root counts) for the sharp small Moore sets, the regular-tournament controls, and the exact Pollock fiber.
