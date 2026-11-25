# A one-edge weighted repair improves the uniform native lower bound

**Conditional on existence of a degree-57 Moore graph.** Both native endpoints over every complement spine have \(H>2^{18927}\); the upper bounds from [the original height window](DEGREE57-HEIGHT.md) stay unchanged. The uniform exclusion of \(E_p\pm1\) now extends through \(p=9464\), for \(E_p=2^{p-1}(2^p-1)\).

## At most one deficient test set per side

Use the [augmented bipartite cut and parity-moment tests](FIBRE-TRANSFER.md#1-uniform-native-bound-both-endpoints-have-more-than-218903-paths) with \(k=1625,\ r=22,\ q=70\) auxiliary sites per half. Of the size pairs with positive demand and both sizes at least 22, the spectral test leaves 246 and the parity-moment test resolves 245. The exception is \(X\) of size 228 opposite the entire half \(Y\) of size 1625; the demand is 3476. To show its deficiency is at most one, test the cap 3474 with \(g=(2A+I)\mathbf1_X\). The minimum square sum is exactly the required value 259236, with equality only at the cap (the strict transfer inequality still holds). Within the half \(L\) of \(X\), equality forces 5 neighbors in \(X\) at every vertex of \(X\), and 6 at all other vertices of \(L\). In \(Y\) it forces 2 at every vertex except a set \(S\) of 224 sites where it forces 3. These are the exact convex minimizers: on \(Y\) use entries 4 and 6 of \(g\); on \(L\) use 11 on \(X\) and 12 outside it.

That equality profile is impossible. It says \(A\mathbf1_X=2\mathbf1+4\mathbf1_L-\mathbf1_X+\mathbf1_S\). Apply \(A+I\) and use \(A^2+A=56I+J\) to get
\[
 4 A\mathbf1_L+A\mathbf1_S=114\mathbf1+56\mathbf1_X.
\]
Thus all neighbor counts into \(S\) are \(2\bmod4\). The vector \((2A+I)\mathbf1_S\) is \(4\bmod8\) outside \(S\), \(5\bmod8\) on \(S\), and sums to 25760. The smallest allowable squared norm starts at the representatives 4 and 5, then adds eight at 1567 of the cheaper outside sites. It gives 254592, but the Moore moment requires 251104. Therefore \(e(X,Y)\ge3475\).

For each fixed full opposing half \(Y\), the deficit \(22(|X|-70)-e(X,Y)\) is a modular set function of \(X\). Two distinct positive-deficit 228-sets on this side would have union and intersection of sizes different from 228; both have nonpositive deficit, a contradiction to modularity. So there is at most one offending subset per side. Place one artificial nonedge with endpoints in the one or two offending subsets, as applicable. This is possible because each of the two candidate target subsets has size 228, above degree 57. All cuts now pass and the augmented graph has a 22-factor.

## Count with a weighted permanent, then discard the fake edge

On the augmented adjacency give genuine and auxiliary edges weight 1, and the artificial edge weight \(z=1/3110\). If no repair is needed omit the edge. In [Gurvits's permanent lower bound, Theorem 2.2](https://arxiv.org/pdf/1106.2844), use the following doubly stochastic test: original-to-original factor edges carry \(1/22\); distribute every original site's remaining fractional deficit uniformly over the 70 opposite auxiliaries. Its deficits sum to 70 on either original side, so each auxiliary row and column sums to one. For entries with fractional weight \(v\), the logarithmic entropy term is \(f(v)=(1-v)\log(1-v)-v\log v\). On \([0,1/70]\) this is concave with \(f(0)=0\), so an original vertex of deficit \(\delta\) contributes \(70 f(\delta/70)\ge70\delta f(1/70)\) on its auxiliary edges. The original contribution is \((21^{21}/22^{20})^{1555}\) apart from a penalty at worst \(z^{1/22}\). Taken together, these give the lower bound
\[
 \operatorname{per}\ \ge\
 (21^{21}/22^{20})^{1555}(69^{69}/70^{68})^{140}z^{1/22}.
\]

All auxiliary-to-original edges of weight one are available in the counted permanent, so **exactly** \((70!)^2\) completions correspond to each original matching of size \(m=1555\), including cases where the original matching uses the artificial edge. Removing the artificial edge in that case gives a distinct genuine matching of size \(m-1\), uniquely recovering the original matching by putting it back. The two size classes are disjoint. Use the [ordered-block switching bound](DEGREE57-FACTORS.md#2-count-compatible-block-orders-by-switching) \(L_m=m! (b/(b+1))^{52929}\) with \(b=m-238=1317\). By Bernoulli,
\[
 \frac{L_{m-1}}{L_m}=\frac{(1-b^{-2})^{52929}}{m}>\frac1{2m}=z.
\]
The first/last unmatched paths are on their respective sides only, so all right-to-left steps across the cut recover the actual matching, including its size, even after dropping the artificial edge. This converts the weighted count to genuine Hamilton paths, for both native charts by working in the mate's converse as usual. Thus
\[
 H\ \ge\ (21^{21}/22^{20})^{1555}(69^{69}/70^{68})^{140}
 \frac{1555!}{(70!)^2}\left(\frac{1317}{1318}\right)^{52929}3110^{-1/22}
 >2^{18927}.
\]
Run `python3 weighted-augmentation.py` for the integer cut checks and the final exact comparison of twenty-second powers. The artificial edge is a weighted proof device, not an edge claimed to exist in a Moore graph or used in any counted Hamilton path.
