# The same apex cut in two gauges

**PROVED.** Let \(Q_T(S)=H(T+x_S)-H(T)\), where the new vertex beats exactly \(S\). For distinct vertices \(u,v\), let \(C_T(u,v)\) count ordered two-slot covers whose first path ends at \(u\) and second starts at \(v\), both nonempty. Set \(W_{uv}=C_T(u,v)+C_T(v,u)\). One form of the cut law is
\[
       Q_T(S)=\sum_{u\to v,\ u\in S,\ v\notin S} W_{uv}.
\]
Indeed the number of insertions of \(x_S\) into a Hamilton path preserving that path after deletion equals one plus the number of \(10\) transitions in its tag tape (append sentinels 0 and 1 and count \(01\) transitions). The extra ordinary insertions at a \(10\) arc \(u\to v\) count \(C_T(u,v)\), by cutting there. Breaking the old path boundary instead contributes \(C_T(v,u)\): the new path reads first slot, \(x_S\), second slot, but the uncovered boundary \(v,u\) points the wrong way. The two cases are disjoint.

The inherited [odd-cycle packing law](../chanpak_minimal_core_grammar_20261122/README.md#the-exact-cost-of-contacts) instead says \(Q_T(S)=2\sum_{u\in S,v\notin S}w_{uv}\), where \(w_{uv}\) sums the even-order directed subset paths from \(u\) to \(v\), each weighted by the complementary path count. Polarizing a cut by
\[
         Q_T(\{u\})+Q_T(\{v\})-Q_T(\{u,v\})
\]
isolates its symmetric weight. Thus \(W_{uv}=2(w_{uv}+w_{vu})\), an even positive integer; the subset \(\{u,v\}\) witnesses positivity from Rédei. The difference of the normalized forward-only gauge \(W_{uv}/2\) on actual arcs and the all-directions gauge \(w_{uv}\) is a skew circulation. On an actual arc \(u\to v\) its coefficient is \(w_{vu}\) in the forward direction and \(-w_{vu}\) in the reverse direction. There is an explicit cycle decomposition of this circulation:
\[
 w_{vu}=\sum_{\substack{D\text{ even directed simple cycle}\\u\to v\in D}}H(T-V(D)).
\]
Close each even-order path \(v\to\cdots\to u\) by \(u\to v\). Conversely cut that particular arc of the cycle. The weights are unchanged. Each directed cycle sends just as much across a cut as it brings back, which proves why this gauge correction is invisible in the scalar apex increments. This is an exact even-cycle correction to the odd-cycle packing observer, not an analogy with a continuous flow.

## What a one-apex exponent step must solve

Given the endpoints of an actual marked unit seam at \(E_p\pm1\), extending both by the same contacts \(S\) lifts it to \(E_{p+1}\pm1\) if and only if both weighted cuts equal \(E_{p+1}-E_p\). The common selection problem is therefore the joint cut value \((Q_{T_-}(S),Q_{T_+}(S))\), not either marginal. The two gauges give equivalent exact inverse tests; in the forward-only gauge the mixed quadratic coefficients are \(-W_{uv}\) and the linear ones are weighted outdegrees.

In particular no single-arc seam survives **every** common apex attachment with its same signed gap. If the two cut polynomials were equal on the Boolean cube their polarized weights would coincide on all pairs. Their single-vertex cut at the tail of the flipped edge would then differ by that positive weight, since all remaining orientations are unchanged. More quantitatively at least one quarter of the contact sets change the signed gap: the difference is a nonzero multilinear polynomial of degree at most two. Over the reals a nonzero Boolean multilinear polynomial of degree \(d\) is nonzero at at least a \(2^{-d}\) fraction of the cube. To see this, take a highest-degree monomial and pin the other variables; its top coefficient survives every pinning and at least one point in every selected \(d\)-cube is nonzero.

For arbitrary two hosts, equality of all labelled apex increments is equivalent to equality of every \(W_{uv}\) and zero weighted divergence of the reoriented arcs (use the old directions and their common \(W\), with zero weight on unchanged arcs). Necessity follows from polarization and the singleton cuts; sufficiency follows from the displayed quadratic coefficients. Thus any two distinct hosts with identical labelled extensions must differ along a positive weighted Eulerian subgraph, not on a lone edge. This classifies invisibility to the entire one-apex cut observer and gives a finite test for a proposed marked extension symmetry.

For a one-arc change, the weight at the flipped pair itself is invariant: its two endpoints belong to different slots in either defining cover, so their mutual arc never occurs. The corresponding oriented term in the cut changes purely by reversal. See [the full cross-component classification](ABSOLUTE-SEAMS.md#an-ordinal-factorization-cannot-hide-a-larger-cross-component-seam) for why the lower rail’s ordinal two-factor construction cannot place the desired bridge.
