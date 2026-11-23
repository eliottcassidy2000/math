# What the two height-nine cores remember after grafting

Take the strong five-vertex cores $A_8$ (the OOOI deletion word) and $A_{10}$ (OIOI), in canonical lexicographic pair masks 8 and 10, respectively. Both have Hamilton height 9. The following gives an exact contextual equivalence, then a sharp reverse decoder for individualized contacts.

## Every homogeneous context is blind

Both cores have precisely three directed cyclic triangles and one directed cyclic 5-cycle. For each integer $j\ge1$, the allowed-empty ordered path-cover read is
$$Q_j(A_8)=Q_j(A_{10})=j^5+6j^3+2j.$$
Indeed apply the tournament odd-cycle packing identity to the tournament induced in each of the $j$ slot colors and sum over all slot colorings. An odd directed cycle must have all vertices assigned the same color and contributes a factor two. On five vertices there is room for at most one cycle, so only the displayed terms survive. The companion replay enumerates the cycles directly.

When either core is a homogeneous module, deleting it from a full path leaves an outside skeleton and some ordered occupied gaps; each surviving skeleton contributes the corresponding nonempty ordered slot count. Inclusion–exclusion expresses these counts from the common allowed-empty reads above. Thus **every fixed outside shell, including arbitrary outside tournament, singleton padding, and every orientation of all marked outside edges, has the same absolute height on these twins**. One cannot solve the branch recognition problem by reading more homogeneous slot registers.

## One individualized apex with unrestricted ordinal padding

Fix a marked apex $v$ and place the core between arbitrary transitive pads $L\to A\to R$. The apex contacts to each pad vertex may be arbitrary, and its out-neighbors within the core form a mask $S$. For a nonempty directed path $P$ within a block $C$, give it weight $H(C-V(P))$, with $H(\varnothing)=1$, and record its parity by using its number of vertices. Define the two parity reads of $K,U,V,W$ within $C$ by counting, respectively, all such paths, paths starting in $S$, paths ending outside $S$, and paths with both endpoint restrictions (use the appropriate contact mask for each block). For the two odd-order cores the unrestricted nonempty read $(K_e,K_o)=(33,42)$ is common: the two complementary parts of a two-slot covering have opposite parity, $Q_2=84$, and the empty first path has weight 9.

**Complete signature.** All Hamilton heights obtained from these two cores by such one-apex, arbitrarily tagged ordinal continuations depend on $(A,S)$ only through
$$\sigma(A,S)=(W_e,U_o,U_e,V_o,V_e).$$
This follows from the odd-cycle packing expansion. Contributions with no cycle through the apex give 9. If its odd cycle is $v\to P\to v$, then $P$ has even vertex count, starts at an apex out-neighbor and ends at a non-out-neighbor. Its contribution is twice the complement height. The path splits into ascending nonempty runs through the three ordinal blocks. If first and last runs are in the same block, use its $W_e$ and multiply by the complement height of the other blocks. Otherwise use the first block's $U$, the last block's $V$, and each intermediate block's unrestricted $K$ including the empty term of weight $H(C)$. Multiply in the parity ring with $s^2=1$, extract the even component, and multiply by the complement heights of any untouched exterior blocks. Summing by first and last blocks counts each path once. This proves sufficiency of the signature for all pad lengths and all their contact words.

**Five explicit inverse probes.** Let $h_0$ be the unpadded apex height. Probe one or two singleton vertices in the left pad, all out-neighbors of the apex, and no right pad; call the heights $l_1,l_2$. Probe one or two on the right pad, all in-neighbors of the apex, and no left pad; use $r_1,r_2$. Then
$$
\begin{aligned}
W_e&=(h_0-9)/2,\\
U_o&=(r_1-h_0)/2,& U_e&=(r_2-h_0)/2-2U_o,\\
V_o&=(l_1-h_0)/2,& V_e&=(l_2-h_0)/2-2V_o.
\end{aligned}
$$
For instance, with two right vertices the nonempty suffix in the pad has two singleton choices and one two-vertex choice; the required even total parity gives $2U_o+U_e$. This also proves necessity for indistinguishability by any ordinal continuation of this kind.

Across the 64 pairs $(A,S)$ there are 60 signatures. The only four collisions, using integer masks in the fixed canonical rows, are
- $(8,0)$ vs $(10,0)$ and $(8,31)$ vs $(10,31)$, the two homogeneous orientations;
- $(10,10)$ vs $(10,20)$ and $(10,11)$ vs $(10,21)$.

Both cores have trivial automorphism group, so the last two collisions are genuinely different rooted contacts. Every nonhomogeneous one-apex signature identifies the core, even when it cannot identify the contact. An unpadded height alone need not identify either. Its full out-neighbor-size/height atlas is printed by the replay; in particular 41 occurs on $A_8$ but not $A_{10}$, while 43 occurs on $A_{10}$ but not $A_8$.

## A sixth probe gives a complete reverse decoder

The last four collisions all separate using the same second apex $w$: on the unpadded configuration, let $w$ dominate core vertex 0 only and lose to the other five vertices, including $v$. The four ordered response pairs in the order above are $(27,31),(17,19),(65,77),(59,65)$. Therefore the five ordinal probes followed, if necessary, by this one individualized probe decode the full pair $(A,S)$ uniquely. The decoder treats the original core labeling as its fixed canonical representative and the first apex as designated. It does not infer a core from a single order or a single output height.

Run `python3 twin-continuations.py` to reproduce the atlas. For a shorter decoder in a different probe family, see [the two-output decoder](README.md#the-shortest-decoder-in-the-one-extra-apex-probe-family).
