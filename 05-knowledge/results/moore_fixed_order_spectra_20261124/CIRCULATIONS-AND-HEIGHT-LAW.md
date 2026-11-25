# Frozen-cell circulations and the height law of a guarded +2 family

The two frontiers continue to distinguish a scalar reading from the support that can produce it. For Moore, arbitrary frozen cells force differences into a small spectral subspace. For Chan–Pak, the same digit recurrence controls a continuum-scale height and any fixed congruence, but only a tiny fraction of the exact integers can come from one chart.

## 1. Two reference-free filters for arbitrary frozen cells

Assume the [degree-57 edge-normalized deck](TWO-COORDINATE-RIGIDITY.md) exists. Mark an arbitrary set \(I\) of cells allowed to change; every other vertex is frozen in the sense that its **whole incidence row** is known. The branch–label bipartite graph \(F\) has an edge for each allowed cell, oriented from branch to label; let \(C\) be its incidence matrix. If two completions differ by a symmetric adjacency difference \(D\) supported on \(I\times I\), each column has zero sum along every row and column (one contact per other line of each kind). Thus \(CD=0\).

- In a fundamental cycle basis \(Z\) of \(\ker C\), a zero row forces a fixed cell. Identical rows of \(Z\) force identical rows of \(D\). Two distinct identical nonzero rows of \(D\) would share at least two deleted neighbors, contrary to square-freeness (at least two deletions by the [Latin-row argument](ENDPOINT-AND-SIGNED-DECK.md)). So delete both cells and recompute. In particular, bridges and simple-cycle blocks vanish: on the latter the circulation is an alternating vector with each sign repeated at least twice. Any cactus residual pattern is unique if completable.
- For \(u=(i,j)\), its difference row circulates on \(F\) **with lines \(i,j\) removed**. Each touched line has one plus and one minus, so the nonzero entries occupy vertex-disjoint alternating cycles of that reduced graph. A difference position \(uv\) therefore requires a cycle through \(v\) after deleting the lines of \(u\), and reciprocally. If the reduction is a forest then \(u\) is fixed (in particular if all allowed cells lie in just two rows or two columns).
- The rank is at most the cycle rank of \(F\), or at most \(\dim(\ker C\cap\ker W)\) with \(W=A[O,I]\) the fixed edges from the frozen vertex set \(O\) to \(I\). Indeed the off-block Moore equation gives \(WD=0\). One can also impose the same-line and per-pair zero restrictions on \(D=ZSZ^{\mathsf T}\) as linear equations in a symmetric form \(S\). If the equations force \(S=0\), there is no switch. A reference completion is not needed to compute these obstructions.

### Pin the two eigenspaces separately, using only the induced frozen subgraph

On \(\mathbf1^\perp\), write \(P=(A+8I)/15\), \(Q=(A'+8I)/15\). These are orthogonal projections, and \(D\) is \(15(Q-P)\) on that space and zero on constants, so \(\|D\|\le15\). The image of \(D\) is invariant under both projections: the kernel of \(Q-P\) is invariant because if \(Pv=Qv\), then \(Q(Pv)=Q(Qv)=Qv=P(Pv)\). Orthogonality gives invariance of the image. The restrictions to the image have equal rank (the global ranks agree and they coincide on the kernel). There is no common image and no common kernel on this image space because \(Q-P\) is invertible there. Consequently each restriction has rank **exactly half** the dimension of the image. Thus \(\operatorname{rank}(D)\) is even and the image contains equal numbers of localized \(7\)- and \(-8\)-eigenvectors of either completion. With symmetric edge distance \(M\), \(\|D\|_F^2=2M\). The [distance estimates](ENDPOINT-AND-SIGNED-DECK.md) imply
\[
\begin{array}{c|c|c}
\text{assumption}& \operatorname{rank}D\text{ at least}&
\text{localized dimension needed separately for each sign}\\ \hline
\text{any nonzero difference}&14&7\\
\text{every changed vertex loses at least three edges}&26&13
\end{array}
\]
This proves uniqueness for residual cycle rank at most 13 or 25, respectively.

The separate dimensions are computable even without the missing entries. The projectors on the **whole** space are
\[
 E_+=\frac{A+8I-(65/3250)J}{15},\qquad
 E_-=\frac{-A+7I+(50/3250)J}{15}.
\]
Their ranks are 1729 and 1520 by the adjacency trace. A localized eigenvector vanishes on \(O\); the codimension of that condition within its eigenmodule is the rank of the principal Gram matrix \(E_\pm[O,O]\), because \(E_\pm[O,O]\) is the Gram of the projected frozen coordinate vectors. Therefore an incomplete deck pins the available dimensions to
\[
 1729-\operatorname{rank} E_+[O,O],\qquad
 1520-\operatorname{rank} E_-[O,O].
\]
If **either** is below 7, the completion is unique if it exists, regardless of the sum of the dimensions or the residual cycle rank. After a separate exclusion of two-loss vertices, replace 7 by 13. An induced Gram matrix that is not positive semidefinite or exceeds the total eigenspace rank is already inconsistent. These are exact solver preconditions; integer projectors scaled by a common denominator allow certifying rank lower bounds modulo a prime.

For a Hoffman–Singleton control, freeze all but 22 scattered cells of the \(6\times6\) normalized deck (every line has some unfrozen cells). The residual graph is connected with cycle rank 11. The ranks of the two frozen Gram minors modulo 101 are 26 and 21 (the full multiplicities are 28 and 21). Hence the localized negative module is zero and the completion is pinned, even though the residual graph has nontrivial circulations. The cell list and standard 50-vertex graph are in the replay. The same proof works with both eigenvalues replaced by \(2,-3\). This illustrates separate pinning, not a claim of a degree-57 witness or a fast full search.

### The common transport group has a character constraint

[Smith and Montemanni, 24 November 2026, §§5 and 7](https://www.mdpi.com/2227-7390/14/23/4318) rule out normalized transports in any regular group of order 56. This is an independent external citation for a narrower class than the [repo's earlier block/solvable-group obstruction](https://github.com/eliottcassidy2000/math/commit/1ea98cba). The spectral argument also applies to an arbitrary invariant label subspace, not just the block-constant module.

Indeed if the transport group preserves a mean-zero label subspace \(V\) of dimension \(d\), place one copy on the satellites of each of the 57 branches and zero on root and heads. This is \(A\)-invariant, orthogonal to constants, and \(A\) has trace zero on it because the within-branch blocks vanish. The eigenvalues \(7,-8\) then have multiplicities \(152d/5,133d/5\). Thus every nontrivial constituent dimension of the permutation character must be divisible by five. This screens possible **transport groups**, not automorphism groups. One remaining test is sharp transitivity: the direct transport and the two-hop transports for a fixed nonhub branch pair, including the hub route, form a sharply transitive 56-element set. By [Müller–Nagy's intersection test](https://arxiv.org/abs/1005.1598), subsets \(B,C\) with all \(|B\cap gC|\) divisible by a prime not dividing \(|B||C|\) would rule out a candidate group; no such certificate is claimed here for the Gewirtz action.

## 2. Optimal one-apex gap and exact arithmetic specimens

Use the [two-flank exact inverse](FIBRE-TRANSFER.md). For any fixed positive half-gap \(d\), the unique normalized minimizing pair (up to reflection and common padding) is \((0,2d-1)\), giving \(h(d-1)\) and \(h(d-1)+2d\). Indeed if the smaller flank is \(a>0\), the larger is \(c\ge2d+1\), and \(2ac\) in the lower count already exceeds \(d\). If \(a=0\) the only other possibility is \(c=2d\), whose height \(h(d)+d+(d\bmod2)\ge d+1\) exceeds \(h(d-1)\le d\). Appending a 1 does not change \(h\), proving the stated optimum. The upper \(3d\) is attained exactly when \(d=1\) or \(d-1\) is a power of two: \(h(v)\le v+1\) has equality only for zero and powers of two, by the even/odd height recurrences. The optimum lower endpoint is 1 exactly when \(d\) is a power of two (no inversion in the unpadded word for \(d-1\)).

In the **different** short-guard six-site family, \(B=2684355418\) yields literal one-arc-apart 41-vertex tournaments with Hamilton counts **2203103445089** and **2203103445091**, both prime. A second payload \(B=2684355568\) gives **2225910104699** and **2225910104701**, also both prime. Use masks 177/181 and apex outcontacts \(\{0,1,4\}\) in the six-site block, a one-site guard, and the binary expansion of \(B\) as in [the guarded construction](GUARDS-AND-DECK-DISTANCE.md). The replay checks the exact count formula and trial division through the integer square roots. This is not an infinitude statement.

## 3. Binary reversals give joint real height and modular control

For an \(L\)-bit \(v\) write \(R=\operatorname{rev}_L(v)\) and \(N_d=\sum_i b_i b_{i+d}\) for its bit autocorrelations. The digit statistic from the [short-guard inverse](TWO-COORDINATE-RIGIDITY.md) satisfies
\[
 q(v)=vR/2^{L+1}-N_0/4+3N_1/8-\sum_{d\ge2}N_d/2^{d+2}.
\]
Expand the ordinary-times-reflected digit product. Its upper-triangle weights are already the pair weights in \(q\) except at adjacency, and the diagonal and reflected lower triangle give the corrections. Since \(0\le N_d\le L\), the error in \(q-vR/2^{L+1}\) has absolute size at most \(3L/8\). In particular \(h(v)=v-2q(v)+1-(v\bmod2)\), so for optimal one-apex gap \(2d\) the upper endpoint divided by \(d\) approaches uniform on \([2,3]\) as \(d-1\) runs through growing dyadic slabs.

The short-guard count satisfies \(H(v)/2^L=(v/2^L)(1662-29R/2^L)+O(L2^{-L})\). In a fixed leading prefix, the top digits of \(v/2^L\) and \(R/2^L\) are uniform on their specified intervals and eventually disjoint. Thus the limit is \(XU\) with \(X\) uniform on the leading prefix's dyadic interval and \(U\) independent uniform on \([1633,1662]\). This law holds **within any admissible fixed residue of \(H\bmod m\)**: condition on any finite outer digit blocks and let the middle mix the finite state from the preceding round to uniform. A fixed suffix leaves \(q\bmod m\) uniform and independent of the resulting first coordinate, which is still uniform modulo the odd part of \(m\). Conditional on that coordinate the only possibly fixed residues are parity (always odd) and the factor 29, over which 1662 is a unit. Therefore the height residue is uniform on its admissible classes even given the outer digits, and approximating a real interval by those digits gives the product limit. This controls a macroscopic height interval and a fixed CRT filter jointly.

For example arbitrary finite exact odd-prime valuations on \(H\) and on \(H+2\) are achievable in any positive-mass height interval, as long as the assignment respects local compatibility. At a shared odd prime one cannot make both valuations positive; prescribing both zero is possible because two residues are forbidden and the prime is at least three. CRT and the product law give the result. Conversely actual twin-prime descriptions in the growing slabs have limiting frequency zero: for each fixed prime cutoff the fraction with neither member divisible by an odd prime below it tends to \(\prod_{3\le p\le z}(1-2/p)\), which tends to zero. This gives no uniform-in-cutoff estimate or twin-prime existence beyond the explicit specimens.

All the distribution statements are about **descriptions** and can include collisions. Independently, \(H(v)\ge1604v+427\) caps the proportion of distinct lower endpoints among the odd integers by \(1/802\) even with unrestricted prefixes. Neither full residue reach nor a near-injective scan of a single chart can give absolute coverage. A new family or a more adaptable scalar interface remains necessary.

Run `python3 circulation-height.py` for the frozen Gram rank witness, all four primality tests, the digit error through \(2^{16}-1\), and 65536 small two-flank optimality comparisons.
