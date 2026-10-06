# Poisson patches: retain the ports, forcing, and boundary bill

**PROVED:** exact bounded Schur elimination, associative patch composition, the localized residual bill, and the exact-lumping obstruction below. **FINITE-EXACT:** the declared matrix universes and ten explicitly completed Collatz controls. **OPEN:** independently checked infinite-domain patches supplying a positive readout for every source. No new ROOT coverage is claimed.

[Program](../../04-computation/experiments/collatz_poisson_patch_20261005.py) and [output](collatz_poisson_patch_20261005.out).

## 1. Inheritance and the object being combined

The closest mechanism is P1--P2 of [the source-sensitive Poisson dual](collatz_poisson_source_dual_20261005.md): a bounded global adjoint subsolution gives a source floor; a positive finite-support witness contains an actual ROOT path. The kernel and its killed ROOT boundary come from [Green weights](collatz_green_weight_20261005.md). The canonical hostile is a contracted component disconnected from ROOT. The corrected near miss is to validate local equations while deleting their boundary residuals. The least-used sidecar here is the **signed forcing retained at the interface**.

Two recovered operations inform the audit without supplying a Collatz theorem. [THM-3501, Rule 30 universal-cover Green potential and slack-holonomy seam](../../01-canon/theorems/THM-3501-rule30-universal-cover-green-potential-and-slack-holonomy-seam.md) keeps an exact boundary seam when a local potential is changed. The corrected [THM-1770, first-return renewal](../../01-canon/theorems/THM-1770-the-localisation-lemma-as-first-return-renewal-and-the-pair-only-closure.md) warns that primitive return terms need not vanish separately. Here the kernel is nonnegative, while the source forcing may be signed; its terms are retained and summed, not declared separately positive.

The board is **selected source / infinite interior / retained ports / first-return weights / signed residual / finite observer**. The map is elimination of interior coordinates. It preserves the full bounded inequality and its root readout. It forgets individual interior values only after recording both the return kernel and transported forcing. Those data, labels, and an operator norm are the required sidecars.

This is an elementary operator construction with the complete proof below, not a claim of novelty for Schur elimination.

## 2. Kernel and exact patch theorem

Use the primitive sibling bases and notation of the inherited Green note. At price \(r=1/2\), put

\[
 P=A^*,\qquad
 (P\phi)(c)=\sum_{G(b)=c,\ b\ne1}2^{-k_b-1}\phi(b),
 \qquad \|P\|_{\infty\to\infty}=\kappa=6/7.
 \tag{1}
\]

The ROOT source row of \(A\) is deleted, so the ROOT column of \(P\) is zero. Let \(g=(I-A)^{-1}\delta_1\). For bounded \(q,\phi\), the inherited inequality
\((I-P)\phi\le q\) gives \(\phi(1)\le\langle q,g\rangle\).

Partition the **actual labelled base space** into an interior \(I\) and retained ports \(J\), with \(1\in J\). Either set may be infinite. All blocks act on the corresponding bounded-function spaces. Define

\[
 R=(I-P_{II})^{-1},\quad
 K=P_{JJ}+P_{JI}RP_{IJ},\quad
 b=q_J+P_{JI}Rq_I.
 \tag{2}
\]

The resolvent exists, is nonnegative, and has norm at most \((1-\kappa)^{-1}\).

**Patching theorem.** A bounded port function \(h\) extends to a bounded global subsolution with those port values iff

\[
 (I-K)h\le b.
 \tag{3}
\]

When this holds, a valid extension is

\[
 \phi_J=h,\qquad \phi_I=R(q_I+P_{IJ}h).
 \tag{4}
\]

Indeed any global subsolution satisfies
\((I-P_{II})\phi_I\le q_I+P_{IJ}h\). Positivity of \(R\) gives \(\phi_I\le\widehat\phi_I\), where the latter is (4). Increasing the interior values reduces the port defect, because \(P_{JI}\ge0\). The resulting port inequality is exactly (3). Conversely direct substitution of (4) gives equality on the interior and (3) on the ports. This proves both directions, including signed \(q,h\).

The reduced kernel is contractive:

\[
 K\ge0,\qquad K\mathbf1\le\kappa\mathbf1.
 \tag{5}
\]

For a direct proof, \(u=RP_{IJ}\mathbf1\le\mathbf1\), since \((I-P_{II})\mathbf1\ge P_{IJ}\mathbf1\); in fact \(u\le\kappa\mathbf1\) from \(u=P_{II}u+P_{IJ}\mathbf1\). Thus
\(K\mathbf1=P_{JJ}\mathbf1+P_{JI}u\le P_{JJ}\mathbf1+P_{JI}\mathbf1\le\kappa\mathbf1\).

The entries of \(K\) are sums of weights of paths from one port to the next with all intermediate vertices in \(I\). They are **first-return weights**, not averages of vertices in one equivalence class.

If one replaces \(R\) by \(R_N=\sum_{j=0}^{N-1}P_{II}^j\), with \(R_0=0\), then

\[
 \|K-K_N\|_\infty\le\kappa^{N+2},\qquad
 \|b-b_N\|_\infty\le
 \frac{\kappa^{N+1}}{1-\kappa}\|q_I\|_\infty.
 \tag{6}
\]

For the first bound use the exact remainder \(P_{JI}P_{II}^Nu\) and \(u\le\kappa\). The second uses \(P_{JI}P_{II}^NRq_I\). These are truncation bounds for a specified kernel; they do not certify a guessed return matrix.

## 3. Composition, signed readouts, and localized error

Eliminating two disjoint interiors successively gives exactly the same \((K,b)\) as eliminating their union. One proof groups each first-return path according to its successive visits to the intermediate ports. Nonnegativity justifies the kernel sums; bounded signed forcing is absolutely summable by the contraction. Equivalently, both elimination orders give the unique harmonic extension for every final port value and every forcing. **The forcing must be composed as well as the matrix.**

The maximal bounded reduced solution is

\[
 h_*=(I-K)^{-1}b.
 \tag{7}
\]

Every subsolution of (3) satisfies \(h\le h_*\). Therefore a positive root readout exists iff \(h_*(1)>0\). For finite ports this is an exact finite rational linear-algebra test whenever the supplied coefficients are exact rationals.

Define the nonnegative occupation vector

\[
 z=(I-K^*)^{-1}\delta_1,\qquad
 \|z\|_1\le(1-\kappa)^{-1}.
 \tag{8}
\]

Then \(h_*(1)=\langle z,b\rangle\). Thus \(\langle z,b\rangle\le0\) is an explicit dual certificate that **this exact patched signed-source problem** has no positive root subsolution. It is not a divergence certificate for an unrelated source, nor an obstruction to different forcing.

For an approximate reduced inequality

\[
 (I-K)h\le b+\varepsilon,\qquad \varepsilon\ge0,
\]

the global measurement satisfies

\[
 \boxed{\quad
 \langle q,g\rangle\ge h(1)-\langle z,\varepsilon\rangle.
 \quad}
 \tag{9}
\]

The vector \(z\) is exactly \(g\) restricted to \(J\), by eliminating the interior of the transposed Green equations. If interior and port defects are bounded by \(e_I,e_J\ge0\), their exact transported error is
\(\varepsilon=e_J+P_{JI}Re_I\). This keeps the location of the error; replacing it with a uniform norm can be much more expensive.

For an explicit completed control, the source 7 has the base path

\[
 7\longrightarrow11\longrightarrow17\longrightarrow3\longrightarrow1,
\]

with weights \(1/2,1/2,1/4,1/4\). Retain ports \((1,11,7)\), eliminating \((3,17)\). Their occupation values are \((1,1/32,1/64)\). For signed forcing \(q=\delta_7-\delta_3/32\), the root readout is

\[
 1/64-(1/32)(1/4)=1/128.
\]

An error \(\eta\) located only at port 7 costs \(\eta/64\) in (9). These values come from a replayed path; this example tests signed cancellation and error transport, not new coverage.

When \(q=\delta_s\) and \(s\in J\), the interior forcing vanishes. Then (4) has \(\|\phi_I\|_\infty\le\kappa\|h\|_\infty\), so \(\|\phi\|_\infty=\|h\|_\infty\). An independently established positive port readout \(h(1)\ge\delta\) therefore supplies the inherited path deadline: any integer \(T\) with \(\kappa^T\|h\|_\infty<\delta\) exceeds the length needed by the maximum-principle extraction. Alternatively the selected Green floor \(\delta\) feeds the inherited fixed-price counter deadline. The missing premise is the global validation of (2)--(3), not an assumed ROOT certificate.

## 4. What compression retains in the actual functional graph

Each nonroot base has a unique forward parent \(G(b)\). Consequently **each column of the actual first-return matrix \(K\) has at most one nonzero entry**. Starting at a port \(b\), follow its unique forward parent path until its first later port. If that port is \(c\), then \(K(c,b)\) is the product of the traversed edge weights. If no later port is reached, the column is zero. A return to the same port can occur in a cycle; the root column remains zero.

This statement is a path characterization, not an algorithm deciding zero columns in finite time. Positive entries need a checked path or a proved symbolic family that implies such a path. Merely labelling two ends of an unresolved interior does not establish an entry.

A disconnected two-cycle with edge weights \(1/2\) illustrates the boundary. Eliminating one of its vertices leaves a retained self-return weight \(1/4\), and its forced coordinate has the unique value \(4/3\). The ROOT readout is still zero. Strict contraction and a positive coordinate in another component do not replace root forcing.

## 5. Exact lumping has a different and stricter information bill

Let a partition of the base space be **exactly row-lumpable** for \(P\): for every two vertices in one block and every target block, their total transition weights into that block agree. A bounded function constant on the blocks is then acted on by a genuine quotient kernel, with the same strict row-sum bound.

**Singleton-forcing obstruction.** Such a block-constant function cannot satisfy
\((I-P)\phi\le\delta_s\) with \(\phi(1)>0\) unless the block containing \(s\) is a singleton. Otherwise another vertex of that block has exactly the same left side and right side zero. Every quotient block then obeys the homogeneous inequality \(h\le\overline P h\). Iteration and boundedness give \(h\le0\), a contradiction. This applies even to a countable partition with bounded \(h\).

In the functional Collatz graph, an exactly lumpable partition that isolates a nonroot source \(s\) must also isolate \(G(s)\). Only this parent has a positive transition into the singleton \(\{s\}\). Any other vertex in the parent's block would have transition weight zero, contradicting exact lumpability. Repeat this argument along the actual parent path.

Thus a finite exact quotient that certifies a rooted source must retain **every distinct vertex of its strict base path as a separate state**. If a path were infinite and distinct, no finite exact lumping could isolate it in this way. This is a finite-memory restriction on exact lumping, not on all compressed proof descriptions. Schur elimination is permitted to remove intermediate vertices precisely because it stores their weighted return segment rather than identifying them with another vertex.

The smallest checked averaging hostile is a three-vertex chain with weights \(1/2\). Merging ROOT and the middle vertex fails lumpability: one has zero transition into the target singleton, the other has weight \(1/2\). Averaging those rows invents a transition for one representative. The checker rejects it.

## 6. Infinite patches and the boundedness boundary

The theorem genuinely permits infinite interiors. As a non-Collatz synthetic control, let ROOT move to \(v_0\) with weight \(a=1/2\), and each \(v_i\) move to \(v_{i+1}\) and a target port with weights \(b=c=1/4\). The reduced ROOT-to-target entry is \(ac/(1-b)=1/6\). The bounded extension has value \(c/(1-b)=1/3\) on every interior vertex. Truncating after \(N\) possible exits has exact error \(acb^N/(1-b)\). This kernel has multiple parents for the target, so it is explicitly **not** the Collatz functional graph.

The local checks alone are insufficient: on the three-vertex chain, assigning value one at ROOT and zero elsewhere satisfies the zero interior equation but leaves an unpaid positive ROOT defect. Deleting that port produces a false certificate. Nor can the contraction assumption be dropped; a unit self-loop makes the resolvent singular. Both are exact hostile controls.

If a positive weight \(V\) has a separately proved bound \(PV\le\kappa_VV\), \(\kappa_V<1\), the same patch theorem applies after conjugating to \(P^V=V^{-1}PV\) and dividing forcing and functions by \(V\). The reduced kernel conjugates by the port restriction of \(V\). This is a conditional weighted extension; the bounded proofs here do not infer a suitable \(V\) from finite tests.

The independently checked [weighted Poisson growth note](collatz_weighted_poisson_growth_20261005.md) supplies two such envelopes, \(V_p=((3b+1)/4)^{1/12}\) and \(V_l=(25+\operatorname{bitlength}(3b+1))/28\), both normalized to \(V(1)=1\).

## 7. Reproduction and scope

Run either command from the repository root:

```text
python 04-computation/experiments/collatz_poisson_patch_20261005.py
python -O 04-computation/experiments/collatz_poisson_patch_20261005.py
```

The finite universe is all 512 nonnegative \(3\times3\) kernels with entries zero or \(1/4\), each with all four port sets containing vertex zero; 128 explicit four-state kernels test nested elimination. The checks use signed forcing, harmonic domination, exact dual identities, and contraction. A separate complete small universe consists of all rooted increasing-parent trees on 2 through 5 labelled vertices and all their set partitions; accepted exact lumpings must preserve singleton ancestry. Ten listed actual Collatz sources are resolved literally within a declared 1,000-base-step control cap before their finite parent-closed matrices are formed. The script never infers an infinite-domain inequality from this census.

The useful next obligation is now explicit: supply globally justified return kernels and transported signed charges on a manageable port set, with a positive root readout after the localized error bill. A positive analytic construction could certify a source. No such universal construction, new residual-source coverage, or exclusion of hypothetical Collatz divergence is asserted here.
