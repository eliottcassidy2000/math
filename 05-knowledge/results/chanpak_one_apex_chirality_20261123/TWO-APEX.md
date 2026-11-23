# Two individualized vertices are the exact universal orientation threshold

**Proved.** Permit arbitrary labelled contacts from at most \(k\) exterior vertices to an old tournament \(T\); every other exterior vertex sees the whole old \(T\) uniformly. Uniform triangle wraps of a completed probe are permitted. For \(k=1\), a tournament and its converse are indistinguishable precisely under the all-slot divergence criterion below. For \(k=2\), the complete observer recovers the entire old labelled orientation. This is for absolute Hamilton-count reads; the two-apex argument does not recover the baseline from signed reversal gaps alone.

Write \(q_m(U)\) for ordered \(m\)-slot path covers allowing empty slots, and \(n=|T|\). The homogeneous probe \(C_3[U,I_{j-1},I'_{j-1}]\) with transitive exterior blocks reads a triangular combination of nonempty cover counts up to \(j\): each revisit uses at least one vertex from each exterior block, and alternating single-vertex visits attains \(j\). Its positive diagonal recovers the cover counts and hence \(q_m(U)\) for all positive \(m\). This is the previously established uniform-interface inversion; only finitely many wraps are needed for a finite old block. All such wraps are uniform toward the original \(T\), including when \(U\) already contains the individualized vertices.

## 1. Audit of the one-apex converse criterion

Define \(D_m(v)\) as the number of ordered \(m\)-covers with the distinguished first path starting at \(v\) minus those with that path ending at \(v\). For an apex \(x\) beating \(S\),
\[
 q_m(T+x_S)-q_m(T^{op}+x_S)=
 m\sum_{v\in S}\bigl(D_m(v)+D_{m+1}(v)\bigr).
\]
Indeed fix the slot containing \(x\). An isolated \(x\) does not contribute to the difference. If \(x\) borders exactly one old path in that slot, the two possibilities contribute \(2\sum_{S}D_m\), using \(\sum_vD_m(v)=0\). If it borders two old paths, split the slot into two named nonempty slots. Let \(K(u,v)\) count covers with the first ending at \(u\) and second starting at \(v\); reversing each path and exchanging these two named slots transposes \(K\) in the converse. Its cross-cut antisymmetric sum is its column-minus-row sum on \(S\), equal to \(\sum_SD_{m+1}-D_m\) because the other distinguished slot is required nonempty. Multiply by the \(m\) slot choices.

Each \(D_m(v)\) is polynomial in \(m-1\) of degree at most \(n-2\): split off the distinguished path's vertex set, with the singleton term zero, and evaluate the remaining \(q_{m-1}\) by selecting occupied slots. Singleton apex probes and uniform wraps therefore imply \(D_m(v)+D_{m+1}(v)=0\); when its left side vanishes at \(n-1\) positive integers it vanishes identically, and a nonzero polynomial cannot satisfy \(D(m+1)=-D(m)\) (compare leading coefficients). Conversely, expand every host path by the outside itinerary and its occupied old runs. With one marked external vertex, the only nontrivial reads are the start/end difference of one distinguished nonempty run or the skew cross-cut kernel of its two neighboring runs. The same transpose argument applies; these differences are the coefficients of \(\binom{m-1}{r-1}\) in \(D_m(v)\). Thus zero divergence is sufficient for all hosts with only one individualized vertex, not just for the homogeneous wraps used in necessity. The cases \(n\le1\) are trivial.

In particular zero divergence implies regularity: the leading coefficient of \(D_m(v)\) is the signed outdegree of \(v\), as the order-two distinguished paths give all degree-\((n-2)\) terms. Regularity alone is not sufficient: the explicit seven-vertex tournament with outgoing rows \((74,28,25,112,97,7,38)\) is regular but has \(D_1=(0,-1,1,-1,1,1,-1)\). Vertex-transitivity *is* sufficient because each start/end statistic is separately constant across vertices, with equal totals. The 21-vertex chiral Cayley core from the one-apex account gives a nonisomorphic collision; even the directed triangle gives a labelled converse collision.

## 2. A three-contact Fourier probe with two individualized vertices

For any completed probe \(U\) of order \(N\), the [odd-cycle collection identity](https://arxiv.org/abs/2307.05569) gives
\[
 q_m(U)=\sum_{\mathcal C}2^{|\mathcal C|}\,
 m^{N-\sum_{C\in\mathcal C}(|C|-1)}.
\]
In the complementary-subset convolution, a packed cycle is wholly assigned to one slot, as is each uncovered vertex; the cycle gives weight two in its slot. In particular the coefficient of \(m^{N-2}\) gives twice the triangle count, and the coefficient of \(m^{N-4}\) gives
\[
 L(U)=2c_5(U)+4p_{3,3}(U),
\]
where \(p_{3,3}\) counts vertex-disjoint unordered pairs of directed triangles. Both coefficients are recovered by finite integer interpolation from the uniform probes above.

To recover the signed outdegree \(d_j=\sum_{k\ne j}s_{jk}\), first read the triangle count of \(T\). Then add a single apex beating only \(j\); its new triangles correspond to arcs \(j\to k\) for \(k\in T-j\). This reads the ordinary outdegree and hence \(d_j\).

Fix distinct old \(i,j\). Add two individualized vertices \(x,y\), fixing all their contacts except the three adjacencies of the scaffold
\[
 y\longrightarrow i\longrightarrow x\longrightarrow j .
\]
Use \(z_e=\pm1\) for each of these three arcs measured in the scaffold direction. Denote the signed sum of readouts \(L\) for the eight contact patterns, with multiplier \(\prod z_e\), by \(\mathcal A_{ij}(T)\). A contribution surviving this sum must use every tested adjacency, by cancellation at an unused bit. No pair of disjoint triangles can contain a connected four-vertex path: one triangle cannot contain four vertices, and two disjoint ones cannot split a connected path. Thus only 5-cycles survive. They consist of the scaffold, read forwards or backwards, and one extra old vertex \(k\notin\{i,j\}\) joining its ends. For one such \(k\) the signed contribution is
\[
 2(A_{jk}A_{ky}-A_{kj}A_{yk})=s_{jk}+s_{ky},
\]
where \(A_{uv}=[u\to v]\), \(s_{uv}=2A_{uv}-1\); the minus sign is from the odd number of toggled edges on reversal. Therefore
\[
 \mathcal A_{ij}(T)-\sum_{k\notin\{i,j\}}s_{ky}=d_j-s_{ji}.
\]
The calibration sum uses only known exterior contacts. This recovers \(s_{ji}\) for every old edge, with eight contact patterns per pair plus a finite family of homogeneous height probes. The order-two case already has the degree read.

In particular, for a regular old core, no calibration or absolute degree probe is necessary to distinguish it from the converse: the alternating sum difference is \(-2s_{ji}\ne0\). Every old core whose converse escaped even the strongest one-apex test is regular by §1. Hence it is caught by some two-apex read, including the specified nonisomorphic Cayley pair. The required threshold is precisely two **individualized** vertices, not necessarily two total added vertices; the extraction of polynomial coefficients uses further homogeneous probes.

In the chiral Cayley example, take (i=0,j=1), let (x	o y), and make every other untested old vertex beat both new vertices. In increasing three-bit mask order, the differences in (L) (core minus converse) are (0,0,0,0,0,0,0,2), by direct 3-/5-cycle enumeration. Thus the all-forward scaffold, where (x) beats only old label 1 and (y) beats only old label 0, already has different feature values. The implication from the feature alone requires uniform readout shells. For this particular core, the [separate whole-order DP](direction-height.out) also verifies a difference in the unwrapped counts.

The [replay](two-apex-divergence-check.py) verifies 24 calibrated two-apex probes, 1,008 one-apex convolution identities, and this chiral example.

Quantitatively the displayed powers of (m) all have parity (N), with smallest exponent at least 1. Thus (q_m(U)=m^epsilon P(m^2)), with (epsilon=1) for odd (N) and 2 for even (N), and only (lceil N/2ceil) positive evaluations are needed. The largest triangular calibration therefore uses at most (2(lceil N/2ceil-1)) additional uniform exterior vertices. Eight three-bit patterns per old pair suffice; there is no exponentially long scaffold in this two-apex argument.
