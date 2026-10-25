# Signed flags, effective monoid models, and arbitrarily deep observer cancellation

Status: **PROVED** identities and model constructions below; **FINITE-EXACT** independent recurrence/permutation checks. The Session58 observer and its Walsh basis are inherited, not new. This packet gives new-to-this-recovered-route deductions, without a literature-wide novelty claim. It does not prove an all-size tournament-class separation theorem, an F1 model of the tournament quotient, or a prize conjecture.

## 0. Inheritance and scope

Closest proved mechanism: Session58 `observer-profile-algebra/THEOREM.md`, §§1–3, diagonalizes the ordinary graph/complement observer through the signed-order coefficient

\[
 \alpha_X=\frac1{n!}\sum_{p\in S_n}\epsilon(p,X),\qquad
 \epsilon(p,X)=(-1)^{\#\{ij\in E(X):i<j,\ p(i)>p(j)\}}.
\]

A graph is **Royle-even** here precisely when this epsilon equals +1 for every graph automorphism. This is neither even vertex degree nor, in general, the sign of the induced permutation of the edge set. For a Royle-even graph, its signed Walsh orbit sum is nonzero. If its number of edges is even, this sum is converse-invariant; if additionally alpha is zero, it lies in the linear kernel of the Session58 observer. These inherited facts are the only tournament-completeness input.

Canonical hostile: ordinary graph isomorphism forgets the order needed by the graph/tournament chart, first at n=3. Corrected near miss: loss of linear/mixture information does not prove a collision of two individual tournament profiles. Least-used sidecar: a full flag retains a reference chamber and its inversion length, before setting q=1 or averaging the chamber away.

Live concept board: signed Reynolds coefficients; Schubert cells; extremal insertion; torus-rank filtration; reference-order descent. Portfolio: the assigned F1 bridge is the anchor; high-order cancellation is the niche; effective realization of a virtual class is the wildcard.

## 1. The signed Schubert polynomial and its honest geometry

For a graph X on the ordered vertices 0,...,n−1, set N=binom(n,2), and let ell(p) be the usual inversion length. Define

\[
 D_X(q)=\sum_{p\in S_n}\epsilon(p,X)q^{\ell(p)},\qquad
 Q_n(q)=\sum_{p\in S_n}q^{\ell(p)}=[n]_q!,
 \quad [k]_q=1+q+\cdots+q^{k-1}.
\]

Then

\[
 D_X(1)=n!\alpha_X,\qquad D_X(0)=1,\qquad
 q^N D_X(q^{-1})=(-1)^{|E(X)|}D_X(q). \tag{1}
\]

The last identity follows from p↦w0p, where w0 reverses all ranks. Every pair changes its inversion status, so length becomes N−ell and epsilon changes by (−1)^{|E|}. In particular D is never the zero polynomial. If |E| is even and D(1)=0, the vanishing order is a positive **even** integer. Indeed if D=(q−1)^sR with R(1)≠0, reciprocity gives R(1)=(−1)^sR(1). Odd edge count similarly forces odd vanishing order. This is stronger than merely checking that the first derivative vanishes.

Let

\[
 P_X^\pm(q)=\sum_{\epsilon(p,X)=\pm1}q^{\ell(p)}
             =\frac{Q_n(q)\pm D_X(q)}2. \tag{2}
\]

These are actual nonnegative integral polynomials, not formal fractions needing a positivity conjecture.

**Geometric realization.** The full flag scheme Flag_n over Z has Schubert cells C_p≅A^{ell(p)}. Give the *disjoint union scheme*

\[
 Z_X^\pm=\coprod_{\epsilon(p,X)=\pm1}\mathbb A^{\ell(p)}_{\mathbb Z}
\]

the obvious cellwise immersion into Flag_n. This is not a claim that the image is a single locally closed subscheme; an arbitrary collection of cells need only be constructible. Each Z_X^± is an actual disjoint union of affine spaces, has an affine torification, and has a literal monoid-scheme model as a disjoint union of affine monoid spaces. Its finite-field count is P_X^±(q). In the Grothendieck ring of varieties, with L=[A^1],

\[
 [Z_X^+]-[Z_X^-]=D_X(L). \tag{3}
\]

Thus D is initially a **virtual difference** of genuine geometric objects, not automatically the count of a variety. Equation (1) identifies alpha with a normalized *Euler/count-polynomial specialization* of that difference. It does not manufacture a field of cardinality one.

CITED background: López Peña–Lorscheid, *Torified varieties and their geometries over F1*, arXiv:0903.2173v3, Definition 1.5 and §1.3.5, Example 1.17, supplies the torification convention and C_p≅A^{ell(p)}. Their discussion explicitly warns that a Schubert torification of the full flag variety is not thereby known to be affine. Our disjoint-cell models avoid making that extra assertion. Source: https://arxiv.org/html/0903.2173v3 . The elementary polynomial, insertion, and model arguments here are supplied below.

## 2. Extremal insertion factorization

**Lemma.** If the smallest or largest labelled vertex v of a graph on k vertices is isolated, then

\[
 D_X(q)=[k]_q D_{X-v}(q).
\]

If it is universal (adjacent to every other vertex), then

\[
 D_X(q)=[k]_{-q}D_{X-v}(q). \tag{4}
\]

Proof: write a permutation as the ordered list of its vertices, equivalently p^{-1}; inversion length is unchanged by inversion. After fixing the order of the other k−1 vertices, insertion of an extremal v in the k available slots creates 0,1,...,k−1 ordinary inversions. It creates that same number of edge inversions if v is universal, and zero if v is isolated. Summing the slots proves (4). The extremal-label hypothesis cannot be silently omitted.

This gives a broad exact class: for any graph with a deletion sequence in which each currently extremal vertex is either isolated or universal, its signed polynomial factors into q-integers, with the sign choice recording universal versus isolated deletion. This is a labelled/decorated construction, not a theorem about a label-free polynomial invariant.

### The inherited connected family

Let F_r have vertices a,b,z_1,...,z_r,w and edges ab, az_i, bz_i, aw, in the original order a<b<z_1<...<z_r<w. Remove the universal minimum a, the now isolated maximum w, and the now universal minimum b. The r remaining vertices are independent. Hence for every r≥1,

\[
 \boxed{D_{F_r}(q)=[r]_q!\,[r+2]_q\,[r+1]_{-q}\,[r+3]_{-q}.} \tag{5}
\]

For odd r, both signed factors have even length, so q=1 has **exact multiplicity two**. Since

\[
 [2s]_{-q}=(1-q)[s]_{q^2},\qquad
 \left.\frac{[2s]_{-q}}{q-1}\right|_{q=1}=-s,
\]

we get, with n=r+3,

\[
 [ (q-1)^2 ]D_{F_r}(q)=n!/4,\qquad D_{F_r}''(1)=n!/2,
\]

and

\[
 \frac{D_{F_r}(q)}{Q_n(q)}=\frac14(q-1)^2+O((q-1)^3). \tag{6}
\]

For even r, D(1)=r!(r+2) is nonzero. The r=1 polynomial identity remains true, but F_1 has an odd automorphism; its signed Walsh invariant is zero, so it is not an observer-kernel witness. The inherited valid witness range is odd r≥3.

## 3. Connected Royle-even kernels of arbitrarily high jet order

Construct H_{r,j}, r≥2 and j≥1, by starting with r isolated vertices and successively adding

\[
 U,I,U,I,\ldots,I,U
\]

with j universal additions U and j−1 isolated additions I. An addition is universal or isolated relative to the vertices already present. Give each new vertex the largest label. Then

\[
 n=r+2j-1,\qquad |E|=j(r+j-1),
\]

and repeated use of (4) gives

\[
 \boxed{D_{H_{r,j}}(q)=[r]_q!\,
 \prod_{t=1}^j[r+2t-1]_{-q}\,
 \prod_{t=1}^{j-1}[r+2t]_q.} \tag{7}
\]

**Automorphisms and the inherited invariant.** The initial r vertices all have degree j. The isolated vertex added after the t-th universal vertex has final degree j−t, for 1≤t<j. The t-th universal vertex has degree r+j+t−2, for 1≤t≤j. For r≥2 these degree classes are pairwise distinct except for the initial r vertices. Thus Aut(H)=S_r, acting only on that initial independent block. Such permutations reverse no edge in the construction order, so H is Royle-even under the exact epsilon convention. The final universal vertex makes H connected.

Suppose now r≥3 is odd. Each of the j signed q-integers in (7) has even length, and all other factors are nonzero at 1. Therefore

\[
 \operatorname{ord}_{q=1}D_{H_{r,j}}=j,
 \quad
 [ (q-1)^j ]D_{H_{r,j}}=\frac{(-1)^j n!}{2^j},
\]

\[
 \frac{D_{H_{r,j}}(q)}{Q_n(q)}
   =(-1/2)^j(q-1)^j+O((q-1)^{j+1}). \tag{8}
\]

For even j≥2, the edge count is even, H is Royle-even, and alpha=0. Thus its nonzero signed Walsh invariant is an invisible converse-even mode of the Session58 observer. For every finite jet bound d, choose an even j>d: all derivatives D^{(k)}(1), k≤d, vanish, although D is nonzero. This disproves a uniform finite-jet detection claim **for these individual gauged signed polynomials**. It does not assert a rank theorem for a fully decorated observer matrix, nor identify derivatives with independent observer samples.

The original F_r is graph-isomorphic to H_{r,2}. Formula (5) is separately proved in its original label gauge; no label invariance is being assumed to transfer (7).

## 4. Effective monoid realization of the virtual difference

For odd r and even j, set

\[
 R_{r,j}(q)=[r]_q!\,
 \prod_{t=1}^j\left[\frac{r+2t-1}{2}\right]_{q^2}\,
 \prod_{t=1}^{j-1}[r+2t]_q\in\mathbb N[q].
\]

Then (7) becomes the effective factorization

\[
 \boxed{D_{H_{r,j}}(q)=(q-1)^jR_{r,j}(q).} \tag{9}
\]

Write R=Σ_d c_dq^d. Define the actual Z-scheme

\[
 Y_{r,j}=\coprod_d\coprod_{a=1}^{c_d}
                  (\mathbb G_m^j\times\mathbb A^d).
\]

Its monoid-scheme model is the same disjoint union of affine charts

\[
 \operatorname{Spec}\bigl(\mathbb F_1[x_1^{\pm1},\ldots,x_j^{\pm1},y_1,\ldots,y_d]\bigr),
\]

where this notation means the free pointed commutative monoid on the indicated invertible/noninvertible generators, and base extension is its monoid algebra over Z with the absorbing zero identified with 0. Thus this is a literal model in the pointed-monoid sense, not an inference from a polynomial evaluation. Directly,

\[
 \#Y_{r,j}(\mathbb F_q)=D_{H_{r,j}}(q),\qquad
 [Y_{r,j}]=[Z_H^+]-[Z_H^-]\quad\text{in }K_0(\mathrm{Var}_{\mathbb Z}). \tag{10}
\]

The equality in (10) is a scissor-ring identity, **not** an isomorphism between Y and a difference of flag loci; no geometric cancellation bijection is asserted.

The coordinate torification of Gm^j×A^d has binom(d,k) tori of dimension j+k. Consequently the coefficients of D(1+T) are nonnegative, the smallest torus dimension is exactly j, and the number of j-dimensional tori is

\[
 R_{r,j}(1)=n!/2^j.
\]

Thus the arbitrarily high q=1 zero has an exact geometric meaning: **minimum torus rank**, not absence of a geometric model. In particular its Euler/count-polynomial value at 1 is zero despite an explicit nonempty F1 model.

Crucial convention warning: N(1) is not being identified with the cardinality of Hom(Spec F1,Y) in the pointed-monoid category. Already Gm has N(q)=q−1, but a morphism from the trivial-unit pointed monoid sends its invertible generator to 1, giving one F1-valued point. The Euler specialization and a literal monoid-valued-points functor are distinct constructions. No statement here erases that distinction.

## 5. The reference-order averaging no-go

A genuine model for a signed count does not repair the order-loss in the tournament chart.

Let O be the set of n! linear orders on the labelled vertex set, b a reference order, and d(b,o) the number of pairwise disagreements (Kendall inversion length). For a tournament T and an unlabelled complementary-graph-pair type c, let B(T,o)=c mean that observing T in order o gives c. Define the reference-weighted profile

\[
 W_{T,b}(c;q)=\sum_{o\in O}q^{d(b,o)}\mathbf1_{B(T,o)=c}.
\]

Its total mass is Q_n(q), and at q=1 it is the inherited count profile A_T(c). For each fixed o, the relative orders b run once through S_n, so Σ_b q^{d(b,o)}=Q_n(q). Interchanging finite sums proves

\[
 \boxed{\frac1{n!}\sum_bW_{T,b}(c;q)
   =\frac{Q_n(q)}{n!}A_T(c).} \tag{11}
\]

For any q with Q_n(q)≠0 (in particular every real q>0 or formal q near 1), normalizing by Q_n gives exactly A_T/n! after uniform reference averaging. Hence that canonical reference-erasing operation preserves the entire original linear/mixture kernel at every q. It does not expose a formerly invisible invariant mode. It also does not create same-latent higher moments such as the two-observer repair in Session58.

An equivalent representation-theoretic check: for any finite-group representation rho and Reynolds projection P, the normalized weighted average R_w=Σ_g w_g rho(g)/Σ_g w_g satisfies PR_w=P. Thus a nonuniform weighted output only offers new information before projecting to invariant output. This assertion is independent of the choice w_g=q^{ell(g)}.

**Smallest gauge hostile.** The two isomorphic Royle-even two-edge stars on three ordered vertices give

\[
 D_{\{01,02\}}(q)=1+q^3,
 \qquad D_{\{01,12\}}(q)=1-2q-2q^2+q^3.
\]

They are neither equal nor negatives. At q=1 their values differ by the permissible signed-basis gauge, but that relation does not extend to q-polynomials. D therefore is not an unlabelled graph invariant even up to global sign. The q-refinement keeps extra reference data; it is not a deformation of the original invariant eigenvalue alpha without additional construction.

**Effectivity is gauge-specific, not just the polynomial.** Swapping vertices 1 and 2 in the original six-vertex F_3 (mask 255) gives the isomorphic mask 1599. Exact computation gives D_1599(2)=−1953 and quadratic Taylor coefficient −60. Thus that presentation cannot be a point-count polynomial of any actual variety at every finite field, despite the effective model in the original gauge. There is no contradiction: (9) retained a special construction order.

**Even the jet depth is gauge-specific.** H_{3,4} in its construction order has zero order four. Swapping vertices 0 and 8 gives zero order two, with quadratic Taylor coefficient 5760. The arbitrarily deep examples therefore disprove bounded-jet detection only when those specified gauges are retained; they do not say an optimal choice of gauge cannot reduce the needed depth. Both examples are checked exactly in the companion script.

**Exact scope of the no-go:** equation (11) concerns reference-independent inversion weights followed by *uniform* forgetting of the reference order. It does not rule out tournament-dependent weights, nonlinear correlated observations, intrinsic decorations, other F1 approaches, or retaining the reference order as part of the object.

## 6. Finite-exact certificates and hostile controls

`python3 verify.py` and `python3 -O verify.py` use only the standard library and explicit exceptions. The independent paths are (a) all permutations with rank-sign epsilon and (b) a last-vertex subset recurrence on permutation words. The latter obeys

\[
 D_S(q)=\sum_{v\in S}(-1)^{|\{u\in S:u>v,uv\in E\}|}
 q^{|\{u\in S:u>v\}|}D_{S\setminus\{v\}}(q).
\]

The exhaustive universe is all labelled simple graphs on 1≤n≤5, no connectivity filter. Checks compare both full polynomials, reciprocity, nonnegative parity-cell coefficients, and the even-edge zero-order parity. Targeted checks cover threshold families through n=12, the deep example H_{3,6} on n=14 (zero order six), original F_r gauges through r=8, and literal reference averaging for every labelled tournament through n=4.

For the two inherited order-six invisible masks, with lexicographic pair bits,

| mask | zero order at 1 | D''(1) |
|---|---:|---:|
| 255 | 2 | 360 |
| 2015 | 2 | 408 |

The second polynomial is exactly

\[
 1+q+2q^2+3q^3-q^4+q^5-4q^6-3q^7-3q^8-4q^9
   +q^{10}-q^{11}+3q^{12}+2q^{13}+q^{14}+q^{15},
\]

and factors as

\[
(q-1)^2(q+1)(q^2+1)(q^2-q+1)(q^2+q+1)
(q^6+2q^5+3q^4+5q^3+3q^2+2q+1).
\]

These are finite exact calculations, not a new census. Hostiles explicitly preserve r=1's invalid kernel status in the theorem, verify the n=3 gauge failure, the negative q=2 count and depth drop under relabelling, verify nonzero polynomial/zero Euler value, and verify that reference averaging recovers the old observer rather than a purported repair.

## 7. Map/loss ledger and boundaries

| Source → target | Map/preserved data | Lost information / required sidecar |
|---|---|---|
| Ordered graph X → signed flag-cell pair | p↦(epsilon(p,X),ell(p)); exact D | Does not retain all graph incidence; retain X and order for interpretation |
| Signed cell pair → D(1) | Euler/count specialization | Loses length grading; H has exact arbitrarily high cancellation order |
| Threshold signed difference → Y_{r,j} | Exact class and all finite-field counts | Only a K0 equality, not a flag-locus isomorphism or tournament decoding |
| Weighted profile with reference b → average over b | Equation (11), original invariant law | Loses the only extra reference information; old kernel returns |
| Monoid model → polynomial value at 1 | Euler/torus-rank shadow | Not the literal monoid Hom-point count |

No ordinary graph/tournament quotient bijection follows. No individual tournament-profile collision is exhibited. The unbounded jet result is about the fixed-gauge signed-polynomial family, not about all possible intrinsic q-deformations or arbitrary finite-dimensional nonlinear observations.
