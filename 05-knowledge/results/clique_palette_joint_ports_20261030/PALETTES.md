# All-demand Paley7 palettes and a seven-vertex integer hole

**Status: PROVED, using the exact certificate specified below; computational instances FINITE-EXACT.** Local labels here are packet identifiers, not reserved global theorem IDs.

## Definitions and inherited mechanism

A tournament has a specified orientation on every unordered pair, with no ties. Write \(\alpha(T)\) for its largest transitive induced subtournament and \(\chi(T)\) for its dichromatic number. A demand vector d assigns \(d_v\) distinct colors to vertex v; each color's support must induce a transitive tournament. Let \(\chi_T(d)\) be the minimum number of colors, and \(\chi_k(T)=\chi_T(k\mathbf1)\). Supports may repeat: they are then different colors.

A tournament is transitive iff it has no cyclic triple. This follows by taking a shortest directed cycle: any cycle of length at least four has a chord shortening it. Thus the complete transitive-support system is the independence system of its cyclic-triple hypergraph. It need not be an ordinary graph independence system.

The October 21 palette packet already proves \(\chi_k(P_7)=\lceil7k/3\rceil\) and the uniform substitution rule. We extend the *unequal*, named-port demand problem. The general substitution principle also has prior literature; see SOURCES.md.

**P1 (nonuniform substitution).** If \(T=Q[T_1,\ldots,T_s]\), with every cross-block arc oriented as Q, then
\[
\chi(T)=\chi_Q\bigl(\chi(T_1),\ldots,\chi(T_s)\bigr),\qquad
\alpha(T)=\max_{S\ {\rm transitive\ in}\ Q}\sum_{i\in S}\alpha(T_i).
\]
Empty blocks have both invariants zero. Restrict any coloring to each block: its distinct local colors are at least \(\chi(T_i)\), and their quotient support is transitive. Delete excess local color occurrences if needed. Conversely, choose a quotient palette and inject each block's coloring into its allocated colors. The union of local transitive pieces over a transitive quotient is transitive. The same two directions prove the alpha formula. The identical argument gives the k-fold version, with child demands \(\chi_k(T_i)\). Taking rational multiples gives the fractional version.

## The exact unequal-demand law

Set \(P=P_7\), vertices 0–6, with \(i\to j\) iff \(j-i\bmod7\in\{1,2,4\}\). There are 14 cyclic triples and 50 transitive supports including empty: 1+7+21+21.

**P2.** For every \(d\in\mathbb Z_{\ge0}^7\),
\[
\boxed{\chi_P(d)=
\max\left\{\max_i d_i,\quad
\max_{C\ {\rm cyclic}}\left\lceil{\sum_{i\in C}d_i\over2}\right\rceil,\
\left\lceil{\sum_i d_i\over3}\right\rceil\right\}.}
\]
Each bound is necessary: a color uses a vertex at most once, a cyclic triple at most twice, and P at most three times.

For sufficiency define
\[
Q=\{x\ge0:x_i\le1,\ \sum_{i\in C}x_i\le2\ (14\ C),\
\sum_i x_i\le3\}.
\]
The supplied certificate proves both
\[
Q=\operatorname{conv}\{\mathbf1_S:S\text{ transitive in }P\}
\quad\text{and}\quad
qQ\cap\mathbb Z^7=\underbrace{(Q\cap\mathbb Z^7)+\cdots+(Q\cap\mathbb Z^7)}_{q}
\]
for every integer q≥1. This is the integer decomposition property, not just a fractional relaxation.

### Why the finite geometry proves the all-size statement

The exact checker reconstructs all supports and all 29 inequalities itself. It checks 988 full-dimensional simplices with vertices among those 50 supports, all inside Q. Their facet incidences are:

| Data | Checked count |
|---|---:|
| Unpaired facets, on declared supporting hyperplanes | 2,364 |
| Paired facets, with opposite-side determinant signs | 2,770 |
| Simplices of absolute determinant 1 / 2 / 3 | 855 / 131 / 2 |
| Nonzero fundamental lattice points, counted per simplex | 135 |
| Their heights 2 / 3 / 4 | 102 / 32 / 1 |

Every defining inequality has a full six-dimensional boundary-face witness.

Let U be the finite union of these closed simplices. If U did not cover Q, take a generic segment from the interior of a simplex to a point of the open complement inside Q. Avoid the codimension-two face skeleton. At its first exit from U the segment crosses a simplex facet. An unpaired facet lies on the boundary of Q and cannot be crossed inside Q. A paired facet has a full simplex on its other side, so cannot be an exit. This contradiction proves U=Q. The argument only requires a cover; it does not rely on an unverified claim that interiors are disjoint. All simplex vertices lie in the displayed support hull, so the equality of polytopes follows.

For a simplex let A have columns \((1,\mathbf1_S)\). Express an integer point \((q,d)\) in its cone with nonnegative coefficients. Subtract their floors. The remainder is an integer point in the half-open fundamental parallelepiped. Absolute determinant one leaves only zero. For determinant p=2 or 3, the adjugate identity makes every coefficient a multiple of 1/p; the checker exhausts all residues in \(\{0,\ldots,p-1\}^8\) satisfying \(Ar=0\bmod p\). It verifies exactly p residue classes and a decomposition of every nonzero remainder into legal supports at precisely its integer height. Adding the floor copies proves integer decomposition for every q. For q=0, d=0 is immediate.

This proves P2. It also proves that the fractional weighted cost is the same maximum with all ceilings removed. The certificate is checked with exact Fraction determinants independently of the producer's Bareiss determinants; exploratory floating-point triangulation is not a proof dependency.

**Computational controls.** An independent demand dynamic program agrees for all \(4^7=16,384\) vectors in \(\{0,1,2,3\}^7\). Literal substitution graphs agree with P1 and P2 in all 64 cases obtained by replacing up to three P7 ports by C3 and the remaining ports by singletons.

## A small hole that survives every separate prime test

**P3.** Let H have bit rows
\[
[74,100,25,50,35,5,60],
\]
where bit j in row i means \(i\to j\). Then
\[
\alpha(H)=4,\quad \chi_f(H)=2,\quad \chi(H)=3,\qquad
\chi_k(H)=2k\quad(k\ge2).
\]

Here are explicit witnesses, with a mask denoting a transitive support.

| Demand per vertex | Palette masks | Cost |
|---|---|---:|
| 1 | 1,6,120 | 3 |
| 2 | 11,21,102,120 | 4 |
| 3 | 11,19,37,92,102,120 | 6 |

Every row is checked for exact demand and legality. Exhausting all complementary subsets proves that no two-coloring exists. The dual weights
\[
(1/2,1/4,1/4,1/4,1/4,1/4,1/4)
\]
sum to 2 and have sum at most 1 on every transitive support, proving the lower bound 2k. Supports 37,58,67,92 each with weight 1/2 give a fractional coloring of total 2.

Every k≥2 belongs to the additive semigroup generated by 2 and 3. Repeat the displayed two-fold palette when k is even; use one three-fold palette and two-fold repeats when k is odd. This proves the all-k formula. The finite checker also replays the actual palettes through k=64.

The maximum of the simpler bounds \(|S|/\alpha(H[S])\), over nonempty subsets S, is only 7/4. Thus three distinct quantities occur:
\[
\max_S |S|/\alpha(H[S])=7/4<\chi_f(H)=2<\chi(H)=3.
\]
THM-1850 already records a seven-vertex \(\alpha=4,\chi=3\) counterexample to the naive density formula. We do not claim to discover that phenomenon or a new isomorphism class. The additions here are an explicit audited witness, its full demand ray and its recursion.

**P4 (amplification).** Let \(H^{[r]}\) be the r-fold lexicographic power, r≥1. Then
\[
|H^{[r]}|=7^r,\quad \alpha(H^{[r]})=4^r,\quad
\chi_f(H^{[r]})=2^r,\quad
\chi(H^{[r]})=3\,2^{r-1}.
\]
Moreover \(\chi_k(H^{[r]})=2^rk\) for every k≥2. Apply P1 and P3 inductively; fractional multiplicativity also follows from the limit \(\chi_k/k\). The additive integer/fractional gap is \(2^{r-1}\). In particular no absolute C can make
\[
\chi(T)\le C\,|T|/\alpha(T)
\]
hold for every tournament: these examples have ratio \((3/2)(8/7)^r\).

At fractional budget \(2^r\), the feasible scaling degrees include 2 and 3 but exclude 1. Thus the minimum p-adic exponent among feasible denominators is zero for *every* prime p, yet an integral optimal palette does not exist. This is a native realization of the earlier mixed-2/3-clock warning. Local prime denominators describe saturation data; they do not certify membership in the original nonnegative semigroup. No median/XOR descent theorem can be imported without proving its relation-preservation hypotheses, which have not been established for these palette constraints.

**Boundary.** P2 is a special normality theorem for P7. H disproves its indiscriminate extension to all seven-vertex tournaments. Neither result says that fractional feasibility solves general tournament coloring.
