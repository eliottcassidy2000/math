# Rooted graph costs, even-star tomography, and finite probability patching

Date: 2026-10-20. Status: **PROVED** for the elementary finite statements below;
**FINITE-EXACT** for the declared controls. No new theorem identifier is reserved.
The attached *Blueprint: countable colorings of linear graphs* is an attributed
source blueprint, not a Lean build or an independently certified theorem here.
The attachment attributes the argument to Plasma AI and the prose write-up to
Claude Opus 5.5. Its SHA-256 is
`a546357d37b4f99966592c6b9eff303d35c2171288af5ff8a503cf67ccfd9b0a`.
Its whole text, including closure chains, derived signs, minimax, and the final
Hilbert-space colouring step, was read. The results below do not require its
countable-colouring headline as a dependency.

## 1. Inheritance, novelty boundary, and the concept board

The attachment defines the rooted extension cost: after signs on roots R have
been fixed, minimize the number of monochromatic edges by choosing signs on a
finite interior S. Its invariant J bounds the expectation of each such cost by
C|S|. The final trace energy in the blueprint, in contrast, uses only pairwise
agreement probabilities.

Three recovered repository mechanisms make that difference useful:

- **THM-2221**, `01-canon/theorems/THM-2221-tournament-context-cut-metric-and-pinned-transport-response.md`:
  pairwise cut data are exactly sufficient for pinned exterior Hamming
  transport, but discard higher correlations. Its full core kernel remains a
  sidecar. We do not replace a symmetric boundary cost with a tournament.
- **THM-3396**, `01-canon/theorems/THM-3396-four-bit-pairwise-independent-fourier-cone.md`:
  the four-bit parity laws and their Gram blindness are already established.
  The parity example itself is inherited, not claimed as a new discovery.
- **Boundary palette compression**, `05-knowledge/results/boundary_palette_compression_20261011.md`:
  complete named-port equality relations glue exactly and retain colouring
  lifts. Here their zero-cost relation is upgraded to a full min-plus cost
  table. Positivity of listed lifts alone still does not prove completeness.

The new connection is an exact completion of the pairwise observable: even-leaf
star extension costs form a basis of every globally sign-flip-invariant function
on a fixed finite boundary. Thus they reconstruct precisely the part of a law
that any unpinned rooted bipartization cost can observe. The finite replication
lemma below also isolates a quantitative remnant of the blueprint's closure
argument, with a five-vertex path showing why one copy is insufficient.

Anchor / niche / wildcard: finite J observability; the older pinned-context and
palette interfaces; finite closure-copy budgets. Concept board: root law,
cut matrix, interior elimination, star parity, named ports, replication budget.
Used research move: retain shared labels and products before scalar totals.
No literature-priority claim is made for the general Fourier or minimax tools.

## 2. The exact boundary-cost interface

Let F=(R disjoint union S,E) be a finite simple graph. Write x_v in {−1,+1},
and define

\[
 f_F(x_R)=\min_{x_S}\sum_{uv\in E}\frac{1+x_ux_v}{2}.
 \tag{1}
\]

The attachment restricts every edge to meet S. We permit root-root edges in
this section as well; those contribute known boundary terms. All claims in
Sections 3–5 apply in either convention.

**Proposition 1, exact operations.**

1. Global sign reversal preserves the table: f_F(−x_R)=f_F(x_R).
2. If F and H have exactly the same named root set, disjoint interiors, and
   disjoint edge sets, then f_{F union H}=f_F+f_H pointwise.
3. Forgetting a root r gives f_new(x_{R\{r}})=min_{x_r} f_F(x_R).
4. Adding a root edge uv increments the table by (1+x_ux_v)/2.
5. The zero set of f_F is exactly the proper two-colour boundary relation.

**Proof.** Reversing all interior signs with all roots preserves every edge
summand. Disjoint interiors are independent minimization variables, so the
minimum of the sum is the sum of the minima. Forgetting introduces one more
variable to minimize. The last two statements follow directly from (1). QED.

These operations generalize to a fixed palette q>=1, with edge cost
1[colour(u)=colour(v)]. Quotient the named boundary assignments by permutation
of colour names, retaining their equality partition with at most q blocks.
The same add/min rules hold, and the zero-support table is exactly the existing
Sigma_q interface. Internal geometry and reconstruction witnesses are not
encoded by the cost table.

**Two load-bearing guards.** If an edge lies in both fragments it is counted
once in their graph union. The correct rule subtracts its known boundary cost;
it does not add that cost twice. If an interior vertex is shared, additivity
can fail even without shared edges: with roots u,v of opposite signs, the
fragments u-y and y-v each have cost zero, but their union has cost one. The
missing datum is the shared sign y. Promote it to a named port before gluing.

Copies therefore obey a precise law. Glue m copies of F only along their named
roots, with private interiors and owned edges; then their table is m f_F.
If W is a further shared set subsequently forgotten, the table is

\[
 \min_{x_W}\sum_i m_i f_i(x_M,x_W),
 \tag{2}
\]

not the sum of separately minimized f_i. Equation (2) is the finite interface
behind the blueprint's replicated patterns.

## 3. The exact three-root compression boundary

For A subset R let chi_A(x)=product_{v in A}x_v. Every real boundary function
has the Walsh expansion

\[
 f(x)=\sum_{A\subseteq R}\widehat f(A)\chi_A(x),\qquad
 \widehat f(A)=2^{-|R|}\sum_x f(x)\chi_A(x).
 \tag{3}
\]

**Theorem 2.** Every rooted cost (1) has zero coefficients at odd |A|. Hence:

- with at most three roots, its expectation under any law is determined by
  the pairwise correlations E[x_ux_v], equivalently the pairwise cut matrix;
- more generally this remains true with arbitrarily many roots if each
  connected component of the interior graph F[S] meets at most three roots;
- four roots are the first possible obstruction, already for an induced
  unit-distance star.

**Proof.** Substitute x -> −x in (3); sign invariance kills odd coefficients.
With at most three roots, only degrees zero and two remain. For the component
claim, minimize independently over the components of F[S]. Each summand sees
at most three roots and is therefore quadratic; root-root edge terms are also
quadratic. The four-star example below proves sharpness. QED.

The component condition is sufficient, not necessary. A larger boundary may
still have vanishing higher coefficients because of its actual cost table.
The exact test for pairwise sufficiency of one fixed table is
widehat f(A)=0 for every |A|>=4. It is a statement for all probability laws,
not merely one chosen law: perturbing the uniform law by a small character
with a nonzero coefficient supplies a hostile if the condition fails.

For the four-leaf star with one free centre,

\[
 g_4(x)=\min(\#\{x_i=+1\},\#\{x_i=-1\})
 =\frac54-\frac14\sum_{i<j}x_ix_j+\frac14x_1x_2x_3x_4.
 \tag{4}
\]

Let mu_+ and mu_- be uniform on the eight words of even and odd minus-sign
parity. These are the two quartic half-cube controls already in THM-3396.
They have identical marginals on every proper subset of the four roots,
including identical pairwise disagreement 1/2. But

\[
 E_{\mu_+}g_4=3/2,\qquad E_{\mu_-}g_4=1.
 \tag{5}
\]

Thus a single pattern's J test at, for example, C=5/4 and |S|=1 distinguishes
laws invisible to the entire pairwise cut matrix. This is a statement about
that test, not a claim that either law satisfies all of J for an arbitrary
graph. At fixed pairwise correlations the omitted quartic is exactly

\[
 E[x_1x_2x_3x_4]
 =4E g_4-5+\sum_{i<j}E[x_ix_j].
 \tag{6}
\]

The geometry is literal: centre (0,0), leaves (1,0),(0,1),(−1,0),(0,−1).
All four centre-leaf distances are one; no leaf-leaf distance is one. Hence
this loss boundary occurs inside exact unit-distance graphs, without claiming
any new Hadwiger–Nelson chromatic bound.

## 4. Even stars are a complete boundary-law instrument

For any nonempty even set T subset R, let

\[
 g_T(x)=\min(\#\{v\in T:x_v=+1\},\#\{v\in T:x_v=-1\}).
 \tag{7}
\]

This is the extension cost of one new centre adjacent precisely to T.

**Theorem 3, star basis and tomography.** The functions

\[
 \{1\}\cup\{g_T: \varnothing\ne T\subseteq R,\ |T|\text{ even}\}
 \tag{8}
\]

form a basis of all functions f on {−1,+1}^R satisfying f(−x)=f(x).
Consequently the expectations E_mu g_T determine exactly the symmetrized law

\[
 \mu^{\rm sym}(x)=\tfrac12(\mu(x)+\mu(-x)),
 \tag{9}
\]

and hence the expected cost of every finite unpinned-root bipartization gadget
on R, regardless of its interior size. For r=|R|>=1, precisely 2^(r−1)−1 star
measurements are used in addition to total mass one. They are linearly
independent; no one can be removed from a globally complete linear instrument.
For the empty boundary only the constant remains.

**Proof and effective inverse.** A star on T depends only on T and has only
even Fourier coefficients. If |T|=2m, its coefficient at the top character is

\[
 \widehat g_T(T)
 =2^{-2m}\sum_{k=0}^{2m}(-1)^k\binom{2m}{k}\min(k,2m-k)
 =\frac{(-1)^m}{2^{2m-1}}\binom{2m-2}{m-1}\ne0.
 \tag{10}
\]

For completeness, split the sum at m, use symmetry, and use
k binom(2m,k)=2m binom(2m−1,k−1). The alternating partial-binomial identity
sum_{j=0}^h (−1)^j binom(N,j)=(−1)^h binom(N−1,h)
gives the middle sum as (−1)^m 2 binom(2m−2,m−1), including m=1 directly.
Ordering even subsets by size makes the change from stars to even characters
triangular with nonzero diagonal, proving (8).

Let m_A=E_mu chi_A and m_empty=1. If z_T=E_mu g_T, recursively recover

\[
 m_T=\frac{z_T-
   \sum_{A\subsetneq T,\ |A|\text{ even}}\widehat g_T(A)m_A}
  {\widehat g_T(T)}.
 \tag{11}
\]

Walsh inversion then gives

\[
 \mu^{\rm sym}(x)=2^{-r}\sum_{A\subseteq R,\ |A|\text{ even}}m_A\chi_A(x).
 \tag{12}
\]

It cannot recover odd moments: every probe has exactly the same value at x
and −x. Conversely (12) recovers all remaining information, so the loss is
exactly global sign orientation. QED.

This is not a claim that a given ambient graph supplies all required star
centres, or that arbitrary unit-distance multi-star gluings are realizable.
The theorem identifies the required observable; physical availability of each
probe is a separate sidecar. The linear expansion may have negative
coefficients, so upper bounds on star costs alone do not automatically give
upper bounds on every gadget cost. A sharp hostile already has two roots:
the two-edge path has cost g_{a,b}=1[x_a differs from x_b], while the
three-edge path has cost 1−g_{a,b}. At deterministic equal endpoint signs,
every star on this two-root boundary costs zero but the three-edge path
costs one. Complete numerical observations are not the same thing as
one-sided inequalities on those observations. The tomography theorem is
specific to two signs; the q-colour interface above makes no analogous
star-completeness claim for q>=3.

**Reading the old four-bit cone through J.** In THM-3396's unbiased,
pairwise-independent cone, (4) gives E g_4=5/4+d/4. Its four cubic parameters
a_i disappear because unpinned-root costs are globally sign-invariant. Thus
this consumer needs only the quartic d after whitening; the old five-coordinate
sidecar is necessary for the full law but is larger than this consumer needs.
This is an exact operation-specific compression, not a contradiction of the
older minimality theorem.

## 5. A finite replication-to-patching theorem

Fix finite disjoint sets M,W, a probability law mu on signs of M, and r>=1
nonnegative cost tables f_i(x_M,x_W). Let s_i>=0 be the interior sizes (or
specified nonnegative charges), C>=0, and an integer copy scale k>=1. Suppose
for every integer vector

\[
 0\le m_i\le k,\qquad \sum_i m_i\le k+r-1,
 \tag{13}
\]

we have the finite aggregate test

\[
 E_\mu\min_{x_W}\sum_i m_i f_i(x_M,x_W)
 \le C\left(|W|+\sum_i m_i s_i\right).
 \tag{14}
\]

**Theorem 4, quantitative finite successor.** There is a joint law lambda on
M union W preserving mu exactly and satisfying every bound

\[
 E_\lambda f_i\le Cs_i+\epsilon_k,\qquad
 \epsilon_k=\frac{C}{k}\left(|W|+\sum_i s_i-\min_i s_i\right).
 \tag{15}
\]

The extension can be chosen as a mixture of at most r+1 deterministic kernels
x_W=kappa(x_M). Graph-theoretically, hypothesis (14) follows when the relevant
private pattern copies exist and the old rooted law satisfies J for those
aggregate patterns. No infinite closure or forced-sign assumption is made.

**Proof.** Let rho_i>=0 have sum one and set m_i=ceil(k rho_i). Then (13)
holds. Since the f_i are nonnegative,

\[
 k E_\mu\min_{x_W}\sum_i\rho_i f_i
 \le E_\mu\min_{x_W}\sum_i m_i f_i.
 \tag{16}
\]

Put delta_i=m_i−k rho_i. Each delta_i lies in [0,1), and their sum is an
integer at most r−1. Therefore, for nonnegative s_i,

\[
 \sum_i\delta_i s_i\le\sum_i s_i-\min_i s_i.
 \tag{17}
\]

Combine (14)–(17) and divide by k. For every rho, the minimum expected excess
sum rho_i(f_i−Cs_i), over deterministic kernels, is at most epsilon_k.
The minimum can be selected separately for each finite root word, so this is
an ordinary finite zero-sum game. Finite minimax, or separation of its cost
polytope from the lower orthant, gives a mixture with every excess at most
epsilon_k. Every kernel leaves the law of M unchanged. Caratheodory's theorem
in the r-dimensional cost space reduces the mixture to at most r+1 kernels.
QED.

The theorem gives a finite stopping criterion: for a desired tolerance epsilon,
a copy scale k>=C(|W|+sum s_i−min s_i)/epsilon suffices if all its specified
aggregate tests and copy witnesses have actually been verified. It supplies
neither those copies nor a finite chromatic bound automatically.

### Why one copy does not give an exact successor

Take the five-vertex path

\[
 a\;—\;s\;—\;w\;—\;t\;—\;b,
 \tag{18}
\]

with M={a,b}, signs x_a=+1,x_b=−1 fixed, and C=1/3. Every pattern over M
has expected extension cost <=C|S|: a nonzero cost is possible only when all
four path edges and all three interior vertices are present, in which case
it is exactly one. This also holds if one or both endpoint roots are omitted.

Now set N=M union {w}. Its two rooted length-two path patterns each have one
private vertex. Their costs are complementary indicators: one is
1[x_w differs from x_a], the other 1[x_w differs from x_b]. Their sum is
identically one, so no law extending the fixed endpoint signs can make both
expected costs at most 1/3.

The one-copy aggregate tests nevertheless pass:
min(m_1,m_2)<= (1+m_1+m_2)/3 for m_i in {0,1}. Two copies of each would expose
the failure: 2>5/3. The exact failed implication is “J over M persists over an
arbitrary intermediate N.” The blueprint's closure/replication requirement
is substantive, and (15) records an honest finite replacement with slack.

## 6. The connection contracts and next targets

**Cut contexts -> rooted costs.** Send a sign law to its pairwise cut matrix.
This preserves THM-2221's exterior Hamming response and all three-port-component
cost expectations. It destroys higher even moments needed after minimizing
an interior sign. Even-star expectations restore exactly the lost observable
law; the four-star is the cheapest decisive test.

**Palette relations -> quantitative gluing.** Send a complete colour boundary
relation to the zero-support of its min-plus cost table. The full table retains
bipartization defects, exact copying, and forgetting. Passing back to zero
support loses positive costs and multiplicities. Named shared ports, edge
ownership, geometric realization, and lifts must be retained.

**Closure copies -> finite probability patching.** Send copy counts to the
weighted finite cost game (14). It preserves all old-root marginals and buys
simultaneous inequalities with tariff (15). It loses the infinite-trace forcing
and transfinite coverage needed by the source theorem. The five-path hostile
prevents claiming exact patching from bounded capacity alone.

Two concrete **OPEN** research targets remain worthwhile:

1. For an actual geometric gadget bank, classify which even-star measurements
   can be physically implemented while retaining distinct vertices and all
   cross-unit edges. Abstract completeness does not settle this realization
   problem. A useful first target is the four-port quartic channel on an
   already verified bank, before higher arities.
2. Replace the generic tariff in (15) by the exact support function of the
   available copy-count polytope. This is a finite optimization problem once
   the bank and its geometric copy witnesses are supplied. A theorem about
   that bank would be a quantitative assembly tool, not a global plane bound.

No inference here settles LRC(14), Hadwiger–Nelson, tournament H>=disc, or the
source's countable-colouring theorem. In particular, on a finite graph,
“countably colourable” by itself has no quantitative finite-colour content.

## 7. Reproduction and exact scope

Run `python -B check_linear_graph_tools.py` and
`python -B -O check_linear_graph_tools.py` from this directory. All checks use
exact integer/rational arithmetic and remain active under optimization.
The saved output is `check_linear_graph_tools.out`.

The controls cover:

- top star coefficient by independent direct-cube and binomial routes for
  even sizes 2 through 16;
- reconstruction of all even moments and symmetrized atoms for explicit
  nonsymmetric rational laws on 1 through 8 roots;
- every edge subset touching the interior, every prefix-root split, on
  0 through 5 vertices; relabelling covers arbitrary choices of named roots;
- complement symmetry, zero-support against independent proper-colouring
  search, the component-three-port criterion, and reconstruction of every
  tested graph table in the star basis;
- all 64 pairs of three-root stars with disjoint private centres;
- all 2,025 two-constraint/two-strategy minimax problems with payoffs in
  {0,1,2} and half-integer targets from 0 through 2;
- every rooted edge subpattern of the five-path under the stated law;
- the sharp round-up tariff on the declared finite rational simplex grid;
- explicit shared-interior, shared-edge, parity-law, one-sided-star-bound,
  and finite-copy hostiles.

The files contain no Lean theorem or claimed kernel audit. The finite exact
controls audit the finite results, not the blueprint's infinitary proof.
