# Occupancy powers, tail-marked trees, and exact path growth

**Status: PROVED.** The square-zero path resolvent, determinant identities,
directed Matrix–Tree theorem, tournament parity mechanism, and quadratic
directed-cut birth field are known or recovered here. The tail-marked
cofactor and its Boolean-lattice extraction are a self-contained application
of those mechanisms; no claim of literature priority is made. The checks in
[receipt.json](receipt.json) are **FINITE-EXACT**, with their universe stated
in §8. None of these formulas proves that every odd Hamilton-path total except
7 and 21 occurs.

## 1. Inheritance and the missing coordinate

The relevant inherited mechanisms are:

* [THM-1460, arborescences as the determinantal relaxation of Hamilton paths](../../../../01-canon/theorems/THM-1460-arborescences-are-the-determinantal-relaxation-of-H.md): a rooted Hamilton path is an out-arborescence with every out-degree at most one. Only this general structural statement and the standard Matrix–Tree theorem are used here; the finite extremal observations in that file are not promoted to theorems.
* [THM-1745, the leaf-graded arborescence filtration](../../../../01-canon/theorems/THM-1745-leaf-graded-arborescence-filtration-721-shadow.md): the one-leaf layer is exactly the Hamilton-path layer. Historical spectrum prose in that file is not a dependency.
* [Squarefree cycle bridges, SQ1–SQ2](../../squarefree_cycle_inventory_20261019/SQUAREFREE_CYCLE_BRIDGES.md): retain the identity of each consumed resource; vertex labels eliminate repeated vertices, whereas labels on some other objects eliminate a different kind of repetition. Naive anticommuting variables can cause unwanted cancellation.
* [Entry-resource machines, M2](../../entry_resources_20261113/MACHINES.md): the square-zero Hamilton-path resolvent already appears explicitly, with its exponential state cost.
* [Multiport defect grammar, §§1–2](../../multiport_defect_grammar_20261103/construction/RESULT.md): endpoint inventories multiply as disjoint path covers, and a birth splits at its new vertex. The empty-cover unit must be retained.
* [THM-4097, the exposed-slot ear kernel](../../../../01-canon/theorems/THM-4097-order-nine-strong-ear-spectrum-solid-interval-and-lane-extension.md), [THM-4114, §6](../../../../01-canon/theorems/THM-4114-ocf-mobius-positivity-tropical-defect-layer-and-opposite-ear-cut-curvature.md), and [THM-4128, §1](../../../../01-canon/theorems/THM-4128-johnson-slice-support-envelope-and-exposure-centrality-criterion.md): the quadratic birth field, submodularity, and the stronger exact directed-cut representation already exist in canon. [THM-4162, §1](../../../../01-canon/theorems/THM-4162-rooted-pair-mixed-two-ear-tensor-and-enumeration-free-johnson-cosets.md) records their coefficient-minimal endpoint interface. Section 7 recovers this old mechanism in the occupancy notation; it is not a new result.

The live comparison is between five objects: unrestricted walks, occupancy
states, endpoint inventories, branching trees, and Boolean-lattice
coefficients. Their common missing coordinate in an ordinary adjacency power
is **which vertices have already been used**.

Throughout, a graph is a finite loopless directed graph, with possible digons
unless explicitly called a tournament. Its adjacency convention is

\[
A_{uv}=1\quad\Longleftrightarrow\quad u\longrightarrow v.
\]

Let \(H(G)\) count directed Hamilton paths, with vertex order distinguished,
and let \(H_{r,t}(G)\) count those starting at \(r\) and ending at \(t\).
Singleton paths are allowed. Set \(H(\varnothing)=1\) only as the empty-cover
convention used in generating functions.

An ordinary power \((A^k)_{uv}\) counts walks with \(k\) edges. A Hamilton
path on \(n\) vertices has \(n-1\) edges, but fixing that exponent does not
remove repeated vertices. Two distinct failures are:

* \(TT_n\) has \(H=1\), while \(A^n=0\).
* For \(C_3\mathbin{\triangleright}TT_1\), \(H=3\) but
  \(\mathbf1^\top A^3\mathbf1=6\). For \(C_3\) plus two isolated vertices,
  \(\mathbf1^\top A^4\mathbf1=3\) although \(H=0\).

Thus the exponent needs an occupancy carrier, not just a length correction.

## 2. An exact matrix-power correction

Work in the commutative ring

\[
\mathcal R_V=\mathbb Z[x_v:v\in V]/(x_v^2:v\in V),\qquad
x_S=\prod_{v\in S}x_v.
\]

Its basis is indexed by subsets of \(V\), and
\(x_Sx_T=0\) exactly when \(S\cap T\ne\varnothing\). Put
\(X=\operatorname{diag}(x_v)\). Then

\[
\boxed{\quad
P_G(x)=1+\mathbf1^\top X(I-AX)^{-1}\mathbf1
      =\sum_{S\subseteq V}H(G[S])x_S.
\quad}                                                   \tag{2.1}
\]

Indeed, a term of \(\mathbf1^\top X(AX)^k\mathbf1\) is an oriented walk
weighted by all its visited vertices, including its starting vertex. It
vanishes if any vertex repeats; otherwise it is exactly a simple path.
The inverse is a finite geometric series because the positive-degree ideal
is nilpotent. In particular,

\[
H_{r,t}(G)=[x_V]\,x_r\big((AX)^{n-1}\big)_{r,t}.           \tag{2.2}
\]

The same known occupancy device also counts **directed Hamilton cycles**.
For \(n\ge2\), let \(\operatorname{HC}(G)\) count such cycles modulo
cyclic rotation, without identifying a cycle with its reverse. Then

\[
\boxed{\quad
\operatorname{HC}(G)=\frac1n[x_V]\operatorname{tr}((AX)^n).
\quad}                                                   \tag{2.3}
\]

A trace term is a closed walk with a specified starting vertex. Its \(n\)
arrival variables consume every vertex exactly once precisely when it is
a directed Hamilton cycle. Each such cycle has \(n\) possible starting
vertices, giving the divisor \(n\). Ordinary numerical
\(\operatorname{tr}(A^n)\) still includes closed walks with repeated
vertices; the squarefree coefficient is what removes them.

The path formula (2.2) is equally an ordinary integer matrix computation on the enlarged states
\((S,v)\), where \(\varnothing\ne S\subseteq V\) and \(v\in S\). Its
edges are

\[
(S,u)\longrightarrow(S\cup\{v\},v)
\quad\text{if }v\notin S\text{ and }A_{uv}=1.
\]

Starting at singleton states and taking \(n-1\) steps to full-set states
counts Hamilton paths. The lift is acyclic and has \(n2^{n-1}\) states.
It makes a “perfect fill” exact: every named vertex is consumed once, with
no added geometric assumptions.

The variables commute. If an exterior-algebra realization is desired, use
even pairs \(x_i=\theta_i\eta_i\), not single odd generators. Also,
setting all \(x_i=1\) is **not** a homomorphism from \(\mathcal R_V\) to
\(\mathbb Z\): it would send \(x_i^2=0\) to \(1=0\). Coefficient extraction
or the subset lift, rather than that substitution, performs the count.

## 3. Why the tournament hypothesis gives oddness

The matrix determinant lemma applied to
\(I-AX+\mathbf1(\mathbf1^\top X)\) gives, for every digraph,

\[
\boxed{\quad P_G(x)=
\frac{\det(I+(J-A)X)}{\det(I-AX)}.\quad}                    \tag{3.1}
\]

The denominator has constant term one, so is a unit in \(\mathcal R_V\).
For a tournament,

\[
A+A^\top=J-I,\qquad J-A=I+A^\top.
\]

Because \(X^2=0\),
\((I+A^\top X)(I+X)=I+X+A^\top X\). Transposition and
\(\det(I+UV)=\det(I+VU)\) give
\(\det(I+A^\top X)=\det(I+AX)\). Hence

\[
\boxed{\quad
P_G(x)=\prod_{v\in V}(1+x_v)
\frac{\det(I+AX)}{\det(I-AX)}.
\quad}                                                   \tag{3.2}
\]

Modulo two the numerator and denominator of the ratio agree. Therefore

\[
P_G(x)\equiv\prod_v(1+x_v)\pmod2.                          \tag{3.3}
\]

Every induced tournament has an odd number of Hamilton paths, including
the full tournament. The hypothesis is used at the exact step
\(J-A=I+A^\top\). This is a sufficient structure for a uniform oddness
theorem, **not** an iff characterization of individual odd-count digraphs:
a directed three-vertex chain has \(H=1\), a bidirected edge has \(H=2\),
and the empty graph on two vertices has \(H=0\).

Equations (2.1)–(3.3) belong to the established matrix-algebra route to
Hamilton paths; the proof is included to expose its retained variables and
hypotheses. Specifically, Irving–Omar, [*Revisiting The Rédei-Berge Symmetric
Functions via Matrix Algebra*, §2.2, Lemma 1](https://arxiv.org/html/2412.10572v3#S2.SS2)
gives the determinant quotient, and the same paper's Theorem 5 proof works
in \(\mathbb F_2[x]/(x_i^2)\) and reduces the path enumerator to
\(\det(I+X)\). Their convention uses \(XA\), equivalent by transposition
or the row/column version of the path resolvent. Primary text consulted
2026-11-13. These are recovered identities, not new parity proofs.

## 4. Tail-marked Matrix–Tree gives the exact endpoint inventory

Here use ordinary commuting polynomial variables \(z_v\) first. Weight
the edge \(u\to v\) by **the tail variable** \(z_u\), and put

\[
L(z)=\operatorname{diag}(A^\top z)-ZA,
\qquad Z=\operatorname{diag}(z_v),
\qquad \tau_r(z)=\det L(z)^{(r)}.                           \tag{4.1}
\]

The superscript deletes row and column \(r\); the columns of \(L\) sum
to zero. The directed Matrix–Tree theorem gives

\[
\tau_r(z)=\sum_{T\text{ out-arborescence rooted }r}
                 \prod_vz_v^{\deg_T^+(v)}.                \tag{4.2}
\]

An out-arborescence is directed away from its root. Project (4.2) to the
commuting square-zero ring. A tree survives exactly when every out-degree
is at most one. A connected rooted tree with that property is a single
directed spanning path. Its sole vertex of out-degree zero is its endpoint.
Consequently,

\[
\boxed{\quad
\tau_r(z)\bmod(z_v^2:v\in V)
   =\sum_{t\in V}H_{r,t}(G)\,z_{V\setminus\{t\}}.
\quad}                                                   \tag{4.3}
\]

This holds for every loopless digraph. At \(n=1\), the empty cofactor is
one and records the singleton path. It also directly refines the inherited
leaf filtration, since every rooted tree satisfies

\[
\#\text{leaves}=1+\sum_v\max(0,\deg_T^+(v)-1).             \tag{4.4}
\]

The squarefree projection is its zero-branching-excess, one-leaf layer.

**The choice of tail is essential.** If edge \(u\to v\) is instead marked
by its head \(z_v\), every nonroot vertex already has in-degree exactly
one in every out-arborescence. The head-marked cofactor is simply

\[
a_r(G)\prod_{v\ne r}z_v,                                 \tag{4.5}
\]

where \(a_r\) counts all out-arborescences. Squarefree projection then
discards nothing. For \(TT_3\), rooted at its source 0,

\[
\tau_0^{\rm tail}=z_0^2+z_0z_1
       \longmapsto z_0z_1,
\qquad
\tau_0^{\rm head}=2z_1z_2.                                \tag{4.6}
\]

Thus “mark each visited vertex” suffices for walks but not for trees:
every spanning tree already visits each vertex. In a tree, mark the
**branching resource** that the desired path must use at most once.
The ordinary unweighted cofactor counts a relaxation:
\(H(G)\le\sum_ra_r(G)\), often strictly. Even its parity need not match
Hamilton-path parity: \(TT_3\) has \(H=1\) and total tree count two.

## 5. Boolean-lattice extraction using ordinary determinants

The cofactor \(\tau_r\) in (4.1) is homogeneous of degree \(n-1\).
For an endpoint \(t\), ordinary zero-one evaluations recover its exact
Hamilton count:

\[
\boxed{\quad
H_{r,t}(G)=
\sum_{S\subseteq V\setminus\{t\}}
(-1)^{n-1-|S|}\tau_r(\mathbf1_S).
\quad}                                                   \tag{5.1}
\]

To prove this, consider a monomial of \(\tau_r\). A monomial involving
\(z_t\) contributes zero to every evaluation. For any other monomial,
inclusion–exclusion vanishes unless its support is all
\(V\setminus\{t\}\). That support has size \(n-1\), equal to the total
degree; thus each exponent must be one. This selects precisely the
coefficient in (4.3).

This is an exact numerical determinant procedure with an exponential
support cost: (5.1) uses \(2^{n-1}\) cofactors per ordered endpoint pair.
It does not turn Hamilton-path counting into a polynomial-time determinant
calculation. The square-zero ring, subset lift, and Boolean Möbius sum are
three representations of the same retained occupancy information.

## 6. Exact deletion, birth, and arc-change kernels

Return to the square-zero ring, and let

\[
E_u(x)=\sum_{\text{simple paths ending at }u}x_{V(P)},
\qquad
S_v(x)=\sum_{\text{simple paths starting at }v}x_{V(P)}.     \tag{6.1}
\]

These have no constant term; singleton paths contribute \(x_u\) and
\(x_v\). A vertex deletion is exact specialization:
\(P_{G-w}=P_G|_{x_w=0}\).

Add a new vertex \(z\), with \(a_u=1\) if \(u\to z\) and \(b_v=1\)
if \(z\to v\). In an arbitrary digraph these masks are independent;
in a tournament \(b_v=1-a_v\). Then

\[
\boxed{\quad
P_{G+z}=P_G+x_z
 \left(1+\sum_ua_uE_u\right)
 \left(1+\sum_vb_vS_v\right).
\quad}                                                   \tag{6.2}
\]

Every new simple path either avoids \(z\), or splits uniquely at \(z\)
into a prefix and a suffix. The two units allow either side to be empty,
and their product counts the path consisting of \(z\) alone. Square-zero
multiplication enforces disjointness of the two sides. This proves (6.2)
at every coefficient, not only at the Hamilton coefficient.

For nonempty old \(V\), define

\[
e_u=[x_V]E_u,\quad s_v=[x_V]S_v,\quad
C_{uv}=[x_V]E_uS_v.                                      \tag{6.3}
\]

Then \(C_{uu}=0\) and

\[
H(G+z)=\sum_ua_ue_u+\sum_vb_vs_v+\sum_{u,v}a_ub_vC_{uv}. \tag{6.4}
\]

The matrix \(C\) counts ordered two-path covers with the indicated exit
and entry ports. It is the sidecar needed in addition to one-path endpoint
totals. For an empty old graph, (6.2) instead gives \(H(\{z\})=1\).

Changing a single arc weight \(A_{uv}\) by \(\delta\) gives the exact
full-inventory identity

\[
\boxed{\quad P_{G'}-P_G=\delta E_uS_v.\quad}                \tag{6.5}
\]

A simple path uses that arc at most once, and splitting there gives the
product. The old endpoint inventories may be used: any factor already
using the changed arc repeats one of its endpoints across the product and
vanishes. In particular, reversing a tournament arc \(u\to v\) gives

\[
\boxed{\quad P_{G'}-P_G=E_vS_u-E_uS_v.\quad}               \tag{6.6}
\]

A simple path cannot use both opposite arcs, so there is no mixed term.
Taking top coefficients in (6.5) or (6.6) gives literal count change,
while keeping all coefficients preserves future deletion and growth data.

## 7. Recovered ear field: an exact even directed cut

For a tournament, substitute \(b=1-a\) into (6.4), use
\(\sum_vs_v=H(G)\), and collect the Boolean monomials:

\[
H(G+z)=H(G)+\sum_u\lambda_ua_u
                 -\sum_{u<v}\gamma_{uv}a_ua_v,            \tag{7.1}
\]

where

\[
\lambda_u=e_u-s_u+\sum_vC_{uv},\qquad
\gamma_{uv}=C_{uv}+C_{vu}.                               \tag{7.2}
\]

Both \(\lambda_u\) and \(\gamma_{uv}\) are even integers. Indeed, setting
\(a=0\) gives the old odd count (the new vertex is a source), and setting
\(a\) to each singleton shows \(\lambda_u\) even. Setting it to each pair
then shows \(\gamma_{uv}\) even. This uses the all-order tournament
parity theorem, not a finite observation. Thus

\[
\boxed{\quad
\frac{H(G+z)-H(G)}2
 =\sum_u\frac{\lambda_u}2a_u
       -\sum_{u<v}\frac{\gamma_{uv}}2a_ua_v.
\quad}                                                   \tag{7.3}
\]

Moreover \(\gamma_{uv}\ge0\), because it is a path-cover count. Therefore
the birth-count function on incoming-neighbour subsets is submodular:
adding another incoming neighbour can only decrease the marginal change
of making a fixed vertex incoming.

The older directed-exposure identity is stronger: no auxiliary orientation
or flow assignment is needed. Count good gaps in old Hamilton paths at
each vertex. Their incoming-minus-outgoing count is \(e_u-s_u\), so

\[
e_u-s_u=\sum_{v\to u}C_{vu}-\sum_{u\to v}C_{uv},\qquad
\lambda_u=\sum_{v\to u}\gamma_{uv}.                       \tag{7.4}
\]

Here a two-path cover with an available joining arc is exactly a Hamilton
path with that arc distinguished. Consequently, for the **incoming** set
\(I=\{u:u\to z\}\),

\[
\boxed{\quad
H(G+z)=H(G)+\sum_{v\notin I,\ u\in I,\ v\to u}
                       \gamma_{uv}.
\quad}                                                   \tag{7.5}
\]

Thus half the birth increment is the incoming-cut capacity in the original
tournament with integer capacities \(\gamma_{uv}/2\). THM-4128 and
THM-4162 instead use the new vertex's **outgoing** set \(S=V\setminus I\),
so their same formula is an outgoing cut. This orientation convention must
be retained when reusing the kernel. In their notation \(Q_{uv}=C_{uv}\)
and \(c_{uv}=2w_{uv}=\gamma_{uv}\).

In particular \(H(G+z)\ge H(G)\), with equality exactly when this incoming
cut has zero capacity. Source and sink attachments always have zero
capacity. The scalar birth field remembers \(H\), the original orientation,
and the symmetric capacities; it forgets individual path witnesses and the
rest of the induced-subset endpoint inventory. This is an exact structural
restriction on possible birth profiles, not a guarantee that every desired
value is attained. Finding orientations and changing parents so their
weighted-cut values cover prescribed totals remains a separate problem.

## 8. Exact controls and scope

Run from the worktree root:

```bash
python3 -B 05-knowledge/results/spanning_tree_fill_20261113/matrix/check.py
```

[check.py](check.py) uses only integer arithmetic and the Python standard
library. Its universe is **all 1099 labelled tournaments of orders 1–5**
and **all 69 loopless directed graphs of orders 1–3**, including digons in
the latter class. It compares permutations, subset dynamic programming,
occupancy matrix powers, determinant expansions, direct parent-function
tree enumeration, and ordinary numerical cofactors.

The saved receipt includes 5606 rooted tree comparisons, 1166 full growth
identities, 802 arc-change identities, 674 Boolean Möbius endpoint
extractions, and 1099 even birth-coefficient controls. It checks
induced-subset coefficients, not just full-set
totals. Its named hostile and positive controls include:

| Graph | Hamilton paths | Sum of \(A^{n-1}\) entries | Total out-arborescences |
|---|---:|---:|---:|
| \(TT_4\) | 1 | 1 | 6 |
| \(C_3\triangleright TT_1\) | 3 | 6 | 9 |
| \(C_3\) plus two isolates | 0 | 3 | 0 |
| Three-vertex out-star | 0 | 0 | 1 |
| Three-vertex directed chain | 1 | 1 | 1 |
| Bidirected edge | 2 | 2 | 2 |
| Empty two-vertex graph | 0 | 0 | 0 |
| \(TT_3\) | 1 | 1 | 2 |

The distinctions are structural. Walk powers forget occupancy; ordinary
tree cofactors forget branching; head marking still forgets branching;
tail marking with squarefree extraction restores the exact path layer.
All three exact representations retain exponential subset information.
They give deletion and growth laws for actual graphs, while a full
Chan–Pak target-coverage construction is still not supplied here.

## 9. Known one-missing-pair witnesses for 7 and 21

**KNOWN / FINITE-EXACT, recovered rather than newly discovered.**
[THM-5071, (4c)](../../../../01-canon/theorems/THM-5071-weighted-cycle-formula-for-digraph-path-totals-and-one-reciprocal-pair-exclusions.md)
and its [literal witness output](../../kps_s18_ternary_chains_digraph_weights_20261109/out/one_defect_witness.out)
give the following loopless digraphs. Rows are binary adjacency rows in
vertex order \(0,1,\ldots,n-1\). Every unordered pair except the displayed
one has exactly one arc; the displayed pair has neither arc.

| \(H\) | Adjacency rows | Missing pair |
|---|---|---|
| 7 | `00001 10000 01000 11100 01110` | \(\{0,2\}\) |
| 21 | `000011 100000 110000 111000 011100 001110` | \(\{1,5\}\) |

The order-five 7 witness was already present in Session234's census;
THM-5071 supplies the literal rows above and the order-six 21 witness.
These are a sharp scope comparison with tournaments: allowing just one
missing pair permits both totals that tournaments exclude. The general
one-missing-pair attainment conjecture
[HYP-9640 (D)](../../../hypotheses/HYP-9640-chan-pak-as-a-matrix-semigroup-orbit-and-one-defect-digraph-spectra.md)
remains **OPEN**; these two witnesses do not prove that conjecture.

The arc-change kernel (6.5) makes their tournament completions explicit:

| Original \(H\) | Added arc | Added paths \([x_V]E_uS_v\) | Tournament \(H\) |
|---:|---|---:|---:|
| 7 | \(0\to2\) | 2 | 9 |
| 7 | \(2\to0\) | 2 | 9 |
| 21 | \(1\to5\) | 20 | 41 |
| 21 | \(5\to1\) | 2 | 23 |

The checker independently counts these six literal graphs by permutations
and subset dynamic programming, checks exactly one missing pair and no
reciprocal pair in each starting graph, and checks all four arc increments.
They are named controls, not a new census or a new attainment theorem.
