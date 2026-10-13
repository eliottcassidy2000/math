# Which finite abelian Cayley graphs can actually be planar unit-distance graphs?

2026-10-13. **PROVED below:** an exact realization criterion, source-bound
closing-relation certificates, and a sharp direction-family deletion repair.
**FINITE-EXACT:** the declared small groups, finite fields and algebraic
coordinate controls. **CITED / supplied accepted papers:** the methodological
interfaces in section 1; their headline theorems are not re-audited or used as
premises of the elementary proofs. **OPEN:** the current six-colour plane
construction and the unresolved exterior-colouring problem in
$\mathbb Q(\zeta_{55})$. No finite-field chromatic lower bound is transferred
to the real plane here. No literature priority claim is made.

## 1. Inheritance and the missing information

The closest proved geometric mechanism is the unit-rhombus identity in
[unit_rhombus_interfaces_20261011c](unit_rhombus_interfaces_20261011c.md):
four **distinct** planar points forming a unit four-cycle satisfy the
parallelogram relation. The recent direction work uses such identities in
[EI19/fish](ei19_fish_placement_q_zeta55_20261012.md). This packet retains every
commutative square of a finite Cayley source and its closing group relations,
rather than only retaining the lengths of the directions.

The finite-field source is the explicitly marked graph of
[THM-4624, finite unit-circle Cayley rigidity](../../01-canon/theorems/THM-4624-finite-beckman-quarles-for-unit-circle-cayley-graphs-affine-rigidity-except-clebsch-hamming-and-ternary-frame.md).
Its colour numbers serve as targets of residue-colouring maps. Such a map is
a graph homomorphism towards a finite graph, not a faithful Euclidean
realization of that finite graph. The canonical small hostile is $K_4=G_2$;
the stronger hostile below is already three-colourable.

The corrected near miss is the current
[EI21/$\kappa(q)$ note](ei21_kappa_q_20261012.md): the in-repository result
$\kappa(13)\ge5$ is a lower bound, not an equality. The separately cited
$\kappa(13)=6$ has the provenance and replication limitations recorded in
[the current logic synthesis, section 10](logic_lenses_jc2_collatz_hn_20261012.md).
Our nonrealizability result does not depend on either chromatic assertion.

The live concept board is: complete bilinear/incidence data; commutative
squares; finite closing relations; marked source embeddings; and deletion
with an explicit cost. The least-used sidecar is the **closing relation**:
locally consistent unit directions can still force a repeated source vertex
to occupy the same point as a distinct one.

The following are precise paper interfaces, not deductions from analogies.

| Supplied source | Retained mechanism | Transfer here and boundary |
|---|---|---|
| [Reconstruction of Function Fields from Mod-ell Milnor K-Theory](<C:/Users/Eliott/Downloads/mod-ell-bogomolov-pop.pdf>), OpenAI, October 5, 2026, Theorem 1.1, pp. 2–4; common incidence and tuple values, pp. 19–24 and 26–28 | Full degree-one/degree-two groups with the bilinear product; simultaneous tuple incidence | Retain all source squares and the same closing relation. A norm or separate direction count cannot replace them. No field-reconstruction theorem is being applied to a graph. |
| [The Deligne–Drinfeld conjecture](<C:/Users/Eliott/Downloads/paper-28.pdf>), OpenAI, September 23, 2026, sections 4–5, pp. 15–16 | A leading projection requires a kernel argument; deletion operators retain endpoint terms | Our mixed difference vanishes on every distinct-direction square, but its finite-period kernel must still be tested. These are different operators on different objects. |
| [Ordinary two-point correlations of multiplicative functions](<C:/Users/Eliott/Downloads/final.pdf>), OpenAI, September 24, 2026, section 4.3 and Lemma 4.6, pp. 18–19 | Deleting exceptional vertices/edges has a declared weighted-average cost | Section 5 gives an exact whole-direction deletion budget. It changes the graph and makes no claim to preserve its original chromatic obstruction. |

The PDF theorem page was visually inspected; the cached full text supplied
the cited section context. The accepted paper conclusions are not new results
of this packet.

## 2. Complete realization criterion

Let $A$ be a finite abelian group and $S\subset A\setminus\{0\}$ a symmetric
set. The simple Cayley graph joins $x$ to $x+s$ for every $s\in S$. Choose one
representative from each inverse pair $\{s,-s\}$, with order-two directions
counted once, and call the resulting set $R$. Let $U=\langle S\rangle$.

A **weak unit realization** means an injective map $z:A\to\mathbb R^2$ taking
every source edge to distance one; additional unit distances between nonedges
are allowed. An **induced unit realization** also forbids those additional
distances. Both notions retain all source vertices distinctly; edge crossings
are permitted, so neither is a planar graph embedding.

**Theorem 1.** The following are equivalent:

1. $\operatorname{Cay}(A,S)$ has a weak unit realization.
2. The natural sum map $\prod_{s\in R}\langle s\rangle\to U$ is an isomorphism.
3. Every connected component is the marked Cartesian product of the cycles
   $C_{\operatorname{ord}(s)}$, interpreting the order-two factor as $K_2$.
4. The graph has an induced unit realization with real algebraic coordinates.

The word “marked” in (3) refers to the specified generators. There is no
permission to change the source edge relation during the implication.

**Necessity.** Fix $s,t\in S$ with $t\ne\pm s$. The four vertices
$x,x+s,x+s+t,x+t$ are distinct. Their unit-cycle images obey

$$z(x)+z(x+s+t)=z(x+s)+z(x+t).$$

Consequently

$$\Delta_s z(x+t)=\Delta_s z(x),\qquad
\Delta_s z(x)=z(x+s)-z(x). \tag{1}$$

Let $B_s$ be the subgroup generated by the other inverse pairs. Equation (1)
makes $\Delta_s z$ invariant under translation by $B_s$.
Put $N=\operatorname{ord}(s)$ and $d=\operatorname{ord}(s+B_s)$ in $A/B_s$.
Then $d\mid N$ and $ds\in B_s$. Along the $s$-cycle the edge increments have
period $d$. Closing that cycle gives

$$0=z(x+Ns)-z(x)=\frac Nd\bigl(z(x+ds)-z(x)\bigr). \tag{2}$$

The target is a real vector space, so $z(x+ds)=z(x)$. Injectivity forces
$ds=0$ and therefore $d=N$. Thus $\langle s\rangle\cap B_s=\{0\}$ for
every representative $s$. In any relation $\sum_s a_s s=0$, its individual
term $a_s s$ lies in this intersection, so every term vanishes. The sum map
is injective and, by definition of $U$, surjective.

**Sufficiency, including nonedges.** For a cyclic factor of order $n\ge2$,
choose a regular unit polygon with vertices

$$b_{n,j}=1+\zeta_n+\cdots+\zeta_n^{j-1},\qquad
0\le j<n,\quad b_{n,0}=0.$$

For $n=2$ this is simply $\{0,1\}$. Its distinct vertices have distance one
exactly when their indices are adjacent modulo $n$. For factors
$n_1,\ldots,n_r$, choose unit complex numbers $\theta_i$ and set

$$z(j_1,\ldots,j_r)=\sum_{i=1}^r\theta_i b_{n_i,j_i}. \tag{3}$$

Every product edge is unit. A difference involving exactly one factor is
nonzero and is unit exactly for a source edge. A difference involving at
least two factors is a sum of at least two nonzero independently rotated
vectors. Its squared norm is not a constant function of the rotation
parameters: the pairwise relative-angle Fourier coefficients are nonzero.
Thus neither collision nor an unwanted unit distance holds identically on
the rotation torus. There are finitely many such conditions. Their proper
algebraic zero sets have empty interior, so choices avoiding them exist.

One can choose the rotations from the dense rational parametrization
$\theta(t)=(1-t^2+2it)/(1+t^2)$ with $t\in\mathbb Q$. Coordinates in (3)
are then algebraic. Different components can be separated by sufficiently
large integer translations, avoiding both coincidences and cross-component
unit distances. This proves (4); the remaining implications are immediate.

The distinctness hypothesis is essential: a noninjective colouring map to a
triangle can exist even when (2) forces collisions. Extra Euclidean unit
edges cannot repair the contradiction, because necessity only used the
original edges.

## 3. Constructive certificates and coordinate search

`classify(source)` works in the original finite group. For each marked $s$ it
enumerates $B_s$ together with explicit generator words. A failure returns
a `Collapse` packet containing $s,N,d$ with $1\le d<N$, $d\mid N$, and a
word in the other directions equal to $ds$. Its checker verifies the same
source, the proper period, the closing relation, and the four-distinct-vertex
condition for every square direction it uses. Equations (1)–(2) then prove
nonrealizability. No proposed coordinates or floating tolerances are inputs.

If there is no failure, a `Product` packet retains the complete bijection
from cyclic coordinates onto the generated subgroup. The verifier checks
its source, completeness and injectivity. Disconnected cosets remain
separate, rather than being silently identified.

There is also an explicit terminating algebraic coordinate search behind
the sufficiency proof. Let $v=|U|$. Take the collision polynomial for each
pair and the unwanted-unit polynomial for each nonedge, after substituting
the rational rotations and clearing the positive denominator
$\prod_i(1+t_i^2)^2$. Each has degree at most four in each variable and is
nonzero by the proof. If their number is $B$, their product has individual
degree at most $4B$. A nonzero polynomial over a field cannot vanish on the
whole grid $\{0,\ldots,4B\}^r$; this follows by induction from the univariate
root bound. Searching this finite grid therefore succeeds. This is an
existence bound for a potentially expensive search, not an efficient
coordinate-height claim.

The implementation uses the cheaper sequential version: add a factor and
test integer rational-rotation parameters through an explicit cap. It
returns `None` on cap exhaustion, which is **inconclusive**. Exact arithmetic
in a cyclotomic quotient checks every pair, so both vertex distinctness and
the complete induced edge set are authenticated. Identical rotations for
$C_3\square C_3$ are a collision hostile; a different rational rotation gives
an induced nine-point realization.

## 4. All finite norm circles: exactly one realizable exception

Let $q=p^e$ be a prime power and let

$$G_q=\operatorname{Cay}\bigl((\mathbb F_{q^2},+),
\{u:u^{q+1}=1\}\bigr).$$

This is the finite-field unit graph of the inherited canon. It is not
defined using Euclidean distances. The additive group is an elementary
abelian $p$-group of dimension $2e$; every nonzero direction has order $p$.
There are $(q+1)/2$ inverse pairs when $p$ is odd, and $q+1$ when $p=2$.

The norm-one directions span the whole additive group. Indeed, choose a
generator $\zeta$ of the cyclic norm-one group, of order $q+1$. Its degree
over $\mathbb F_p$ is $\operatorname{ord}_{q+1}(p)=2e$: this order divides
$2e$, whereas any proper divisor is at most $e$ and would imply
$q+1\mid p^d-1<q+1$. Hence $1,\zeta,\ldots,\zeta^{2e-1}$ is a basis made
entirely of norm-one directions.

**Corollary 2.** Among all prime powers, $G_q$ has a weak planar unit
realization if and only if $q=3$. In that case
$G_3=C_3\square C_3$ and has an induced algebraic realization.

For odd $p$, Theorem 1 would require $(q+1)/2$ independent one-dimensional
subgroups in dimension $2e$. For $e=1$ this is possible only at $p=3$.
For $e\ge2$, $(3^e+1)/2>2e$ (start at $e=2$ and induct). In characteristic
two, $2^e+1>2e$ for every $e\ge1$. At $q=3$ the two inverse pairs are
independent because they span, so the positive direction follows from
Theorem 1. No finite numerical census is needed for this all-$q$ statement.

The smallest useful odd-characteristic closing receipt is in
$\mathbb F_{25}=\mathbb F_5[\alpha]/(\alpha^2-2)$:

$$s=(1,0),\quad t=(2,2),\quad u=(2,3),\qquad s+t+u=0.$$

All three are norm-one directions. Thus $s=-t-u\in B_s$, so $d=1<N=5$;
the marked $s$-increment is constant around a five-cycle and must vanish.

A three-colourable hostile separates geometry from chromatic number even
more sharply. In $\mathbb F_3^4$, use directions
$e_1,e_2,e_3,e_4,-e_1-e_2-e_3-e_4$ and their negatives. The map
$x\mapsto x_1+x_2+x_3+x_4\pmod3$ is a proper three-colouring. Nevertheless
the five generator subgroups cannot be independent in dimension four, so
the graph has no weak unit realization. Its identification with $G_9$ is
inherited from THM-4624; the explicit three-colouring and fold obstruction
are checked here independently of that identification.

## 5. A sharp positive repair with a declared deletion model

Permit deletion only of **whole direction families**: deleting $\{s,-s\}$
means removing all edges $\{x,x+s\}$ for every source vertex $x$. Vertices
are retained, and the resulting Cayley graph is realized afresh. This does
not preserve the removed constraints or require extending an old placement.

For $G_{p^e}$, Theorem 1 says the retained direction representatives must be
linearly independent. At most $2e$ families can remain. Conversely, the norm
circle contains a basis, so $2e$ can remain and the graph becomes
$C_p^{\square 2e}$ (use $K_2^{\square 2e}$ at $p=2$). The exact minimum
number of deleted families is consequently

$$\begin{cases}
(q+1)/2-2e,&p\text{ odd},\\
q+1-2e,&p=2.
\end{cases} \tag{4}$$

For odd $p$ this deletes $q^2$ edges per family; for $p=2$, $q^2/2$ edges.
The independent basis product is three-colourable for odd $p$ and
two-colourable for $p=2$: sum proper cycle colours modulo three, or sum bits
modulo two. Thus this repair explicitly loses any original higher chromatic
obstruction. It does not provide a bound for arbitrary individual-edge
deletions or for a graph outside this Cayley category.

At $q=13$, for example, seven inverse-pair families must be reduced to at
most two, so five whole families must be deleted. This is a concrete reason
why a finite-field colour obstruction cannot simply be drawn in the plane.
The next constructive HN obligation is to retain enough genuinely Euclidean
constraints in a different source graph; increasing the finite field alone
does not meet it.

## 6. Exact controls, interfaces and scope

The full declared small-group universe is every inverse-pair subset in
$\mathbb Z_n$ for $n=2,3,4,5,6,8$, and in
$\mathbb Z_2^2,\mathbb Z_2\times\mathbb Z_3,\mathbb Z_3^2,
\mathbb Z_2\times\mathbb Z_4$: **100 sources**, including disconnected and
edgeless cases. Independent product-size checks find 56 product sources and
44 fold sources. The finite norm-circle controls use $q=2,3,4,5,7,9,13,17$
with directly checked irreducible finite-field moduli. Exact algebraic
induced realizations cover $K_2,C_3,C_5,K_2^{\square3},C_3\square C_3,
C_3\square C_4,C_5\square C_5$. Boolean/float packet aliases, wrong sources,
wrong closing words and reducible field moduli are rejected.

| Object | Retained map/predicate | Information that cannot be discarded |
|---|---|---|
| Finite-field unit graph | Norm-one edge relation on the actual finite source | Finite addition and all closing relations |
| Source Cayley squares | Injective planar unit realization implies (1) | Vertex distinctness; source direction labels |
| Quotient by the other generators | Marked period $d$ versus actual order $N$ | The nonzero source displacement $ds$ |
| Direct-product certificate | Complete cyclic coordinate map | Original-source vertex names and connected cosets |
| Algebraic coordinate realization | Every intended edge and every nonedge | All pair distances, not a sampled direction norm |
| Direction-family deletion | Exact minimum (4) and a realizable basis graph | The removed edge constraints and their old colour obstruction |

There is no new six-chromatic point configuration, no resolution of the
plane six-versus-seven question, and no contradiction with finite-field
residue colourings. No tournament is introduced: the intrinsic relation
here is symmetric unit adjacency, and arbitrarily orienting it would lose
the closing-relation obligation.

[Script](../../04-computation/experiments/unit_cayley_realization_20261013a.py)
and [saved output](unit_cayley_realization_20261013a.out):

```text
python -B -X utf8 04-computation/experiments/unit_cayley_realization_20261013a.py
python -B -O -X utf8 04-computation/experiments/unit_cayley_realization_20261013a.py
```

Normal and optimized runs match the saved output at **1,147 exact checks**;
all checks remain active under `-O`. LF SHA256:
`0d8bca8fd5deb3c2c449f939495aca49c2ff5a214d5075532cfc6008ef16c4bd`.
