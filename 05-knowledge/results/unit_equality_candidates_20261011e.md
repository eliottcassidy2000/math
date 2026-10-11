# Exact unit-edge equality candidates in the Moser field

**Status.** PROVED elementary certificate implementation of an inherited field obstruction; FINITE-EXACT candidate and lower-palette controls. No six-colour forced-equality fragment is produced, and no conclusion about the whole plane being six- or seven-colourable is drawn.

The useful search result is a stopping rule: arbitrary finite constructions inside the common rational Moser field have an explicit four-colouring. Consequently an omitted pair that is still exactly one unit apart can always be separated with four, hence six, colours. Increasing the number of these rotations, sums, or shell points cannot overcome this obstruction.

## 1. Inheritance, target, and retained information

The closest current result is **THM-4558**, [residue colourings of field planes](../../01-canon/theorems/THM-4558-residue-colourings-of-field-planes-the-heegner-tower-is-at-most-four-chromatic.md), especially Lemma R and its known Moser-field value. That note explicitly repairs the historical rotation-count/Heegner roadmap (MISTAKE-576 and MISTAKE-579). The present field upper bound is **not a new chromatic-number theorem**. The new work here is an exact source-owned certificate implementation, its elementary denominator proof, and an authenticated candidate ledger.

Other inherited interfaces are [finite algebraic witnesses](finite_algebraic_witness_20261007.md), [unit-rhombus interfaces](unit_rhombus_interfaces_20261011c.md), and [boundary context weights](boundary_context_weights_20261011d.md). In particular, all eleven Moser-spindle edge deletions are already known positive controls at palette three. Merely blocking one-colour swaps, or obtaining many equal-port colourings, does not establish forced equality.

A candidate consists of distinct exact plane points, a retained subset of their unit edges, and two marked vertices $u,v$ with $|u-v|=1$ whose edge is omitted. It is **nonvacuously forced equal at palette $k$** when the retained graph has a proper $k$-colouring and every such colouring gives $u,v$ the same colour. Equivalently, the retained graph is $k$-colourable but adding the omitted edge is not.

This is deliberately a possibly non-induced fragment. The omitted pair remains a genuine geometric unit edge, not a nonunit pair that happens to be forced equal. The implementation retains coordinates, complete unit incidence, fragment incidence, both ports, and palette. Relabelling or translating a different source does not silently authorize replacing the source in a certificate.

The concept board is: exact algebraic coordinates; local residue frames; component origins; omitted-edge equality; and palette-dependent support. The decisive map is from source offsets to a four-element residue field. It preserves inequality on unit edges, loses metric and coordinate information, and therefore needs the original coordinates and unit-incidence sidecar.

## 2. A complete elementary four-colour certificate

Put
\[
 a=(1+i\sqrt3)/2,\qquad b=(5+i\sqrt{11})/6,\qquad
 L=\mathbf Q(a,b)=\mathbf Q(\sqrt{-3},\sqrt{-11}).
\]
The two quadratic fields are distinct, so $1,a,b,ab$ are a $\mathbf Q$-basis. Their relations are
\[
 a^2=a-1,\qquad b^2=(5/3)b-1.
\]
Let $\mathbf Z_{(2)}$ denote rationals with odd denominator. The ring
\[
 R=\mathbf Z_{(2)}[a,b]
   =\mathbf Z_{(2)}+\mathbf Z_{(2)}a+
     \mathbf Z_{(2)}b+\mathbf Z_{(2)}ab
\]
has this exact four-coefficient description: reduction by the two monic relations stays in that span, and independence proves uniqueness.

Let $\mathbf F_4=\mathbf F_2[\alpha]/(\alpha^2+\alpha+1)$. There are two ring homomorphisms
\[
 \rho_+(a)=\alpha,\quad \rho_+(b)=\alpha,\qquad
 \rho_-(a)=\alpha,\quad \rho_-(b)=\alpha^2.
 \tag{1}
\]
Both defining quadratics reduce to $X^2+X+1$. Complex conjugation reduces to squaring under either map, since $\bar a=1-a$ and $\bar b=5/3-b$.

The paired map
\[
 R/2R\longrightarrow \mathbf F_4\times\mathbf F_4
 \tag{2}
\]
is bijective. For example, the images of coefficient parities $(u,v,w,t)$ are
\[
 (u+t)+(v+w+t)\alpha,\qquad
 (u+w+t)+(v+w)\alpha .
\]
If both vanish, the four coefficients vanish. Both sides have sixteen elements. The script also exhausts the sixteen states and checks addition, multiplication, and conjugation.

**Unit-denominator lemma.** If $z\in L$ satisfies $z\bar z=1$, then $z\in R$ and both residues in (1) are nonzero.

Indeed choose the least $d\ge0$ such that $w=2^dz\in R$. If $d>0$, then $w\bar w=2^{2d}$ reduces to zero in both fields. Thus $\rho_\pm(w)^3=0$, so both residues vanish. By (2), all four coefficients of $w$ are even in $\mathbf Z_{(2)}$, and $w/2\in R$, contrary to minimality. Hence $d=0$. Reducing $z\bar z=1$ now gives $\rho_\pm(z)^3=1$.

This proves the assertion for every rational denominator, not just the finite coefficient box used below. One residue map alone does not prove the denominator lemma: $a-b$ is nonzero modulo $2R$ but $\rho_+(a-b)=0$. The second residue detects its leading parity.

**Colouring theorem.** Every finite unit-distance graph with its points in $L$ has an explicit proper four-colouring.

First compute the full exact unit graph of the point set. In each connected component choose an actual source vertex $z_0$. Every offset $z-z_0$ is a sum of unit edge vectors, hence lies in $R$. Set
\[
 c_\pm(z)=\rho_\pm(z-z_0).
 \tag{3}
\]
An edge difference is a unit vector, with nonzero residue, so its endpoints receive different colours. The two frames give two valid colourings, not eight colours jointly needed by the graph.

Component origins are essential. The spindle translated by $1/2$ has coordinates not all in $R$, so reducing absolute coordinates is invalid; subtracting a retained component origin repairs the problem. The finite algorithm chooses the smallest original vertex label in each full-unit-graph component, making the certificate deterministic.

It follows immediately that **every** candidate of Section 1 whose coordinates lie in $L$ has a separating colouring at every palette $k\ge4$. Use (3) on its full unit graph, including the omitted edge, then restrict to the fragment. Extra missing edges cannot invalidate this colouring.

This excludes all finite repeated Moser rotations, rational combinations, and Minkowski sums that remain in this one common field. A globally isometric copy inherits the abstract colouring, but arbitrary mixtures of outside-field rotations are not covered by this source coordinate implementation. No field-independent upper bound follows.

## 3. Exact finite candidate search

The script uses the basis above and rational arithmetic throughout. Distinct construction paths that give the same algebraic point are identified before forming vertices. The complete unit graph is then computed; every one of its edges is tried as an omitted marked edge. No approximate distance tolerance is used.

| Construction | Vertices | Unit edges | Components |
|---|---:|---:|---:|
| Triangular shell of radius 1 | 7 | 12 | 1 |
| Triangular shell of radius 2 | 19 | 42 | 1 |
| Moser spindle | 7 | 11 | 1 |
| Two shells, one rotated by $b$ | 13 | 24 | 1 |
| Shell plus rotated triangle, by Minkowski sum | 21 | 57 | 1 |
| Shell plus rotated shell, by Minkowski sum | 49 | 180 | 1 |
| Union of six $b$-rotated spindles | 22 | 41 | 1 |
| Four-coefficient box $\{-1,0,1\}^4$ | 81 | 312 | 1 |
| Spindle translated by $1/2$ | 7 | 11 | 1 |
| Spindle scaled by $1/2$ | 7 | 0 | 7 |

There are **690 labelled unit-edge omission instances**, each supplied with a verified separating six-colouring from each residue frame. They are a declared construction catalogue, not an isomorphism census. The finite search finds no six-equality candidate, and the theorem explains this entire catalogue at every size.

The 168 rational-parameter rotations
\[
 \frac{m+n\sqrt{-11}}{m-n\sqrt{-11}},
 \qquad -6\le m,n\le6,\quad(m,n)\ne(0,0),
\]
are separate unit-denominator controls, counted with their parameter labels. Their canonical four coefficients all have odd denominators, including instances whose raw formula initially has an even denominator.

For the positive control, the spindle points are
\[
 0,\ 1,\ a,\ 1+a,\ b,\ ab,\ b+ab.
\]
They have exactly eleven unit edges. The full graph has a checked three-colour refutation with 139 proof nodes. Deleting any one edge gives exactly **12 labelled proper three-colourings**, all identifying the deleted endpoints. Each claim is checked both by a positive colouring plus full-edge refutation, and by exhaustive enumeration of all $3^7$ assignments. Every same fragment has a separating four-colouring and six-colouring. Equality forcing therefore does not persist merely by enlarging the palette.

## 4. Certificate API and hostile boundaries

The companion [exact script](../../04-computation/experiments/unit_equality_candidates_20261011e.py) provides:

- **FieldColour / field_colour / verify_field_colour:** original exact points, chosen residue frame, component origins, and a reauthenticated full-unit-graph four-colouring.
- **Fragment / verify_fragment:** original points, retained actual unit edges, an omitted actual unit pair, and an exact palette.
- **separate:** consumes the same source's field certificate and exports a proper colouring separating the marked pair.
- **Equality / verify_equality:** consumes a proper fragment colouring and a complete checked refutation after adding the marked edge. The positive colouring is required to rule out vacuous equality.

Coordinates and incidence are checked before cached graph lookup or equality comparisons. Boolean and floating aliases, coincident points, malformed ports, an included rather than omitted edge, and a nonunit marked pair are rejected. A translated source certificate is rejected by **separate** even when its abstract colouring happens to agree; a new correctly marked source certificate can instead be produced.

The key hostiles are:

1. A single residue frame kills $a-b$ but does not authorize division by two.
2. Absolute reduction fails on a half-translated spindle; component origins restore it.
3. A three-equality fragment ceases to force equality with four or six colours.
4. Scaling the spindle by $1/2$ changes the unit graph completely; algebraic point provenance does not preserve unit incidence.
5. An isomorphic or translated graph cannot replace the source of an already issued receipt.

No five-chromatic source bank was imported or asserted here, and no large six-colour SAT search is claimed. A useful next candidate must first escape this authenticated four-colour field screen (and any stronger applicable inherited screen), while retaining exact geometry. Escaping a screen is necessary for this route, not evidence that a forced-equality fragment exists.

## 5. Reproduction and scope

Run from the repository root:
~~~
python -B -X utf8 04-computation/experiments/unit_equality_candidates_20261011e.py
python -B -O -X utf8 04-computation/experiments/unit_equality_candidates_20261011e.py
~~~

All gates are active under optimization. The declared finite universe comprises sixteen parity states, their operation tables, 168 unit rotations, the ten point sets above, 690 omissions in both frames, eleven complete lower-palette controls, and malformed/source-replacement hostiles. The separate all-size theorem rests on Section 2, not extrapolation from those controls.

Normal and optimized execution agree on **2,701 always-active checks**. The saved [output](unit_equality_candidates_20261011e.out) has LF-normalized SHA256 **e44cc9ec7bda51ec4f6e64b69a2521e35373cfba79848f7c4f29defa61e8ad96**. This package constructs exact separating receipts and small genuine equality receipts; it constructs **no six-colour equality gadget and no six- or seven-chromatic plane witness**.
