# Boolean realization energy and marked tournament transport

**Status:** PROVED for the lemmas below; FINITE-EXACT for the displayed
tournaments and complete specified censuses. Same-agent proof audit with
independent computational paths. The constructions specialize standard
Boolean algebras, probability and conditional expectation; no global
priority claim is made.

## The object that a weighted orientation actually represents

Let \(\Omega_n\) be the \(2^{\binom n2}\) labelled tournaments, with one bit
\(e_{ij}\) for each \(i<j\); bit one means \(i\to j\). The Boolean algebra
of real-valued functions is
\[
 {\cal A}_n=\mathbb R[e_{ij}]/(e_{ij}^2-e_{ij})
           \cong\mathbb R^{\Omega_n}.
\]
Every class has a unique multilinear representative. Write \(B[f]\) for
this reduction. The Hamiltonian-path function \(H\) is the sum of the
orientation indicators of all vertex permutations. Each path uses each
unordered pair at most once, so its ordinary path polynomial is already
multilinear.

A weighted point \(x\in[0,1]^{\binom n2}\) specifies independent Bernoulli
edges. Evaluation of the multilinear representative is expectation:
\[
 L_x(f)=B[f](x)=\sum_{T\in\Omega_n}\mu_x(T)f(T).
\]
It is a positive unital linear functional. It is multiplicative for every
pair of functions exactly when every edge probability is zero or one:
\(L_x(e_i^2)-L_x(e_i)^2=x_i(1-x_i)\).
Thus a polynomial value at a fractional point is a mean, not a Boolean
orientation. A nonzero gradient does not change this type distinction.

More generally, a positive unital linear map
\(\Phi:\mathbb R^Y\to\mathbb R^X\), for finite sets \(X,Y\), is a stochastic
kernel: \(\Phi f(x)=\sum_y K(x,y)f(y)\). It is an algebra homomorphism iff
each row of \(K\) is a point mass. Proof: apply multiplicativity to the
idempotent singleton indicators; their images are zero or one and sum to
one. The converse follows by direct substitution.

This gives a precise constructor test. A deterministic Boolean source
map induces an algebra homomorphism. A weighted average induces only a
positive linear map. The multiplicativity defect
\[
 D_\mu(f,g)=L_\mu(fg)-L_\mu(f)L_\mu(g)
\]
is covariance. Vanishing defect for \(H\) alone says that \(H\) is constant
on the support; it does not make the orientation itself deterministic.

## Target energy and a constructive rounding lemma

**Lemma E1, PROVED.** For an integer target \(h\), define the Boolean
realization energy
\[
 {\cal E}_h(x)=B[(H-h)^2](x).
\]
Then
\[
 {\cal E}_h(x)=\mathbb E_{\mu_x}(H-h)^2
 =\operatorname{Var}_{\mu_x}(H)+(H(x)-h)^2.
\]
It is nonnegative and vanishes iff every positive-probability completion
has count \(h\). Its support is exactly the Boolean face obtained by
fixing coordinates whose probabilities are zero or one. Consequently
\[
 \min_{x\in[0,1]^N}{\cal E}_h(x)
 =\min_{T\in\Omega_n}(H(T)-h)^2.
\]
Proof: the energy is a convex combination of its Boolean values, and every
Boolean vertex is allowed. The same assertions hold on any fixed face,
in particular a fixed-Hamiltonian-path chart.

**Lemma E2, PROVED.** Given exact evaluations of \({\cal E}_h\), a weighted
point on a face can be rounded to a Boolean vertex with no greater energy,
using at most one clamp per free coordinate. For a coordinate with
probability \(p\),
\[
 {\cal E}_h(x)=(1-p){\cal E}_h(x_{i=0})+
                    p{\cal E}_h(x_{i=1}).
\]
Choose an endpoint no larger than this average and repeat. Fixed arcs
remain fixed. If \(h\) is odd and \({\cal E}_h(x)<4\), the rounded tournament
has exactly \(h\) Hamiltonian paths, since all such counts are odd and
distinct odd integers differ by at least two.

This is a sufficient constructive certificate, not an efficient method
for finding a point below the threshold. Exact conditional energy queries
can be expensive. The finite demonstrations for \(h=17,71\) deliberately
start near known solutions; they do not discover new spectrum values.

The same clamping proof works for any real function on a finite cube.
If its values belong to an integer lattice and the target is an integer,
energy below one suffices. The threshold four above uses ordinary
tournament parity specifically; guarded path counts need not be odd.

There are three different inverse requests:
a **mean inverse** satisfies \(\mathbb E H=h\);
a **support inverse** is a positive law supported on \(H^{-1}(h)\);
a **point inverse** supplies an actual tournament in \(H^{-1}(h)\).
Energy zero upgrades a mean inverse to a support inverse, and selecting a
positive atom gives a point inverse. The strict inequality in E2 is a
stronger rounding convenience for product laws, not a necessary condition
on every weighted point around a feasible target.

## The sharp gap at seven and twenty one

**Lemma E3, PROVED using the inherited spectrum exclusions.** Every
positive probability law on tournaments of any fixed order satisfies
\[
 \mathbb E(H-h)^2\ge4\qquad(h=7\text{ or }21).
\]
The input is the all-order exclusion of these two counts plus Rédei
parity. Each support value is a distinct odd integer, so its squared
distance from \(h\) is at least four. Equality holds iff the law is
supported on \(h-2,h+2\). At mean \(h\), equality requires equal masses.

More generally, if every actual count lies outside the interval
\((a,b)\), then \((H-a)(H-b)\ge0\). At mean \(\mu\in[a,b]\),
\[
 \operatorname{Var}(H)\ge(\mu-a)(b-\mu).
\]
This is a gap certificate for the positive second-moment model.

The constant four is sharp even among independent-edge laws. Full arc
masks use lexicographic pair order starting at bit zero:

| Order | Target | Boolean masks | Counts | Fractional edge |
| --- | --- | --- | --- | --- |
| 5 | 7 | 4, 260 | 5, 9 | (2,4) |
| 6 | 21 | 20, 28 | 19, 23 | (0,4) |

These exact one-edge pairs are also recorded in the incoming
[Session142 bridge, Proposition 2](../session142_tournament_odd_spectrum_20261102/jacobian_bridge/PROOFS.md).
They were independently reproduced here; the pairs and the expectation
interpretation are recovered inputs, not claimed as new discoveries.
Our additions are the target energy, its sharpness and rounding bound,
the interior certificates, and the explicit local-search obstruction.

Giving that edge probability \(p\), with the other arcs fixed, gives
\(H(p)=h-2+4p\). At \(p=1/2\), the mean is \(h\), the derivative is four,
and the Boolean energy is four. The naive squared polynomial
\((H(x)-h)^2\) is zero. Both pairs share a directed Hamiltonian path;
the certificate gives its normalization to a single edge of the
fixed-path chart.

These examples also extend to regular points in the **full open cube**.
Move every other probability to \(\epsilon\) or \(1-\epsilon\) according
to its Boolean endpoint, with \(\epsilon=1/1024\). Let \(A,B\) be the two
conditional means at the distinguished edge. Exact rational calculations
give \(A<h<B\). Set \(p=(h-A)/(B-A)\); then all probabilities lie strictly
between zero and one, the mean is exactly \(h\), and the derivative is
\(B-A>0\). Independent weighted-path sums verify both means. The energy
is now strictly greater than four because every tournament has positive
weight. Here “regular” concerns the polynomial differential, not equal
vertex scores.

Positivity is indispensable: signed weights \(1,-3,3\) on actual counts
\(1,3,5\) have total mass one, mean seven, and zero “energy.”
They are not a probability law. Mean alone cannot exclude seven or
twenty one because each is inside the convex hull of attained counts.
One additional second-moment scalar suffices to test support concentration
for a known positive law; it does not supply all information needed for
future attachments.

The exclusions are inherited, not new: see
[THM-1370 / h-spectrum-omits-7-21-all-n](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md).
Busch's [primary paper](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/)
proves the sharp strong-tournament floor and recalls parity and component
multiplicativity. The [original 2016 discussion](https://mathoverflow.net/questions/232751/the-number-of-hamiltonian-paths-in-a-tournament)
already records the two exclusions. Completeness of the unrestricted
positive odd spectrum outside these values remains OPEN.

## Exact energy still has strict local traps

**FINITE-EXACT witness, with a PROVED local consequence.** On six vertices,
mask 72 has outgoing-neighbour bit rows
\[
 (16,9,3,5,14,31)
\]
and \(H=15\). For target 17 its energy is four. Flipping the fifteen edges
in lexicographic order gives the counts
\[
 13,13,13,3,33,9,9,13,31,9,13,31,13,31,27.
\]
Every squared error is strictly larger than four. This is a strict local
minimum even in the continuous product cube: the derivative into the
cube in every coordinate is positive, and all remaining terms are of
second or higher order in distances from the vertex. There is a small
neighbourhood in which their sum cannot cancel the positive linear term.

Yet count 17 is attained on six vertices. The full labelled sublevel
graph \((H-17)^2\le4\), with single-edge flip adjacency, has 769 connected
components, of which 768 are isolated non-solutions. Thus even neutral
steps cannot rescue monotone single-edge descent from these vertices.
This does not exclude multi-edge moves, nonmonotone search, or another
certificate. The failed implication is “exact global objective implies
a locally convergent inverse,” not the rounding lemma.

## Shared sources change the composition law

**Lemma E4, PROVED.** For an ordered join \(A\Rightarrow B\), all arcs
point from \(A\) to \(B\), so every Hamiltonian path visits \(A\) first and
then \(B\). Hence \(H(A\Rightarrow B)=H(A)H(B)\).
For independent random inputs, first and second moments multiply, and
\[
 \operatorname{Var}(XY)=
 \operatorname{Var}(X)\operatorname{Var}(Y)
 +(\mathbb EX)^2\operatorname{Var}(Y)
 +(\mathbb EY)^2\operatorname{Var}(X).
\]
For two copies of the same source \(T\), the mean is instead
\(\mathbb E H(T)^2=(\mathbb EH(T))^2+\operatorname{Var}(H(T))\).

In the count 5/9 law above, independent copies yield join counts
\(25,45,81\) with probabilities \(1/4,1/2,1/4\), of mean 49.
Copying the same source twice yields \(25,81\) with probabilities
\(1/2,1/2\), of mean 53. All outcomes are genuine ten-vertex tournaments.
The discrepancy four is exactly the covariance forgotten by replacing
a shared source with independent marginals. The verifier constructs and
counts all four joins directly.

## Actual metagraph reconstruction requires transport

The fixed-path chart \(X_n=\{0,1\}^{\binom{n-1}{2}}\) fixes
\(0\to1\to\cdots\to n-1\). Its free coordinates \((i,j)\), \(j\ge i+2\),
are ordered lexicographically, with bit one meaning \(i\to j\).
Let \(C\) flip every free bit. The merged metagraph is the span
\[
  M_n\ \longleftarrow\ X_n\ \longrightarrow\ X_n/\langle C\rangle .
\]
The first map forgets the path marking and takes the tournament class
modulo converse. The second forms literal complement lines. Lines can
be parallel and loops retain two half-edges. This is the inherited
[THM-781 / hamiltonian-path-inverse-metagraph-fibre](../../../01-canon/theorems/THM-781-hamiltonian-path-inverse-metagraph-fibre.md)
and [THM-796 / three-sorted-recursive-tiling-line-node-incidence](../../../01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md).

Use **transport coherence** for the requirement that every actual legal
Hamiltonian-path re-marking has the correct source and destination
fibre, with its permutation and possible converse retained as a witness.
Boolean realization alone does not imply this property.

The inherited order-five hostile exchanges the class labels only at
chart tiles 4 and 26. It is still a genuine partition of the entire
Boolean cube with all correct fibre masses. Thus its indicator functions
satisfy every Boolean idempotence/convolution identity. However the legal
Hamiltonian path \((0,4,1,2,3)\) in tile 4 re-marks it as tile 9.
The true merged labels of 4 and 9 are both 12; the fake labels are
41 and 12. The direct verifier replays this exact obstruction. The
[original realizability certificate](../merged_metagraph_realizability_20261020/REALIZABILITY.txt),
sections 7–8, proves the stronger colour and Fourier invariants also survive.
This packet does not re-audit every one of those stronger invariants.

Accordingly, a reconstruction must discharge separate obligations:
real Boolean orientations; the correct isomorphism/converse fibres;
legal marked operations; and the joint source law when weighted copies
share ports. These are independently testable conditions, not synonymous
uses of the word “realizable.”

## The literal seventeen and seventy one connection

**FINITE-EXACT, complete seven-vertex chart.** Among all 32,768 fixed-path
tiles, exactly twelve complement lines join path counts 17 and 71.
The count-17 tiles have 21 distinct complement counts:
\[
 27,29,31,33,35,41,43,45,49,51,53,55,59,65,67,69,71,73,77,85,93.
\]
They occupy two merged tournament nodes, each of mass 34.

More strongly, tiles 24511 and 31 are the **same ordinary tournament**:
re-mark tile 24511 along the legal path
\[
 (0,5,6,3,4,1,2)
\]
to obtain tile 31, without converse. Both have count 17. Their complements
are tiles 8256 and 32736 with counts 71 and 33. Thus even the actual
unmarked isomorphism class does not determine its complement destination.
No improvement of a scalar formula can fix this: the missing coordinate
is the selected path port.

For the deliberately restricted command language consisting only of
repeated \(C\), with output \(H\), the pair
\[
  S(x)=(H(x),H(Cx))
\]
is the smallest exact deterministic observer. Indeed \(C^2=1\), so this
pair determines every future output; differing pairs are separated by
either the empty command or one \(C\). At order seven it has exactly
3,008 states. This is not a minimum carrier for re-marking or attachments;
those commands require a fresh congruence test.

The all-order blueprint is therefore an operation-relative partition:
two sources can be identified only if every allowed future word has the
same legality and observable result. This is standard behavioural
equivalence. It clarifies the target without claiming an all-order
classification of that partition.

