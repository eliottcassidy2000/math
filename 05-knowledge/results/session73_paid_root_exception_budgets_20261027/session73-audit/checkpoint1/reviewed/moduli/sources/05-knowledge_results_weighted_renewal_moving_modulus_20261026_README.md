# Weighted renewal, moving precision, and a relocation-closed port forest

2026-10-26. **PROVED scoped**, with self-contained proofs and **FINITE-EXACT**
controls; classical ingredients are credited. Internal adversarial review only,
not an independent-agent audit. Collatz convergence, native ROOT positivity,
and a compact grammar closed under all complement-line operations remain **OPEN**.

## Inheritance and concept board

The anchor is source-specific renewal energy; the niche is relocation of a
marked tournament; the wildcard is adaptive arithmetic precision. The live
concepts are **cut energy, named source, boundary carry, exception forest,
height/precision, and integer clock**.

Closest proved mechanisms:

- [Native-return packet N3–N6](../native_return_necklaces_20261026/README.md):
  the separable-permutation/two-port grammar, exact reset-segment inverses,
  discount costs, and affine Collatz cylinders. We extend its grammar under
  relocation and replace selected discount weights by a test of supplied weights.
- [THM-4807](../../../01-canon/theorems/THM-4807-untraced-arcs-are-cycle-closing-chords-elementary-legal-path-reconstruction-endpoint-decoder-one-guard-separation.md):
  all Hamilton guards reconstruct the labelled tournament. Thus an unrestricted
  guard-faithful state quotient cannot identify distinct labelled tournaments.
- [THM-4808](../../../01-canon/theorems/THM-4808-berggren-tree-is-non-backtracking-branch-character-pole-frame-every-modulus-obstruction-cusp-pairing.md):
  every fixed modulus has opposite-character endpoint collisions; the pole frame
  restores the missing determinant. Its near-Ramanujan prime census is finite
  evidence, not a uniform branch-parity theorem.
- [Collatz renewal R1, R4, R6](../collatz_renewal_source_20261011.md): finite
  excursion laws, the noninteger repeat limit 19/5, and rational positive atomic
  source weights with explicit residue-alias tails. R6 is the least-used useful
  sidecar here: changing precision can pay for growing response energy.

Canonical hostiles: the N3 order-four relocation outside the old grammar;
one reset followed by a long integration block; arbitrary finite all-one
Collatz valuation words; and the THM-4808 root-residue collision. Corrected
near misses: Haar renewal laws do not prove integer-source convergence;
complement can change labelled modules without changing the unlabelled tree;
refinement alone does not supply a uniform energy bound.

Incoming commit `31002464ae`, [Session72 synthesis](../session72_even_graph_unifying_carriers_20261026/session72-root/UNIFYING_THEOREMS.md),
adds two useful distinctions: state separation versus response rank, and
pushforward of masses versus conditional averaging of observables. Its
test-return results count whole-permutation tests, not elementary move clocks.
The constructions below preserve these distinctions.

## R1. A cut criterion for supplied weighted renewal gains

**PROVED; classical discrete weighted Hardy mechanism.** On a finite segment
with zero initial state, consider

\[
 x_1=f_1,\qquad x_j=a_jx_{j-1}+f_j\quad(2\le j\le L),
\]

where the nonzero gains are real and the supplied input/output weights
\(v_i,u_j\) are strictly positive. Set \(p_1=1\),
\(p_j=a_2\cdots a_j\), and

\[
 B=\max_{1\le k\le L}
 \left(\sum_{j=k}^L u_jp_j^2\right)
 \left(\sum_{i=1}^k\frac1{v_ip_i^2}\right).
 \tag{1}
\]

For the operator \(R:f\mapsto x\), with the indicated weighted norms,

\[
 \boxed{B\le\|R\|^2_{\ell^2(v)\to\ell^2(u)}\le4B.} \tag{2}
\]

Zero gains split the recurrence into independent reset segments. In that case
take the maximum B over segments. A collection of source-labelled recurrences
with the same positive source mass \(\mu_s\) multiplying both norms has norm
given by the supremum over sources and segments. Source masses cancel in each
individual B; a small average B cannot replace this supremum. No source is
rerandomized and no gain is changed by choosing the test.
This criterion controls every forcing vector. Failure of a uniform B bound
does not rule out cancellation in the particular native signed current; that
more selective estimate remains a possible, weaker route.

**Proof.** The substitution \(x_j=p_j y_j, f_i=p_i g_i\) reduces the problem
to prefix summation with weights \(U_j=u_jp_j^2,V_i=v_ip_i^2\); signs are an
isometric diagonal gauge. Write \(W_k=\sum_{i\le k}V_i^{-1}\) and
\(T_k=\sum_{j\ge k}U_j\). The test \(g_i=V_i^{-1}1_{i\le k}\) proves
\(\|R\|^2\ge T_kW_k\).

For completeness, a Schur proof supplies the constant4. The normalized
matrix is \(K_{ji}=\sqrt{U_j/V_i}\,1_{i\le j}\). Take
\(h_i=(V_iW_i)^{-1/2}\) and \(l_j=(U_jW_j)^{1/2}\). Then

\[
 \sum_iK_{ji}h_i\le2l_j,
 \qquad \sum_jK_{ji}l_j\le2B h_i.
\]

The first estimate uses
\(\sum_{i\le j}(W_i-W_{i-1})/\sqrt{W_i}\le2\sqrt{W_j}\).
For the second, summation by parts and \(T_j\le B/W_j\) give

\[
 \sum_{j\ge i}U_j\sqrt{W_j}
 =T_i\sqrt{W_i}+\sum_{j>i}T_j(\sqrt{W_j}-\sqrt{W_{j-1}})
 \le\frac{2B}{\sqrt{W_i}}.
\]

The last sum is bounded by the integral of \(t^{-2}\) from
\(\sqrt{W_i}\) to infinity. Weighted Cauchy–Schwarz applied to the two
Schur inequalities gives \(\|K\|^2\le4B\). Direct sums prove the remaining
claims. The argument extends to infinitely many coordinates by finite
truncation precisely when the cut supremum stays finite. QED.

The underlying criterion is classical, not a new Hardy inequality: see
[Bui–López-García–Tran, introduction](https://www.cpp.edu/faculty/fal/paper_bui_lopezgarcia_tran.pdf)
for its discrete formulation and
[Muckenhoupt, 1972](https://www.impan.pl/en/publishing-house/journals-and-series/studia-mathematica/all/44/1/98433/hardy-s-inequality-with-weights)
for the foundational weighted inequality. The all-gain/source-labelled
formulation and boundary accounting are the application here.

### What a block cut actually costs

Cut after a. The cross-block matrix has entries \(p_j/p_i\), i≤a<j;
it is rank one, with squared weighted norm exactly

\[
 C_a=\left(\sum_{j>a}u_jp_j^2\right)
       \left(\sum_{i\le a}(v_ip_i^2)^{-1}\right). \tag{3}
\]

It carries the scalar boundary state \(x_a/p_a\). Artificially declaring a
reset deletes this port. A one-dimensional interface therefore can have
unbounded energy: with unit weights and gains, C_a=a(L−a). Compact boundary
dimension is not a bounded boundary norm. This is the exact relationship to
the tournament port problem, not an assertion that the two operators coincide.

## R2. Refinement must preserve clock and the weighted embedding

**PROVED.** Split each coarse atom i into named children with positive weights
whose sum is its old weight. Let J copy an observable to its children and
E take its conditional weighted mean. Then J is an isometry, E=J*, and EJ=I.
For input and output refinements define

\[
 R_{\rm ref}=J_{\rm out}R E_{\rm in}. \tag{4}
\]

Its norm equals the original norm exactly, and
\(R_{\rm ref}J_{\rm in}=J_{\rm out}R\). Both equalities follow from the
isometries and contractions; the reverse norm inequality tests copied inputs.
All children retain the parent's integer timestamp. This is refinement of an
observation, not insertion of additional dynamical steps. Mass pushforward
is summation, whereas E averages observable values; these are different maps.

For an actual refined operator A write D=A−R_ref. Then
\(\|A\|\le\|R\|+\|D\|\). Agreement on copied inputs only says DJ_in=0;
it does not bound D on zero-mean fibre modes. For identical input/output
fibres, \(JRE+M(I-JE)\) already exhibits arbitrarily large hidden norm as
M grows, despite perfect agreement on all copied coarse inputs.

**Minimal hostile to mass-only refinement.** Split one mass1 atom into two
mass1/2 atoms and substitute the two-step prefix-sum matrix. On input (1,1)
the energy ratio is5/2 instead of1. For r equal children it is
\((r+1)(2r+1)/6\). Refining descriptions without inserting clock steps
avoids this error; inserting genuine steps must pay their carry.

## R3. A rational passive two-port identity on every positive Collatz branch

**PROVED scoped.** For an actual positive odd source n, let
\(k=v_2(3n+1)\), \(n'=(3n+1)/2^k\). The relative linear multiplier is

\[
 b=\frac{(3/2^k)n}{n'}=\frac{3n}{3n+1},\qquad
 d=1-b^2=\frac{6n+1}{(3n+1)^2}>0.
\]

For arbitrary real or rational relative error e and forcing f, retain a
second output port:

\[
 e'=be+f,\qquad r=de-bf.
\]

Direct expansion gives the lossless identity

\[
 \boxed{(e')^2+\frac{r^2}{d}=e^2+\frac{f^2}{d}.} \tag{5}
\]

This is the elementary two-port orthogonal dilation of a scalar contraction,
written without square roots. For zero forcing it describes relative
differences along the same affine valuation branch; finite differences need
the same legal branch word. It is not an assertion that nearby integers have
the same valuation word.

Along any finite actual orbit prefix, (5) telescopes:

\[
 e_L^2+\sum_{j=1}^L\frac{r_j^2}{d_j}
 =e_0^2+\sum_{j=1}^L\frac{f_j^2}{d_j}. \tag{6}
\]

Multiplication by each original source mass and summation preserves the
identity. The physical elapsed clock is \(\sum_j(1+k_j)\), with each term
attached to that actual transition; the identity does not replace it by L.
Thus we have a universally defined, exactly energy-controlled *weaker*
carrier for positive sources, at all finite depths.

Its price is explicit:

\[
 \frac1d=\frac{(3n+1)^2}{6n+1}\sim\frac{3n}{2}.
\]

Neither this forcing budget nor a positive ROOT observation has been supplied
by (6). Forgetting r also forgets the cancellation responsible for equality.
For negative integer n, b>1 and d<0; the same algebra is an indefinite identity,
not a positive energy estimate. This is a concrete sign boundary.

**Why relative contraction is still insufficient.** The first L valuation
digits of \(n_0=2^{L+1}-1\) are all1, with
\(n_j=3^j2^{L+1-j}-1\). On this entire prefix, each b_j is so close to1
that every contiguous product exceeds1/2: the sum of their defects is at most
\(L/(3n_0+1)<1/2\). Consequently the relative renewal operator with unit
input/output weights has squared norm at least
\((L+1)(2L+1)/24\), by testing all-ones forcing. This refutes a uniform
unpaid all-prefix energy bound even after relative normalization. It does not
refute a paid first-exit or signed ROOT-current bound.

## R4. Changing modulus can pay for growing energy, without a uniform norm

**PROVED conditional extraction lemma; no ROOT conclusion assumed.** Recover
the existing R6 prior on the named odd source 2m+1:

\[
 \mu_m=(m+1)^{-r}-(m+2)^{-r},\quad r\ge1.
\]

All masses are positive rational numbers and sum to1. For q>m, the residue
cell C_q(m) contains m and its larger aliases. Its off-source mass obeys the
inherited bound \(E_q(m)\le(r+1)q^{-r-1}\). Let a response f_t at depth t
have a *proved*, possibly growing energy budget
\(\sum_s\mu_s|f_t(s)|^2\le G_t\). If
\(S_{q,t}=\sum_{s\in C_q(m)}\mu_s f_t(s)\), then

\[
 \left|\frac{S_{q,t}}{\mu_m}-f_t(m)\right|
 \le \frac{\sqrt{(r+1)G_t}}{\mu_m q^{(r+1)/2}}. \tag{7}
\]

This is Cauchy–Schwarz on the actual off-source atoms, not a uniform law on
residue lifts. It preserves the selected m and its positive mass.
For other supplied atomic masses the same proof uses their actual off-source
mass E_q(m), with error sqrt(E_q(m) G_t)/mu_m. Substituting the convenient
prior above for a different native ROOT measure would change the problem;
its tail rate must be proved for the original measure instead.

Given any certified finite G_t and any desired epsilon_t>0, choose a nested
dyadic q_t, divisible by q_(t−1), with q_t>m and

\[
 q_t^{r+1}\ge\frac{(r+1)G_t}{\epsilon_t^2\mu_m^2}. \tag{8}
\]

Then the extraction error is at most epsilon_t even if G_t grows without
bound. Taking epsilon_t=2^(-t) makes these errors summable. For example G_t=9^t
requires log2(q_t) at most t log2(36)/(r+1)+O_(r,m)(1), apart from any larger
initial precision. Increasing the modulus changes only observation precision;
the native operator and its integer clock stay fixed. Operator norm growth
can therefore be exchanged for address precision **for this aliasing error**.

There is a second, independent precision budget. A fixed Collatz valuation
word of length h and total valuation S has native source modulus 2^(S+1).
To also determine b endpoint bits in the odd-index coordinate, retain source
index precision q divisible by 2^(S+b), as follows from N6's exact affine
cylinder. The implemented modulus must meet **both** this branch/endpoint
condition and (8), as well as nest the preceding modulus. The dynamics clock
is h+S. Refining only to (8) does not by itself certify a native branch word.

The original sharper tail asymptotic is r*zeta(r+1)/q^(r+1). The odd-zeta
connection when r is even is inherited and exact, with no irrationality
assumption. The rational majorant in (8) suffices for the algorithm.

**Positive certificate with explicit remaining obligation.** If an independently
proved lower bound satisfies \(S_{q,t}>\sqrt{(r+1)G_t/q^{r+1}}\), then
f_t(m)>0. Neither a global signed energy budget G_t nor that positive cell
lower bound for the native ROOT observable follows from the new coordinate
choice. Equation (8) removes a fixed-precision obstruction; it cannot create
the directional sign that a trace or averaging operation discarded. The
universally controlled output in (6) is not automatically a ROOT score.

## T1. A compact relocation-closed extension of the two-port grammar

**PROVED for all orders.** A forest frame consists of:

1. named vertices and a bijective total rank r;
2. an undirected forest D on those vertices;
3. the actual active Hamilton path P and every marked external/shared port;
4. source labels, supplied source weights, and the existing operation clock.

The tournament is decoded exactly by

\[
 x\to y\quad\Longleftrightarrow\quad
 (r(x)<r(y))\ \mathbin{\mathrm{XOR}}\ (\{x,y\}\in D). \tag{9}
\]

It needs O(n log n) combinatorial bits for ranks, forest and path, excluding
the representation cost of externally supplied arbitrary weights and clocks.
It is a restricted exact representation, not a quotient of all tournaments.

**Old grammar inclusion.** A separable permutation supplies the total rank.
The old rule forces its consecutive presentation edges forward; precisely
the disagreeing ones enter D. Thus D is a subgraph of the presentation path,
and hence a forest. The former order-four escape remains represented after
re-marking because neither r nor D has to be rebuilt.

**Relocation closure.** Replace P by any Hamilton path of the same labelled
tournament, keeping r,D and every named mark unchanged. All arc queries and
future guards are still exact. If a presentation is reindexed, transport
every mark, rank and forest endpoint through the same bijection. Nothing
averages source weights or merges two shared ports.

For a single cut-and-insert relocation, at most three genuinely new arcs need
checking: the removed vertex's old gap and its two new neighbours. This is a
guard-evaluation bound, **not** a claim that an arbitrary re-marking costs one
elementary native move. The supplied operation clock is preserved separately.

**Binary closure.** Take disjoint framed A,B, with their active paths. For J,
put all A ranks before all B ranks and use D_A union D_B. All cross arcs are
forward. For star, put B ranks before A ranks and add the single defect edge
from the last active vertex of A to the first active vertex of B. Its cross
arcs are all backward except that hinge. The union of two forests and one
cross edge remains a forest. Both paths concatenate legally and all internal
arcs are unchanged. This recovers the literal star operation from N1, even
after prior relocations. Full converse reverses total rank and path while
retaining D and named ports. All statements follow directly from (9).

The four readout maps retain their earlier source-archive convention. I and
the singleton J/star reset constructions stay in this carrier; the old free-
edge complement F needs the extension below. This is not a closed four-map
compiler with a uniformly bounded exception budget.

### Why the extension is a port forest, not always two current endpoints

After re-marking, a future hinge can attach to a vertex that was interior in
the old expression. Repeated attachments can give an old subtree many boundary
ports. The global defect graph stays a forest and has at most n−1 edges, but
its local interface width need not stay two. Preserve the exact named endpoint
of every exceptional edge; do not replace it by the current interval endpoints.
Port incidence is retained, while the obsolete fixed-interval two-port claim
is intentionally not carried across these operations. No capacity-one routing
or unrestricted metagraph connectivity theorem is asserted.

### Complement boundary and its explicit budget

Let tau_P reverse every arc except those of P. The exact frame update is

\[
 r' = n-1-r,\qquad D'=D\mathbin\triangle E(P). \tag{10}
\]

Thus a general k-exception carrier is closed under this operation with
\(|D'|=|D|+(n-1)-2|D\cap E(P)|\). Relocations do not increase k. The same
formula allows signed cancellation of exceptions without inventing new ports.

The forest condition on this **specified frame update** fails already for
r=(0,1,2), D={{0,2}}, P=(0,1,2): D' is the triangle. This is not a proof that
the resulting tournament lacks *some other* forest frame. Reframing may repair
individual cases; the stronger obstruction below closes forest repair in
general. A uniform linear exception budget after arbitrary complement words
is **OPEN**. THM-4807 prevents using guard-indistinguishability to erase
genuine labelled states, but does not forbid compact restricted carriers.

## T2. Paley7 forces an actual cycle, under every rank choice

**PROVED.** The forest-frame class itself is not closed under legal path
complement, even allowing arbitrary changes of rank afterward.

Let P7 have x→y iff y−x lies in {1,2,4} modulo7, and use the Hamilton path
P=(0,1,2,3,4,5,6). Its predecessor S=tau_P(P7) has the forest frame

\[
 \text{rank order }(0,3,1,6,4,2,5),\qquad
 D=\{04,02,23,35,15,56\}.
\]

These six edges form a tree. Direct application of (9) gives S, and P is legal
in both S and P7 because the path edges were retained. The next tau_P sends
S to P7.

P7 has fourteen directed triangles and each edge belongs to two. This is the
old [two-Fano/torus-face mechanism](../camion_busch_gaps_polyhedra_collatz_20261002.md),
section6.4, now used as a rank obstruction; the [Paley ranking work](../paley_ranking_preregistration_20261017.md)
already recorded the natural order's optimum at7. Here is a self-contained
short proof of the required bound. The triangles are
{i,i+1,i+3} and {i,i+2,i+3}, i modulo7; inspection of the three nonzero edge
differences verifies the incidences. Every directed triangle has a backward
edge under every total order. Double counting gives

\[
 |D|\ge14/2=7>6.
\]

Thus no rank order gives P7 a forest exception graph. Natural rank order
attains7 with

\[
 D_{\rm nat}=\{03,14,25,36,05,16,06\}.
\]

It is connected and unicyclic: delete06 and it becomes a tree. Hence the
minimum exception-cycle rank of P7 is exactly1. One genuine cycle is the
smallest repair for this witness. This cycle rank is not identified with
the genus of the separate torus embedding.

**Arbitrarily many necessary cycles in one complement.** For every r≥1,
take the ordered join of r disjoint copies of the predecessor S and concatenate
their named paths. It still has a forest frame. One path complement gives
the iterated star of r copies of P7. Every rank order restricts on each copy
to a rank order of P7, so the backward-edge graph contains a cycle in each
of r vertex-disjoint blocks. These cycles are independent in its cycle space,
forcing beta≥r. Conversely, use natural ranks within each P7, reverse the
block ranks, and retain the r−1 actual hinge edges. The r unicyclic graphs
joined by bridges have beta=r. Thus the minimum exception-cycle rank of the
output is **exactly r on 7r vertices**, after a single legal complement of
a forest-encoded input. No constant cycle budget can close this carrier.

Nor can adding fresh helper vertices repair the forest claim while retaining
the old induced tournament: restriction of any forest frame is a forest frame,
contradicting the retained P7 block. This obstruction preserves the named old
vertices; it is stronger than a statement about a graph minor or path count.

**Reroute.** Allow an exception graph with a recorded cycle budget
beta(D)=|E(D)|−n+components(D). Relocation and full converse preserve beta;
J and star add the two inputs' budgets, since the latter adds a bridge.
For complement, (10) remains exact and

\[
 \beta(D\triangle E(P))\le\beta(D)+|E(P)\setminus D|.
\]

Deleting edges never increases cycle rank and adding one increases it by at
most1. This gives an honest compositional debt, not a uniform bound. The
seven-vertex case and its joins prove these units cannot always be removed by
reframing. The compact all-word target is now a linear bound on the *minimal paid
exception budget*, with rank changes and all named ports retained.

## M1. What a changing modulus can and cannot recover in Berggren dynamics

**PROVED corollary of THM-4808.** Its hostile word
\(AD^{M-1}UD\) sends (2,1) to

\[
 (m,n)=(16M+2,4M+1),\quad\chi=-1,
 \quad c=m^2+n^2=272M^2+72M+5. \tag{11}
\]

It has the root's endpoint residues modulo M but the opposite branch character.
Therefore an endpoint-residue-only character decoder for **all** pairs of
height c≤H fails whenever \(272M^2+72M+5\le H\), even if M was chosen from
H. Conversely M>sqrt(H) determines both positive coordinates exactly, so it
determines their unique Berggren word and its character. These bounds give
order-sqrt(H) necessary/sufficient modulus scale up to constants for this
uniform endpoint-only problem. They do not apply to projective slope alone.

A nested dyadic choice M_t>m_t gives an exact positive-coordinate lift along
the tree, since m_t strictly increases. It uses a height certificate; arbitrary
compatible residues alone do not provide one. A variable modulus with no
height information has not escaped the finite-stage ambiguity in (11).

The smaller alternative is inherited: retain the pole column (alpha,gamma)
alongside (m,n). Then chi=alpha*n−gamma*m; the full frame modulo3 already
distinguishes +1 and −1. Thus greater endpoint precision and retaining a
history coordinate are two distinct repairs. Neither changes the physical
branch clock or produces an estimate for a signed source-weighted sum.

## Connection contracts and next decisive tests

| Source → target | Exact map / preserved predicate | Loss or remaining price | Cheapest next test |
|---|---|---|---|
| Reset segments → weighted renewal | gain conjugation, cut product (1) | arbitrary weights may have unbounded B | evaluate B with the actual native ROOT weights, retaining source index |
| Coarse atoms → fine atoms | copy J and conditional mean E | zero-mean fine modes need D budget | measure A−JRE for an actual refinement |
| Collatz affine branch → passive two-port | (b,d), identity (5) | forcing cost ~3n/2; ROOT sign absent | identify the native current as f and account for its full cost |
| Atomic source prior → adaptive residues | (7)–(8), same m and clock | needs certified G_t and positive cell signal | use a genuine signed ROOT readout, not an unsigned surrogate |
| Old grammar → forest frame | exact edge decoder (9) | local two-port width may grow | search minimal exception frames after tau, retaining mark transport |
| Two Fano triangle systems → exception budget | each cycle needs a backward edge | triangle counts do not determine the full named path language | the explicit forest→Paley7 transition forces exactly one cycle |
| Berggren endpoint → bounded-height lift | (11) and M>sqrt(H) | precision grows with height | compare endpoint bits versus pole-frame sidecar on identical words |

The proposed reusable move is **retain the boundary carry, then pay it with
weights or precision**. Here it has exact but different realizations: a
rank-one renewal state, a named exception edge, and a source-alias tail. None
of these identifications makes the three major open problems equivalent.

## Reproduction and scope

Run `python3 verify.py > receipt.json` and
`python3 -O verify.py > receipt.optimized.json` from this directory. The checker
uses exact integers/Fractions and explicit exceptions, never removable assert.
See [AUDIT.md](AUDIT.md) for universes, hostile controls, and remaining gaps.
The manifest pins the artifacts, not external web pages or live mainline.
