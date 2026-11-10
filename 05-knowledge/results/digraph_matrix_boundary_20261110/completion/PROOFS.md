# Why exact tournament responses require general digraphs

**PROVED:** the named contraction and lifting laws below. **FINITE-EXACT:**
the stated exhaustive small controls and named examples. The full tournament
Hamiltonian-path spectrum remains **OPEN**. Local labels CT1–CT3 belong only
to this file. No external priority claim is made.

## Inheritance and the retained sidecar

The closest proved mechanisms are Session234
[DG6/DG7](../../session234_digraphs_20261109/PROOFS.md): weighted tournament
resolution and the endpoint response to an added arc. Session235
[PB2/PB3](../../session235_arc_pair_models_20261109/theory/PROOFS.md) already
defines the full pair-use polynomial, its deletion derivatives, and the
intersection/union of all tournament completions. The tile identity in
[THM-5041, matrix-tiling-pair-dictionary](../../../../01-canon/theorems/THM-5041-matrix-tiling-pair-dictionary.md)
is also inherited. Their scalar-overlap warning is not a new claim here.

The additional concrete objects are an order-four collision of the entire
resolution class list, an explicit contraction formula for each mixed
response, and its converse realization by a tournament with twice as many
vertices. These are elementary consequences of retaining the ordered ends
of every selected path block. Unlike ordinary contraction of a tournament
module, these blocks need not have uniform exterior incidences.

All paths below are ordered, spanning, and vertex-simple. Loops never enter;
H(empty)=1. For a digraph D and a set S of unordered vertex pairs, define
N_D(S) to be the number of Hamiltonian paths using every pair of S
consecutively, in whichever direction is present. Thus N_D(empty)=H(D).
This is Session235's mixed pair-use response, not a new path statistic.

## CT1. Endpoint contraction computes every mixed pair response

If the undirected graph (V,S) is not a linear forest, then N_D(S)=0.
Indeed, the consecutive unordered pairs of a vertex-simple Hamiltonian path
form one undirected path, so they contain neither a cycle nor a vertex of
degree at least three.

Now suppose (V,S) is a linear forest. Include its isolated vertices as
singleton blocks. Choose a direction along each nontrivial component; there
are two choices per component, not two choices independently per edge.
Discard a choice if one of its directed internal edges is absent from D.
For each surviving choice alpha, let its ordered blocks be

    B_i=(b_i,1,...,b_i,k_i).

Define a loopless quotient digraph Q_alpha by

    Q_alpha[i,j]=D[b_i,k_i, b_j,1],  i != j.                (CT1)

Then

    N_D(S) = sum_alpha H(Q_alpha).                        (CT2)

**Proof.** A Hamiltonian path counted by N_D(S) traverses each selected
component as one consecutive block, in one of its two directions. Its block
order is a Hamiltonian path of Q_alpha by the definition of its interblock
arcs. Conversely, expanding the ordered blocks of a quotient path gives
exactly one Hamiltonian path of D using S. The two operations are inverse.
Every orientation is counted once, including singleton components, which
have only one choice. QED.

Either a missing pair or a reciprocal pair can appear in the quotient even
when D is a tournament. Each phenomenon already occurs at order three:

- In TT3 with order0→1→2, contract the selected edge0→2. The blocks(0,2)
  and(1) have no interblock arc in either direction. The constrained count
  and quotient H are both0.
- In C3 with0→1→2→0, contract0→1. The blocks(0,1) and(2) form a digon:
  1→2 and2→0. The two constrained paths are012 and201, so the quotient H is2.

The second example explains why declaring11 to mean an empty edge cannot
be used inside this contraction rule: it changes the correct answer2 to0.
The symbol interpretation must be specified before counting.

For the inherited polynomial P_O, its response
N(S)=(partial_S P_O)(1) is therefore an explicit sum of smaller ordinary
digraph counts. The dependence on whole selected sets is essential; a
cycle or branching set has zero response even when each separate pair has
positive response. This makes the larger digraph category a closure class
for exact tournament path constraints.

## CT2. Every digraph is such a quotient of a tournament

Let D be any loopless digraph on m named vertices. Replace vertex i by the
two named vertices i_in=2i and i_out=2i+1. Construct a tournament T_D by:

1. Orient i_in→i_out for every i.
2. For i!=j, orient i_out→j_in exactly when D[i,j]=1; otherwise orient
   j_in→i_out.
3. Orient all remaining in/in and out/out pairs arbitrarily. The explicit
   checker uses increasing physical names.

These prescriptions are compatible. The pairs {i_out,j_in} and
{j_out,i_in} are distinct physical pairs, so D[i,j] and D[j,i] can be chosen
independently. Put F={i_in→i_out : 0<=i<m}. Then

    #{P in Ham(T_D): every arc of F occurs in P} = H(D).  (CT3)

The same is true of N_T_D applied to the underlying matching pairs because
each has only its prescribed orientation in T_D.

**Proof.** If every selected arc occurs, each(i_in,i_out) is consecutive.
Contract these m directed blocks. The quotient entry from block i to block j
is exactly D[i,j]. Expansion and contraction are inverse on Hamiltonian
words, as in CT1. QED.

If m>=2 and D is strong, then T_D is strong: directed walks in D lift through
the prescribed block arcs, and every vertex can enter and leave its block
along such lifted walks. This assertion does not require a Hamiltonian cycle
in D and does not assert the converse.

**Count universality has a restricted target.**
[THM-5071(4a), weighted-cycle-formula-for-digraph-path-totals-and-one-reciprocal-pair-exclusions](../../../../01-canon/theorems/THM-5071-weighted-cycle-formula-for-digraph-path-totals-and-one-reciprocal-pair-exclusions.md)
already realizes every h>=1 by a semicomplete digraph. For h>=2 its sink
back-arc construction uses m=ceil(log2 h)+1 vertices: the binary digits of
h-1 select the back arcs. The highest digit supplies the back arc from the
sink to0, so this particular D is strong. Combining it with CT2 realizes h
as a selected-matching response of a strong tournament on
2ceil(log2 h)+2 vertices. This bound is not claimed minimal. The value1 has
the two-vertex tournament with its sole arc selected; zero has the lift of
an edgeless pair. Thus every nonnegative integer is a selected-matching
response of some tournament.

This does **not** realize every integer as the full tournament count.
In the two literal inherited counterexamples the distinction is:

| Input digraph H | Tournament order | Selected matching response | Full tournament H |
|---:|---:|---:|---:|
|7|8|7|79|
|21|10|21|2621|

The receipt gives the complete named tournament matrices. The matching is
(0,1),(2,3),... in each case. Arc-deletion inclusion–exclusion independently
checks the constrained totals. The full counts are odd and are neither7
nor21, as required by the inherited tournament exclusion theorem.

## CT3. Resolution isomorphism classes lose the needed overlap

Consider two strong semicomplete digraphs, written as row bitmasks:

    D6=(6,9,10,1),       D8=(6,5,8,3).

Both have exactly one reciprocal pair, {0,1}. Each of their two tournament
resolutions has H=5. More strongly, all four resolved tournaments are
isomorphic to the unique strong order-four tournament. Nevertheless,

    H(D6)=6,  H(D8)=8.

Deleting both arcs of {0,1} leaves4 and2 Hamiltonian paths respectively.
These are exactly the intersection sizes of the two resolved path sets.
Consequently the inherited inclusion–exclusion identity gives

    5+5-4=6,       5+5-2=8.

The receipt supplies every Hamiltonian word and explicit canonical relabeling
certificates for the four resolutions. Thus even an indexed list of the
resolved tournament *isomorphism classes*, as well as their counts, can
lose the cross-resolution identification of paths. Keeping the complete
labelled matrices would retain that identification and recover the digraph
by their union; that stronger representation is not contradicted.

Order four is minimal for this one-reciprocal-pair scalar collision.
Fixing that pair as{0,1} loses no cases by relabeling. Exhausting the one
order-two and four order-three possibilities shows that each indexed
resolution-count vector there has only one possible union count. There
is no order-one graph with a reciprocal pair. Since full resolution classes
determine their H-values, this also proves minimality for the class-list
collision under the same one-pair hypothesis.

For arbitrary a,b>=0, the ordered sums TT_a→D6→TT_b and TT_a→D8→TT_b keep
these H-values and matching resolution-class lists. Every padding vertex
is forced into its ordered singleton position; the core isomorphisms extend
by fixing the padding. This proves such overlap losses at every order>=4.
The original two cores are strong; the padded graphs generally are not.

## Exact scope, controls and transfer contract

`check.py` uses two independent counting methods: literal permutation words
and endpoint subset DP. It checks CT1 on every loopless labelled digraph of
orders0 through4 and every selected unordered-pair set:4,166 graphs and
262,666 responses, including122,944 cycle/branching-support zeros. These
are unrestricted four-state matrices, not only tournaments or semicomplete
ones. It checks CT2 by literal constrained words on all70 digraphs of
orders0 through3, and separately by arc-deletion inclusion–exclusion.
The order8/order10 target lifts use that latter independent check. Eighteen
padding controls, legal-orientation controls, and cycle/degree-three hostiles
are included. Normal and optimized receipts are byte-identical.

Reproduce from the repository root:

    python3 -B 05-knowledge/results/digraph_matrix_boundary_20261110/completion/check.py > /tmp/completion-normal.json
    python3 -B -O 05-knowledge/results/digraph_matrix_boundary_20261110/completion/check.py > /tmp/completion-optimized.json
    cmp /tmp/completion-normal.json /tmp/completion-optimized.json

The connection contract is exact:

| Source → target | Preserved predicate | Lost information | Needed sidecar |
|---|---|---|---|
| D with selected path forest → endpoint quotient(s) | paths containing all selected pairs, by bijection | unselected paths; other incidences of internal block vertices | ordered block names, their directions, and selected forest |
| arbitrary D → T_D with selected matching | H(D) equals the constrained matching response | equality with full H(T_D) | the entire selected matching |
| semicomplete D → resolved tournament class list | each resolved class and its individual H | shared identity/intersection of Hamiltonian words | labelled resolutions or joint pair-use response |

The first two rows explain why arbitrary digraphs can help evaluate exact
tournament responses. The third explains why counting their individual
isomorphism classes, or matching sequence values, does not by itself solve
the tournament spectrum problem.
