# Three strong seed states and the full named observer

**PROVED:** the scalar-to-named-response implications and the sharp five-run
bound below. **FINITE-EXACT:** the six retained complete response polynomials.
No spectrum completeness, full-response collision, or prime-core exclusion
is claimed.

## Incoming results and the actual connection

The seed theorem and its independent audit are
[Session185 alternative](../session185_proof_frontiers_20261104/alternative/PROOF.md)
and [alternative audit](../session185_proof_frontiers_20261104/alternative_audit/AUDIT.md).
Use physical names0..5 and outgoing masks

    Q0=(16,5,33,39,46,3),
    Q1=(16,5,49,39,42,3),
    Q2=(18,4,49,39,42,3).

The only changes are reversal of4->2 and then1->0. Each contains the cycle
0->4->3->1->2->5->0. Replace5 by TT_r on names5..r+4,
in increasing order, and call the result T_j(r). Then, for all r>=1,

    H(T_j(r))=4*3^r+9*2^r-3+2j, j=0,1,2.          (8)

These are three consecutive odd values, not three consecutive integers.
The nonzero three-run coefficient4 and two-run coefficient21 come from
the complete 120-exterior-order gap enumeration in the incoming proof:

    H(Q_j[A])=(27+2j)H(A)+21q_2(A)+4q_3(A).         (9)

For TT_r, q_1=1, q_2=2^r-2, q_3=3^r-3*2^r+3.
Substitution in (9) gives (8), retaining both exponential terms.

Reversing each path in each fixed ordered slot proves
q_k(A)=q_k(A^op). Thus all three scalar states in (9) are blind to module
converse. The same is true after **any number of uniformly attached ports**,
as proved by incoming **G2** in
[Session186 gate](../session186_proof_frontiers_20261104/gate/PROOF.md).
That statement concerns the restricted uniform query language.

If A has a homogeneous pair, **P3** in [the direct-query proof](PROOF.md)
separates the two actual full ambient responses. In particular, for every
r>=2 the adjacent names5,6 in TT_r are homogeneous, so their two calibrated
queries separate T_j(r) from Q_j[TT_r^op] and from its whole converse.
The scalar three-point interval therefore coexists with exact named
orientation detection at every r>=2.

Incoming **G1** makes the spectral connection quantitative with another
measurement. For r>=3, let a=5,b=6,w=7. Then

    C_bw-C_aw+C_wa-C_wb
      =2H(Q_j[TT_(r-2)])
      =8*3^(r-2)+18*2^(r-2)-6+4j >0.              (10)

Both whole and internal converses give its negative. Contract the first
pair to z, apply C2, then contract the remaining consecutive pair z,w;
each of the two cover terms becomes the Hamilton count in (10).
This is inherited from G1, not a new formula claimed here.

Session186 also contributes **G3**'s stronger scalar hostile: modules
A=(4,1,2,7,15), B=(10,0,3,6,15) have the same entire ordered-cover
profile (3,42,162,240,120), although only A has a homogeneous pair.
Every scalar uniform substitution context agrees, but a named two-port
query separates them, including inside all three strong seeds.
This prevents upgrading scalar profile equality to twin detection or
full response equality. B is pair-free, not prime.

## Two nonuniform ports raise the run cap from three to five

**Theorem I1.** In any extension of Q_j[A] by two arbitrarily attached
ports x,y, a Hamilton word has at most five maximal A-runs. This bound is
attained in all three seed states.

**Proof.** Each nonterminal A-run must be followed by a distinct vertex
outside A which it beats. Among old exterior vertices only0 and1 are
available. The ports x,y supply at most two further followers. Thus there
are at most four nonterminal runs and one terminal run.

For sharpness take A=TT_5 on5..9, x=10,y=11, x->y, and attachment masks
a=32=e_5, c=64=e_6. The word

    5,x,6,y,7,0,4,3,8,1,2,9                        (11)

is Hamilton in every state. Its five A-runs are the singletons
(5),(6),(7),(8),(9). The port edges are exactly those prescribed; every
other displayed edge is common to the seeds. QED.

The integer counts by maximal A-run number under this same query are:

| State | 1 run | 2 runs | 3 runs | 4 runs | 5 runs |
|---|---:|---:|---:|---:|---:|
| 0 |27|1046|5236|4896|576|
| 1 |29|1034|5128|4824|576|
| 2 |31|1046|5140|4824|576|

These are counts of actual words with nonuniform contacts. They are not
universal scalar coefficients multiplying q_k(A). The five-run word
does not contradict the usual at-most-three *old paths* after deleting
two ports: old paths can themselves pass through exterior vertices and
visit A repeatedly. These are different decompositions.

More generally m ports give the upper bound m+3 for this seed family,
by exactly the same follower argument. No general sharpness assertion
for every m or every module is required here.

## Full named cubic contrast, with the full responses retained

[response_tables.json](response_tables.json) retains the **entire**
Boolean multilinear F and F_op for j=0,1,2 at r=2 and r=7: six fixed
cores of orders7 and12. It does not quotient old names by any automorphism.

For each core the ordered mask list is zero, all singleton masks, then
all unordered pair masks. Raw values are row-major in (a,c).
The same indexing stores coefficients of
product_(u in a) a_u * product_(v in c) c_v, including zero coefficients.
The artifact stores F, F_op, and their exact difference separately.
Because each port has at most two incidences in a Hamilton word,
bidegree(2,2) is a complete polynomial representation of each full response.
Shared labels across the two ports are retained as distinct variables.

At order12 this is6241 values and coefficients per response. Its whole
converse contrast has no quartic or pure quadratic terms and is exactly
reconstructed from the1753-query battery of **T4**,
[twelve principles](../two_port_twelve_principles_20261103/PROOF.md).
The nonzero contrast coefficients by degree are:

| r | State | Linear | Mixed quadratic | Cubic |
|---:|---:|---:|---:|---:|
|2|0|14|42|286|
|2|1|14|42|294|
|2|2|14|42|290|
|7|0|24|132|1568|
|7|1|24|132|1576|
|7|2|24|132|1576|

For example at r=7,state0 the two calibrated queries at names5,6 are
(41763,53317); for the whole converse they are (148509,110585);
for the internal module converse they are (220857,162937).
The Hamilton count is9897 in all three cases. Thus the full cubic has
an actual named separation while preserving the incoming scalar clock.

Do not transfer this cubic reconstruction to a partial module converse.
Incoming [Session186 structure](../session186_proof_frontiers_20261104/structure/PROOF.md)
and its [independent review](../session186_proof_frontiers_20261104/structure_review/REVIEW.md)
give an exact contact formula and a sharp five-vertex hostile:
T=(24,17,19,6,8), T'=(8,17,19,6,9), reversed pair0,4, with

    [a_1 a_3 c_1 c_4](F_T-F_T')=2.                 (12)

The present counter replays its actual16-query certificate. Its quartic
coefficient is the deletion-cover skew described in that incoming proof.
Low-degree measurements can distinguish the pair while quartic terms
remain necessary to reconstruct its entire ambient response.

## Exact reproduction and limits

Run `python response.py` or

    python -O response.py --output response.optimized.json --tables response_tables.optimized.json

The C++ predecessor-subset DP evaluates42,668 actual attachments, including
all complete interpolation batteries for both orientations, six high-weight
controls per response,18 internal-reversal pair controls, and the quartic
hostile. Python forward subset DP independently checks24 queries on the
six cores, literal permutations check three order-nine extensions, and
a run-count DP verifies (11) and the entire run table. The scalar law is
separately checked at all r=1..7 and j=0,1,2.

These are finite controls for the cited all-order results. Neither six
full response polynomials nor the scalar interval enumerate prime
twelve-vertex cores. The prime hidden-converse branch remains OPEN.

