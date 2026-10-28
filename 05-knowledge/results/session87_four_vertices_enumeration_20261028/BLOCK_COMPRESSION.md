# One-active-block compression is a partial-spread problem

**Session 87 — 2026-10-28. Status:** elementary criteria and query bounds
**PROVED** below; partial-line-spread capacity **CLASSICAL / CITED**, with a
self-contained binary proof; the accompanying finite controls are reported
separately. The partial-spread, block-distance and field-slope mechanisms are
classical; this is their explicit port-query application and verification,
not a claim of a new partial-spread or coding theorem.

## Inheritance and the question

Session 86, `05-knowledge/results/session86_relations_attachment_carriers_20261028/NORM_PORTS.md`
and `norm_ports.py`, identifies its F64 norm-one carrier with the nine nonzero
vectors in three coordinate planes. In its declared observation models the
costs are six fixed normalized-partner queries, five adaptive such queries,
and four fixed unrestricted linear queries. Its 651-plane census finds 36
optimal kernels. We generalize the *declared observation problem*, not the
multiplicative norm law, arithmetic labels, or tournament edges.

The live concept board is: one-active-block states; query vocabulary and
adaptivity; partial spreads; kernel block-distance; and operation closure.
The hostile is a pair of individually injective blocks with the same image.
The least-used sidecar is the full difference set, rather than only U itself.
The inherited near miss is the assumption that a counting bound is attainable
by a linear map. The relevant meta-patterns are “search the statement before
the method” and “correct the object before sharpening the technique.”

## 1. Exact criterion, including the exception

Let r,k be positive integers, let V=V_1 direct-sum ... direct-sum V_r over
F2, with dim(V_i)=k, and set

    U = union_i (V_i minus {0}).

A fixed s-query binary linear observation is a linear map L:V -> F2^s;
redundant output coordinates may be deleted. Define block-weight(x) to be
the number of its nonzero coordinate blocks.

**Theorem 1 (k>=2).** The following are equivalent:

1. L is injective on U.
2. Each restriction L|V_i is injective and the k-dimensional image spaces
   W_i=L(V_i) meet pairwise only in zero.
3. Every nonzero element of ker(L) has block-weight at least three.
4. L is injective on U union {0}.

**Proof.** If 0!=w in V_i is killed, choose a in V_i outside {0,w}; then
a and a+w are distinct nonzero states with equal outputs. This uses k>=2.
Thus (1) forces injectivity on each block. A nonzero vector in W_i intersect
W_j would then be the common output of nonzero states in different blocks.
Conversely, injectivity within and between the blocks proves (1). For k>=2
the set of differences of distinct U states is exactly the set of vectors
with block-weight one or two: a weight-one w=a+(a+w) as just above, and a
weight-two vector is the sum of its two active blocks. Avoidance of this set
is precisely (3). Blockwise injectivity also makes every U output nonzero,
so adding the zero state makes no change. QED.

Thus the images form a **partial k-spread** in the binary output space, or
equivalently a constant-dimension subspace code with minimum subspace distance
2k. Conversely, any such family, together with chosen isomorphisms V_i->W_i,
defines a valid L by summing the component images. This is an exact reduction,
not merely a necessary condition.

**k=1 exception.** U consists of the basis vectors e_1,...,e_r. Its outputs
need only be distinct columns, and one column may be zero. Therefore the
minimum s is ceil(log2 r), including s=0 when r=1. With the zero state included,
the columns must instead be distinct and nonzero, and the minimum becomes
ceil(log2(r+1)). Without zero, the kernel must avoid weight-two vectors, but
may contain a weight-one vector. Do not apply Theorem 1 to this exception.

For r=1 and k>=2 the minimum is k. No undefined “minimum distance of the
zero code” convention is needed: throughout, the kernel condition is a
universal condition on its nonzero elements, hence vacuous for kernel {0}.

## 2. General lower bounds and a sharp range

Theorem 1 immediately yields, for k>=2,

    r(2^k-1) <= 2^s-1.

If r>=2, any two image spaces have a 2k-dimensional direct sum, so the
additional bound s>=2k holds. Counting alone can be much weaker: two
3-dimensional blocks contain only 14 allowed states, but need six fixed
linear queries rather than four.

Identify F2^k with the field E=F_(2^k). In E^2 the spaces

    W_a = {(x,a*x): x in E}, for a in E,
    W_infinity = {(0,x): x in E}

have dimension k over F2 and meet pairwise only in zero. They partition all
nonzero output vectors. Choosing any r of them proves the exact minimum

    s_min(r,k)=2k,  whenever 2 <= r <= 2^k+1 and k>=2.

For larger r the problem remains the classical partial-spread existence
problem; no universal sufficiency of the counting inequality is claimed.

**Three-block kernel classification.** For r=3 the optimal map has rank 2k,
so its kernel K has dimension k. Its nonzero vectors have all three blocks
nonzero. Each coordinate projection K->V_i is therefore an isomorphism.
Identifying K with V_1 gives exactly

    K = {(x,A*x,B*x): x in F2^k},   A,B in GL(k,2).

Every such graph works. Consequently the number of optimal kernels for a
fixed labelled block decomposition is |GL(k,2)|^2. At k=2 this is 6^2=36,
explaining the Session 86 census. These are kernels, not ordered probe lists
and not isomorphism classes of kernels. A convenient representative is

    L(x,y,z)=(x+y,x+z),   K={(a,a,a):a in F2^k}.

At k=2, every good kernel has a four-dimensional annihilator in the dual
space. Its ordered bases number |GL(4,2)|=15*14*12*8=20,160. Consequently
there are 36*20,160=725,760 ordered optimal four-query lists, or 30,240
unordered four-query sets. A nondegenerate pairing identifies these forms
with arbitrary probe partners, so the same counts apply to Session 86's
unrestricted F64-partner model. These counts do not apply to its restricted
nine-partner model.

## 3. Exact answer for two-bit blocks, in every output dimension

Let M(s) be the maximum number of pairwise trivially intersecting
two-dimensional spaces in F2^s, for s>=2. Then

    M(s) = (2^s-1)/3,  s even;
           (2^s-5)/3,  s odd.

This is the classical partial-line-spread formula, recorded for general q
in Section 2.1, equation (2), of Honold--Kiermaier--Kurz,
[*Optimal Binary Subspace Codes of Length 6, Constant Dimension 3 and Minimum
Subspace Distance 4*](https://arxiv.org/pdf/1311.0464), v2, 2014. The following
independent binary proof makes the obstruction and construction explicit.

**Upper bound.** Every two-dimensional space contains three nonzero vectors.
When s is even the counting bound gives the displayed answer. When s is odd,
the counting bound would allow (2^s-2)/3 spaces, leaving exactly one nonzero
vector uncovered. But the XOR of all nonzero vectors of F2^s is zero for
s>=2, and the XOR of the three nonzero vectors of every plane is zero. The
XOR of the uncovered set must be zero. A singleton nonzero hole is impossible,
so at least four holes remain, giving (2^s-5)/3.

**Construction.** The bases s=2 and s=3 each admit one plane. For s>=4 put
m=s-2 and write the output space as F2^m direct-sum F2^2. Choose a binary
linear transformation T on F2^m for which both T and T+I are invertible.
Such a T exists for every m>=2: use direct sums of the companion blocks

    dimension 2: T(e0)=e1, T(e1)=e0+e1;
    dimension 3: T(e0)=e1, T(e1)=e2, T(e2)=e0+e1.

Every m>=2 is a sum of 2s and at most one 3, and each block and its translate
by I is invertible (their characteristic polynomials have neither 0 nor 1
as a root). For every a in F2^m take

    G_a = span{(a,1,0),(T(a),0,1)}.

These 2^m graph-planes are mutually disjoint outside zero: in a collision
with common nonzero final coordinate (u,v), their first-coordinate
difference is (u I+v T)(a+b), forcing a=b. The only three possible nonzero
operators are I,T,I+T, all invertible. The graph-planes also avoid the first
summand outside zero. Adjoin an optimal partial spread inside that summand.
This constructs 2^(s-2)+M(s-2) planes, exactly the displayed formula. QED.

Thus the minimum fixed linear query count is the least s>=2 with r<=M(s):

| blocks r | states 3r | minimum fixed unrestricted linear queries |
|---:|---:|---:|
| 1 | 3 | 2 |
| 2 | 6 | 4 |
| 3 | 9 | 4 |
| 4 | 12 | 4 |
| 5 | 15 | 4 |
| 6 | 18 | 5 |
| 7 | 21 | 5 |
| 8 | 24 | 5 |
| 9 | 27 | 5 |
| 10 | 30 | 6 |

Counting is insufficient at r=2,10,42,170,..., namely
r=(2^s-2)/3 for odd s>=3. These are precisely the one-hole cases.

### Explicit five-bit witness for six through nine blocks

In integer bit masks for F2^5, the nine triples of nonzero output vectors are

    {8,16,24}, {9,18,27}, {10,20,30}, {11,22,29},
    {12,19,31}, {13,17,28}, {14,23,25}, {15,21,26}, {1,2,3}.

Each triple has XOR zero and no two share a point. The holes are {4,5,6,7}.
Use any r triples, choose two vectors from each as the images of that block's
basis, and extend linearly. This explicitly supplies all r=6,...,9 cases.

## 4. The observation models must stay separate

A **block-local query** is any binary linear functional supported on one
input block. This definition works for all k; it does not posit a nondegenerate
alternating form on odd-dimensional blocks. When each block carries a
nondegenerate pairing, all such functionals can be realized by block-local
partners. Session 86's orthogonal binary planes are a special case.

For k>=2, the exact fixed block-local cost is rk: every block restriction
must have rank k, and local queries contribute rank only to their own block.

The exact adaptive block-local worst-case cost is

    k,      r=1;
    rk-1,   r>=2.

**Lower bound.** Follow the all-zero answer branch. If the queries made in
block i have rank t_i, that block retains 2^(k-t_i)-1 nonzero candidates.
Before this branch is a leaf, it has at least two candidates and its zero
child cannot be empty. At a leaf,

    sum_i (2^(k-t_i)-1) <= 1,

forcing sum_i t_i >= rk-1, hence at least rk-1 queries. For r=1 the binary
counting bound on 2^k-1 states instead gives k.

**Upper bound for r>=2.** Ask k-1 independent basis queries in each block,
in order. If any answer is one, the active block is identified; finish its
remaining basis queries. This uses at most (r-1)(k-1)+k=r(k-1)+1<=rk-1.
If all answers are zero, exactly one candidate remains in each block; test
at most r-1 of these candidates with the omitted local functional. Total:
r(k-1)+(r-1)=rk-1. QED.

For r=3,k=2 this is precisely fixed-local 6, adaptive-local 5, and
fixed-unrestricted-linear 4. These are two independent changes: permitting
adaptivity and enlarging the query vocabulary. There is no universal ordering
between adaptive-local and fixed-unrestricted costs: at r=2,k=2 they are
three and four respectively. Nor does the fixed-linear answer quantify the
optimal adaptive strategy using unrestricted queries.

## 5. Operation closure is the load-bearing boundary

**Four bits identify the nine declared states; they do not losslessly encode
all six input bits.** The ambient map has a two-dimensional kernel. U is not
additively closed. The map preserves addition of outputs, but after adding
inputs the result can leave U, where injectivity was proved.

For example, for any nonzero a,

    L(a,a,0) = L(0,0,a)

under L(x,y,z)=(x+y,x+z). One state is two-block-supported and the other is
one-block-supported. This is a concrete first failed implication, not just a
warning about possible information loss.

More generally, let C_t={x: block-weight(x)<=t}, including zero. A linear
map is injective on C_t **iff** every nonzero kernel vector has block-weight
at least 2t+1. Indeed the difference set C_t+C_t is exactly C_(min(2t,r)), by
partitioning the support into two parts of size at most t. This is the usual
syndrome-decoding distance criterion in the block metric. In particular,
full additive closure is V and requires all rk bits. A closed nonlinear
operation on U may be transported through the bijection L|U, but it must be
defined by decoding, applying that operation, and re-encoding; it is not
automatically XOR, multiplication in the output field, or a graph operation.

**Connection ledger.** Source: U with its labelled block decomposition.
Target: the disjoint nonzero parts of output k-spaces. Map: L as above.
Preserved predicate: equality/identity of declared U states, and the explicitly
transported operations if supplied. Destroyed information: ambient kernel
cosets outside U. Needed sidecar: the source restriction and chosen block
isomorphisms, plus a lawful native operation if composition is requested.
Cheapest hostile: the displayed one-block/two-block collision.

At five two-bit blocks the 15 outputs exhaust F2^4 minus {0}; this is an
honest connection to the 15-point binary projective carrier. It is **not** an
identification with five norm-one F16 ports, seven Fano points, the
Hamming[7,4,3] kernel, or the Hamming[15,11,3] kernel merely because a number
matches. Each such transfer needs its own source, target, and operation map.
The companion K4 lane supplies a separate rank-four incidence/parity
realization only after choosing the appropriate adjacent-edge pairing
scaffold; opposite-edge matchings are a hostile, not a default realization.

## 6. Reproduction and audit contract

**Run status:** the normal and optimized-Python commands below passed with
byte-identical JSON receipts. An independently coded JavaScript calculation
also matched the subspace counts, all 36 diagonal kernels, spread sizes
through dimension 10, the explicit nine-line witness, and all listed adaptive
minimax values. The general theorem/proof and closure boundary received an
independent adversarial reasoning check; they do not rely on extrapolation.

Run from this directory:

    python3 check_block_compression.py > block-compression.json
    python3 -O check_block_compression.py > block-compression-optimized.json
    cmp block-compression.json block-compression-optimized.json
    node check_block_compression_independent.js > independent-javascript-receipt.json

The checker declares its full finite universes. It compares the three
criteria on every binary matrix in five small parameter sets; enumerates
every F2^6 subspace of dimensions 0,1,2,3 with no orbit filter; compares the
36 good planes against all twisted-diagonal GL(2,2)^2 graphs; verifies the
recursive maximal-spread witnesses for output dimensions 2 through 10 and
all r<=10 codes; independently solves exact candidate-subset minimax for
small block-local adaptive models; and checks the k=1, cross-block collision,
and additive-closure hostiles. In the latter, the 37 at-most-two-block states
in F2^6 have only 16 outputs, and the sparse-difference identity is checked
for t=1,2,3. Explicit exceptions remain active under python -O.

The general proofs, not extrapolation from this finite run, establish the
all-parameter results. No exhaustive classification of all maximal partial
spreads or all output maps at the larger parameters is claimed.
