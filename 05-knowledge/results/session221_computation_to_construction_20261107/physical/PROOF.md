# Every even single-arc increment: a physical exclusion-tape construction

Status: PROVED elementary identities and constructions; FINITE-EXACT replay.
Independent audit: PASS for the main physical theorem; see
[physical audit](../audit/PHYSICAL_AUDIT.md). No literature-priority claim.
Chan–Pak remains OPEN.
This packet produces prescribed **differences** of Hamiltonian path counts,
not a tournament with an arbitrary prescribed absolute count.

## Main result

For every integer N>=1 there are strongly connected tournaments U_N,V_N
of the same order

    5 + ceil(log_2 N)

such that V_N is obtained from U_N by reversing **one arc**, and

    H(V_N)-H(U_N)=2N.

For N>=2 their exact counts are

    H(U_N)=12N-9-f(N-1),   H(V_N)=14N-9-f(N-1),                 (A)

where

    f(0)=0,  f(2u)=f(u),  f(2u+1)=f(u)+u+(u mod 2).

For N=1 the counts are 9 and 11 on five vertices. The displayed order is
optimal **within the adjacent-chain one-apex flip construction defined below**,
not among all tournament arc flips. Without requiring strong connectivity,
N>=2 has order 4+ceil(log_2 N) and counts

    4N-3-f(N-1),   6N-3-f(N-1).                               (B)

For N=1 there is a three-vertex pair with counts 1 and 3.

## Inheritance and concept board

Closest mechanism: the October 20 packet
`05-knowledge/results/feedback_word_hamilton_spectrum_20261020/proofs.md`,
Sections 2--4. It already proves the weighted-inversion count, the exact
adjacent-exchange response, strict rearrangement, strong-interior exclusion
of increment 2 from order 6, and the sorted-word extremal count. Thus
(1), (3), the monotonicity consequence, and (15) below are **KNOWN**, not
new discoveries. We repeat their proofs only to make the constructor
self-contained. The new parallel formula (13) is not needed for the
all-even selector. Another canonical route is
THM-4876, `01-canon/theorems/THM-4876-deletion-and-flip-increments-of-the-hamiltonian-path-count-an-increment-of-two-forces-at-most-three-paths-the-increment-ladder-thm-013-for-all-orders-and-what-extension-spectra-forbid-for-chan-pak.md`,
part (6). Its general arc-flip decomposition, part (4), is also inherited.
The exact digital evaluator is already in
`05-knowledge/results/session133_fractran_spectrum_refinement_20261101/path_counter/THEORY.md`:
its bit convention is opposite to this packet's. We do not claim a new
one-apex transfer, a new general arc-flip formula, or a new evaluator.
The binary f register is inherited from Session216 three-core numeration.

Canonical hostile: equal H does not determine continuation. Corrected near
miss: a machine evaluating H is not a constructor selecting a graph for a
target. Least-used sidecar: the two physical endpoints of an adjacent tape
exchange and their independent left/right binary charges.

Concept board: one-arc surgery; exclusion cellular automata; dyadic response
registers; exact fixed-order difference coverage; absolute-count sidecar.
The proposed new scoped content is the complete marked-gap interval, its explicit
all-even inverse, the optimal order within this encoding, the strong lift,
and the parallel interaction identity. THM-4204's spatial inverse monoid and
THM-4210's marked Rule30 current do not imply these constructions; no Rule30
universality or finiteness assertion is used. THM-4857/4858 are different
vertex-automaton and TM-simulation carriers. The operation here is the actual
number-conserving radius-one exclusion rule, not an abstract matrix evaluator.

Input reference: main commit 38fdf9f234159aaedb2d0ba1e0f7bcdefb49e817.
Shared repository was read only; all artifacts are owned session files.

## 1. Physical tape and the inherited counting proof

For a word x=(x_0,...,x_{n-1}) define T(x) on apex a and chain vertices
v_0,...,v_{n-1}. Every chain arc is v_i -> v_j for i<j. Set

    x_i=1 iff v_i -> a;  x_i=0 iff a -> v_i.

This completely specifies a tournament. Define

    w(1)=1,  w(d)=2^(d-2) for d>=2,
    C(x)=sum_{i<j}(1-x_i)x_j w(j-i).

The inherited count is

    H(T(x))=1+2C(x).                                          (1)

For completeness, here is a direct proof, not an appeal to a cycle formula.
A Hamiltonian path is uniquely P,a,Q, where P and Q are increasing lists
partitioning the chain. If both lists are nonempty, let p=max P,q=min Q.
For p>q, legality requires x_p=1,x_q=0 and the p-q-1 intervening vertices
may independently join either list, giving 2^(p-q-1) choices. For p<q,
necessarily q=p+1 and the legal cuts have pattern 10. The two extreme cuts
contribute 1-x_0 and x_{n-1}. Thus the total cut contribution is

    1-x_0+x_{n-1}+#10 = 1+#01.

The p>q contribution at distance1 supplies one additional #01, and at
larger distances supplies 2w(d). This proves (1). At n=0 the single apex
has one path directly.

For n>=2, T(x) is strong iff x_0=0 and x_{n-1}=1. These bits give the
Hamilton cycle a,v_0,...,v_{n-1},a. Conversely x_0=1 makes v_0 a source,
and x_{n-1}=0 makes v_{n-1} a sink.

## 2. Local exchange is a single physical arc flip

Suppose x_i x_{i+1}=10 and let x' replace this pair by 01. The literal
T(x),T(x') differ in two apex arcs. However, reverse only the chain arc
v_i -> v_{i+1} in T(x), obtaining F_i(T(x)). Relabel v_i and v_{i+1}.
Because the chain positions are adjacent, all other chain incidences are
unchanged; the reversed chain arc becomes forward and the apex bits swap.
Therefore

    F_i(T(x)) is isomorphic to T(x').                         (2)

This is why the main theorem is about **one**, not two, reversed arcs.
No graph editing is hidden in an arithmetic state transition.

The inherited exact response is

    H(T(x'))-H(T(x))=2(1+L_i(x)+R_i(x)),                      (3)

    L_i(x)=sum_{j=0}^{i-2}(1-x_j)2^(i-j-2),
    R_i(x)=sum_{j=i+3}^{n-1}x_j2^(j-i-3).

Empty sums are zero. The immediate outside neighbors x_(i-1),x_(i+2)
do not occur. To prove (3), the swapped pair gains one in C. Each earlier
zero changes its weight by w(d+1)-w(d); each later one changes similarly.
This difference is zero for d=1 and 2^(d-2) for d>=2, proving the formula.
The inverse exchange 01->10 has the negative increment evaluated in the
corresponding 10 configuration. Thus an individual exchange is reversible.

### Direct injection and the unmatched paths

There is an explicit injection from Hamilton paths of T(x) into those of
F_i(T(x)). Keep every path avoiding the reversed arc unchanged. In a path
using the consecutive pair v_i,v_(i+1), swap those two vertices. The apex
cannot be immediately before that old pair, since a does not beat v_i;
it cannot be immediately after it, since v_(i+1) does not beat a. All other
vertices are chain vertices and have identical incidence with the adjacent
pair. The modified path is therefore legal. Images of these modified paths
use the reversed arc, while unchanged images do not, proving injectivity.

A new path using the reversed pair is in the image exactly when the apex
is not adjacent to it. Thus the unmatched paths are exactly the two block
languages

    a,v_(i+1),v_i       or       v_(i+1),v_i,a.               (3a)

Their respective counts are

    [i=0] + [i>0]x_(i-1)
       + sum_{j=i+2}^{n-1} x_j 2^(j-i-2),

    [i=n-2] + [i<n-2](1-x_(i+2))
       + sum_{j=0}^{i-1}(1-x_j)2^(i-j-1).                   (3b)

For the first block, all earlier chain vertices must lie before the apex;
the later ones may choose either side, with the last before-apex vertex
required to beat it. The second block is the reversed argument. The two
immediate-neighbor contributions cancel their complementary boundary terms
when (3b) is added, giving exactly 2+2L_i+2R_i. This independently explains
the positive physical surplus counted by the inherited response formula.

## 3. Complete marked-gap spectrum and optimal order

Fix n>=2 and a legal site 0<=i<=n-2. Put

    a_i=max(i-1,0),   b_i=max(n-i-3,0).

As the outside tape varies, L_i and R_i independently range over
[0,2^a_i-1] and [0,2^b_i-1]. Consequently the exact set of half-increments
at this marked site is the full integer interval

    1 <= N <= 2^a_i+2^b_i-1.                                 (4)

There is even an exact multiplicity. Write u=n-2-a_i-b_i for the number
of response-null outside bits. The number of tapes with a 10 at i and
half-increment N is

    2^u max(0, min(2^a_i-1,N-1)-max(0,N-2^b_i)+1).             (5)

This counts the pairs L+R=N-1 and the free null bits.

For n>=3 the largest value in (4), over i, is 2^(n-3), attained at either
end. Those end intervals already attain every smaller positive N. Hence
with v=n+1 vertices the unrestricted exact positive gap set is

    {2,4,...,2^(v-3)}  for v>=4.                              (6)

At v=3 its positive gap set is {2}. It follows that for N>=2 the smallest
order in this class is exactly 4+ceil(log_2 N).

For strong tapes, the exchange cannot touch the first or last position, so
1<=i<=n-3, n>=4. Now a_i+b_i=n-4, and the maximum of (4) is at most
2^(n-4). The bits enforcing strongness give the more precise fixed-site
interval

    1+[a_i>0]2^(a_i-1)+[b_i>0]2^(b_i-1)
       <= N <= 2^a_i+2^b_i-1.                               (7)

Each is a full interval, but their union need not start at 1: n=5 permits
only N=2, and n=6 permits only N=3,4. We do not assert fixed-order strong
coverage below the maximum. Equation (7) gives the lower bound
v>=5+ceil(log_2 N), with v>=5 also at N=1. The construction next attains it.

## 4. Explicit all-even selector and the strong lift

Let M=N-1 and l=bitlength(M), with bitlength(0)=0. Write the l binary digits
of M most significant first, complement them, and append one 1. Call this
word p. Use the pair

    x=p10,   x'=p01.                                         (8)

The last pair's left register is exactly M, its right register is zero;
(3) already proves the desired increment 2N, at order l+4.

To identify absolute counts, for any binary word representing M,

    f(M)=2 sum_{i<j} b_i b_j w(j-i).                          (9)

The stated recurrence proves (9) by appending one bit; f is even. For the
complemented digit word followed by 1, summing the geometric weights from
each zero gives

    C(p)=M-f(M)/2.

Appending10 adds M, while appending 01 adds 2M+1. Therefore

    H(T(p10))=4M+1-f(M),  H(T(p01))=6M+3-f(M),                (10)

which is (B). This includes N=1 at order4; the smaller pair10/01 gives
counts 1/3 on three vertices.

For N>=2, p begins with 0. Append a common 1 to (8):

    y=p101,  y'=p011.                                        (11)

These words begin 0 and end 1, hence both tournaments are strong. The added
bit is the immediate right neighbor of the exchange, so (3) keeps the
increment unchanged. Directly, for any nonempty word z of length k,

    H(T(z1))-H(T(z))=comp_k(z)+(1-z_(k-1)),                   (12)

where comp_k is the integer represented by the complemented k-bit word.
For both words in (8), (12) equals 8M+2. Adding this to (10) proves (A).
The order is l+5=5+ceil(log_2 N), as required. For N=1 use 0101/0011;
their energies are 4/5, giving 9/11 on five vertices.

Finally apply the single chain-arc flip (2) rather than changing apex arcs.
Strongness and H are preserved under the explicit endpoint relabeling,
so the resulting literal graphs prove the main theorem.

## 5. Parallel Rule184: exact interactions and a Lyapunov law

Let S be any set of legal 10 sites in the original tape. They are necessarily
disjoint. Simultaneously swap them, obtaining x^S. Quadratic expansion of
(1) gives

    H(T(x^S))-H(T(x))
      =2[sum_{i in S}(1+L_i+R_i)
           +sum_{i<j in S}w(j-i-1)].                         (13)

For two sites at distance d>=2, the mixed second difference of C is
w(d-1)+w(d+1)-2w(d)=w(d-1). There are no higher interactions because C is
quadratic. This proves (13), including the smallest distance 2. The strict
rearrangement sign itself was already proved in the October 20 packet.

Taking S to be every legal 10 implements radius-one Rule184 on a finite
tape with closed boundaries x_(-1)=0,x_n=1:

    x_i' = (x_(i-1) AND NOT x_i) OR (x_i AND x_(i+1)).         (14)

The rule conserves the number r of ones. It is not claimed globally
reversible; only the individual marked exchange is reversible. Equation
(13) makes H strictly increase unless x=0^(n-r)1^r. Each exchange decreases
the ordinary number of 10 inversions by one, so the process terminates in
at most r(n-r)<=floor(n^2/4) nontrivial rounds. At its endpoint,

    H=2^n-2^r-2^(n-r)+3       if 0<r<n,
    H=1                       if r=0 or r=n.                 (15)

Indeed summing the geometric weights across the unique zero/one cut gives
C=((2^r-1)(2^(n-r)-1)+1)/2. Thus Rule184 supplies an exact physical
monotone graph update, but not independently additive digit moves.

Hostile control:1010 ->0101 has H=3->9. The two individual increments in
the original tape sum to 4; their interaction adds another 2. Forgetting this
pair term gives a false parallel constructor. At each parallel step, the
pairwise disjoint adjacent chain-arc flips followed by disjoint endpoint
transpositions give the literal physical graph implementation of (14).

## 6. Why this is not yet an absolute-count constructor

The programmable increment is target-complete, but the base counts in
(A)/(B) retain f(N-1). Subtracting that base algebraically does not produce
a tournament with the difference as its number of paths. In particular,
the even difference is not itself a valid Hamilton-path count. A physical
operation turning a paired difference into an absolute odd count, with its
needed endpoint coordinates, is still missing.

Neither H alone nor a bounded local neighborhood closes the update.
The equal-H tapes 00010 and 01010 both have H=9 and 10 at site 3, but their
increments are 8 and 6. The discarded register changes the outcome.

There is also a sharp equal-length context-free rewrite obstruction.
For every word z of length k and integer g>=1,

    H(T(z 0^g 1))-H(T(z))=2^g(comp_k(z)+1).                  (16)

The suffix zeros add no energy, and the last 1 proves the formula by a
geometric sum. If distinct equal-length blocks u,v had a constant additive
H-change under replacement u->v in every outside context, compare the
empty right context with 0^g1. Equation (16) would force comp(u)=comp(v),
so u=v, a contradiction. Thus a nontrivial fixed tape rewrite cannot be
a context-independent additive gate in this encoding. This is not an
obstruction to target-dependent rewrites, restricted languages, other
encodings, or the Chan–Pak conjecture.

## 7. Reproduction and validity gate

Run from this directory:

    python3 verify.py > VERIFY.json
    python3 -O verify.py > VERIFY_OPT.json
    cmp VERIFY.json VERIFY_OPT.json

The checker uses explicit `require` tests, so -O disables no checks. Its
literal Hamilton-path subset DP is separate from both the weighted energy
and the P,a,Q partition count. The graph constructor is explicit; the
single-flip target is verified by its literal adjacency matrix, its exact
one-arc difference, the endpoint transposition, and strong reachability.

The recorded finite universe is:
- all 1023 tapes of lengths 0..9: literal subset DP, partition count, energy;
- all 20481 marked 10 exchanges in lengths 0..12;
- all 321 marked exchanges in lengths 2..7: literal path injection and
  exact unmatched-path language, not only count agreement;
- all 66 fixed-position gap histograms, including exact multiplicities,
  and all 45 strong fixed-position intervals;
- all 8191 parallel updates in lengths 0..12;
- all 4059 switch subsets in lengths 0..9;
- all 2047 complete Rule184 trajectories in lengths 0..10;
- N=1..4096 for both explicit difference families and all literal arc flips;
- N=1..64 independently recounted on literal graphs for each family;
- all 6141 observer suffix tests at lengths 0..10 and gaps 1,2,5;
- a 1001-bit selected increment, verifying its exact value on 1005 vertices.

The large-target test uses N=2^1000+2^601+137, on 1005 vertices in the
unrestricted family; it checks the exact difference, not an empirical growth
fit. Finite tests support the displayed proofs and do not supply unproved
infinite coverage. The script emits several literal graph witnesses, the
same-H hostile, the parallel-interaction hostile, and all scope labels.
