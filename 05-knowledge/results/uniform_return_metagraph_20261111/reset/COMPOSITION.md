# The reset gadget belongs to a completely explicit two-parameter family

**PROVED in this note.** The operation, Hamilton count, retained-root
response, residual count, recurrence, and rooted family enumeration are
known for every parameter. This is one family inside the metagraph, not
a classification of every tournament or a recursive Chan–Pak proof.

## 1. Actual composition and its five-dimensional carrier

Start with BA_k from the old ready pentagon, k>=1. Mark its01 tail vertex
B, keeping the original root r and ports a,b. The cell operation G appends
two successive global sinks x,y, then z beating exactly B and x. Repeat
G t times, always using the SAME old B. Call the result V_(k,t).
The new uniform reset lands at V_(k,1).

Each new cell is x→y→z→x. Every earlier vertex except B dominates every
vertex of the cell; B beats x,y and loses to z. Removing B therefore
leaves the ordinal join

    X_k ▷ C3 ▷ ... ▷ C3,                              (1)

where X_k=BA_k−B and the t triangles have identical apex cut{x,y}.
This is a literal graph decomposition, not a recurrence inferred from
initial terms. G preserves r and its two in-ports and keeps the graph
strong: the old graph reaches x, and the new cycle reaches z→B.

Use the inherited [five-dimensional interface algebra](../../session252_fivefold_bridges_20261111/algebra/PROOF.md),
T_C(z)=[[d,O,J],[0,d+K,I],[0,0,d]], with ordinal concatenation as matrix
multiplication. For the coloured triangle the two exact matrices are

    M+ = [[3, 6,2], [0,12, 3], [0,0,3]],
    M− = [[3,−2,0], [0, 0,−1], [0,0,3]].                       (2)

They retain the signed central channel as well as the positive one.
For t>=1 their relevant entries are

    I+_t=(12^t−3^t)/3,       J+_t=2(12^t−3^t)/9,
    I−_t=−3^(t−1),          J−_t=2(t−1)3^(t−2),                (3)

where J−_1=0. These follow by multiplying upper-triangular matrices:
the two-step off-diagonal sum is
(12^t−3^t−9t3^(t−1))/81. Substitution of the entries in(2) cancels
the positive t3^(t−1) term, while the negative middle diagonal vanishes.
The surviving negative central channel is exactly the polynomial factor
in the final3^t term.

For D=2^k, the initial coloured block X_k has

    d=6D+1,       O+=15D+2,       O−=−3D,
    H(BA_k)=14D+3.                                             (4)

The count d is the old A_k family. The two remaining identities have an
all-k positive-run certificate: the first-two-bin-starting-in-cut count
has run coefficients(32,15,0,0,...), and the endpoint count has
(6,3,0,0,...). The ordinary initial counts have coefficients(13,6)
and(31,14). Replacing the free prefix by TT_k evaluates each pair as
A1+A2(2^k−2). A direct enumeration of all marked placeholder words proves
the displayed coefficients; the second endpoint identity uses the inherited
alternating endpoint theorem O(−1)=−E_out. Thus(4) is a fixed finite
combinatorial calculation that proves every k, not two fitted samples.

Multiplication T_X M^t and H=d+J(1)+J(−1) give

    H(V_(k,t))=(14D+3)3^t+d(J+_t+J−_t)+O+ I+_t+O− I−_t.          (5)

This derives the complete response below.

## 2. Closed formulas for all structural parameters

For t>=1 define

    A_t=[19·12^t+(4t+22)3^t]/3,
    B_t=[8·12^t+(2t+17)3^t]/9,
    q_t=[24·12^t+(6t+21)3^t]/9,
    C_t=[16·12^t+4(t+1)3^t]/9,
    R_t=[5·12^t+(2t+5)3^t]/9.

Then for every k>=1,

    order(V_(k,t)) = k+5+3t,
    H(V_(k,t)) = A_t·2^k+B_t,
    s_r=q_r=q_t,          F_r−H=C_t,
    H(V_(k,t)−{r,a,b})=R_t.                                  (6)

At t0 use(A,B,q,C,R)=(14,3,5,2,1). The boundary is stated separately:
the zero eigenvalue in M− has a genuine t0 transient.

Formula H=A_tD+B_t follows from(3)–(5). To obtain q_t, delete the root
and the free prefix; deleting B from what remains leaves t+1 identically
coloured C3 blocks. Applying(2) to their product gives

    H(apex over r such blocks)=[2·12^r+(2r+5)3^r]/9.

Set r=t+1. Prefix dominance makes s_r=q_r=q_t. The general two-bin
prefix formula H=A_tD+q_t−C_t gives C_t=q_t−B_t. For R_t, remove both
ports and r: deleting B leaves one singleton followed by t triangles;
B beats that singleton. Its initial matrix has d1,O+1,O−−1,J±0,
which gives the stated R_t through the same product.

The first values illustrate exactly where the five-core bank enters and
where it ceases to be the entire residual:

|t|A_t|B_t|q_t|C_t|R_t|
|---:|---:|---:|---:|---:|---:|
|0|14|3|5|2|1|
|1|102|17|41|24|9|
|2|1002|149|417|268|89|
|3|11250|1605|4725|3120|993|
|4|132354|18657|55701|37044|11637|

At t1 the residual is exactly TT_k→S5a. At later t it is still a
specified apex over an ordinal sequence of triangles; its growth law
and all response coordinates in(6) remain known.

## 3. Composition laws and exact rooted enumeration

The preparation operation P prepends a new11 prefix vertex. It commutes
with G as a physical graph operation, up to the evident newborn labels:
the new prefix vertex beats every gadget vertex, and G always anchors at
the same B. Hence(P,G) supplies an actual N² family of commuting growth
operations, retaining r,a,b,B.

For fixed t, preparation gives H_next=2H−B_t and preserves q_t,C_t.
For fixed k, all t>=1 satisfy

    H_(t+3)=18H_(t+2)−81H_(t+1)+108H_t.                       (7)

Its characteristic polynomial is(X−12)(X−3)^2. The12^t coefficient and
the t3^t coefficient in(5) are positive, so the minimal scalar recurrence
has order3 for each k. This is a specialization of the older fixed-block
recurrence mechanisms, including [the sequence grammar](../../sequence_gf_20261110/grammar/PROOF.md);
no claim is made that arbitrary substitutions are C-finite.

The parameter pairs give distinct ROOTED classes: C_t is strictly
increasing in t, and after t is fixed the order identifies k. Therefore
their exact rooted order-generating function is

    x^6 / ((1−x)(1−x^3)),

or x^9/((1−x)(1−x^3)) if only the new t>=1 targets are included.
This counts this specified family, not all rooted tournaments with the
same response. No unrooted isomorphism or automorphism claim is needed.

## 4. Every possible ready root is classified, in both orientations

For every k>=1,t>=1, V_(k,t) has NO1-ready or2-ready root, including
arbitrary re-marking and taking the converse. This is an all-parameter
statement, not just the retained root's C_t profile.

The inherited positive defect bound gives gap>=d(d−1) for a root of
indegree d. Only vertices of indegree1 or2 need consideration.
In the forward graph these are exactly the following candidates:

* The retained r has gap C_t>=24.
* The first prefix vertex has sole in-neighbour r. Its gap equals the
  number of Hamilton paths starting at r after deleting that prefix vertex.
  If k>=2 this is q_t>=41. If k1, the t1 deletion has21 such paths; later
  cells cannot decrease that number, since r beats every added vertex
  and the ordinary insertion injection can keep r first. Thus its gap>=21.
* The second prefix vertex, when present, has two in-neighbours(r,u1).
  Deleting all three leaves an ordinal prefix of the fixed kernel W_t,
  of H=q_t. The joint-defect bound gives gap>=2q_t>=82.
* Port a has indegree2 only when k1. Deleting it and its two in-neighbours
  leaves the directed triangle B→b→r→B as an induced subgraph, so its
  joint defect is at least6. At larger k its indegree is already>=3.

Every other vertex has indegree>=3. These cases follow directly from the
cell construction: old vertices dominate every new cell except B, which
acquires an additional in-neighbour z; each new cell vertex itself has
at least six in-neighbours.

In the converse, only the LAST cell's x,y,z have indegree<=2, equivalently
outdegree<=2 in the original. Earlier cell vertices dominate all three
vertices of a later cell; all old BA vertices have outdegree>=4 once t>=1.
Write S=V_(k,t−1).

* Original y has only out-neighbour z. Its converse gap is the number of
  paths ending at z in V−y. But x is a global sink there, so this number
  is zero. This is a real zero-gap root, not a ready root.
* Original x has only out-neighbour y. In V−x, each Hamilton path P of
  S−B gives the distinct path P,z,B,y ending at y. Hence its converse gap
  is at least H(S−B)=(6·2^k+1)3^(t−1)>=13.
* Original z has the two out-neighbours B,x. Its converse joint defect is
  at least2H(V−{z,B,x})=2H(S−B)>=26.

This exhausts every degree-eligible root in both orientations. It explains
why the next recursive return requires a physical source change rather
than selecting a different root in the new output.

## 5. Exact exclusion from the near-square prime sequence

The entire family V_(k,t), k>=1,t>=0, is disjoint from A002496, the
primes of the form n²+1. This is an exact infinite-family filter, not a
claim about the unknown infinitude of that sequence.

* At t0, H−1=14·2^k+2 is2 modulo4, hence cannot be a square.
* At t1, H=17(6·2^k+1) is composite.
* At t>=3, (5) factors as3^(t−2) times a positive integer, and H>3.
  Therefore H is composite.
* At t2, H=1002·2^k+149. If H−1=n², then n=2m and
  m²=501·2^(k−1)+37. For k1,2 the right side is respectively2,7
  modulo8; for k>=4 it is5 modulo8. At the only remaining k3 it is
  2041, strictly between45² and46². There are no solutions.

Thus the repeated gadget supplies completely known tournament families
and exact enumeration filters, while its only potentially prime layer
beyond BA is t2 and even that layer contains no square-plus-one value.
The higher layers' forced factor3 also removes them from the prime/7p
frontier used by the repository's Chan–Pak reduction.

## 6. The exact recursive-closure boundary

G has a known law, but it is not the next ternary branch. For example

    V_(k,1): H=102D+17,
    V_(k,2): H=1002D+149,
    3H(V_(k,1))+2=306D+53.

The response changes from C24 to C268 rather than returning to a1/2-ready
interface. The all-future O/Z/J comparison recovered in the incoming
[H13 dynamics note](../../h13_paired_descent_20261111/dynamics/PROOF.md)
can compare rooted costs after a chosen re-marking, but equality of those
cost bounds does not supply the missing exact target. A different surgery,
chart switch, or controlled variable-length cell word is still required
for an iteratable ternary proof.

Reproduce with `python -B composition.py`. It verifies the marked initial
run coefficients by an independent word enumeration, matrix products
through t20, literal kernel path counts through t4, actual graphs and
root profiles over k1..4,t0..2, and physical commutation of P and G.
The proof of every parameter is the decomposition and matrix identities,
not those finite controls.
