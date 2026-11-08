# A full-rank rooted tournament graft

**Status: PROVED, with freshly reproduced FINITE-EXACT controls.**
This is a reconstruction, not restoration of the unavailable original Session225 files.
No literature-priority, smallest-gadget, scalar-coverage, or completeness claim is made.

## 1. Objects, conventions and the operation

All graphs are finite tournaments. A directed path has a specified vertex order.
Set H(empty)=1. When a path-cover bin is allowed to be empty, there is exactly
one empty path in that bin. Bins are ordered; no factorial division is used.

For a nonempty rooted tournament (A,a), define the column vector

    x(A,a) = (D,E,U,d)^T,

where D=H(A), E counts Hamilton paths ending at a, d=H(A-a), and U counts
ordered pairs (P,R) of vertex-disjoint directed paths covering A, with R
nonempty and starting at a and P possibly empty.

For a second nonempty rooted tournament (B,r), form G_B(A,a) by retaining
the two induced tournaments, directing

    every vertex of B -> a,
    every vertex of A-a -> every vertex of B,

and marking r as the new root. These rules specify every cross-edge.

The operation preserves a four-coordinate linear description, not the scalar
Hamilton-path count alone. In particular, the singleton source and the
source-rooted directed edge both have D=1, but the gadget below produces
D'=5 and D'=37 respectively.

## 2. Gadget statistics and the general closure theorem

For (B,r), let:

- H=H(B), and Q be the number of ordered two-bin spanning directed-path covers,
  with either bin allowed to be empty.
- s and e count Hamilton paths starting and ending at r, respectively.
- alpha counts ordered two-bin spanning covers with the first bin nonempty
  and starting at r; the other bin may be empty.
- delta is the corresponding count with the first bin ending at r.
- T counts ordered three-bin spanning covers with the first bin nonempty and
  starting at r; the other two bins may be empty.
- q=H(B-r), and b be the ordered two-bin cover count of B-r, allowing empty bins.

Swapping the two ordered bins shows that alpha and delta can equivalently
designate the second bin. In particular alpha is the gadget's U coordinate.

For |B|>=2,

    x(G_B(A,a),r) = K_B x(A,a),

with rows

    [ H,     -H,       Q-H,         2H-Q       ]
    [ e,     -e,       delta-e,     e-delta    ]
    [ alpha, s-alpha,  2T-alpha,    3alpha-s-2T]
    [ q,     -q,       b-q,         2q-b       ].

If B is a singleton, use the same first three rows but replace the fourth
row by [1,0,0,0]. The complete singleton matrix is

    J = [1 -1 1 0]
        [1 -1 0 0]
        [1  0 1 0]
        [1  0 0 0].

### Proof of the first row

The only possible transition from B back into A is into a. A Hamilton path
therefore has at most one such transition.

With no B-to-A transition, all A precedes all B. The A path must not end at a,
because a loses to B. This gives (D-E)H paths.

If the transition occurs, deleting the B runs leaves an arbitrary path P in
A-a and a path starting at a, jointly covering A. There are U such splits.
Exactly d have the latter path equal to the singleton a. For these, the
whole B path occurs immediately before a, contributing dH.

In the remaining U-d splits, the a-starting path has at least two vertices.
A nonempty B run occurs before a and an optional B run after that old suffix.
There are Q-H assignments of B to these ordered runs: subtract the H covers
whose first run is empty. This contributes (U-d)(Q-H).

The decomposition is unique in both directions, so

    D' = (D-E)H + dH + (U-d)(Q-H),

which expands to the stated row.

### Proof of the second row

If the Hamilton path ends at r, its final B run is nonempty. With no B-to-A
transition, the contribution is (D-E)e. With a transition, the old a-starting
suffix must be nontrivial; the two B bins have a nonempty first bin and a
nonempty final bin ending at r. Their count is delta-e. Hence

    E' = (D-E)e + (U-d)(delta-e),

giving the second row.

### Proof of the third row

Count covers (P,R) of the graft, with R nonempty and starting at r.
Across both disjoint paths there is still at most one B-to-A transition,
since every such transition must enter the unique vertex a.

1. No transition: all of A lies in P. If its Hamilton path does not end at a,
   an optional B suffix may follow it, giving (D-E)alpha. If it ends at a,
   that suffix is forbidden and R must be a Hamilton path of B starting at r,
   giving Es.
2. The transition lies in R: its initial B run starts at r. A nontrivial old
   suffix from a gives (U-d)T. The three ordered B bins are the initial R run,
   terminal R run, and terminal P run. If the old suffix is just a, a terminal
   R run is forbidden, giving d alpha.
3. The transition lies in P: R lies entirely in B and starts at r.
   A nontrivial old suffix gives (U-d)(T-alpha), because the B run immediately
   before a must be nonempty. Subtracting alpha removes exactly the assignments
   with that run empty. If the old suffix is just a, P has no terminal B run,
   giving d(alpha-s).

These cases are disjoint and exhaustive. Ordinary path bins may be empty,
including when A is a singleton. Summing gives

    U' = (D-E)alpha + Es
         + (U-d)T + d alpha
         + (U-d)(T-alpha) + d(alpha-s),

which expands to the third row.

### Proof of the fourth row and boundaries

Deleting the new root r leaves the same unrooted graft with gadget B-r.
When B-r is nonempty, the first-row proof applies with H,Q replaced by q,b.
When B is a singleton, deleting it leaves A, so d'=D instead.

For singleton B, H=s=e=alpha=delta=T=q=b=1 and Q=2, giving J above.
For singleton A, x=(1,1,1,1), and the formulas give D'=H, E'=0, U'=2alpha.
This is exactly the graph B followed by the old sink a. No empty rooted
source or empty rooted gadget is part of the theorem.

## 3. The explicit four-vertex gadget

Use vertices 0,1,2,3 and root 0. Outgoing-neighbor bit masks are

    B = [8,1,3,6].

Equivalently, its arcs are 0->3, 1->0, 2->0, 2->1, 3->1 and 3->2.
The spanning cycle 0->3->2->1->0 proves that B is strong.

Direct finite path and path-cover enumeration gives

    (H,Q,s,e,alpha,delta,T,q,b) = (5,32,1,1,5,7,17,1,8).

These tiny counts are independently checked by subset DP and permutation/split
enumeration in verify.py. Substitution in the proved general formula gives

    K = [5 -5 27 -22]
        [1 -1  6  -6]
        [5 -4 29 -20]
        [1 -1  7  -6],

and exact elimination gives det(K)=-8. Negative displayed coefficients are
a change-of-coordinate effect: the preceding combinatorial case counts are
nonnegative counts of actual paths.

## 4. Infinite strong family and exact recurrence

Start with the directed edge 0->1 rooted at 0, whose outgoing rows are [2,0]
and whose state is x_0=(1,0,2,1)^T. Repeatedly apply the gadget, each time
marking vertex 0 of the newly added copy. After j applications the graph
has 4j+2 vertices and state x_j=K^j x_0.

If |A|>=2, a reaches every vertex of A, and B is strong, the graft is strong:
every old vertex other than a goes directly into B; a can first reach another
old vertex and then B. Every vertex of B reaches a and hence every old vertex.
Strong connectivity inside B supplies the remaining connections.

The initial source-rooted edge meets this reachability hypothesis. Thus every
iterate with j>=1 is strong, and the hypothesis persists at its new root.

Exact matrix multiplication verifies

    K^4 - 27K^3 - 37K^2 + 31K - 8I = 0.

Consequently h_j=H(A_j) satisfies, for every j>=0,

    h_(j+4) = 27 h_(j+3) + 37 h_(j+2) - 31 h_(j+1) + 8 h_j,

with

    h_0=1, h_1=37, h_2=1113, h_3=31349.

The next values are 886465, 25060261, 708463337 and 20028510133.
The first three nontrivial states are

    x_1=(37,7,43,9),
    x_2=(1113,234,1224,277),
    x_3=(31349,6561,34585,7785).

Literal graph DP checks these at orders 6,10,14. Later displayed values are
exact matrix consequences of the proved closure, not claimed literal
enumerations of their larger graphs.

## 5. Separation from one fixed universal O/Z/J word

Here the operations have these explicit meanings:

- O adjoins a global sink and keeps the old root.
- Z adjoins a vertex beaten by every old nonroot vertex and beating the old
  root; the old root is retained.
- J uses the same graph as Z but marks the newly added vertex.

Their four-coordinate matrices are

    O = [1  0 0  0]       Z = [1 -1 1  0]
        [0  0 0  0]           [0  0 0  1]
        [0  0 2  0]           [0  0 2 -1]
        [0  0 0  1]           [0  0 0  1],

together with J from Section 2.

For O, every Hamilton path ends at the sink, so E'=0; in a rooted two-bin
cover the sink can be appended to either bin, yielding U'=2U. The other
coordinates are D'=D and d'=d.

For Z, the new vertex either occurs before the old root or ends the whole
path, yielding D'=U+D-E. A path ending at the old root must have the new
vertex immediately before it, so E'=d. In a root-starting two-bin cover the
new vertex is appended to either bin, except when that bin is the singleton
old root. Exactly d such forbidden placements exist, so U'=2U-d.
Deleting the old root leaves a global sink on A-a, giving d'=d.
The J formula was proved as the singleton-gadget case.

All three matrices are singular. Every nonempty fixed word in them is
therefore singular. The empty word is the identity.

The physical state vectors supplied by the rooted singleton, source-rooted
edge, rooted directed triangle [2,4,1], and source-rooted transitive triangle
[6,4,0] are respectively

    (1,1,1,1), (1,0,2,1), (3,1,3,1), (1,0,4,1).

The matrix with these vectors as rows has determinant -4, so physical states
span the full four-dimensional rational space. Thus agreement of two linear
operators on every physical rooted tournament forces matrix equality.

Since det(K)=-8, K is neither a nonempty O/Z/J-word matrix nor the identity.
No single fixed O/Z/J word therefore simulates this graft's four-coordinate
action on every rooted tournament.

This does not exclude input-dependent words, individual matching states or
Hamilton-path counts, other graph operations, or broader contexts. It is
not a scalar-spectrum separation.

## 6. Reproduction and scope

Run, in this packet's directory,

    python3 verify.py > replay.json
    cmp replay.json CHECK.json
    python3 -O verify.py > replay.optimized.json
    cmp replay.optimized.json CHECK.optimized.json

Only the Python standard library is required. Explicit exceptions, rather
than assertions, enforce every check. The program generates all labeled
tournaments of orders 1 through 4 with root fixed at label 0: 75 sources
and 75 gadgets, with all 5625 pairs checked by direct subset DP. It does
not use a missing original gadget bank or an unrecorded filter.

Additional checks comprise 75 direct permutation gadget-stat comparisons,
76 permutation profile comparisons (including the first six-vertex graft),
225 literal O/Z/J controls, 1161 applicable strongness controls, the four
physical basis states, exact determinants, the matrix recurrence, and the
three larger literal iterates.

The proof supplies the all-order theorem; finite checks are controls, not
an inference of universal validity. This packet makes no claim to an
improved unrestricted finite coverage range, an eventual cutoff, or the
Chan-Pak all-odd completeness conjecture, which remains OPEN.
