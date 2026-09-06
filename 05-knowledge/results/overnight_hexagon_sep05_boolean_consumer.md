# Boolean cycle adjacency: two losses, a complete Fourier obstruction, and the legitimate comparison

Status: **PROVED ELEMENTARY STATEMENTS, WITH INDEPENDENT PROOF AUDIT**;
small full-matrix controls are **FINITE-EXACT**. The exact all-order Boolean
spectrum remains **OPEN**. No new theorem identifier is requested.

## 1. The actual consumer and three different spectral questions

The native object is the simple orbit graph E_n^(D) of Eulerian graphs
under toggling a simple cycle of length 3,...,D+1. The original current
consumer is [THM-4073, exact cycle distance](../../01-canon/theorems/THM-4073-even-graph-diameter-layer-exact-cycle-distance.md),
§5: its full unweighted spectrum, diameter, clique number and chromatic
number remain undetermined. It does not state a conjectural exact Boolean
gap formula. We do not manufacture one from the weighted theorem.

Let M be the cumulative integer multiplicity matrix and let B be its
off-diagonal Boolean support. Write d_i=sum_j B_ij and N=sum_k N_(n,k).
The following are different objects:

```text
adjacency:                 B;
combinatorial Laplacian:   L_B=diag(d_i)-B;
simple random walk:        P_B=diag(d_i)^(-1)B.
```

B is symmetric but generally irregular. The continuous-time generator
-L_B has uniform stationary measure; P_B has degree-proportional
stationary measure. In contrast M/N has orbit-size-proportional stationary
measure. The bounds below concern L_B, not an unspecified adjacency gap
or the absolute mixing rate of the bipartite triangle walk.

The closest mechanism is
[THM-4078, signed-dual spectrum and Boolean noncommutation](../../01-canon/theorems/THM-4078-even-graph-triangle-quotient-spectrum-and-boolean-noncommutation.md).
Its Fourier orbit sums diagonalize M, not B. The canonical hostile is its
n=4 noncommuting Boolean triangle/four-cycle pair. The corrected near miss
is [MISTAKE-495](../../01-canon/MISTAKES.md): a canonical vertex quotient did
not make a fixed spanning-tree edge image canonical. The least-used
sidecar here is the native symmetric edge conductance, including the
orbit volume rather than only the positive-support graph.

The live board is: orbit volume; directed generator multiplicity;
cross-length support overlap; degree-adjusted loop deletion; native
Dirichlet energy; actual sampling clock. No tournament or LRC transfer is
asserted. The old parity-sensitive graph/two-graph identification is not
substituted for the current Fourier duality.

## 2. Two different multiplicities are lost

Inside one cycle length, M_k(i,j) counts the labelled generators carrying
a fixed representative of i to orbit j. Booleanization forgets that
number, retaining only whether it is positive. Even after symmetrization
using orbit sizes, the edge weights need not be equal. At n=4 the triangle
matrix in the order empty,C3,C4 is

```text
M3=[[0,4,0],[1,0,3],[0,4,0]].
```

The two off-diagonal products M_ij M_ji are 4 and 12. These products are
invariant under diagonal similarity, so neither an orbit-size diagonal
change nor an additional common scalar can turn M3 into the unweighted
path adjacency. Correct stationary weighting does not remove edge weights.

Across lengths, the cumulative simple graph is a **union**, not the sum
of the simple layer graphs:

```text
B=Bool_offdiag(sum_k M_k),
B need not equal sum_k Bool_offdiag(M_k).             (1)
```

The first possible overlap is at n=5. Both a triangle and a five-cycle
toggle can join [C3] to [C4]. With F=(12,23,31), the triangle (12,24,41)
gives a four-cycle after XOR; the five-cycle (12,23,34,45,51) also does,
now giving (31,34,45,51). Isolated vertices extend this to every n>=5.
It is minimal: at n<=4 the only two possible lengths, 3 and 4, have
opposite parity, while a fixed pair of endpoint classes determines the
parity of the toggled cycle.

There is an important distinction concerning loops. Deleting them changes
the adjacency algebra and can change the row-normalized walk. But if the
degree is adjusted consistently, it leaves the weighted combinatorial
Laplacian **exactly unchanged**:

```text
diag(N-M_ii)-(M-diag(M_ii))=NI-M.                    (2)
```

Thus loop deletion alone is not an explanation for loss of a Laplacian
gap. Flattening the off-diagonal multiplicities and changing the vertex
measure are the relevant changes for that question.

## 3. The first signed Fourier mode fails at every possible order

For n>=4, the single-negative-edge Fourier orbit sum on Eulerian graphs is

```text
psi(F)=sum_(e in E(K_n)) (-1)^[e in F]
      =binom(n,2)-2|E(F)|.
```

No affine shift of psi is an eigenvector of the Boolean **triangle**
Laplacian for any n>=4.

For n>=6, the triangle class has exactly four neighbors: empty (the same
triangle is toggled), C4 (one common edge), two triangles sharing one
vertex, and two disjoint triangles. These exhaust intersections of two
triangles. Consequently

```text
L_B psi(empty)=6,
L_B psi(C3)=-6+2+6+6=8.
```

If L_B psi=lambda(psi-a), subtracting these equations gives
lambda=(6-8)/6=-1/3, impossible for a positive-semidefinite Laplacian.

At n=5 the disjoint-triangle neighbor is absent, so the triangle value is
2 and the same two equations would give lambda=2/3. But K5 is Eulerian
and has exactly one triangle-toggle orbit neighbor, K5 minus a triangle;
therefore L_B psi(K5)=-6. Comparing empty and K5 instead gives
lambda=12/20=3/5, a contradiction. This endpoint argument was supplied
by the independent `lrc_components` referee.

At n=4 the Boolean triangle graph is the path empty--C3--C4, and

```text
psi=(6,0,-2), L_B psi=(6,-4,-2).
```

The three columns 1,psi,L_B psi have rank three, which excludes an affine
eigen-equation. At n=3 there are only two states and this obstruction
does not apply. The all-order failure is an analytic statement, not an
extrapolation from the small matrix census.

## 4. The valid conductance comparison

Let s_i be the number of labelled Eulerian graphs in orbit i. Detailed
balance from THM-4073 gives

```text
c_ij=s_i M_ij=s_j M_ji>0 on each Boolean edge.
```

Let s_min,s_max be the extreme orbit sizes and c_min,c_max the extreme
positive c_ij. For real functions f on the common orbit set, define

```text
Q_B(f)=sum_(i<j, B_ij=1)(f_i-f_j)^2,
Q_w(f)=sum_(i<j, B_ij=1)c_ij(f_i-f_j)^2,
V_1(f)=min_a sum_i(f_i-a)^2,
V_s(f)=min_a sum_i s_i(f_i-a)^2.
```

The convention counts each unordered edge once in both forms. Directly,

```text
c_min Q_B<=Q_w<=c_max Q_B,
s_min V_1<=V_s<=s_max V_1.                           (3)
```

The weighted gap lambda_w is the infimum of Q_w/V_s over nonconstant
functions, and the Boolean gap lambda_B is the infimum of Q_B/V_1.
Applying (3) to all functions for the lower bound, and to a weighted gap
eigenfunction for the upper bound, proves

```text
(s_min/c_max)lambda_w <=lambda_B
 <=(s_max/c_min)lambda_w.                            (4)
```

This elementary variational proof needs no imported comparison theorem.
It is not an equality transfer and may be very nonsharp.

The cumulative weighted gap is now proved to be lambda_w=2E by
[THM-4427, all cumulative gaps](../../01-canon/theorems/THM-4427-all-cumulative-signed-cycle-gaps-by-transposition-rigidity.md)
and its inherited smaller-layer theorems, with E=sum_k (n-2)!/(n-k)!.
For D=2 use THM-4078's triangle theorem. The empty orbit has size one, so
s_min=1, while every s_i<=n! and every M_ij<=N. Hence the explicit
all-order survivor, for n>=3 and 2<=D<=n-1, is

```text
gap(L_B)>=2E/(n! N)
         =4 kbar/[n(n-1)n!],
kbar=(sum_(k=3)^(D+1) k N_(n,k))/N.                 (5)
```

In particular the triangle Boolean gap is at least
`12/[n(n-1)n!]`. This is a sufficient quantitative bound, not a claimed
sharp scale. Conversely the indicator of the empty orbit, after centering
in uniform measure, gives the native upper bound

```text
gap(L_B)<=(D-1)q/(q-1), q=number of Eulerian orbits. (6)
```

It uses THM-4073's exact empty degree D-1. Even this one-vertex test shows
why the weighted gap's factorial growth cannot be assigned to the Boolean
Laplacian.

## 5. An exact sampler that retains the missing multiplicity

Propose a uniformly random labelled cycle among the N generators, giving
orbit proposal probability M_ij/N. If the target is a different orbit j,
accept with probability 1/M_ij; otherwise stay. The off-diagonal accepted
probability is exactly 1/N on every Boolean edge. The resulting kernel is

```text
P_thin=I-L_B/N.                                      (7)
```

It is stochastic since d_i<=N, symmetric, and reversible in the uniform
orbit measure. Its constant holding clock is not the distinct-neighbor
walk P_B; conditioning on a move would give P_B and change the stationary
measure to degree weighting. Formula (7) is a native transition compiler,
not a way of recovering Boolean eigenvalues from the weighted eigenvalues.
The directed target multiplicity is the required extra input.

## 6. Exact controls and remaining work

```bash
python3 -B 04-computation/overnight_hexagon_sep05_boolean_consumer.py
python3 -B -O 04-computation/overnight_hexagon_sep05_boolean_consumer.py
```

The [script](../../04-computation/overnight_hexagon_sep05_boolean_consumer.py)
constructs every Eulerian state, orbit and cycle prefix for n=3,...,6,
using the inherited direct primal constructor. Independent literal edge
subset recognizers check the complete Eulerian universe and every
connected 2-regular cycle set through n=5. It retains orbit sizes, checks
detailed balance and weighted commutation, verifies the known weighted
gaps, tests the affine eigenvector obstruction and cross-length collisions,
and checks the exact loop identity, thinning probabilities, energy and
centered-variance comparisons. Characteristic-polynomial digests refer to
L_B explicitly, not to adjacency or a normalized walk.

The n=5 triangle gap is `2-sqrt(2)`, whereas its all-cycle Boolean gap is
`(9-sqrt(17))/2`. The first cross-length support overlap consists of six
edges at n=5,D=4; there are 23 at n=6,D=4 and 44 at n=6,D=5.
These are finite diagnostics, not formulas at larger orders. The
[stored output](overnight_hexagon_sep05_boolean_consumer.out) records the
full ten-row universe and 965 explicit failure gates. Normal and optimized
runs agree; the semantic digest is
`3ec28083073c6221af814f2c36fc29ea47a951dfe1c3a42b10ef1ccc4388ca25`.

The independent `lrc_components` referee audited the all-order Fourier
obstruction, comparison normalization and the degree-adjusted loop identity.
Root independently reviewed the complete written proof, both comparison
directions, the explicit lower/upper bounds and thinning: **PASS**.
Normal and optimized outputs byte-match. Raw-LF SHA256:

```text
source efe4c34392382ef85e045b5db9d0eadb056b60a42336136fa12ae83d7abefa0e
output 4c42829e8237c342522312fcfe4903df0929c6ab0fa8b81dc90e3ba13eba8c0c
```

The next useful question is the
Boolean graph's own conductance profile or a sharper native comparison,
not another signed-cycle weight census. Source-to-target map: retain the
common orbit vertices and Boolean support, preserve reachability and path
lifting, lose directed and cross-length multiplicities, and recover only
the statements licensed by the retained conductance and measure sidecars.
