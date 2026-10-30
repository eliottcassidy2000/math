# Two native complements can hide an almost full-rank grid

2026-10-30. **PROVED elementary, self-audited; FINITE-EXACT controls.**
No independent-agent audit or new A000568 term is claimed. The statements
concern labelled tournaments and actual legal Hamilton paths. They do not
assert that either move is black rather than blue.

## Inherited mechanism and the new question

[THM-4815, dipole principle](../../../01-canon/theorems/THM-4815-dipole-principle-native-move-distance-bounds-endpoint-relation-cut-or-colour-exactness-and-the-p7-parity-correction.md),
A1/A3, already proves the native law K'=-K+2A_P, the two-move XOR law,
converse covariance, and score dependence only on the endpoints.
[THM-4813, spectral exception budgets](../../../01-canon/theorems/THM-4813-collatz-numerator-coordinate-forcing-cocycle-gauss-sum-intertwiner-and-sharp-spectral-exception-budgets.md)
already supplies the operator-norm tariff below4 per move, crediting Session75.
None of these is new here. The question is whether a bounded number of
moves leaves a correction of bounded rank or bounded treewidth, especially
when the named endpoint ports and the scores are unchanged.

The answer is no, including along a path of prime tournaments. The positive
replacement is an exact labelled edge-incidence sidecar; its size is bounded
by the number of moves times N, rather than by a fixed number of coordinates.

## G1. An explicit legal two-move construction

Let m>=3 be odd, N=m², and label the vertices (r,c), 0<=r,c<m.
Let P traverse rows successively, left-to-right on even rows and right-to-left
on odd rows. Let Q traverse columns successively, top-to-bottom on even
columns and bottom-to-top on odd columns. Both go from a=(0,0) to
b=(m-1,m-1). Here top-to-bottom means increasing row index.

Their common edges are precisely the m-1 row connectors of P and the m-1
column connectors of Q. All lie on the outer boundary, and their directions
agree: both boundary row indices and both boundary column indices are even.

Prescribe P forward in T0 and every Q-only edge backward. These prescriptions
are consistent. Complete all other edges arbitrarily. Then P is Hamiltonian
in T0, Q is Hamiltonian in T1=tau_P(T0), and T2=tau_Q(T1) exists legally.
Writing H=E(P) symmetric-difference E(Q), one has

    K2-K0=2(A_Q-A_P),
    |H|=2(m-1)²,       maximum_degree(H)<=4,
    K2*1=K0*1,         ||K2-K0||op<8.                 (G1)

H is Eulerian: the odd-degree sets of the two paths are the same {a,b}.
Every edge of the interior (m-2) by (m-2) square grid belongs to H. Thus its
treewidth is unbounded. For a direct lower bound, the row-column crosses in
an s by s grid form pairwise intersecting connected sets; fewer than s
vertices miss both a row and a column and hence miss their cross. In any tree
decomposition the corresponding bag subtrees have a common point (the Helly
property for subtrees), so some bag has at least s vertices. Taking s=m-2
gives treewidth at least m-3. No graph-minor preservation claim is needed.

## G2. Rank and norm separate sharply

Restrict (K2-K0)/2 to the interior. Conjugating by the diagonal signs
(-1)^(rc) produces

    J_s tensor I_s - I_s tensor J_s,       s=m-2,    (G2)

where J_s is the skew matrix of the directed path, with superdiagonal1 and
subdiagonal-1. Indeed the vertical Q edges originally have sign(-1)^c,
and the horizontal -P edges have sign-(-1)^r. The gauge removes both signs.

The s distinct complex eigenvalues of J_s are 2i cos(j*pi/(s+1)), 1<=j<=s.
Consequently exactly s of the s² differences in (G2) vanish. Thus

    rank(K2-K0)>=(m-2)(m-3)=N-5sqrt(N)+6.           (G3)

This rank computation has a second algebraic proof. J_s is cyclic: the
vectors e1,J_s e1,...,J_s^(s-1)e1 are triangular with nonzero leading entries.
Its commuting algebra has dimension s over any field; an intertwiner for
J_s and its transpose has the same dimension because a matrix is similar to
its transpose (here the alternating diagonal gauge already supplies the
similarity to -J_s). Vectorizing (G2) gives such an intertwiner equation.
Hence its rank is s²-s over any field as well. The checker uses fields3 and5.

The bound<8 in (G1) and the lower bound approaching N in (G3) hold
simultaneously. Small operator norm, unchanged scores, and fixed endpoint
ports do not imply a bounded-rank correction. Two moves already rule out
an all-order carrier whose exact edge correction has bounded treewidth.
This is a scope obstruction to these particular carriers, not an information
lower bound against every imaginable symbolic representation: P and Q
themselves have short descriptions.

## G3. All three tournaments can be prime and have small norm

For every odd m>=9 the unspecified edges can be chosen so that T0,T1,T2 are
all prime and

    max_(i=0,1,2) ||Ki||op<=10sqrt(N)+4.            (G4)

This is an existence statement by elementary random completion, not a
constructive rounding algorithm and not the ramified-bank theorem.

Proof of primality. More generally prescribe any edges with incidence degree
at most4, and orient all other edges independently and uniformly. A fixed
proper subset S of size s>=2 is a module only when, for each outside vertex,
all its cross edges have one direction. There are at most 2^(N-s) such cross
patterns. At least s(N-s)-4min(s,N-s) cross bits remain free, so

    Pr(S is a module)<=2^(-(s-1)(N-s)+4min(s,N-s)).

For s<=N/2, use binom(N,s)<=N^s to dominate the sum by
sum_(s>=2) q1^s with q1=N*2^(4-N/4). For s>N/2 put t=N-s and use
sum_(t>=1) q2^t with q2=N*2^(5-N/2). For N>=64, N<=2^(N/8), hence
q1<=1/16 and q2<=2^-19. The module failure probability is at most

    (1/16)²/(1-1/16) + 2^-19/(1-2^-19).

The prescribed support E(P) union E(Q) has degree at most4. Each T_i has
this same prescribed support and an independent uniform completion on the
other coordinates: the fixed flips preserve their joint uniform law. The
three events need not be independent; their union bound is less than1/64.

Proof of norm. Write K0=C+R, with C containing only prescribed entries.
Then ||C||op<=4. For real unit vectors u,v,

    u^T R v = sum_(free i<j) epsilon_ij (u_i v_j-u_j v_i),
    sum_(i<j)(u_i v_j-u_j v_i)²=1-(u.v)²<=1.

The elementary Rademacher exponential bound gives
Pr(|u^T Rv|>=t)<=2 exp(-t²/2). A1/4-net of the real unit sphere has at
most9^N elements (the disjoint-ball volume argument). Two such nets and
||R||op<=2 max_net |u^T Rv| give

    Pr(||R||op>10sqrt(N)) <= 2*9^(2N)*exp(-25N/2)<2exp(-8N).

This plus the module failure budget is less than1. In every state the free
part is R or -R, while the prescribed part still has degree at most4.
Thus the SAME norm event proves all three bounds, without paying a new
tariff per move. This proves (G4) jointly with primality, fixed ports,
unchanged scores and (G3).

### G3a. A legal path-program compiler with prime completions

The argument extends beyond this grid. Supply r>=1 labelled permutations
P1,...,Pr of N vertices, as proposed native paths. Start with epsilon=0 and
H empty. At step t, each consecutive directed pair u->v of Pt forces the
initial base entry to be

    B_uv=(-1)^(epsilon+[uv in H]).

Then update epsilon and H by the parity rule. The program is realizable
if and only if these forced entries agree whenever they name the same
unordered pair (with orientation reversed when the order is reversed).
Necessity follows from the decoder. For sufficiency complete unforced
entries arbitrarily; induction authenticates every supplied path.
This test retains the named vertices and ports and costs O(rN) dictionary
operations. A conflict has a two-prescription witness on one actual edge.

**PROVED.** If the test passes and N>=max(64,32r), there is a completion
for which all r+1 states are prime and all satisfy

    ||Kt||op<=10sqrt(N)+2r.                         (G5)

The prescribed union has degree Delta<=2r<=N/16. In the module estimate
above replace4 by Delta. With N<=2^(N/8), the two geometric ratios are
at most2^(-N/16) and2^(1-5N/16). The union over r+1<=N/32+1 states is
less than1/64: at N=64 it is bounded by the displayed three-state
budget, and (N/32+1)2^(-N/8), as well as (N/32+1)2^(1-5N/16), decreases
thereafter; their geometric denominators only improve. The free matrix
is (-1)^t R at time t, so one norm event controls the whole program.

This yields prime realizations of any CONSISTENT bounded-length port-labelled
path program at sufficiently large N. It is not a completion of an already
fully specified source tournament, and not a claim that every abstract path
word is consistent. P followed immediately by its reverse is a minimal
conflict: the first move keeps every P edge, so the reverse is illegal.
Nor is it a new ordinary isomorphism-count recurrence. It is a reusable
native-composition existence theorem with an explicit legality test.

## G4. What carries over from the new ramified theorem

The incoming [Session110 logarithm-free proof](../session110_overnight_20261030/LOGFREE_PRIME_TOURNAMENTS.md)
establishes prime balanced lifts with norm<=200sqrt(pN), uniformly in p,
at N=p^(2k), k>=2, under its simple-ramification assumptions. Its residue
charges force norm>=sqrt(pN) for every state of that bank. Its actual
new-row estimate is below191sqrt(pN);200 closes the induction with slack.

Starting from any such labelled base B, overwrite only the prescribed grid
edges to construct T0. The overwrite changes at most four entries in every
row, each by at most2, so ||K0-B||op<=8. Therefore the three norms are at
most200sqrt(pN)+8,+12,+16, respectively. The overwrite may destroy both
primality and the ramified multiplier action. We do not identify this
construction with (G4), where primality is obtained differently. In
particular10sqrt(N) outside the bank is not an improved constant for a bank
whose lower bound is sqrt(pN).

## G5. An exact closed extension, and the remaining quotient

Fix the labelled base B. Store (epsilon,H), where epsilon is move parity
and H is the set of unordered pairs traversed an odd number of times.
Decode by

    K_uv=(-1)^epsilon B_uv (-1)^[{u,v} in H].

On an actual Hamilton path Q update epsilon by1 and H by symmetric difference
with E(Q). This is an exact complement-closed marked carrier. After r moves,
|H|<=r(N-1) and its incidence degree is at most2r. Re-marking by a permutation
acts simultaneously on B, H, the current paths and every named port.
Converse is compatible with reversing each legal path. It is a repackaging
of the inherited XOR law, with a clear accounting of the retained data.

For A000568, taking orbits under Aut(B) only counts objects with the base
marking retained. Different marked states can decode to isomorphic ordinary
tournaments via permutations not fixing B. A target-isomorphism transporter
or canonical target decoder is still required. The existing
[future-response recurrence through order200](../tournament_future_boundary_20261020/README.md)
is already stronger than another small ordinary census. No new unrestricted
recurrence or new A000568 coefficient follows from this packet.

The useful next problem is a joint **edge-incidence plus transporter**
compiler, allowing bounded-degree corrections rather than demanding bounded
rank. The grid shows what that compiler must accept. The incoming
[module enumeration](../session110_overnight_20261030/MODULE_ENUMERATION.md)
likewise retains a single global gauge and a joint tuple of ports; separate
component orbit counts cannot replace either joint object.
