# A one-arc factor-one shield with a programmable mirrored exterior

Use a transitive spine and two apices A,B. At each spine site i,
write the letter 2a_i+b_i, with a_i=1 if A beats that site and
b_i=1 if B beats it. Let Δ(w)=H(w,A→B)−H(w,B→A).

## Contract the arc instead of recomputing all Hamiltonian paths

Every path not using the mutual apex arc survives unchanged.
For paths using A→B consecutively, contract these two vertices
together and partition the spine into the before and after runs.
The former is increasing and must end at a vertex beating A,
and the latter is increasing and must begin at a vertex beaten by B.
A run may be empty with no corresponding endpoint condition.

If the before run ends at i and the after run begins at j, then
either j=i+1 (one partition), or j<i (the i−j−1 interior
sites can be assigned to either run). The single-empty-run
terms are 1−a_last and b_first. Thus, writing

    F(j,i) = (1−a_i)b_j−(1−b_i)a_j,

the exact signed gap is

    Δ(w)=b_first−a_first+b_last−a_last
        + Σ_{j<i} 2^(i−j−1) F(j,i)
        + Σ_{i=0}^{n−2} F(i+1,i).                         (1)

This uses no cycle-packing identity.

For a word v, let L_a(v), L_b(v) be the binary addresses of its
rows read left to right, and let R_a(v), R_b(v) be the addresses
of the complement rows read right to left. The weighted cross
term between end-to-end consecutive words u and v is therefore

    L_b(u) R_a(v)−L_a(u) R_b(v).

This is a determinant of marked addresses, not a comparison
between ordinary counts of distinct tournaments.

## A factor-one seam of arbitrary height

In a word u c^t v with both u,v nonempty, let (a,b) be the bits
of c. For t≥1, (1) gives Δ=K2^t+C, where C is independent of
t and

    K=(L_b(u)+b)(R_a(v)+1−a)−(L_a(u)+a)(R_b(v)+1−b).     (2)

Indeed the u→v term crosses a guard of t sites and picks up 2^t.
The other two cross terms use the sums 2^t−1 for the guard.
Within the guard both the weighted and the adjacent local terms
have coefficient b−a; their geometric and linear parts add to
(b−a)(2^t−2). The u/c and c/v boundaries in the last sum of (1)
do not depend on t.

For u=12,v=2,c=0, the pairs (L_a,L_b)=(1,2) and
(R_a,R_b)=(0,1) give K=0. Direct calculation at t=1 gives Δ=−2,
hence flipping A→B to B→A adds 2 for every t≥1 without editing
the guard. Both orientations are strong: B beats the first
spine vertex; A beats vertex 1 and the last; each spine vertex
reaches A or B, and the middle zeroes beat both apices. For example, A→y→B→x→A holds, and y→each zero→A.
Adding the last site, which is reached from both A and the
zeroes, preserves strongness because it points to B. The first
spine vertex and a chosen guard end give a one-apex cycle
through B with arbitrary parity-valid guard sublists; thus
the count grows with the guard height.
The guard is an effective shield, not merely transitive
source/sink padding.

## A live mirrored exterior retaining the very same +2 move

Let rc(u) mean reverse u and replace every letter z with 3−z.
Choose any word u of length k≥2 over {0,3} which starts and
ends with 3. Define

    W(u) = u 0120 rc(u).

**Theorem.** H(W(u),A→B)=H(W(u),B→A)+2 for every such u.
Both orientations are strong. There are 2^(k−2) literal
exterior choices at exactly the same order 2k+6, and only
one arc differs within each pair.

To prove it, let p=0120. Under rc the two orientations swap
by taking the converse tournament and reversing the spine,
so the matched blocks u and rc(u) give cancelling internal
terms in (1). Both outside endpoint terms are zero because
both apex rows agree at each exterior site. The weighted
cross term directly from u to rc(u) vanishes, since the
right reversed complement addresses equal the original
left addresses after reversing and complementing. The
sum of the weighted cross terms through the middle is

    L_b(u)[R_a(p)−L_a(p)]−L_a(u)[R_b(p)−L_b(p)].

Here both defect brackets equal 9, while L_a(u)=L_b(u);
thus the sum vanishes identically. The left adjacent
boundary is from 3 into 0, and the right from 0 into 0, so
both local terms of (1) vanish. The patch's own first and
last bits agree, so its two (zero) endpoint terms may be
included without a correction. Therefore Δ(W(u))=Δ(0120)=2.

For strongness, both apices beat the first spine vertex,
and the last spine vertex beats both of them. Every spine
vertex is reachable from the first and reaches the last.
Indeed the first/last contact arcs give a directed odd cycle
through A and any even subset of the spine containing
those two sites. By the positive odd-cycle-packing formula
this yields H≥1+2^(2k+2), which grows with k.

The mirrored edge is a single exact local repair at an
exponentially programmable set of exterior words. It is
not a proof that every large target number has a compatible
anchor. The exterior bits at the two ends compensate and
are not independent of one another. The exact-count
multiplicities have been checked at widths 2 through 13 with no
within-width collisions, but
literal address multiplicity alone proves neither that
all anchors are distinct nor that their spacing is uniformly
bounded. This is the next prescribed-baseline problem.

## A longer atom with both exterior banks independent

The mirrored shield ties its right exterior to its left. There is also
a composable shield whose two exteriors vary independently. The atom

    P = 133220032220,  Q = rc(P) = 311103311002

has marked addresses

    (L_a(P), L_b(P)) = (1950,3600) = 150 (13,24),
    (R_a(P), R_b(P)) = (2145,3960) = 165 (13,24).

For Q the two address vectors interchange. The literal subset DP
gives H(P,A→B)=186,399 and H(P,B→A)=186,397, and (1) gives Δ(P)=2.
Taking the converse gives Δ(Q)=−2.

For nonempty words U,V the exact concatenation identity from (1) is

    Δ(UV)=Δ(U)+Δ(V)+D(L(U),R(V))+λ(last(U),first(V)),
    D(L,R)= L_b R_a − L_a R_b,
    λ(x,y)=a_x−b_x+b_x a_y−a_x b_y.                    (3)

The last correction includes the loss of the two old interior
endpoint terms. Both address vectors remain proportional to (13,24)
under arbitrary concatenation of P and Q, by the ordinary shift
laws; hence D vanishes at each join. The possible outgoing boundary
letters are 0 and 2, and the possible incoming ones are 1 and 3.
For all four ordered pairs λ vanishes: if x=0 both terms are zero;
if x=2 then a_x−b_x=1 and b_y=1. Therefore

    Δ(any P/Q word) = 2 (#P − #Q).                      (4)

For example put a single P between arbitrary neutral pairs PQ or
QP to its left and to its right. Every such construction has exactly
the same physical +2 edit, a single reversal of the shared mutual arc.
If there are k neutral pairs split between the two sides, there are
2^k literal choices, completely independent across the split.
More generally for a given n and d with |d|≤n and d≡n (mod 2),
there are binom(n,(n+d)/2) literal words of n atoms with prescribed
signed edit 2d. Each atom has both orientations strong: B beats its
first site, its last beats B, and A interacts in both directions with
its spine; the first reaches the last along the spine. Concatenating
them preserves strongness because their strong induced subgraphs
share both apices.

The factor one here is exact at the **port edit**, and the auxiliary
choices are not copies of a scalar H. It is not a theorem that the
strong tournament has H=1 or that these words cover all prescribed
baselines. Also, the gap range at fixed atom length is only linear
in the number of atoms.

### Why simply repeating an aligned pump cannot widen that range exponentially

Fix nonempty U,G,V and consider U G^t V for t≥1. Write q=2^|G|.
From the binary shift laws, the leading coefficients of the two
address vectors as affine functions of q^t are positive scalars
times

    L(U)+L(G)/(q−1),   R(V)+R(G)/(q−1)

respectively. Summing (3) over the guard copies, the coefficient of
q^t in Δ is D of exactly these two vectors. For example the internal
geometric cross coefficient is D(L(G),R(G))/(q−1)^2, from

    Σ_{i<j among t copies} q^(j−i−1)
       = (q^t−1)/(q−1)^2 − t/(q−1).

The internal local corrections and within-copy gaps grow only
linearly in t. If both address vectors of U G^t V remain on a
fixed common marked ray for infinitely many heights, their leading
vectors are on that ray too. Therefore the coefficient of q^t is
zero, and Δ is affine in t. An exponential widening must change
the ray, use a different seam observation, or avoid fixed repeated
pump words. This is not an obstruction to arbitrary mixed and
adaptively changing constructions.

### An all-width exclusion for power-of-two marked ratios

Suppose both marked vectors of a tape are positive and lie on one
ray of slope B/A=2^j with j>0. Write A,B for its two forward binary
addresses, n for its width, and C=2^n−1. Then B=2^j A and, since the
shift fits in width n, rev(B)=rev(A)/2^j. The second ray condition
is C−rev(B)=2^j(C−rev(A)), which forces

    rev(A)=2^j C/(2^j+1).

The order of 2 modulo 2^j+1 is exactly 2j: it divides 2j since
2^j=−1, while no positive exponent ≤j can return 1 by size.
Consequently 2j divides n, and the explicit geometric series gives
A=(0^j 1^j) repeated and B=(1^j 0^j) repeated as binary words.
The tape itself is therefore (1^j 2^j) repeated, invariant under
reverse-complement. The mutual reversal has gap zero. Slope
2^(−j) is the same after swapping apices, and slope 1 forces
the rows identical. This supplies a precise reason to look for a
mixed marked ray such as 24/13 rather than a power-of-two shield.
