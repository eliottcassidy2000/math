# Integer gauge packets: local changes that cancel on tournament-word cycles

**Status: PROVED IN NOTE; finite controls provided separately.** The graph-coboundary criterion is standard elementary graph cohomology, not a new general theorem. The explicit equal-order tournament lift, exact class/fiber formulas and their scope are the contribution here. Chan–Pak remains OPEN.

## 1. A surprising source and a precise transfer

The historical inspiration was `07-reflections/the-drift-belongs-inside-the-integral.md`: in a lonely-runner averaging problem, interior signed errors cancel while endpoint terms survive. It suggested keeping the *signed valuation change* of a local tournament replacement until the entire word closes. This is a transfer of a proof mechanism, not a reduction of lonely runner to tournaments.

The closest exact tournament carrier is Session268, `automata/PROOF.md`, Theorems1–4: four strong rigid ten-vertex cells G_ij, with

    W = [[215,179],[245,201]],    H(G_ij)=W_ij.

They are built from two five-vertex tournaments B_0,B_1 (lexicographic arc masks8,10, with a1 meaning smaller-to-larger) by reversing the cross arc0→0 in B_i⊕B_j. The four distinct H values distinguish the cells intrinsically. This lane inherits their construction and rigidity. The full path-cover congruence in Session268 `modular/PROOF.md` is stronger than scalar H preservation. The particular physical packet below in fact preserves that full invariant too, using the inherited order-join transform of THM-3166; the general weighted-graph criterion alone does not imply this.

**Hostile:** an independent operation can increase one cell count and decrease another, yet fail to cancel if the two occur different numbers of times. **Sidecar:** the actual directed adjacency word, including its initial and final states. This is the same kind of lost boundary coordinate as the compatibility objects in `07-reflections/four-frontier-compatibility-objects-and-jet-budget-synthesis-codex-20260803.md`, but these objects are not identified mathematically.

Live concepts: graph cocycles; prime-valuation transport; tournament SCC words; exact isomorphism multiplicity; boundary memory; full-cover versus scalar observation.

## 2. Finite certificate for all-word conservation

Let D be a finite strongly connected directed graph. Give each edge e a positive rational old weight w_e and new weight w'_e. For a directed walk P, write W(P)=product of its edge weights, retaining repetitions.

**Theorem A (closed-word gauge criterion).** These conditions are equivalent:

1. W'(P)=W(P) for every nonempty directed closed walk P.
2. The same holds for every directed simple cycle (loops included).
3. There are positive rational potentials lambda_v such that

       w'_uv / w_uv = lambda_v / lambda_u

   for every edge u→v. Parallel edges have the same ratio.

Then an arbitrary open walk u→v has W'(P)/W(P)=lambda_v/lambda_u.

**Proof.** A closed walk splits into directed simple cycles, giving2⇒1. Pick a root r and for each vertex v a directed path r→v. Define lambda_v as the product of edge ratios along that path. Two paths to v have the same ratio because a fixed path v→r closes both, and1 says both closed products are1. An edge u→v appended to a path to u therefore proves3. Potentials telescope on a walk, proving3⇒1 and the boundary formula. The zero-length path gives lambda_r=1. ∎

This is a finite cycle-space certificate: after choosing one path per vertex and its associated potentials, checking the |E| displayed edge equations proves conservation for every word length. Equivalently, take each prime valuation. The signed change vector is an integral graph coboundary whenever old/new weights are positive integers. In a connected graph the space of rational coboundaries has dimension |V|−1; the independent cycle constraints have dimension |E|−|V|+1. No logarithmic numerical test is needed.

For a two-state complete directed graph, diagonal ratios must be1 and the two off-diagonal ratios reciprocal. At a fixed prime, if the old off-diagonal exponents are b,c, the integral exponent redistributions are exactly

       (b+t,c−t),   −b≤t≤c.

Thus there are exactly b+c+1 nonnegative integral distributions, and their closed-word products agree. This count is a bounded lattice interval, not a heuristic similarity of decimal sequences.

## 3. A physical equal-order tournament lift

Fix integers K≥1 and0≤k≤K. Let C3 be the directed triangle and TT_s the transitive s-vertex tournament; TT_0 means no vertices. Use ⊕ for ordered sum: every earlier vertex points to every later vertex. Set

    F_k(0,1)=C3^{⊕k} ⊕ TT_(3(K−k)),
    F_k(1,0)=C3^{⊕(K−k)} ⊕ TT_(3k),
    F_k(0,0)=F_k(1,1)=TT_(3K),
    E_k(i,j)=G_ij ⊕ F_k(i,j).

Every E has exactly L=10+3K vertices. Their H weights are

    W_k = [[215, 179·3^k], [245·3^(K−k), 201]].

These matrices are diagonally gauge equivalent: relative to k=0 take lambda_0=1, lambda_1=3^k. In particular each local01 count is multiplied by3^k and each10 count divided by3^k. Both endpoints are actual integral tournament counts; no formal subtraction or fractional tournament is used.

For an indexed binary circular word x=(x_0,...,x_(N−1)), N≥1, define

    T_(K,k)(x)=E_k(x_0,x_1) ⊕ ... ⊕ E_k(x_(N−1),x_0).

This is an ordinary unmarked tournament, not a marked class. The circular indexing is the imposed compatibility condition; the SCC order of the tournament remains linear, so rotations of x are NOT quotiented out.

**Theorem B (all-size physical gauge packet).** For fixed K,N:

(a) H(T_(K,k)(x)) is independent of k for every x.
(b) For fixed k, x maps injectively to ordinary tournament isomorphism classes.
(c) For a nonconstant x the K+1 classes obtained as k varies are distinct. For either constant x all k give one common class.
(d) The union over x,k therefore has exactly

       2+(K+1)(2^N−2)

ordinary classes at order N(10+3K).

**Proof.** In a circular binary word the number r of01 transitions equals the number of10 transitions. The added triangle count is kr+(K−k)r=Kr, so the product H gains exactly3^(Kr), independent of k. Ordered sums multiply H, proving(a).

The unique ordered SCC decomposition identifies each ten-vertex G block: all other SCCs have order1 or3. The H values215,179,245,201 identify its pair label. Hence the full indexed word is recovered, proving(b). For a nonconstant word there is a01 block. The number of three-vertex SCCs after that G block and before the next G block is exactly k. Thus k is recovered too. Constant words have only diagonal blocks and all-transitive fillers, independent of k. This proves(c),(d). The final filler is visible after the final ten-vertex SCC, so the last cell is not an exception. ∎

## 4. Complete spectrum and exact multiplicities of the packet

Let n0 be the number of zeros and r≥1 the number of zero-runs in the circular word. For a nonconstant word,

    H_K(N;n0,r)=215^(n0−r) 201^(N−n0−r) (179·245·3^K)^r,

where1≤n0≤N−1 and1≤r≤min(n0,N−n0). Constants give215^N and201^N.

**Theorem C (exact observable and automorphism census).** For fixed K,N these values are all distinct, and the packet has exactly floor(N²/4)+2 distinct H values. Its nonconstant H fiber has ordinary isomorphism multiplicity

    (K+1)·(N/r)·binom(n0−1,r−1)·binom(N−n0−1,r−1).

Each constant fiber has multiplicity1. Every nonconstant class has automorphism group of order3^(Kr), so its labelled fiber on a fixed N(10+3K)-set has size

    (N(10+3K))! · (K+1)·(N/r)
      ·binom(n0−1,r−1)·binom(N−n0−1,r−1) / 3^(Kr).

**Proof.** The prime valuations at43,67,179 recover respectively n0−r,N−n0−r,r, including the constant cases; hence no values collide. Distinguish one zero-run start. There are N choices of its indexed position and positive compositions of n0,N−n0 into r parts; each word was counted r times. Multiply by K+1 using TheoremB. Summing min(n0,N−n0) gives floor(N²/4). The rigid G cells and singleton SCCs contribute no automorphisms; each of the Kr triangle SCCs contributes its order3 rotation group, independently. SCC order prevents their permutation. Orbit-stabilizer proves the labelled formula. ∎

Consequently this gives an explicit lower bound for A000568 at the relevant orders, but does not approximate the total census, which is far larger. More importantly, it is an exact *same-order change of ordinary multiplicity without changing support*. It separates the support problem Chan–Pak asks from the fiber/enumeration problem A000568 asks.

## 5. Generating functions and endpoint observability

The Hamilton mass over one k-family is trace(W_k^N), and all its positive integer moments are preserved too: replace each entry by its j-th power. Diagonal similarity remains valid. Therefore the entire scalar H distribution is unchanged, not merely its mean. The finite matrix identity

    sum_(N≥1) trace(W_k^N) z^N = −z d/dz log det(I−z W_k)

is independent of k. The same holds for formal prime-valuation variables in place of the scalar weights, with Laurent-monomial diagonal gauges. This is an exact partition-function / transfer-matrix bridge, not a claim about Yang–Mills or a continuum limit.

On OPEN binary words the distinction reappears: the output ratio between k and0 is3^(k(x_end−x_start)). Thus a scalar observer on closed words cannot recover the hidden transport potential; a single endpoint-marked experiment can. A finite family of closed tests is complete for the cocycle condition, but no number of them determines the vertex potential normalized at one vertex: all potentials are invisible to closed products.

## 6. Full response conservation, and a seam that exposes the hidden location

Let Q_T(t)=sum_d p_d(T)(t)_d be the falling-factorial path-colour polynomial. THM-3166, `01-canon/theorems/THM-3166-falling-factorial-order-join-path-colour-transform.md`, proves Q_(A⊕B)=Q_A Q_B. In particular Q_vertex=t and Q_C3=t³+2t.

**Theorem D (whole-profile gauge invariance).** For the nonconstant word coordinates n0,r above,

    Q_(T_(K,k)(x))(t)
      =Q_G00(t)^(n0−r) Q_G11(t)^(N−n0−r)
       ·(Q_G01(t)Q_G10(t))^r
       ·(t³+2t)^(Kr) · t^(3KN−3Kr).

It is independent of k and of the word beyond its edge counts. All gauge variants are indistinguishable inside EVERY iterated uniform tournament-substitution context.

**Proof.** There are Kr triangle SCCs and3KN−3Kr singleton SCCs, independently of k. Multiply their Q polynomials and those of the ten-vertex SCCs. The falling-factorial basis is linearly independent, so Q equality is full path-cover equality. Apply Session268's contextual congruence. ∎

A nonuniform seam can nevertheless reveal the stored position. Define S(A,a) by adjoining one vertex b, initially a sink, and reversing a→b to b→a. For a literal sink a, S just inserts b before a in the final transitive part and H is unchanged. For A=P⊕C3 and a in its final triangle,

    H(S(A,a))=5 H(P),       H(A)=3 H(P).

The last four vertices are the familiar strong four-vertex H5 tournament; every vertex of P still beats them. One can also derive5 from the inherited one-seam formula H−e+L, since a triangle has (H,e,L)=(3,1,3).

More generally, whenever x_0=0 and x_(N−1)=1, the k=0 presentation ends in a triangle while the k=K presentation ends in TT_(3K). The same operation gives5H/3 versus H on these equal-profile inputs, at every K≥1,N≥2. This is a controlled exposure of a finite stock of triangle factors, not an indefinitely renewable digit step.

Take K=1,N=2,x=(0,1). For k=0 the final filler is C3; for k=1 it is TT3. TheoremD says the two26-vertex tournaments have identical complete path-cover profiles, and both have H=3·179·245=131565. Choose a to be the final labelled vertex in each construction and apply S. Their27-vertex outputs have respectively

    5·179·245=219275,        3·179·245=131565.

This is **not a newly discovered3→5 operation**. The point is that the gauge moved a physical terminal resource while preserving the entire unmarked substitution response, then a single endpoint-sensitive seam exposed it. Thus we have an explicit all-order transport architecture, an exact invisibility theorem for one operation category, and a physical exit from that category. Marking the terminal vertex is indispensable; these are not well-defined functions of an unrooted Q-profile alone.

## 7. What this does and does not advance for Chan–Pak

This supplies actual all-size families of operations which do not preserve H locally but preserve it after compatible composition. It can transport constructions across fibers and grow multiplicities in a controlled, isomorphism-faithful way. It also provides a finite exact test for which proposed local multiplicative repairs can work on every closed word.

It does NOT create arbitrary new odd counts: this particular packet uses only the fixed prime support {3,5,7,43,67,179}. For example13 is absent from its H spectrum at every K,N. No gauge transformation can change the closed-word value. A gauge conservation theorem is therefore an *architecture for retaining successful outputs while rearranging internal resource*, not an existence theorem for missing outputs.

The concrete next target is a **reusable gauge-compatible seam** that goes beyond the terminal3→5 move: after exposing the endpoint resource, it must restore a supply of rooted coordinates suitable for another independently selected target step. Session268's one-seam law identifies the needed statistics (e,s,L,R). TheoremD supplies a transport and one controlled exposure, but not a repeatable arbitrary-output interface. Required sidecar: an exact transformation law and renewal invariant for those four rooted coordinates. Uniform additive coverage still needs a separate target-bearing theorem; no arbitrary subtraction or division of a tournament count is licensed.
