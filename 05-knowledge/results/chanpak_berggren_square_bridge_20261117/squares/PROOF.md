# Scaled three-square chains: exact cycle births and the hidden rotation word

## Provenance and scope

The primitive chain is inherited, not new: [THM-4510 / reflection-orbit source](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/procgen_sumgraph_20260926_reflection_orbits.md), Corollary A3''', gives the three-square Hamilton path from every primitive triple, including the size-15 chain from 8 to 9 using sums 9,16,25. The same archive uses two reflections to explain the rotation. The new derivations below retain the nonprimitive residue fibers, identify exactly when their cycles are born, recover their colored carry word, and use that readout to separate infinite families that have identical underlying graph.

Compare the separate [Berggren norm-one address/absorber theorem](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/square_sum_berggren_absorber_bijection_20261019.md) and [Berggren guard transport](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/collatz_mod6_20260921_berggren_edge_transport.md). The former is a branch-preserving correspondence to an entirely different shell; the latter retains an actual Collatz edge guard and does not automatically intertwine all branches. The [Price tree](https://arxiv.org/abs/0809.4324) organizes the same primitive triples in a different genealogy. Both trees can supply these graphs, but no graph-local action realizing either set of three child maps is asserted.

Only the three named sums are allowed in our graphs. If other square sums are permitted, additional arcs can destroy the component classification.

## 1. Exact component structure at the natural cutoff

Let A<B be coprime positive integers with A+B=C odd, and let d≥1. Form the simple graph on vertices 1 through dB−1 with an edge when its two endpoints sum to dA, dB, or dC. For Pythagorean triples take A=a², B=b², C=c² and d=s².

Use residues x in 0 through M−1, M=dB, and the two modular reflections

    t(x) = −x (mod M),   r(x) = dA−x (mod M).

The ordinary t edges have sum dB. The ordinary r edges have sum dA when x≤dA and dC when x>dA. Remove all loops and the virtual vertex zero. The only ordinary edge lost in removing zero is 0–dA. Thus this is exactly the stated sum graph. The product t∘r is rotation by −dA, which is transitive on each residue class modulo d by gcd(A,B)=1. Both reflections swap residue r with −r modulo d. Therefore their orbits are the union of each pair of residue classes ±r, or one class when it is self-opposite.

* For 1≤r<d/2, the orbit has 2B vertices. Both reflection matchings are fixed-point-free and distinct on it (t(x)=r(x) would imply B divides A), giving a connected simple cycle of length 2B.
* The zero residue orbit before deleting zero is a B-vertex path. Indeed the rotation is transitive, all degrees are ≤2, and the only two missing matching incidences are at zero and at the midpoint of the even one of A,B (multiply by d). Removing zero, an endpoint, gives the B−1-vertex path, with endpoints dA and that midpoint. In Pythagorean applications B≥16; in the abstract B=2 one gets a singleton and the two formal endpoints coincide.
* If d is even there is one other self-opposite orbit, residue d/2, a B-vertex path. Its two endpoints are dC/2 and d times half of the odd one of A,B. To see this, inspect the fixed points of the two modular reflections: t fixes x=M/2 in this orbit just when B is odd; r fixes x=dA/2 just when A is odd, and it always fixes x=dC/2 in this orbit.

Thus there are exactly floor(d/2)+1 connected components. For (a,b,c)=(8,15,17) and s=3, the sums are 576,2025,2601. Through 2024 they give one 224-vertex path with endpoints 288,576 and four disjoint 450-cycles. The decomposition does not require inspecting 2024 individual vertices.

## 2. Cycle births and all-order component counting

Once x≥dB, no dA- or dB-edge is possible; its only possible neighbor is dC−x≤dA. Hence the new vertices above the natural cutoff are leaves attached to existing components until dC−1, then isolated vertices. This leaves the exact cycle structure unchanged beyond the cutoff. For smaller cutoffs, the graph is induced from the previous paths and cycles; the only possible cycles are entire generic fibers. The max vertex of generic fiber ±r is dB−r. Thus the first Betti number is, for every positive integer N,

    β1(N) = #{1≤r≤floor((d−1)/2) : N≥dB−r}.

The edge count for sum T at cutoff N is

    e(T,N)=max(0, min(N,floor((T−1)/2))−max(1,T−N)+1),

by selecting the smaller endpoint. Hence the component count for all N is

    c(N)=N−e(dA,N)−e(dB,N)−e(dC,N)+β1(N).

There is at most one cycle in each component for every N. In the earlier 2024 control the four cycle births are exactly 2021,2022,2023,2024, one immediately after another.

## 3. Colored cycles encode a balanced carry word

For a generic orbit and initial x₀=r∈(0,d/2), traverse alternately by r(·) and t(·). The vertices in the original residue class are

    x_j = (r−jdA) mod dB.

A step with 0<x_j<dA uses the dA-edge, otherwise it uses the dC-edge. Put u=r/d and label the event dA by 1. Its exact bit is

    b_j=floor((u−jA)/B)−floor((u−(j+1)A)/B).

In particular every cyclic block of length ℓ has weight differing from ℓA/B by strictly less than 1, since the floor differences telescope. Over a period there are A dA-edges, B−A dC-edges and B dB-edges. The colored word therefore reads the integer A, not just a phase or a near-equal real ratio.

The slope family A/B from primitive triples is dense in (0,1): rational Euclidean parameters in (0,1) give dense ratios between the two legs, and scaling a rational triple down to primitive leaves its leg ratio unchanged. This supplies exact finite carry words approaching any desired slope, **not** exact finite realization of a prescribed irrational slope or of a tournament target.

## 4. Infinite indistinguishability of uncolored graphs

The cutoff graph type depends only on d and B, whereas its marked edge counts recover A as soon as d≥3. There are infinitely many pairs of primitive triples with the same larger leg b and different smaller leg a. To construct them, put x=6t≥6, b=(x+1)(x+3)(x+5)(x+7). The four factors are pairwise coprime: they are odd, pair differences are 2,4,6, and neither x+1 nor x+7 is divisible by 3. Use two factorizations (b=UV) with

    (U,V)=((x+7)(x+3),(x+5)(x+1))
    (U',V')=((x+7)(x+5),(x+3)(x+1)).

Both ratios are >1 and at most ((13·11)/(9·7)=143/63<1+sqrt(2)) by monotonicity in x of each larger factor divided by the corresponding smaller factor. Thus both ((a,b,c)=((U²−V²)/2,UV,(U²+V²)/2)) are primitive with a<b: U,V are odd coprime and (a/b=(U/V−V/U)/2<1). The two ratios differ, giving different a. At fixed square scaling d=s²≥4 the cutoff graphs are isomorphic uncolored (same paths and cycles) and distinguishable by the actual marked cycles (different small-edge count). This is a concrete example of exact phase data hidden by an uncolored quotient, but it is not a Hamiltonian-path-count construction for tournaments.

## 5. The modular inverse is visible in the split just before closure

At (N=M-d), (M=dB), the only two missing vertices of a generic fiber (pm r) are (M-r) and (M-d+r). In the alternating cycle order
(x_{2j}=(r-jdA) mod  M, x_{2j+1}=r(x_{2j})), the vertex (M-r) is at index −1 (it is joined to vertex r by the B-reflection), and (M-d+r) is at index (2j), where (j) is the inverse (A^{-1}pmod B) chosen from 1 through B−1. The two intervening open pieces consequently have (2j) and (2B-2j-2) vertices, ignoring a zero-length piece.

The first top vertex (M-d+r) arrives for ascending r and joins the pieces into one path; if a piece was empty it merely extends the remaining path. The second top vertex (M-r) arrives for descending r and closes the cycle. Other actions in the self-opposite half residue when d is even are independent. Thus the (unlabeled) two-piece length pair reads the modular inverse up to replacing (j) by (B-1-j), and the labeled alternating order reads it exactly. In the control ((A,B,d)=(64,225,9)), (j=109), and at (2016) the four generic fibers each split as (218+230), besides the 224-vertex zero-residue path. Then vertices 2017–2020 perform four separate path joins and 2021–2024 perform four separate cycle closures.
