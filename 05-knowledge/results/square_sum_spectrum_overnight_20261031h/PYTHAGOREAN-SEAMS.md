# Two square ports give an exact, reusable Hamiltonian product

**Status: PROVED scoped composition theorem; exact finite controls.** A carrier is a specified square-sum subgraph on a finite set of distinct positive integers. It need not be induced, and its vertex set need not be an initial interval. Hamiltonian cycles here are undirected edge sets, with rotation and reversal identified; this differs from the raw directed path counts on the tournament side.

## Port type and join

Suppose carriers G and J have the same designated edge e={a^2,b^2}, where a^2+b^2=c^2. Choose positive integers r,k,l with r^2+k^2=l^2 and such that the square-scaled vertex sets r^2 V(G) and k^2 V(J) are disjoint. Remove the two scaled copies of e, and insert exactly

    {r^2 a^2, k^2 a^2},    {r^2 b^2, k^2 b^2}.

Their sums are(a l)^2 and(b l)^2. All other edges remain square sums under square scaling. Here is a universally available choice with controlled label size. Let M be the maximum input label, and let m be the smallest even integer at least floor(M^(1/4))+2. Set

    r=m^2-1, k=2m, l=m^2+1.

Then gcd(r,k)=1 and r^2>M. A collision r^2 u=k^2 v would force r^2 to divide v, impossible for1<=v<=M. This is an arithmetic disjointness proof; the two sets may interleave numerically. It is stronger than requiring a large empty interval between the source copies. Since m<=4M^(1/4), every output label is at most256M^2.

**Exact product theorem.** Hamiltonian cycles of the joined graph are in bijection with pairs of Hamiltonian cycles of G and J that contain their designated edge e. The bijection cuts e in each input and joins the resulting paths by the two seams.

**Proof of the reverse direction.** The joined graph has exactly two edges across its displayed source partition. Every Hamiltonian cycle must cross this nonempty cut a positive even number of times, hence use both. Inside either side it is a spanning path with the designated endpoints. Restore the deleted e to obtain a Hamiltonian cycle of the original carrier. These operations are inverse and retain the original edge sets and source sides. Merely preserving a minor or a path count would not suffice.

If the inputs have maximum degree3, so does the output: a cut port loses one edge and gains one, and all other incidences are unchanged. The theorem permits variable input edges, but their presence must be compiled jointly with all other local constraints. An invariant input edge costs no branch bit.

## Why the port reproduces

Each new seam has two square endpoints and square sum. For example its endpoints are(ra)^2,(ka)^2, with hypotenuse la. It is used by every output Hamiltonian cycle, so either seam is a new invariant port. Thus the operation can be iterated on two copies of its own output. If an input has exactly2^b cycles and an invariant square port, depth d of binary duplication gives exactly2^(b*2^d) cycles on2^d copies of its vertex set, with a new usable square port. Counts follow from the bijection at every depth, not from finite tests.

More generally, at each level one may join different carriers sharing the same normalized square port; using the same r,k at that level preserves the common output port type. The full source tree, scaling factors and choices give a reversible address. Original numerical labels change by the stated injections; they are not silently treated as unchanged native ports.

**Label-bit bound.** In binary duplication, M_(d+1)<=256 M_d^2. Therefore

    log2 M_d <= 2^d (log2 M_0+8)-8.

The vertex count is2^d times the original count, so the number of bits per label grows at most linearly with the output vertex count. The very large integer labels do not require superpolynomial-size explicit certificates. This bound is not claimed optimal.

## An infinite authenticated input family

For every N=50t+12 with t>=56, the inherited A terminal table has quotient m=2t, residue12, and contains adjacent atoms `single(9), head(9)`. The second atom has value25-9=16. Consequently9--16 is an invariant edge for every branch choice, regardless of the deeper input words or refinements. The pair is checked directly in the complete symbolic table; the all-size inference uses its fixed atom semantics.

Use the primitive refinement from the inherited [closed-carrier theorem](../square_sum_spectrum_overnight_20261031f/README.md) to obtain the exact graph language. At N2812 and2862 there are25 independent bits in each input and all their Hamiltonian cycles use9--16. The join therefore has exactly2^50 Hamiltonian cycles on5674 vertices, maximum degree3, and a new invariant square port. A second self-join has exactly2^100 cycles on11348 vertices.

The test program walks all64 choices in a six-bit product fibre, restores both input cut edges, and recovers the exact original source edge sets and bits. It also checks four full100-bit choices after the second join, every square edge, disjoint ownership, degree two, connected traversal, and that precisely the two declared seams cross the partition. The upper count comes from the reverse-direction theorem and inherited exact input languages.

## Hostile boundary and next problem

An ordinary square-sum edge need not be this port type. The invariant1--3 edge fails the square-endpoint condition. Under the same-coordinate scaling seam, its two required sums would be r^2+k^2 and3(r^2+k^2); they cannot both be nonzero integer squares. That failed implication is why square endpoints are retained in the type.

Coprimality is also essential to the size-only separation guard: the nonprimitive legs6,8 satisfy6^2>16, but6^2*16=8^2*9=576. The primitive-leg condition prevents that erased source-side ambiguity. The actual implementation separately checks disjoint sets as a consequence audit.

The low port is not present in every A terminal state. The first failed index in the tested initial50-state block is2815. The repaired theorem uses the exact quotient parity0/residue12 state, hence N=50t+12, rather than extrapolating the port from a few adjacent examples.

This construction lives on a sparse union of scaled sets. Additional arithmetic edges may exist in the full induced square-sum graph, so the exact count is for the selected carrier. The next useful strengthening is to place such a two-seam join inside an interval constructor without losing its exact ownership and phase conditions. Large scaling by itself does not fill the gaps between labels.

## Stronger port classification: a common square class is exactly enough

The square-endpoint case is part of a complete type criterion for this same-coordinate seam. Let the edge endpoints be A,B>0 with A+B a square. There are positive integer scales r,k making both A(r^2+k^2) and B(r^2+k^2) squares, with arbitrarily large primitive r,k and hence disjoint copies, **if and only if A/B is a rational square**.

Necessity follows by dividing the two new square sums. For sufficiency, write

    A=d a^2, B=d b^2,

where d is squarefree. Divide a,b by their gcd to obtain coprime u,v. Since A+B is a square, there is an integer w>0 with u^2+v^2=d w^2. Fix this primitive norm seed once for the class. For t a sufficiently large positive multiple of2dw put

    r=(t^2-1)u-2tv,
    k=(t^2-1)v+2tu,
    l=w(t^2+1).

Then r+ik=(u+iv)(t+i)^2 in the Gaussian integers, so r^2+k^2=d l^2. The new seam sums are A d l^2=(d a l)^2 and B d l^2=(d b l)^2.

These legs are primitive. If a prime p divided both, the matrix with diagonal t^2-1 and off-diagonals -2t,2t shows that p must divide t^2+1, since u,v are coprime. Thus p is odd, and the two equations imply p divides u^2+v^2=d w^2. But then p divides t, contradicting p dividing t^2+1. Taking t>=max(4v+4,2(floor(M^(1/4))+1)) guarantees r>0 and r^2>M. The same coprimality separation proof applies.

For a fixed norm seed, t=O(M^(1/4)), the maximum output label is O(M^2), and the label-bit recurrence remains linear in the number of binary leaves. The code retains the fixed seed through repeated joins rather than recomputing increasingly expensive seeds from the newly scaled port.

Both seam endpoints remain in square class d, so the type reproduces. This is an exact use of Gaussian norm multiplication and a congruence guard, not an analogy between graphs and number fields. It also pinpoints the obstruction at1--3: its two endpoints have different square classes. Square class7 cannot occur as a positive same-class square-sum port, since u^2+v^2=7w^2 has no primitive solution (reduction modulo7 forces7 to divide both u,v). No claim transfers that obstruction to norm forms in other quadratic rings.

### Actual nontrivial classes and a paid port

The inherited carrier at N2812 has invariant port20--80, in class5 with primitive norm seed(1,2,1). At N2828,578--1922 is an invariant class2 port with seed(17,31,25). Both give exact2^50 first joins and2^100 second joins.

The port160--1440 at N2812 has class10 and seed(1,3,1), but it requires source bit22=1. Conditioning both inputs pays one bit each, giving exact2^48 first joins and2^96 second joins. The retained pin/free-bit sidecar maps each compacted source code back to the original25-bit code. One cannot silently report2^50 for this case.

The conditioned carrier keeps its chosen base edges and only the unfixed switches. Its Hamiltonian cycles are exactly the required-edge face of the original closed carrier: every surviving cycle was an original cycle, and the excluded alternative matching has a removed edge. Thus the graph-complete product theorem applies after conditioning, with the full shared-port constraints retained.

## Synchronizing different port shapes: an explicit duplication cost

Equal square class alone does not make two original ports identical. There is nevertheless an exact repair with a declared cost. Let the two invariant ports be(A,B) and(C,D), with all four labels in one square class d and each pair a square-sum edge. Write A=d a^2 and C=d c^2. Choose one primitive norm frame r^2+k^2=d l^2 large enough for both source carriers.

1. Join two copies of the first carrier with that frame. Its new port is A{r^2,k^2}.
2. Join two copies of the second with the same frame. Its new port is C{r^2,k^2}.
3. Square-scale the first output by c^2 and the second by a^2. Their ports are now exactly equal: d a^2 c^2{r^2,k^2}.
4. Apply the ordinary two-seam join.

If the original carriers have M and N Hamiltonian cycles through their ports, the result has exactly M^2 N^2. Every source occurs twice, and both histories remain independently recoverable. This does not assert the unpaid M N join for arbitrary different port shapes.

The actual example joins the3:4 port9--16 of the2812-vertex carrier with the7:24 port196--2304 of the2861-vertex carrier. A direct unaligned join is correctly rejected. Synchronization then produces exactly2^100 cycles on11346 vertices; four full branch codes are checked by actual square-edge traversal. This turns equality of square classes into a constructive compatibility theorem after a precisely charged duplication, instead of silently forgetting the port shape.
