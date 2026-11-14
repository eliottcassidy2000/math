# Reader gauges and a physical source/reader swap

Status: the identities, swaps, and displayed all-height families are **PROVED**
given the locally inherited Hamiltonian odd-cycle packing identity. The small
reader maximum is **FINITE-EXACT** only in the stated universe. No cofiniteness
or Chan–Pak proof follows.

## Two old templates, and what actually transfers

* [Planar Jacobian target-shear gauge, THM-2230](../../../01-canon/theorems/THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient.md) separated a *response fiber* from the unobservable gauges in its coordinate presentation. Here a physical marked reader is compressed by projecting onto its response to a specified source family. With two unbounded heights the observables are the ordinary path count h, the slope f=h+E+O, and the correction E. The exact switch proof below changes actual marked tournaments inside an (h,f) fiber; merely subtracting two formal readout rows does not produce a reader.
* [The Lonely Runner two-clock tomography, THM-2624](../../../01-canon/theorems/THM-2624-two-clock-root-tomography-and-disjoint-carrier-holotopy-boundary.md) exposes the difference between seeing coordinates by signed recovery and physically joining two carrier charts. Our shifted factor equality below is a genuine physical overlap: the factors agree as integers for every admissible pair of clock heights, with literal source and reader on both sides. Its sparse range does not repair the old Lonely Runner descent or solve absolute-target selection here.
* The neighboring [THM-2680 dilation-reversed fiber product](../../../01-canon/theorems/THM-2680-dilation-reversed-two-edge-clock-fibre-products-and-source-drift-boundary.md) insists on retaining source labels through a clock reversal. The tournament analogue keeps the literal normalized pair (h, oddpart(Q)) on the source side, (h, oddpart(f)) on the reader side, and matches both before swapping. The offset is a sidecar, not automatically transported. The related [Christoffel smooth ladders, THM-4716](../../../01-canon/theorems/THM-4716-balanced-collatz-cycles-are-three-smooth-ladders-and-are-free.md) motivated looking for an ordered unit defect step rather than for saturation of residues; no Collatz identity is used in the proof.

## Exact seam and endpoint coordinates

For a marked tournament (R,q), let h=H(R). Let E and O be the sums over
directed paths ending at q of respectively even and odd vertex count, with
each summand weighted by the Hamilton path count of the deleted complement
(the empty tournament has count 1). Put f=h+E+O and σ=O−E.

Join the source T before the reader R, and reverse just the seam arc
p→q to q→p, where p is the active source root. Then the exact count is
$
 hH(T)+2O E_T+2E O_T
$
where E_T,O_T are the weighted **starting**-path coordinates at p. In the
odd-cycle packing formula the inherited packings stay inside the two
blocks. Every new crossing odd cycle uses p and q and the reversed arc.
It consists of one rooted path starting at p in T and one path ending at q
in R, of opposite parities. The forward connection cannot be the deleted
arc because opposite parities exclude two singleton paths. Once the cycle
is removed the remaining graph is ordinal. The weighted formulas follow.
A successive reversed-seam reader chain, whose port moves to the new block
each time, is one composite reader marked at its first entry; thus the
formula is not restricted to small blocks.

The signed path coordinate has a useful interpretation: σ is exactly the
ordinary number of Hamilton paths **starting** at q. For proof of the dual
starting-path version, let B=R−q, its adjacency A, X=diag(x_i), and
calculate in the commutative ring with x_i²=0. Set M=I+AX and
K=I−AᵗX−X, so M=K+1 1ᵗX by the tournament adjacency relation.
Since K=(I−AᵗX)(I−X), multiplication on the left by X gives
XK⁻¹=X(I−AᵗX)⁻¹. The unrestricted Hamilton-path subset series is
G=1+1ᵗXK⁻¹1. The rank-one inverse formula yields GM⁻¹1=K⁻¹1.
If a is the row of incidences q→B, the starting signed sum is the
coefficient of the full squarefree monomial in
G(1−a X M⁻¹1)=G−a X(I−AᵗX)⁻¹1.
The second term counts B-paths ending in an outneighbor of q; precisely
the paths ending in an inneighbor survive and extend to q. Converse gives
the finishing-path version used above. In particular σ is nonnegative.

## Two free heights and the exact physical swap

Use the source C3[p,A,B] with p→A→B→p, B=TT_m, and D=Q_A−h_A for
h_A=H(A) and Q_A the ordered two-path cover count allowing empty paths.
By signed color splitting of odd-cycle packing blocks, the unrestricted
weighted two-way partition in nonempty B has vanishing signed sum.
A path starting at p alone, or proceeding through a nonempty A-piece then
an optional B-piece, therefore has
(E_T=D2^{m-1}), (O_T=h_A+D2^{m-1}).
For H(T), the odd cycle through p uses nonempty pieces from both A and B
with equal vertex parities; the two signed convolutions over nonempty
pieces each have value −h_A or −1. Adding packings avoiding p gives
(H(T)=2h_A+D(2^m-1)).
Consequently the output is
$
 N=D(f2^m-h)+2(h+E)h_A.
$
For A a transitive pad of k vertices followed by a fixed nonempty module U,
(h_A=h_U), (Q_A=Q_U2^k). Hence
$
 N=(Q_U2^k-h_U)(f2^m-h)+2(h+E)h_U.
$
This yields an explicit membership test for a named pair of reader and
module: factor N−2(h+E)h_U into the two prescribed shifted dyadic factors.
A positive factorization is not enough unless it is of the specified forms.

If two physical source/reader pairs have exchanged factor types
((h_U,\operatorname{oddpart}Q_U)) and
((h,\operatorname{oddpart}f)), swap clock variables and add dyadic
shifts to match them. Then the product term agrees absolutely. The
difference of the offsets is (2(Eh_U-E'h_R)), up to which side is
called first; it cannot be read off from factor matching alone.
For example, for k≥2,m≥1,
$
 3(2^k-1)(2^{m+1}-1)+{6,8}
$
comes from a true source/reader swap:
* for +8, take A=TT_k, B=TT_m, R=C3 pointed anywhere;
* for +6, take A=TT_{m-1}→C3, B=TT_{k-1}, R=TT1.
All graphs are literal reversed-seam joins. The h=3 pair can also be
realized without swapping via the E=0/1 fiber; the claim here is the
factor-exchange mechanism, not exclusivity.

Positive height bounds matter. A finite bank of such fixed-module,
fixed-reader ordered factor choices has only O(log² X) distinct outputs
up to X. A proof of coverage would need target-dependent reader or source
supply, not only more points inside any finite bank.

## Literal wider banks

The four readers of order eight below have the same h=73 and f=146.
The entry q is 6, with zero-based rows written as outgoing bit masks.

| E | Rows |
|--:|--|
| 27 | 8,5,9,82,39,15,183,63 |
| 28 | 2,8,3,53,39,71,159,63 |
| 29 | 2,4,9,51,39,71,159,63 |
| 30 | 2,12,9,81,7,31,183,63 |

For every k,m≥1, A=TT_k and B=TT_m give the four consecutive odd outputs
$
 (2^k-1)(146\cdot2^m-73)+{200,202,204,206}.
$
The first window is 419–425, odd values only. Here the shared factor
also equals (73(2^k-1)(2^{m+1}-1)).
The fixtures additionally give:
* 208·2^m+189,191,…,207 for every m≥1 (ten odds);
* 27·2^k+{−5,−3,…,11} for k≥6 (nine odds);
* 135·2^k+{−1,1,3,5,7,9,11} for k≥4 (seven odds); this
  actually changes A among TT2, C3, and a four-vertex (h,Q)=(5,32)
  source, and changes reader and B height.

Run `python literal_windows.py` for live witnesses and 76 complete graph
counts checked by independent Hamilton subset DP. This checks the
polynomial coefficients as well as two or more small height settings.
The all-height conclusion uses the proved equations, not extrapolation.

For the finite maximum at fixed (h, oddpart f), run
`g++ -O3 canonical_readers.cpp -o /tmp/reader-census` and
`/tmp/reader-census 8 > reader-census.txt`, then
`python reader_census.py reader-census.txt`. The exact class counts are
1,1,2,4,12,56,456,6880 at orders 1–8. Every isomorphism class has at
least one deletion onto an order-(n−1) class, and every incident mask is
tried on that canonical parent. Degree-sorted relabelings with every
within-degree permutation give a canonical adjacency code. The
enumeration agrees with the separately exhaustively labeled universe of
2,732 distinct reader triples through order seven. The normalized fiber
maximum increases from three through order seven to four through order
eight. This assertion concerns one marked reader in the stated size
universe, not arbitrary composites, all possible source changes, or a
global maximum.

## Absolute small-E rigidity

E=0 iff q has no incoming arc; then f=2h. E=1 forces h in {1,3}:
the only predecessor u of q has at most one predecessor y among the other
vertices, since any arc x→y among them gives an additional even path
x→y→u→q. When y exists it is their source and the initial block is a directed
triangle. For E=1 the remaining graph in either case has H=1 and is
transitive. Thus the pointed reader is transitive or a directed triangle
followed by a transitive tail.

E=2 forces h in {1,5}. If q has one predecessor u, the residual after
deleting u and q contributes an odd number at least one and hence must
contribute one, as otherwise it alone exceeds E=2. It is transitive.
For each predecessor y of u in that transitive spine, the extra even paths
to q are in bijection with odd subsets of vertices earlier than y.
This forces y at the second position, optionally with the first position
also a predecessor; that four-vertex core has exactly two cyclic triangles
and H=5, and all later vertices are a transitive ordinal tail.
If q has two predecessors u→v, all 2-vertex-deletion complements are
transitive and no even path of order four exists. No other x beats u, as
x→u→v→q would be such a path, and none beats v because u→x→v→q would
be such a path. All other vertices then form the transitive tail of a
transitive reader. These all-order bounds show why arbitrarily large
reader order cannot be used as a license to assume low correction digits.
