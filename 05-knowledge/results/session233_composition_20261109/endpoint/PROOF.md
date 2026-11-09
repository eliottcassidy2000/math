# Constructive endpoint universality from a binary near-transitive boundary

Status: PROVED, independently audited ([report](../endpoint_audit/REPORT.md)); no global Chan–Pak completeness claim.

## Inheritance

THM-4988-endpoint-721-boundaries-and-faithful-automorphism-graphs.md establishes that the endpoint orbit fibre has cardinality e_T(x)/|Aut(T)_x|, exhibits local 7 and 21, and expressly leaves endpoint-spectrum completeness open in its source packet. THM-250-flip-formula.md supplies the near-transitive single-reversed-edge family and its total H=1+2^(n-2). We use this family, not a new global H construction. Its historical corollary erroneously lists 2 in its set despite the theorem excluding it; that typo is not used here.

Concept board: endpoint binary digits; complete-path parity; rooted stabilizer destruction; one-edge join response; coarse packing mass. The key change is from isolated endpoint examples to a complete boundary selector while retaining the distinction from complete counts.

## Theorem

For every integer N>=1 there is a strong tournament U with a distinguished vertex x such that

- e_U(x)=N;
- Aut(U)_x is trivial, so the endpoint orbit fibre has exactly N elements;
- |V(U)| <= max(5, ceil(log2(N+1))+3).

For each n>=4, a single rigid strong deletion tournament B_n supports all these endpoint values N=1,...,2^(n-2)-1 by changing only the incoming-neighbor set of the added x.

## Proof

Let B_n have vertices 0,...,n-1, all arcs i->j for i<j except n-1->0. It is strong because it contains the spanning cycle 0->1->...->n-1->0.

Its terminal-count vector is

    (1, 1, 2, 4, ..., 2^(n-3), 1).

Indeed the only Hamilton path avoiding n-1->0 is the increasing word, ending at n-1. Every other Hamilton path is an increasing prefix ending at n-1 followed by an increasing suffix starting at 0. If its final vertex is i with 1<=i<=n-2, all vertices greater than i must be in the prefix, i must end the suffix, and each of 1,...,i-1 can independently go into either part. There are 2^(i-1) choices. Ending at 0 has exactly one choice.

The outdegrees of B_n are (n-2,n-2,n-3,...,2,1,1). For n>=4 each score class is a singleton or a pair. An automorphism preserves every score class and cannot transpose a pair because that would reverse its internal arc. Consequently B_n is rigid.

Choose n=max(4,ceil(log2(N+1))+2). Write N=sum_{i=1}^{n-2} epsilon_i 2^(i-1). Adjoin x with i->x exactly when i is an interior vertex and epsilon_i=1; set x->i for all other vertices, including 0,n-1. Deleting the last vertex of an ordered Hamilton path ending at x yields a Hamilton path of B_n ending at one of its incoming neighbors, and this operation is a bijection. Therefore

    e_U(x)=sum_{i:epsilon_i=1} e_Bn(i)=N.

As N>0, x has at least one incoming edge, and x has outgoing edges to 0 and n-1. A vertex with both directions of adjacency to a strong subtournament extends it to a strong tournament: reach its incoming neighbor to reach x, and leave x through an outgoing neighbor to reach the old tournament. Thus U is strong.

Every automorphism fixing x restricts to an automorphism of B_n=U-x. That restriction is identity, hence Aut(U)_x=1. THM-4988 then gives the claimed endpoint quotient. The order is n+1, as asserted.

## What this does and does not imply

This gives every positive boundary count, including the global forbidden numbers 7 and 21. It does not give H(U)=N. The complete count is the sum of all terminal counts, and all the other coordinates change when x is attached. Exact witnesses from this selector have (N,H)=(7,23),(21,187),(37,461),(101,1573),(257,9899). Neither local digits, nor trivial root stabilizer, nor a boundary-orbit equality supplies global arithmetic target placement.

Likewise, adjoining ordered C3 components before U multiplies e and the root stabilizer order by the same power of 3 and preserves the endpoint quotient N, so every positive quotient has towers with nontrivial stabilizers. This follows by ordered SCC concatenation and does not turn N into a complete path quotient.

The one-reversed-edge join specialization was investigated and recovered, not claimed new: THM-4865(f) credits HYP-9457(c) with the two-base response formula. Product root stabilizers extend through this join, so rooted P7 contributions retain a factor 3; local 7/21 counts cannot bypass the stabilizer or response sidecars.

## Exact checks and reproduction

Run `python3 verify_endpoint.py` in this folder. `probe.py` is a subset DP for terminal counts. Independent literal permutation enumeration checks the entire terminal vector for every constructed tournament from bases n=4,5,6 (25 selected incoming-neighbor subsets total), including all roots needed for those values. It also computes actual root stabilizers, rather than assuming rigidity. The DP checks all 119 selectors from bases n=4,...,8 and verifies strong connectivity by reachability. Base terminal vectors and rigidity are checked independently through n=7. Larger named witnesses above are DP-only controls. These finite checks validate the implementation; the bijection and score-class proof establish all orders.


## Root-response sidecar: exact linear selection, mismatched digits

Use the response B_v(A)=sum_{S containing v} s_v(A[S]) H(A-S), with singleton paths allowed and H(empty)=1, from THM-4865. For any one-vertex extension U=A+x,

    B_x(U)=H(A)+sum_{v in N_U^+(x)} B_v(A).

The root path in an ordered two-path cover is either the singleton x, leaving an arbitrary spanning path of A; or it starts x,v, and deleting x gives exactly a rooted two-path cover of A at v. This proves the formula by disjoint cases, independently of OCF.

For our digit selector,

    B_x(U_N)=H(B_n)+B_0(B_n)+B_(n-1)(B_n)
                 +sum_{i=1}^{n-2}(1-epsilon_i) B_i(B_n).

Thus the endpoint coordinate uses binary weights 2^(i-1), while the response coordinate uses the non-binary B_i weights. This is an exact description of the obstruction to treating endpoint universality as a freely selectable term in the one-exit or two-base join law.

`response.py` checks this identity by a separate start-path subset DP for every selector at base orders 4,...,9. At n=6 the base responses are (31,69,44,28,18,17). As N runs from 1 to 15, the child H values are (27,41,33,59,51,65,39,85,71,91,59,109,77,97,47), and B_x values are (155,180,111,196,127,152,83,206,137,162,93,178,109,134,65). These nonmonotone profiles are finite-exact, not a theorem that no other composition can eliminate the sidecar.

Verification robustness: all mathematical checks raise explicit exceptions rather than using removable assertions. Both executable verification scripts have a `--fail` sentinel. Normal and `python3 -O` runs produce byte-identical outputs; the deliberate failure exits nonzero in both modes. See `verification_robustness.out` for receipts. Independent audit PASS is recorded in ../endpoint_audit/REPORT.md.


Session233 close-out note: THM-250's erroneous corollary entry 2 and its adjacent-flip proof sentence have now been repaired in canon; the formula used above is unchanged. See MISTAKE-744.
