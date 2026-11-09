# Exact three-state transfer for bridge-linked transitive blocks

Status: PROVED, with independent exact tests and [separate audit PASS](../transducer_audit/AUDIT.md). No claim of Chan–Pak completeness. All artifacts are self-contained. Reproduction: `python3 bridge_blocks.py` and `python3 check_polynomial.py`.

## Inheritance and target

THM-4951 (`01-canon/theorems/THM-4951-packing-number-lines-tight-boundary-staircases.md`) realizes all admissible tight packing-number splits by disjoint staircase windows. THM-4952 (`01-canon/theorems/THM-4952*`) transfers general ribbon odd-cycle packing polynomials. THM-4881 (`01-canon/theorems/THM-4881*`) warns that exact tripling does not provide additive digits. We retain the lost coordinate: the entire odd-cycle packing polynomial, rather than only its degree. MISTAKE-740 warns that finite specialization tests do not alone prove identities; MISTAKE-742 warns that run-window conclusions require their exact run hypotheses.

Concept board: packing degree, actual count, window seams, open-path parity, scalar digit systems. Anchor: Chan–Pak composition. Niche: compress the opposite run-window endpoint. Wildcard: use a unique forward bridge across each cut instead of bounded back-arc span.

## Family and theorem

For a composition m=(m_1,...,m_k), define B(m) on consecutive blocks of these sizes. Inside each block all arcs point forward. Between different blocks all arcs point backward except each edge from the last vertex of a block to the first vertex of the next block, which points forward. Thus the complete base path is present.

Let P_T(x) be the odd-cycle packing polynomial, cycles weighted by x. Let v=(h,e,o) initially (1,0,0). Process blocks by v <- A_m(x)v. Then P_B(m)(x)=h at the end. Coordinates e,o count configurations consisting of completed odd cycles plus one unfinished directed path ending at the last processed vertex, with even/odd vertex cardinality, respectively. Its initial vertex is arbitrary in the processed prefix. The unfinished path is not given an x weight until it closes.

Here A_1(x)=[[1,x,0],[0,0,1],[1,1,0]]. For m>=2 put a=2^(m-2), d=m-2,

s_0=s_1=a,
t_q=(3^d+(-1)^q*(-1)^d)/2,
N(q,r)=(3^d+(-1)^(q-1)+(-1)^(r-1)+(-1)^(q+r-2+d))/4.

Indices q,r,p are in {0,1}. Matrix entries are:

- A[0,0]=1;
- A[1+q,0]=s_q;
- A[0,1+p]=x*s_(1-p);
- A[1+r,1+p]=t_((r-p) mod 2)+x*N(1-p,r).

All these expressions are integers; the combinatorial proof supplies their nonnegativity.

At the Hamiltonian-path base x=2 this simplifies remarkably:

A_m(2) = [[1,2^(m-1),2^(m-1)],
          [2^(m-2),3^(m-2),3^(m-2)-1],
          [2^(m-2),3^(m-2)+1,3^(m-2)]], m>=2.

Thus unbounded transitive block sizes require only THREE state coordinates, unlike a general ribbon family. This is a concrete alternative composition transducer, not a universal integer decoder.

## Proof

Across every block boundary there is exactly one forward arc. Any directed cycle crossing the boundary uses that arc, and exactly one backward crossing: the number of forward and backward crossings agrees. Thus at most one cycle in a packing crosses that cut, and its processed intersection is one directed path ending at the last processed vertex. This proves sufficiency of the stated three states.

If no path enters a new block, either it remains unused or a new unfinished path is an increasing subset containing its last vertex. The counts of each parity are s_q. If an old path closes in this block, choose an increasing subset containing the first vertex, with opposite parity to the old path, and close by the backward arc. This gives x*s_(1-p). Forbidden two-cycles cannot occur because closure requires odd total length at least three.

For a path that stays open, the old path can be extended both at its tail and at its head. The tail extension is an increasing subset containing the new block's first and last vertices. The head extension is a disjoint increasing subset of the middle vertices, followed by a backward arc to the previous head. Each middle vertex is unused, in the tail extension, or in the head extension. The total number of selected middle vertices determines parity. Their even/odd counts are (3^d +/- (-1)^d)/2. These are t_q. Importantly, omitting the optional head extension is WRONG: it misses the cycle 0->1->3->4->2->0 for block sizes (1,3,1).

Finally the old path can close while a distinct new unfinished path starts. The old-cycle subset contains the first vertex, the new-path subset contains the last vertex, and each middle vertex has three choices (old, new, unused). The character sum for specified subset parities is N(q,r). This gives x*N(1-p,r). These cases are disjoint and exhaustive, establishing the transfer for all x and all block sizes.

## Connection to run-window lines

If every run in w has length at least two, the staircase windows [p,q+2] corresponding to its 1-runs are disjoint. The constructed endpoint t_w is an ordinal sum of these staircases and any uncovered singletons. Consequently

H(t_w)=product over 1-runs of f(length+2),

where f(0)=f(1)=f(2)=1 and f(n)=f(n-1)+f(n-2)+f(n-3) is the staircase/tower tribonacci sequence, NOT Busch's minimum. Nonempty factors begin 5,9,17 and increase thereafter. Thus no factor is3 or11, and the primality of3 and11 excludes both values for this endpoint.

The complementary endpoint t_w+1 is exactly B(m), where its blocks are the staircase windows and uncovered singleton vertices. Arbitrary compositions enlarge this particular family, so an impossibility for B(m) is also an impossibility for the complement of every such canonical run-window construction.

## Strongness and exact small-value obstruction

Every B(m) with at least two blocks and n>=3 is strong: the base Hamiltonian path together with the backward edge from the final to the initial vertex is a Hamiltonian cycle. For n=2 it is transitive, as is every one-block B((n)).

Use Busch's lower bound (CITED via THM-4744, `01-canon/theorems/THM-4744*`): g(3q)=3*5^(q-1), g(3q+1)=5^q, g(3q+2)=9*5^(q-1). The notation g deliberately separates this from the tribonacci f above. Since g(6)=15, a realizer of11 in B would have n<=5. Exhaustive composition transfer through n5 gives values {1,3,5,9,13}, so11 is universally excluded. Since g(7)=25, 19 would require n<=6; it is likewise absent. Thus the union of both canonical endpoints misses11, and ordinal-sum multiplication of any number of these endpoints also misses11, since11 is prime. Tight packing-number universality is insufficient even after closing this family under ordinal sums.

Since g(19)=15,625>10,000, every multi-block B with H<=10,000 has n<=18. The exact pruned state computation covers through n22 and finds precisely1,374 distinct values <=10,000. This is a complete bounded-value census, not a fixed-point or universal-completeness claim. Pruning h>10,000 is sound since the construction extends by induced vertices and Hamiltonian path count cannot decrease (insertion lemma). First gaps:7,11,19,21,35,37,39,43,47,55,61,63,...

## Verification and limits

`bridge_blocks.py` independently verifies x=2 by subset Hamilton-path DP on all511 compositions of orders1..9. `check_polynomial.py` independently enumerates odd directed cycles and disjoint packings for every composition of orders1..7 at x=0,1,2,3,5 (635 comparisons). Receipts give full matrices and small-value enumeration counts. Exact tests support the proof but do not replace it.

The initial discarded transfer omitted head extension and failed at(1,3,1):11 predicted versus13 actual. The repaired transfer explicitly retains this option and passes both independent checks. This is a new local near-miss mechanism worth recording if the theorem is adopted.

The residual necessary enrichment is precise: any composition approach built from these canonical line endpoints and ordinal sums must add a new seam or a new seed already realizing11. Adding11 alone is NOT claimed sufficient. General ribbons, noncanonical tilings with the same tight triple word, and broader modular substitutions are not excluded by this obstruction.
