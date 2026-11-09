# Independent audit: bare intervals and exterior rerooting

**PASS, with the exact scope below.** No counterexample was found to the candidate PROOF.md. This audit did not edit or import the producer's scripts or repository. `check.py` and `physical.cpp` are standalone. The additional unordered-set exception below is a warning against strengthening the candidate's cautious graft statement.

## All-order bare-family theorem

The graphs B_(k,j)+r, with B_(k,j)=X^j ⇒ C ⇒ Y^(k−j), are strong on 12k+9 vertices for all integers k≥0 and 0≤j≤k. Their Hamilton counts are Z_k−2j, with exactly the stated numerator and denominator. The k=0 graph has count577 and rooted state(577,31,923,45). It is a valid one-value boundary case, not a nontrivial interval.

The all-order justification is algebraic, not extrapolation of the finite checks. Literal permutation-and-split enumeration verifies the core seven-state. Ordinal composition verifies the flanked move identity, with only T changing by2. X ends0 and Y starts1, so the same flanked change propagates at every position through h=1 contexts. Recover the bare interface from the original flanked interface using q_L=8q_B, A_L=4q_B+2A_B+45, G_L=4q_B+4G_B+45. These formulas also hold when either word is empty. The prefix10 changes the apex count by precisely2(q_B−G_B); deleting a selected trailing global sink leaves its count unchanged. This gives the proposed baseline and all positions.

An especially convenient exact integrality certificate for the fraction is

    Z_0=577,
    Z_(k+1)=4096 Z_k−264150 k−118165.

It follows directly from the stated numerator. Integrality and oddness hold for every k. No numerical rounding or integer-division assumption is needed.

Strongness has four endpoint cases. The initial selected vertex of X, or a selected vertex inside the initial strong C when j=0, allows r to reach all ordinary vertices. The terminal unselected vertex of Y, or an unselected vertex inside terminal strong C when j=k, allows every ordinary vertex to reach r. This includes k=j=0 because C is strong and the cut and its complement are both nonempty.

## Physical rooted and deletion states

The claimed bare state is correct:

    (D,E,U,d)=(Z_k−2j,45−e_j,A_B+45,45),
    e_j=0 for j<k, e_k=14.

Here U follows by splitting a root-start path into r alone (h choices for the other path) or r followed by a nonempty cut-start path (A choices). E=h−e_j counts deletion paths ending outside the cut. In particular, the final E differs by−14. The bare interval is not a pure-D rooted fiber across all k+1 witnesses.

Independent full-graph subset DP verifies every root of the 9-vertex k=0 bare graph, both 21-vertex k=1 bare graphs, and the original 13-vertex k=0 and both original 25-vertex k=1 graphs. This does not use ordered-block transfer. Literal k=0 deletions additionally verify:

- delete v: (2303,0,7606,45);
- delete v and the trailing sink: (2303,31,3803,45).

Thus deleting a count-neutral global sink does change physical sidecars, exactly as required by the candidate. The original root-v state is(9909,2303,9909,2303), confirming the proposed alternative-root interpretation.

## Additional terminal-graft warning

From N=H(D−E+U)+(Q−2H)(U−d), a common nonempty graft at the bare apex changes each interior indexed step by−2H but the last by+12H. This is universal in the gadget; Q cancels. Consequently the last grafted value lands at index k−7 in the extended arithmetic progression determined by the earlier witnesses.

Do not turn this into a claim that every unordered bare window is destroyed. For k=6 and H=1, the terminal value is N_0+2 and the other six are N_0,N_0−2,...,N_0−10: the unordered set is again seven consecutive odds. For k≥7 the final value duplicates an earlier one. For1≤k≤5 it is separated by a gap; H>1 cannot supply consecutive odds. The candidate only says width preservation cannot be assumed, which is correct.

## Exterior-root-change obstruction

PASS for the stated common graft sequence and protected L_(k,j). Every old protected vertex beats every new gadget vertex, whether the chosen root is r or in R. Hence deleting r always leaves L_(k,j)⇒R. Rooting at r appends a block to R; rooting at an exterior vertex only changes the induced exterior graph and cut. Original v remains outside the cut. The induced R and its cut are common across witnesses when the operations and gadget choices are common.

Composition gives ΔT=2H(R), all other seven coordinates unchanged, hence ΔN=−2H(R). Only for k≥1 is consecutive-odd spacing equivalent to H(R)=1. H(R)=1 iff R is transitive (a nontrivial strong component has a spanning cycle and hence at least three Hamilton paths). Removing terminal selected vertices of transitive R removes global sinks; the last surviving unselected vertex is a closing v and all preceding vertices are ordinary singleton suffix padding. There is no assertion about roots in the protected block, witness-dependent gadget choices, or k=0 forcing H(R)=1.

## Finite named-port obstruction

PASS specifically for the two original k=1 source graphs. Independent physical subset DP confirms all thirteen named-port differences, including ΔJ=8, ΔV=−2 at the trailing port, and ΔJ<−2, ΔV<0 at every leading/core port. Since N=HJ+(Q−2H)V and Q≥2H, the latter inequalities exclude |ΔN|=2 for every nonempty common gadget. The trailing expression12H−2Q is divisible by4 because Q is even (pair each ordered two-cover with its swapped cover), and cannot equal±2. Roots r and v have(ΔJ,ΔV)=(−2,0). These conclusions are universal in the gadget but finite in the source pair. They must not be generalized to every k or arbitrary sequences through protected vertices.

## Reproduction and scope of controls

    g++ -O3 -std=c++17 physical.cpp -o physical
    python3 -I -B check.py > CHECK.json
    python3 -I -B -O check.py > CHECK.optimized.json
    cmp CHECK.json CHECK.optimized.json

Explicit failure gates survive Python optimization. C++ endpoint counts use uint32 storage with an explicit overflow check; all totals and two-cover products use uint64 with uint128 checked product and accumulation gates. Endpoint recurrence sums have at most25 uint32 summands, so their uint64 accumulators cannot overflow. Normal and optimized receipts remained byte-identical to the receipts before adding these overflow gates. The literal full-subset runs include25 vertices and need roughly3.7GB RAM. The source restricts graph order to25. Core states are counted by literal permutations and splits, without the producer's composition or transfer code. Finite composition checks extend through k=30 and physical reachability through k=7. The all-order conclusion is supported by the finite flanked identity plus the displayed algebraic/context argument, not by those ranges. No all-target coverage, minimality, or literature-priority claim is established.
