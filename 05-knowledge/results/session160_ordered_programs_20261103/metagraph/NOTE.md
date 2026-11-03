# A seven-state ordered engine from actual metagraph lines

Status: elementary all-order hub construction and explicit finite example. The priority order is externally supplied controller data; it is not an intrinsic orientation or order of the metagraph. Nothing here gives an unbounded universal machine or settles the tournament path-count spectrum. Compare the original line convention and caution in S151 and the reflection construction in S158. No claim is made of novelty in the literature.

## General cycle hub

Let W be an n-vertex labelled tournament, n>=3, with a directed Hamilton cycle C=(v_0,...,v_{n-1},v_0). Let Z be obtained by reversing every arc not on C; its C arcs remain forward. For every cycle edge e_g=v_g v_{g+1}, let T_g be W with e_g reversed, and P_g=(v_{g+1},...,v_g) the directed path obtained by deleting e_g from C (indices modulo n). All P_g are directed paths in Z also. If K_P reverses every arc except the n-1 arcs of P, then

    K_{P_g}(T_g)=Z, and K_{P_g}(Z)=T_g.                 (1)

Indeed, the non-path pairs comprise e_g and the non-cycle edges, and exactly these differ between T_g and Z. Thus T_g and T_h are joined by an actual two-line walk through Z. On the original vertex labels, the resulting composition reverses exactly the two cycle edges e_g,e_h (and is the identity if g=h). No claim is made that the two-line walk is geodesic, that its two lines are different when g=h, or that all nodes are distinct after quotienting.

In the converse-merged line multigraph, these lines can collapse, identify nodes, or form loops. Distinct T_g and center classes, and individual line identities, must be checked before calling the quotient a star. In particular, a line color is a marked-path property and cannot automatically be read from its endpoint class.

## Blue obstruction to a full quotient star

If P_g is blue, reversal of P_g is an anti-automorphism of Z. It preserves C as an undirected cycle and acts on its edges by a reflection fixing e_g. The same reflection is an anti-automorphism of W: reversing all non-cycle arcs commutes with a permutation preserving C. Thus it sends T_h to a converse copy of T_{reflection(h)}. For n>=3, a reflection has at most two fixed cycle edges and some nonfixed edge. Hence at least two leaves become the same converse-merged class. Consequently, if all n leaves survive as distinct merged nodes, all n exhibited lines are black. For a blue spoke the number of distinct merged leaf classes is at most floor(n/2)+1. If n is even, the reflection also fixes the opposite cycle edge, so that opposite spoke is blue too. This uses a witness reflection and does not assert that all leaves within a reflection orbit account for all possible identifications.

## Seven distinct states

On vertices 0,...,6, let W orient i→j for i<j exactly for

    01, 03, 04, 05, 12, 23, 26, 34, 36, 45, 46, 56;

all remaining pairs point from larger to smaller. In particular 0→1→2→3→4→5→6→0 is a directed Hamilton cycle. Define Z and T_g as above. The normalized score of a tournament is the lexicographically smaller of its sorted outdegree sequence and the sorted outdegree sequence of its converse. These eight normalized scores are distinct, so even isomorphism plus converse cannot identify the center or any two leaves:

    Z    (1,2,3,3,4,4,4)
    T_0  (2,2,3,3,3,4,4)
    T_1  (0,2,3,4,4,4,4)
    T_2  (1,2,2,4,4,4,4)
    T_3  (1,2,2,3,4,4,5)
    T_4  (1,2,3,3,3,4,5)
    T_5  (1,3,3,3,3,4,4)
    T_6  (1,1,3,3,4,4,5)

All seven exhibited spokes are black by the reflection argument (and direct finite check). The black-bundle theorem implies there is another distinct black line joining the same merged pair for each spoke. We select one particular spoke per g, rather than incorrectly treating a bundled pair as one line.

As an independent finite control, endpoint subset-DP gives H(Z)=103 and H(T_g)=129,37,75,59,73,111,49. The eight counts also happen to be different; this is not used for the proof or claimed as a general property.

## Priority program

Use the seven leaf classes as actual machine states, with Z an intermediate phase. Provide the following ordered list of guarded macro-rules:

    g in {0,1,2} -> 3
    g in {2,3,4} -> 5
    g in {5}     -> 6
    g in {6}     -> 0

The first applicable rule wins, so g=2 goes to 3 rather than 5. At state g, traverse the selected spoke labelled P_g to Z, then the selected spoke labelled P_f to T_f, where f is the winning rule's target. The resulting deterministic map is

    0,1,2 -> 3 -> 5 -> 6 -> 0,    4 -> 5.

Every arrow denotes two real line traversals. The controller must remember its intermediate phase and desired destination when visiting the common hub. This program is finite-state with seven states and has no integer multiplication, divisibility semantics, unbounded registers, or universality. It is analogous to FRACTRAN solely in having an externally given ordered list and first-applicable semantics; it does not realize FRACTRAN as such.

## Reproducibility

`python3 check.py` and `python3 -O check.py` emit byte-identical output. The checker builds all tournaments explicitly, checks path validity and both K operations, colors, converse-normalized degrees, actual support of every macro, priorities, and separately computes Hamilton path counts by endpoint subset-DP and direct permutation enumeration. The proof of (1) and the reflection obstruction are general, not consequences of finite enumeration.
