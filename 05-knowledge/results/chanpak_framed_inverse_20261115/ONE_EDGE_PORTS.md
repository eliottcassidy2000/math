# Ordered placement of a single tile: exact mixed-base port loads

This supplements the
[one-edge frame example](FRAMED_SIMPLE_GRAPH.md) by allowing the
root itself to move. All graphs below have the same *one-edge*
backward tile shape at a fixed order, but there are different
path loads from their different ports.

Take a strong remainder W on 0,...,j, all arcs forward except
j→0, where j>=3. It has H=1+2^(j−1) by counting odd cycles
through the backward arc. For an interior root k put d=j−k−1.
Its weighted even path count is exactly
```
E_k = 2^d + 2^(k−2)(2*3^d−1)   (2<=k<=j−1);
E_1 = 2^(j−2) + 3^(j−2)−1;
E_0 = 2^(j−1)−1;
E_j = 2^(j−2).
```
Each root is admissible: W is strong. These are
actual H+2E outputs after the new outlet is added.

For an interior root, first take increasing paths. An odd
selection of vertices greater than k gives the required
even vertex order, so the baseline contribution is 2^d.
When j is omitted the complement retains the backward arc.
A selected odd subset S of the d higher interior vertices
then contributes an extra 2^(j−2−|S|) if any interior vertex
remains in the complement. For k>=2 this is automatic (vertex
1 remains), and the binomial sum is 2^(k−2)(3^d−1).
For k=1 subtract the exceptional all-selected term when d is odd;
the weight for the remaining reversed pair was already included
in the baseline.

If the actual rooted path takes the arc j→0, it has the three
mandatory vertices k,j,0. Each of the k−1 lower interior
vertices is optional after 0, and each of the d higher interior
vertices can occur before j, after 0, or nowhere. Each complement
is transitive. The parity condition selects half the assignments
for k>=2 by toggling one lower vertex. For k=1 the odd-selection
count is (3^d−(-1)^d)/2, giving the stated sum when combined
with the increasing family. At root 0 the only excluded even
subset in natural order is {0,j}; at root j every contributing
path starts j→0 and then selects an even subset of the interior.
These are specifically *labelled* boundary responses, not
undirected path counts. They involve 2 and 3 because the two
subpaths give higher interior vertices three actual choices.

Only O(log² X) distinct outputs can occur through X in this
strong one-tile lane: the guard count already bounds j by
O(log X), and there are j+1 ports. It is useful as an explicitly
solvable inverse family, not universal coverage.

The failure of unordered graph counts is not restricted to
graphs with zero Hamilton paths. Take the two labelled P4 graphs
on 0,1,2,3 with edge sets {01,12,23} and {02,12,13} as backward
tile graphs. Each has one undirected Hamilton path up to reversal.
The associated five-vertex remainders are strong (cycles
0→3→1→4→2→0 and 0→4→1→2→3→0), have H=11 and 9,
and with root 0 give H+2E=37 and 31, respectively. The only
change is the placement of the same unlabelled P4 geometry.
