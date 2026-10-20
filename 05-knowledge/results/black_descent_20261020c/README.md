# All-order black descent, SCC basins, and cyclic-five Petersen reductions

20 October 2026, continuation c. Inherits baseline 7eea99e6eb.

**PROVED:** strict black triangle descent from every strong black-active
tournament, at every order; exact SCC terminal states and a uniform
four-black interchange; a cyclic-five-preserving pole substitution and
an infinite nonflower family reducing to Petersen.
**FINITE-EXACT:** explicit double-star and Loupekine certificates,
complete specified seven-port interface censuses, and the ranks and
spectra of the matching-incidence maps.
**OPEN:** black-core connectivity at arbitrary order, exclusion or
construction of actual non-leaf bridges beyond the checked range,
and an unavoidable admissible reduction in every cyclically
five-connected snark.

## 1. Cycle-changing descent is now an all-order theorem

For a directed Hamiltonian path P=(a,...,b), reverse its arcs. Only the
endpoint scores change, and

    c3(T^P)-c3(T) = s(a)-s(b)-1.

The [new proof](cycles.md) shows that every strong black-active T has a
**black closed** path with s(b)>=s(a). Thus the triangle count strictly
decreases. The earlier fixed-cycle failures survive; the proof changes
the cycle when necessary.

At odd orders the mechanism groups prime-order automorphism orbits into
consecutive blocks, then solves their boundary compatibility by a
sumset argument. At even orders the remaining hypothetical obstruction
has two opposite blue cuts. A terminal-block insertion using indegree
counts produces a black descending cut and excludes that obstruction.

The [SCC consequences](scc.md) are exact at all orders n>=3:

* The strict triangle minima of the full graph are exactly nonstrong
  tournaments.
* On the black-active subgraph, the minima are exactly its nonstrong
  nodes. Every black component reaches this terminal set.
* Contracting any descending forest gives a terminal-basin multigraph
  with the same connected components.
* Its terminal Laplacian has rank m_n minus the number of black
  components. The issue is now irreducibility on SCC boundary states,
  rather than the existence of a strong-node descent.

There is a concrete infinite family of basin connections. For any
tournament B, with A=1->1->B,

    A -> 1 -> C3  --four black path reversals-->  A -> C3 -> 1.

The proof retains all labels, restores every prefix arc, and supplies
the global paths. Removing one of the two leading singletons can turn
the last move blue; this guard cannot be discarded.

## 2. Petersen reductions inside cyclically five-connected snarks

The [snark packet](snarks.md) constructs an actual rooted reduction of
the 30-vertex double-star to Petersen. At an explicitly retained
seven-port boundary,

    |S_old|=192,   |S_new|=120,   S_new strictly contained in S_old.

This direction preserves uncolourability in every unchanged exterior.
Exact equality is impossible for every double-star seven-pole obtained
by deleting three vertices inducing one edge: all 1,080 interfaces have
larger support than a Petersen pole. The successful reduction retains
the exterior {0,1,3} and gives ten explicit edge deletions, twenty
suppressions, and a labelled Petersen isomorphism.

The seven Kempe coordinates add a stronger test than nonnegativity.
Although f_old>=4 f_new entrywise, the nonnegative remainder has a
positive Kempe defect. It is therefore not a physical cubic-pole tensor,
nor a nonnegative combination of such tensors. A physical remainder
would require the coefficient of f_new to be at most two.
This is a realizability obstruction occurring inside a genuine
cyclic-five reduction.

A second move uses the classical NN five-pole with

    f_NN = 3 f_P3.

Both Loupekine snarks reduce to Petersen by a certified rooted instance
of this substitution. A new cut lemma, combined with an exhaustive
fixed-gadget check, proves that NN substitution preserves cyclic
five-connectivity in every suitable host. Replacing a path in J_k gives
an infinite family F_k with

    |V(F_k)|=4k+12,   k odd and k>=5,

of cyclically five-connected nonflower snarks, each reducing first to
J_k and then to Petersen. Whole-graph checks cover 32,40,48,56 vertices
in addition to the all-order argument. The known NN colouring identity
and the cited flower connectivity theorem are explicitly credited.

## 3. What 105 and 17 now connect

The [inheritance and algebra packet](105_17.md) reviews the earlier
matching, Hamming, defect-template, Brown-configuration, Fano and
tribonacci work before constructing two precise bridges.

Marking one of eight ports identifies its 105 perfect matchings with
the size-three matchings of K7. The 105 seven-port junction closures
are the size-two matchings of K7. Their containment matrix R is
3-regular and has rank 91. Its coefficient kernel is exactly the
14-dimensional determinant kernel of the eight-port colouring probes.
This gives a genuine linear quotient correspondence.

The loss after composition matters: if A is the seven-port junction
measurement matrix, rank A=84 but rank RA=70. The complete blindspace
consists of the seven Kempe coordinates and fourteen additional
row-sum-zero edge-weight coordinates. Thus this natural 105-to-105 map
does not preserve the full obstruction test without additional data.

Separately, every Lagrangian spread in the eight-dimensional binary
symplectic space of ten-label even subsets has 17 members. The code
dictionary forces exactly seven matching codes and ten Hamming pieces.
Their distinguished pairs form a complementary 2-factor on ten labels.
Explicit examples yield C10 and C4+C6. The first gives colourable cubic
graphs when a matching is added, so the spread does not force Petersen.

The tribonacci values 17 and 105 also belong to the actual inherited
path-reversed tournament family. A recent inversion theorem studies
that same family with a different operation; the packet records the
precise identification and the operation mismatch.

## 4. The all-size classification ledger

The graph throughout is the **converse-merged complement-line
multigraph**, retaining loops, colours and line multiplicities.
Single-arc-reversal graphs have different laws.

| Aspect | Current all-size result | Remaining evaluation or structural obligation |
|---|---|---|
| Vertices | Ordinary tournament types modulo converse; exact Burnside/SCC descriptions | No small closed list of arbitrary strong types |
| Pure-blue locus | Classified SCC alphabet and palindrome structure; rational counting series | Already inherited, with the affine-only proof gap repaired |
| Degrees and D_n | Exact fibre histogram and finite all-order coefficient algorithm | A simple closed evaluation of extrema is unavailable |
| Maximum M_n | Exact extremal characterization; M_n/n!=2^(-n+O(log n)) | Leading factor and exact extremizers at arbitrary n |
| Strong/nonstrong incidence | Independent nonstrong periphery and exact colour-resolved cut masses | Does not determine individual neighbours |
| Black descent | Now PROVED at every order | None for existence of strict descent |
| Black components | Every component meets the explicitly counted SCC terminal set | Connectivity of its terminal-basin quotient |
| Bridges | No black bridge; degree-three pure-blue nodes attach triply to one black component; all additional bridges must be mixed-to-mixed blue lines | No all-order proof that genuine bridges must meet leaves |
| Abstract wiring | Exact sufficient-and-necessary involution realization conditions | These do not guarantee actual tournament fibres |
| Actual wiring | Literal compatible tournament/path data suffice; score and Walsh packets give necessary constraints | A useful compressed sufficiency theorem remains open |
| Snark reduction | Rooted support inclusion preserves obstruction; named cyclic-five and infinite-family constructions | Unavoidable reductions for arbitrary cyclic-five snarks |

For the degree polynomial convention retained from the
[earlier metagraph packet, Section 6](../boundary_support_20261020_metagraph.md),
put h(T)=H(T)/|Aut(T)|. Then

    D_n(z)=sum_(self-converse T) z^h(T)
           +(1/2)sum_(non-self-converse T) z^(2h(T)),
    M_n=deg D_n.

That packet gives every coefficient through n8 and a finite algorithm
at every n; it does not evaluate M_n in closed form. The
[asymptotic proof](../boundary_support_20261020_asymptotics.md)
retains the precise scope. The n9 connected-core and five-leaf-bridge
classification remains [FINITE-EXACT](../blue_quotient_20261020b/census.md).
The present theorem does not turn it into an all-order classification.

The concurrent [black-detour packet](../defect_transitions_20261020/blue.md)
was read and independently audited during integration. Its boundary
triple `(first,second,last)` gives an exact common labelled target for
two-step transport from nonstrong self-converse sources. This proves
the all-order triple attachment in the table and three-black detours
for the mixed-centre family A->C->A^op with |A|>=2.
Its R4->Z->R4^op hostile shows that prescribed endpoints can force
every available path blue even at a mixed node. These results identify
some basins that our descent theorem alone does not identify.
The [concurrent pole proof](../defect_transitions_20261020/petersen.md)
independently supplies the four-pole mechanism and forced partial Kempe
linkages; simultaneous disjointness remains its missing coordinate.
Our rooted NN carrier and cut certificate supply actual instances,
while the double-star strict inclusion goes beyond exact equivalence.

## 5. Pythagorean data that survive the reductions

For nonzero nonnegative boundary tensors f,g, put U=f tensor g and
V=g tensor f. Their norms agree; U+V and U-V are perpendicular, with
squared leg lengths

    2(M+C^2), 2(M-C^2),
    M=||f||^2||g||^2, C=<f,g>.

Disjoint support is exactly C=0, the isosceles case. The NN equality
scales the triangle and preserves its normalized shape. The double-star
support inclusion preserves the obstructed isosceles condition, although
it need not preserve other angles. This is the precise geometric
connection being used; primitive integer legs and complete/bipartite
graph families do not follow automatically.

The separate formula

    (q-1)^2+(2 sqrt(q))^2=(q+1)^2,  q=4^h,

connects the Euclid triples 3,4,5 and 15,8,17 with the q+1 members of
a trace-symplectic spread. Its explicit parameter link is useful, but
its limiting shape is flat, so it cannot replace the boundary-obstruction
angle.

## 6. Proofs, certificates and reproducibility

* [cycles.md](cycles.md), [cycles.py](cycles.py): all-order descent,
  additive block certificates and complete strong self-converse tests
  at n8 and n9.
* [scc.md](scc.md), [scc_normalization.py](scc_normalization.py):
  all-order minima, basin reduction, exact terminal count and 7,888
  literal contextual macro audits.
* [snarks.md](snarks.md),
  [snark checker](../../../04-computation/black_descent_20261020c_snarks.py):
  rooted minors, full boundary tensors, independent closed-graph counts,
  cyclic cut tests and fixed-gadget connectivity certificate.
* [105_17.md](105_17.md), [105_17_probe.py](105_17_probe.py):
  exact kernels, spectra, characters, Fano recoding and spread witnesses.
* [BOARD.md](BOARD.md): inheritance and the live research portfolio.

Each checker uses explicit failures that survive Python optimization.
The notes distinguish mathematical proofs, fixed finite gadget proofs,
complete finite universes, cited antecedents and open generalizations.
