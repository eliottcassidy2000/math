# Ramsey deficit budgets and native port compatibility

**PROVED reduction, using the CITED completeness of the local graph banks:**
R(3,4,4)<=73 is equivalent to exclusion of646 specified rooted cases below.
**FINITE-EXACT:** every supplied local graph passes its forbidden-clique tests;
an explicit aggregate-inventory model survives the first count relaxation.
**OPEN:** all four requested upper targets73,137,86,198. No new numerical
upper bound is claimed. The fourth tuple is R(3,3,6), as clarified by the user.

## Source check and inheritance

The current primary inputs include
[Boza, arXiv:2603.10851](https://arxiv.org/html/2603.10851v1):
R344<=77, R345<=157, R336<=91 and R444<=229, together with the standard
small two-color values and R334=30. The proof of
[R334=30 by Codish–Frank–Itzhakov–Miller](https://eprints.gla.ac.uk/119080/)
uses local abstractions followed by compatible edge-coloring completion.
The [McKay graph bank](https://users.cecs.anu.edu.au/~bdm/data/ramsey.html)
supplies the complete R44 graph lists on14 through17 vertices. Completeness
is CITED; this packet independently checks each supplied member's predicate,
not a fresh enumeration of every graph on those orders.

The closest repo mechanism is
[ramsey_saturated_ports_20261030/PROOFS.md, P1–P5](../ramsey_saturated_ports_20261030/PROOFS.md).
Its canonical hostile distinguishes exact row sums from column averages on a
27-vertex coloring. Here we retain each edge as one shared variable. Aggregate
counts are used only as necessary conditions. The least-used sidecar is the
complete small R44 port bank, rather than another mean-degree estimate.

## R1. Exact degree-deficit budgets

Let c_i be one less than a proved upper bound for the Ramsey number with
color i's forbidden clique size decreased by one. Every vertex v of a
counterexample has d_i(v)<=c_i. Set a_i(v)=c_i-d_i(v). Then

    a_i(v)>=0, sum_i a_i(v)=sum_i c_i-(n-1).             (5)

These are exact vertexwise equalities; their sum is not a substitute for
their shared-edge realization. The present target budgets are:

| target counterexample order n | color caps | deficit per vertex |
|---|---|---|
| R344,73 | (17,29,29) | 3 |
| R345,137 | (24,56,76) | 20 |
| R336,86 | (17,17,56) | 5 |
| R444,198 | (76,76,76) | 31 |

If R344<=73 were proved first, the R345 deficit becomes16 and the R444
deficit becomes19. It would still not establish either requested target.
Even combining R344<=73 with the hypothetical optimal lower endpoint
R335=45 only makes the basic R345 recurrence142, before any small modular
correction. A joint-neighborhood argument is needed for137. Likewise the
basic R444 recurrence with73 gives218, rather than198.

For R336, the stronger auxiliary assertion R335<=51 would already imply
R336<=85 by the recurrence and its parity improvement. We do not assume or
prove that auxiliary assertion. The86 target is pursued directly via its
five-unit deficit budget.

## R2. A646-case reduction for R344<=73

In a hypothetical73-vertex (3,4,4)-coloring, call the triangle-free color red.
Each red degree is at least14 by (5). They cannot all be14: then every blue
and green degree would be29, impossible for a graph on73 vertices by the
handshake lemma. Choose a vertex of **maximum red degree** r. Thus
r belongs to{15,16,17}, and every global red degree is at most r.

The red neighborhood has no red edges. Its blue/green graph is an R44 graph
on r vertices. There are640,2,1 supplied isomorphism classes at r=15,16,17.
After a global exchange of blue and green, the root's degree triple is one
of

    (15,28,29), (16,27,29), (16,28,28),
    (17,26,29), (17,27,28).                              (6)

Therefore there are640+2*2+1*2=646 root cases. In each, fix the root and
its red-port graph; keep all other edge colors as actual shared variables;
impose global red cap r and blue/green caps29; and forbid every red triangle
and every blue/green K4. The other two root neighborhoods are kept as labelled
ports with the displayed cardinalities. They are not replaced by inventories.

Every actual counterexample can be relabelled into one of these cases.
Conversely any satisfying complete coloring is a counterexample. Thus
excluding all cases is sufficient AND necessary for this target. The
maximum-degree selection and global color swap change the presentation, not
the coloring or its source incidences.

`ramsey_root_cases.json` records all646 native cases. `ramsey_roots.py` gives
a bounded SAT probe using one color variable per actual edge, native degree
constraints, all red triangle clauses, and lazy blue/green K4 clauses.
Every lazy clause is a genuine forbidden clique. An interrupted or conflict-
limited run is UNKNOWN; it excludes no case. A solver-only UNSAT result is
not promoted without a separate proof certificate. This is a concrete search
interface, not a646-case proof hidden behind a count.

## R3. Why just the next inventory is still too weak

The complete supplied bank has edge ranges41..50 at order14,50..55 at15,
exactly60 at16 and68 at17. All16-vertex graphs have eight vertices of each
degree7 and8; the17-vertex graph is8-regular. The actual graph identities
remain available when these profiles coincide.

Here is a consistent **aggregate relaxation**, not a coloring. Give all73
vertices color degrees(16,28,28), and assign the internal edge-count rows
of their red, blue and green neighborhoods to be

    (0,60,60), (96,102,180), (96,180,102).               (7)

Every row sums to its port's binomial edge count. The red row is realized by
the genuine16-vertex local graphs. The other rows satisfy the R334 internal
degree caps, and their diagonal triangle inventories are divisible by3.
They even admit integer internal degree histograms: in the blue port take
20*(7,7,13)+4*(6,8,13)+4*(7,8,12), then swap blue/green for the green port.

Globally these data give4,964 monochromatic triangles,49,056 two-color
triangles and8,176 rainbow triangles, summing to binomial(73,3). The mixed
wedge identity2*T_two+3*T_rainbow=sum_v sum_(i<j)d_i(v)d_j(v) holds exactly.
Thus parity, these degree caps, local edge totals and the global triangle
identity do not exclude73. This witness says nothing about whether all its
putative blue/green ports are realizable or glue together. Those are the
precisely missing constraints. A scalar relaxation cannot recover them from
the same totals without a new theorem.

## R4. Audit of the originally written R366 target

Before the correction to R336, the literal R366<=86 was impossible. A
self-contained check supplies an explicit202-vertex obstruction: take two
copies of the Paley graph on F101, color residues blue and nonresidues green
inside each copy, and all edges between copies red. Red is bipartite; blue
and green cliques remain in one copy. Both internal colors have36,360 K5s
and zero K6s, checked by recursive clique counting and independently by
testing every edge's common neighborhood for K4. Multiplication by a
nonresidue exchanges the two internal colors. Hence R366>=203. This is a
classical-type construction and a typed-target correction, not a new record;
the current survey already records a stronger lower bound338.

## Next obligations

1. Exclude actual shared-edge extensions in the646 cases, beginning with
   the three high-order port graphs. Degree-profile equality is not enough.
2. For R336, retain the17/17/56 degree-budget structure and apply the
   existing H13/R36 port machinery to vertices with small deficits; prove
   compatibility constraints rather than assuming saturation at order86.
3. For R345 and R444, seek local certificate inequalities that pay their
   respective20- and31-unit deficits. First closing R344 would reduce these
   to16 and19, still nonzero and still requiring a proof.

The common lesson with Collatz and LRC is exact: compatible local objects
must use the same source edge, multiplier lift or time label. It is not an
asserted equivalence between the three problems.
