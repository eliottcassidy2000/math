# Independent audit: square-sum ancestry and coloring observations

PASS, with exact scope. This is an independent proof review plus executable
replay, not a new all-size Hamiltonian census or a proof of an island-to-Berggren
correspondence.

## Structural proofs

1. The birth parent formula is exact for positive distinct endpoints. Its only
roots are 1,2,4. Direct ancestor recurrence independently confirms the base
71..99, first triple ancestry22, and last exception70 with ancestors1,4.
For n>=100 the largest square below2n yields a parent p<n and
p>n-2sqrt(2n-1)>71. The latter inequality follows from a positive increasing
quadratic beyond100. This proves the infinite tail and its minimal threshold;
testing through1000 in the independent checker is additional evidence, not
its substitute. The final complete pure-founder islands are also verified
exactly against the recurrence; the infinite tail excludes any later additions.

2. A parent-selector forest is acyclic by strict descent. A hard coloring is
induced precisely when each nonroot has a parent of its own color. The product
of counts of compatible parents is its exact multiplicity. Independent dynamic
aggregation of distinct coloring states yields192 colorings of total weight1024
through25, rather than reproducing the producer's enumeration of all selectors.
The impossible joint event c13=1,c23=4 is verified. The harmonic marginal law
is valid because a vertex's fresh parent choice is independent of all earlier
choices. The statement that the induced colors are independent would be false;
the packet correctly rejects it.

3. The equality-constraint rank proof gives collision probability q^(-r) exactly.
The independent checker computes connected components rather than union-find
for every same-endpoint pair in the supplied Q25..32 path sets. Every recorded
rank histogram and rational union bound agrees. Inclusion-exclusion uses the
rank of the UNION of constraints, correctly retaining shared objects. A separate
tiny independent full-coloring test checks that formula; the producer's stronger
K5 experiment is also replayed in normal and optimized modes.

4. The complete-graph warning is sound: endpoint-fixed vertex paths require
all internal colors distinct, giving q^2(q)_(n-2) for n>=2. Edge words obey
q^(n-1)>=(n-2)! as a necessary counting condition. These rule out a fixed palette
for unrestricted complete-graph families, not for the square-sum family without
an additional path-count lower bound.

5. The first proper two-color failure at15 and first triangle at30 were checked
independently. The three chosen islands have induced edge counts3,4,5, but
vertex counts4,5,6. An all-incidence isomorphism with a single-root Berggren tree
is excluded by roots and merging; the obstruction does not exclude more subtle
quotients, embeddings or separately defined arithmetic maps. The stated
Fibonacci rays and heap-address distinctions are explicitly inherited from
THM-3339/3382, not rederived by this audit.

## Procedural experiment and provenance

Normal and optimized full replays agree with every producer receipt. The
3104 raw formulas,1584 distinct colorings and1549 Q32 decoders refer only to the
specified finite feature grammar. Endpoint information and orientation are
retained. The Q36/Q80 controls are only two supplied pairs, not full test sets.
The high pass rate is correctly interpreted as weak evidence of a unique
arithmetic principle, not proof that these colorings work at all sizes.

The source receipt's SHA256 was verified directly from pinned Git commit
5258d52ac95924d524ee816a7dffa7329ddcf095, and every extracted path list was
compared exactly to that source. Every path is also independently checked for
coverage and square adjacency. Completeness of those path lists is inherited
from Session40 rather than rerun here.

A bibliographic author attribution was corrected during review to Janičková
and Csókási for arXiv2304.05230. No mathematical source was altered by the auditor.
The final source pins in receipt.json identify the reviewed bytes.

## Audit reproduction

The independent checker requires the sibling square-colors directory and the
pinned Git object for source provenance. Its path discovery supports both the local session layout and the published
sibling square-colors layout, locating the containing Git checkout. Replayed producer outputs and source hashes are retained.
No Git mutation was performed by this audit.
