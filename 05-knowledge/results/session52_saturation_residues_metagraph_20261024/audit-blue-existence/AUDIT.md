# Independent audit: internal blue existence translation

24 October 2026. Verdict: PASS for the proved reduction and finite controls.
All-order existence remains CONDITIONAL on the uninspected primary theorem.
No mathematical defect was found. The producer added a proper closed manifest
and pinned local context during this audit; the frozen copy is under `producer/`.

## Source and quantifiers

Independently inspected Farrugia (1999), Section 2.30, printed page 59,
PDF page 75, at https://alastairfarrugia.net/sc-graph/sc-graph-survey.pdf .
The statement explicitly fixes an order-two antimorphism of a regular tournament
of odd order, and covers every self-converse arc and every cycle length except
possibly four and either three or five. Section 2.29 defines self-converse
cycles via reversal by that same antimorphism. Thus the producer has correctly
preserved the specified-involution and specified-arc quantifiers, not merely
existence for some involution. The bibliography identifies Salvi-Zagaglia's
1980 paper as reference 346. This verifies the secondary report only. The
primary paper was not inspected in this audit. N>=7 avoids all stated
exceptions; the short Hamiltonian-corollary sentence is not needed.

## Independent mathematical review

For n=2p, the fixed-point-free condition follows because two fixed vertices
would give an arc mapped to its reverse, and an even involution cannot have
exactly one fixed vertex. Score reversal exchanges the p high and p low
vertices. Adding z with z->high and low->z is forced by regularity, gives
outdegree p everywhere, and extends the given q uniquely while fixing z.

A q-self-converse Hamilton cycle rooted at z must have its ith vertex paired
with its (2p+1-i)th vertex. Its z-successor is high and z-predecessor is low.
Deleting z gives precisely a q-symmetric high-to-low Hamilton path; adding z
is the inverse. The middle arc is the unique arc reversed onto itself by q.
Every pair {v,q(v)} supplies one such arc, and no arc incident with z is fixed.
Hence the reported every-arc theorem implies at least p distinct labelled
paths for that fixed q, not p different quotient words.

For a path word, complementing only nonpath arcs sends its degree to
n-1-d(v)+[v=start]-[v=end]. In the balanced layer the complemented tournament
is balanced exactly when start is high and end is low. Reflection symmetry
of the path is exactly the blue condition. This justifies the internal-blue
translation without an extra structural assumption.

The group argument is sound: tournament automorphism groups have no element
of order two, hence have odd order. If g reverses arcs, g squared has odd
order m; g to the mth power reverses arcs and squares to the identity.
Therefore every self-converse tournament in the stated even-order range has
an anti-involution. This proof is independent of the unrelated classification
claims in the pinned contextual theorem.

## Independent finite verification

`audit.py` enumerates the entire labelled tournament cube, not the producer's
reflection-orbit parametrization. It tests q anti-invariance directly and
uses every vertex permutation, not the producer's signed half-path formula.
At n=2,4,6 the cubes have 2,64,32768 matrices; 2,16,512 pass anti-invariance;
2,8,112 are balanced. Minimum q-symmetric downhill path counts are 1,2,5.
Every paired arc occurs centrally. Aggregate path/cycle counts are 2,16,576.

For every retained matrix the audit checks regularity and anti-invariance
of the extension, literal equality of the path and rooted-cycle sets, and
uniqueness of the fixed middle arc. It also checks the internal endpoint
criterion on every Hamilton path, not merely the symmetric ones. Normal
and optimized execution agree byte-for-byte and use explicit exceptions.
These controls close n=2,4 without relying on the reported exceptional lengths.
No finite test is treated as proof at larger order.

No connectivity, bridge, articulation, minimum-17, or all-order counting
conclusion follows. Regular completion and blue-path conventions are inherited;
no new theorem priority is claimed. Stop at exact translation, source lead,
and finite controls; obtain the primary source before upgrading source status.
