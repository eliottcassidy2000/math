# Independent audit: shared network counts and the small-cell bound

23 October 2026. Auditor: five-cell subtask, separate from the parent
network and seven-cell derivations. **PASS in the declared scopes.** No
claim about A000568 or unnamed-cell quotients is licensed by these counts.

## Network Burnside formula

For a connected named support H with v cells, prime p, internal group C_p
acting regularly at every cell, and fixed edge degrees d_e in1..p-1, write

    A = product_e A_(p,d_e),    B = product_e binom(p,d_e),

where A_(p,d) is the number of labelled p-by-p d-regular binary matrices.
The ordinary orbit count is

    N_H = (A+(p-1)^v B)/p^v.

Audit derivation: for a translation vector a in F_p^v, an edge with exactly
one endpoint translation zero fixes no nonuniform regular block. By
connectedness a nonzero fixed vector must have EVERY coordinate nonzero.
For each of these (p-1)^v vectors, the edge action has p orbits of length p,
each transversal to all rows and columns. A degree-d_e fixed block chooses
exactly d_e such orbits, giving binom(p,d_e) choices independently across
the actual named edges. The zero vector fixes A configurations. Burnside
gives the formula. This uses prime order and connected nonuniform support;
it does not apply unchanged to composite p or disconnected support.

The restriction of a whole-network stabilizer to any root C_p is
injective by the nonuniform block kernel lemma. Hence its order is1 orp.
Counting pairs (configuration, nonidentity stabilizer element) gives
(p-1)^v B. Each symmetric configuration occurs p-1 times, so their number
is (p-1)^(v-1) B. Dividing by the orbit size p^(v-1) yields

    S_H = (p-1)^(v-1) B / p^(v-1),
    R_H = (A-(p-1)^(v-1) B)/p^v.

These count symmetry types within the specified attachment family; they
do not themselves construct a marking scheme for unknown cells. For two
five-cells and degree1, the independent census is N=8,S=4,R=4; for degree2
it is N=88,S=8,R=80. Both agree with the full five-cell block enumeration.

The parent's direct-orbit implementation was inspected: each generator
uses the SAME chosen vertex rotation on all incident edge coordinates,
with inverse domain precomposition at a row endpoint and image translation
at a column endpoint. It does not quotient adjacent attachments separately.
Its orbit-mass and orbit-stabilizer checks are appropriate independent
tests of the formula.

`audit_networks.py` adds a separate path. It does NOT import the parent's
orbit traversal or formula. It enumerates actual fixed binary matrices
for each local translation, then sums their products over all shared
translation vectors. It reproduces the complete matching-network counts
4,16,32,128,13888 and the mixed degree1/2/1 triangle count. The receipt
stores the entire local fixed-count tables and exact Burnside sums.

## Fixed edge-orbits still leave joint attachment information

The parent's further formula for fixed local edge types also passes an
independent orbit-mass derivation. Let beta=e-v+1, let r be the number of
rigid edge orbits, and put a=1 precisely when r=0 and every cycle gain is
balanced; otherwise put a=0. Each rigid edge has p^2 labelled realizations,
and each symmetric edge has p. Fixing all edge orbits therefore gives
p^(e+r) labelled networks. Their global stabilizers ALL have the same
order p^a: this depends only on the presence of rigid edges and the slopes
already fixed by the local edge orbits, not on remaining translation
phases. Dividing by the global orbit size p^(v-a) yields

    joint isomorphism classes = p^(beta+r-1+a).

This counts actual shared vertex rotations. In particular a tree with
two rigid attachments already has p different joint attachment classes,
even though it has no cycle. Each block's own canonical phase must agree
at their shared cell. Thus cycle-space information by itself is not a
complete loss ledger when independently rigid attachments are glued.
The formula is for connected named support and fixed full local edge
orbits, not a signature consisting only of degrees or scalar correlations.

## Nine-cell hostile

The parent's R9=C3[C3,C3,C3] has automorphism group C3^3 semidirect C3,
of order81. The block I3 tensor J3 allows independent internal rotations
on its six triples, but couples the two outer rotations, giving
3^6*3=2187 compatible pairs and a27-element target kernel when a full
source chart is fixed. This matches the complete direct backtracking
enumeration in check_networks.py.

One actual owner in EACH of the six triples is necessary: with fewer
than six, an untouched triple admits its independent rotation. Six are
sufficient: they eliminate the inner rotations and the outer rotations.
Thus the six-owner lower bound is global, not merely a check that the
particular displayed six-point base is irredundant. Naming all nine
vertices on one side still leaves the other side as three C3s. This is
an actual regular nonuniform attachment and a valid sharp boundary for
the arbitrary-odd-cell extrapolation, not a counterexample to the stated
prime-rotation theorem (its internal group is not C9).

## Seven-cell proof audit

Read [seven/PROOF.md](../seven/PROOF.md), independently checking:

1. Automorphism cycle shapes on regular seven-vertex tournaments. Odd
   order is mandatory, the5+1+1 shape violates fixed-vertex degree3, and
   the3+1+1+1+1 shape would require four-vertex internal degrees only0or3.
   Only7 and1+3+3 remain.
2. In a biregular a-by-seven block invariant under1+3+3, the entry at its
   unique fixed column is k modulo3 for every row. It is identically0or1
   (or no invariant row exists), and equal column sums force the whole
   block uniform. Both row and column equitability are used.
3. A proper nondiscrete equitable partition of a regular-seven tournament
   must be1+3+3. The two cyclic triples have matching/complement-matching
   cross block and uniform singleton incidence, so an actual coordinated
   C3 automorphism preserves every partition cell.
4. With a neighboring cell fully named, every target colour cell has a
   common literal column word. A surviving1+3+3 partition would therefore
   supply a forbidden block-kernel automorphism. Discreteness propagates.

These arguments prove two owners per connected seven-cell support
component at every total size. The proof correctly distinguishes this
ordinary-refinement statement from merely having a trivial global
automorphism group; the order31 rigid-but-nondiscrete examples remain
valid hostiles to that false implication.

The sharper one-owner choices depend on the explicitly enumerated
regular-seven isotypes and their point stabilizers. They are scoped as
FINITE-EXACT and do not infer arbitrary seven-by-seven block coverage
from the sampled attachment controls. The all-size kernel and propagation
proofs are independent of that sampled bank.

## Reproduction

    python -B 05-knowledge/results/odd_cell_completion_20261023/five/audit_networks.py
    python -O -B 05-knowledge/results/odd_cell_completion_20261023/five/audit_networks.py

The mathematical audit has no external literature dependency. The scripts
and receipts are exact integer controls, not floating-point evidence.
