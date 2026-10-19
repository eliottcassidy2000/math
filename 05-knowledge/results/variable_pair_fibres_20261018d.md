# Variable ternary labels, pinned shade fibres, and shared source cuts

**2026-10-18. PROVED:** the exact fibre criterion and finite lazy refinement
theorem below. **FINITE-EXACT:** all 570 saved colourings of the old
714-point graph extend to the new 1,482-point graph while every old colour
is fixed. **OPEN:** an equal-colour geometric fragment with no separating
six-colouring, and the plane six/seven question. A bounded unfinished search
is `pending`, never a negative certificate.

## 1. Inheritance and the change of state

The closest proved mechanisms are the
[signed boundary parity compiler](plane_six_reframes_parity_20261011b.md),
[exact shade fibres](snark_boundary_transfer_20261018a.md), and
[joint residue support compiler](coupled_support_branches_20261018c.md).
The last packet's fixed-anchor triangle defeats every inherited two-colour
pair chart, yet admits unrestricted colourings. That is the canonical
hostile. The corrected near miss is to interpret failure of those charts as
failure of all six-colourings. The least-used sidecar is the parity condition
between two *pinned* endpoints of an actual monochromatic path.

The live board is: physical vertices, freely varying labels, source shades,
actual cycle/path provenance, and complete rather than bounded rejection.
The source graph is unchanged. The representation is

\[
\ell_v=c_v\bmod3,\qquad b_v=c_v\bmod2,\qquad
c_v=4\ell_v+3b_v\pmod6.                              \tag{1}
\]

Every vertex now has an independent ternary label. It is not required to
inherit a pair label from an anchor or layer. Formula (1) uses the specific
palette pairs \(\{0,3\},\{1,4\},\{2,5\}\); the three labels are arithmetic
residues, not the three nonzero vectors of an \(\mathbb F_2^2\) flow. The map
preserves all six-colourings because it is a bijection on the six colours.
Interior shades may be eliminated only with their complete fibre criterion.

## 2. Exact fixed-label fibre, including pinned paths

Fix a finite simple graph \(G=(V,E)\), marked vertices \(B\), and literal
source colours \(\phi:B\to\{0,\ldots,5\}\). A candidate label assignment
must satisfy \(\ell_v=\phi(v)\bmod3\) on \(B\). Let \(H_\ell\) retain
precisely the original edges whose endpoints have equal labels. On each
such edge properness is exactly \(b_u\mathbin\oplus b_v=1\); other edges
need no shade condition.

**Fibre theorem (PROVED).** The labels lift to a proper six-colouring
extending \(\phi\) iff:

1. \(H_\ell\) is bipartite; and
2. for every pair of marked vertices in one component, their prescribed
   shade XOR equals the parity of a connecting path.

If these hold, there are exactly \(2^f\) lifts, where \(f\) counts components
of \(H_\ell\) without a marked vertex. Isolated vertices count as components.
A spanning forest supplies the relative shades. Each marked component's
flip is fixed by its first source shade, and the second condition makes all
other pins agree. Every unmarked component has two independent flips. This
also proves necessity and the count.

`shade_fibre` returns either a full colouring with \(f\), actual simple odd
cycles, or actual simple paths with incompatible pinned shades. All paths
and cycles retain physical vertex IDs and are authenticated edge by edge.
Source labels are checked before looking for a shade lift.

## 3. Globally valid lazy cuts

Use one-hot variables \(X_{v,t}\) for \(\ell_v=t\). Retain a Boolean shade
variable only for each marked vertex; under a source pattern it is pinned
to \(\phi(v)\bmod2\).

An actual odd cycle \(C\) yields the three clauses

\[
\bigvee_{v\in C}\neg X_{v,t},\qquad t=0,1,2.          \tag{2}
\]

An actual simple path \(P\) between marked vertices \(u,v\), of length
\(L\), yields the conditional constraint

\[
\bigl[\ell\text{ is constant on }P\bigr]
\ \Longrightarrow\ b_u\mathbin\oplus b_v=L\pmod2.    \tag{3}
\]

For each label \(t\) and each of the two incompatible endpoint-bit pairs,
the code exports one CNF clause: either a path vertex is not labelled \(t\),
or that endpoint-bit pair is not used. Thus (3) has six clauses. Retaining
the endpoint variables is essential: the same path can be inconsistent for
one source pattern and consistent for another. Both (2) and (3) are valid
for *every* full six-colouring of the original graph, not just the assignment
that discovered them. They can therefore be shared across source patterns.

**Complete cut criterion (PROVED).** With the prescribed source labels and
shades, one-hot labels satisfy all odd-cycle cuts (2) and all marked-path
cuts (3) iff the full source colouring extends. One direction follows from
properness. Conversely, a nonbipartite \(H_\ell\) supplies an odd cycle
violating (2); if it is bipartite but its pins disagree, its spanning forest
supplies a path violating (3). Otherwise the fibre theorem constructs the
extension.

This gives a lazy algorithm. Ask for labels satisfying the retained cuts;
if their shade fibre fails, add its actual obstructions and repeat. Every
failed assignment is removed, while no extendible assignment is removed.
For fixed source colours at most \(3^{|V\setminus B|}\) label assignments can
fail, so the mathematical process terminates. No uniform efficient complexity
bound is claimed.

The implementation has explicit iteration and negative-proof budgets.
A positive output contains an independently checked full colouring. A
negative output contains all source-indexed cuts and a complete Boolean
decision tree for the remaining label CNF with the actual source assumptions.
The inherited independent unit-propagation tree checker verifies both
branches and every conflict leaf. If either budget expires before the
required certificate, the result is `pending` and asserts no nonexistence.
The SAT solver is a proposer, not an unverified negative oracle.

## 4. Why the obstruction must be found after choosing labels

In \(K_4\), consider the fixed fundamental cycle basis
\(0120,0130,0230\). Labels \((1,0,0,0)\) make each basis cycle nonmonochromatic,
but the remaining triangle \(1231\) is monochromatic and has no shade lift.
Thus a fixed basis of the original cycle space is not a complete cut bank.
A spanning forest recomputed inside the current \(H_\ell\) does detect the
obstruction; no issue is asserted with that dynamic procedure.

The path \(0-1-2\) gives the complementary hostile. With all labels zero,
source colours \((0,3)\) at endpoints have the wrong shade parity, whereas
\((0,0)\) have the right parity. Sharing an unconditional rejection from
the first source would destroy the second valid source. The conditional
path cut retains both. Allowing the middle label to change also repairs
the first source, as the variable-label solver explicitly demonstrates.

An abstract \(K_7-e\), with endpoints of the missing edge marked, provides
the equality positive/negative control. Equal pins have a verified colouring;
separating pins are rejected after 16 label iterations and 30 actual
odd-cycle cuts, with an independently checkable 23-node negative tree.
No unit embedding is asserted for this control.

## 5. Structured old witnesses and literal extension to the cubic candidate

The old graph and its complete eight-port word bank are supplied by
[coupled support branches](coupled_support_branches_20261018c.md). Three
selected patterns, at canonical indices 0,153,569, are newly solved by the
variable-label procedure on the retained core. They take 10,1,27 iterations,
sharing the cuts from earlier patterns. After restoring all physical vertices,
their fixed full-label fibres have respectively 132,110,118 unpinned
components, and therefore exactly \(2^{132},2^{110},2^{118}\) shade lifts.

In these particular witnesses, 456,479,509 physical vertices have an address
whose label differs from the original anchor's arithmetic label. The
corresponding address counts are 505,525,560. These are measured migrations
in the exported witnesses, not minimal distances over all colourings or
palette relabellings. They explicitly leave the failed anchor-only model.

The new separately authenticated input is
[surviving repair geometry](surviving_repair_geometry_20261018d.md):
1,482 points and all 6,118 actual unit edges. The old physical vertices
remain literally 0 through 713. For each of the 570 saved old full colour
words, **all 714 old colours are pinned**, rather than only the eight original
ports. Relative degree-at-most-five peeling leaves 1,205 vertices and
4,944 edges. This ordinary-colouring support reduction is audited once;
the learned paths/cycles retain actual vertex IDs when lifted back.

**FINITE-EXACT:** every saved old word extends. The shared scan uses 933
label proposals in total, at most 32 for one word; 475 words require no new
cuts. Its final bank has 2,606 odd-cycle cuts and 3,315 conditional pinned-path
cuts. Recorded cycles have lengths between 3 and 45, and paths between 4
and 41; these are observed receipt lengths, not minimal necessary cut lengths.
Each of the 570 full positive words is retained and checked on every actual
edge and every old pin. Their fixed full-label fibres have between 161 and
212 unpinned components, giving at least \(2^{161}\) extensions for each
of these particular old words. This is a lower bound supplied by these
fibres, not a count of all extensions.

The independent geometry analysis subsequently identifies a source-fixing
retraction from the new graph to the old graph: new translated points
\(r+t_s\) map to \(r\), while every old vertex is fixed. Its complete
edge/contact audit is in the geometry packet. That theorem proves that
**every** old colouring extends, a stronger statement than this finite
word bank. The lazy solver receipts are independent verification of the
570 selected instances; they neither establish nor replace the retraction.
The geometry, words, and source hashes are retained separately.

Hence all 123 equal and 447 separating eight-port patterns survive. The
literal pentagon terminals 2,4 remain nonunit-separated, so even a future
equality at those terminals would need an authenticated unit connector.

## 6. Exact next obligation and reproduction

The source-preserving obstruction target is now precise: authenticated new
geometry must make the full label CNF generated by valid cycle and pinned-path
cuts unsatisfiable for every separating source pattern, while one equal
pattern has a full positive lift. Merely forcing label changes on seed
triangles, failing a fixed palette chart, or accumulating many cuts is not
that target. A source-fixing retraction supplies the opposite conclusion.

The controls exhaust all 603 fixed-label fibres arising from all 75 simple
graphs on one through four vertices with the declared pins, compare their
counts with direct shade enumeration, and test learned-cut soundness against
full source colourings. They include the fixed-basis hostile, the two pinned
path sources, complete abstract equality certificates, exact type/source
hostiles, and the zero-budget `pending` boundary.

```text
python -B 04-computation/experiments/variable_pair_fibres_20261018d.py
python -B -O 04-computation/experiments/variable_pair_fibres_20261018d.py
```

The JSON stores the actual graph and source hashes, old-bank hash, all
learned physical walks, full colour words, exact fibre dimensions, and
complete negative-control trees. Normal and optimized runs both pass
**3,668 checks** and reproduce the saved JSON byte-for-byte after LF
normalization. The LF SHA256 is
`e754bac715c83b7503463ca3efa9ba2de2290c7f72de7df53e2ff1fa1982f4f8`.
**Independent proof/API and solver-free receipt audit: PASS.** The reviewer
rechecked all 570 full words against every edge and all 714 old pins, the
source/coordinate/old-bank hashes, all 5,921 actual cut walks, and each
fixed-label component count using an independent union-find implementation.
The relative core and peeling were replayed. A separately reconstructed CNF
and set-based propagation checker validated the complete 23-node negative
tree; the three structured old words and their fibre counts also passed.
