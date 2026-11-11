# Independent audit of named gap transport and seven-state observer loss

**PASS.** The proofs in [GAP_TRANSPORT.md](../fibres/GAP_TRANSPORT.md) and
[ROOT_OBSERVER_LOSS.md](../fibres/ROOT_OBSERVER_LOSS.md) have the required
operation and orientation scopes. No correction was needed.

The four-vertex mixed block is an ordinal module relative to every outside
component throughout each local operation. Thus its H3,H5,H5,H3 values
multiply by the unchanged outside H, even though the intermediate block
is strong rather than a triangle plus singleton. The lower route first
makes the labelled triangle transitive; swapping the singleton across its
three adjacent vertices preserves transitivity at every step. Restoring the
one internal edge restores the original triangle, giving the identical
named endpoint. The high and low routes respectively preserve a lower
floor and an upper ceiling; they must not be conflated.

The normalization count is exactly the inversion count between triangles
and singletons. The audit deliberately crosses the rightmost available
inversion, unlike the producer's leftmost schedule, and recovers the same
named ordinal normal form for all56 stated finite gap vectors. Triangle
order and singleton order never change within their own classes. This
supports the general terminating inversion-count proof.

The complete21-neighbour audit of C3->K1->C3 gives minimum15 among
neighbours with H>=9. Every nontrivial arc-flip path respecting that floor
must take one of these first edges, so15 is the sharp peak under the floor.
The separately constructed path through H3 refutes any unqualified
peak15 assertion. With no floor, the endpoint H9 itself forces peak at
least9, and the low route attains9.

For the observer collision, the common root is a global sink. This makes
E=H, S=0, U=H(deletion), and V=K(deletion). Ordered two-bin covers multiply
over ordinal joins because each bin traverses the component blocks in their
fixed order. The local factors12 and2 are exact. Therefore the complete
seven-coordinate equality is valid for all specified r,g, independent of
gap placement. It is a genuine equality at the same declared root type.

The first-content distinction uses the already proved complete
triangle-chain C1 theorem. Independently of that polynomial implementation,
the audit here constructs all256 one-vertex extensions of EACH displayed
eight-vertex rooted graph and counts their paths by fresh subset DP. Their
gcds are respectively3 and1. Both full seven-states are independently
`(9,9,0,9,288,9,288)`.

The inherited O/Z/J and ordered root-amalgam formulas apply to weak rooted
tournaments, including sink roots, so there is no omitted strongness
hypothesis. The asserted indistinguishability is for common O/Z/J futures
and separately for each fixed two-sided root-amalgam context. No assertion
about arbitrary attachments, arbitrary interleavings with unclassified
operations, or universal metagraph state completeness follows.

Reproduce with `python -B audit_gap.py` and `python -O -B audit_gap.py`.
`gap_audit.json` stores all512 literal extension values, all21 neighbours,
all local flip orders, the alternative-schedule gap universe, larger
seven-state controls, and source SHA256 pins. This supplements the
analytic infinite-family proof rather than replacing it.
