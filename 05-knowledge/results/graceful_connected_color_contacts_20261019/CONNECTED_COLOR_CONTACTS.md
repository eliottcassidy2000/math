# Two-color repair from connected color-contact graphs

## General theorem

Let F be a fixed rainbow forest with r+1 components, and let r free colors be available in an edge-colored graph. A completion uses one edge of every free color and forms a spanning tree containing F. For each free color d, let Q_d be its contact graph on the contracted components of F, ignoring loops and retaining the original edges as realizations of contacts.

If all but at most one Q_d are connected, the completion fiber is connected under moves that change at most two free-color edges simultaneously. Its diameter is at most r. These moves can change the abstract tree shape.

### Proof by monotone target-edge installation

Fix source T and target T*. If there is a color whose contact graph may be disconnected, process that color first; process the remaining colors in any order. At each step, the goal is to install the target edge g of the current color e, preserving all target edges already present. Skip the step if g is already present.

Insert g into T, creating a fundamental cycle C. If the old color-e edge lies in C, remove it; this is a valid single-color move, retaining every target edge.

Otherwise choose an old edge f of C that is not an edge of T*. Such an edge exists because C cannot be contained in the target tree. It is not a frozen edge, and its color q differs from e. Its contact graph Q_q is connected: in the first step the only potentially bad color is e itself; later the bad color's target edge is already installed, so a non-target f cannot have that color.

Remove both the old color-e edge and f, then add g. This leaves a forest with two components. Indeed removing f from the fundamental cycle preserves connectedness after adding g, while the old color-e edge lies outside that cycle and its subsequent deletion creates exactly two components. All F edges remain, so this cut is a union of F components on each side. Connectedness of Q_q supplies an actual color-q edge crossing the cut. Add it. The result is a spanning tree with exactly the original free colors, and only colors e,q changed.

The removed color-e edge was not its target, and f was chosen outside T*, so all previously present target edges remain. Each nontrivial step installs a new target color edge. After at most r moves every free target edge is present, hence the result is T*. This proves the bound.

No vertex buffer is needed. In particular the earlier abstract leaf-buffer/contact theorem is a special case: its hypotheses make every full contact graph Q_d connected. The hostile-family bound improves from2k−1 to k, without asserting that k is an optimal diameter.

## Every upper-tail fiber with three free colors

In the integer-difference graph on0,…,m, freeze any acyclic selection of colors4,…,m and leave1,2,3 free. For every m≥3, the completion fiber is connected under at-most-two-color moves, with diameter at most3.

For m≥5, Q_1 is connected because all difference-one edges form a spanning path. Q_2 is connected because all difference-two edges form the two parity chains and F contains a color-five edge, which joins those chains. Therefore only Q_3 can be disconnected, and the general theorem applies with r=3.

The remaining cases are finite boundaries, not a growing census. For m=3, there is no frozen edge and exactly4 graceful trees; any two differ in at most two colors, giving diameter1. For m=4, the only frozen edge is{0,4}; the12 completions have two-color diameter2. The included exhaustive boundary check enumerates every allowed color-start tuple and verifies these exact statements.

The pre-existing two-free-color disconnected fibers also become connected under width two trivially, since a completion has only two free edges. The new argument treats three free colors uniformly for unbounded m; it does not settle arbitrary numbers of free colors.

## A concrete failure of the all-connected hypothesis

For k=3,m=7, take F={(0,4),(1,6),(0,6),(0,7)}, with colors4,5,6,7. Its components are H={0,1,4,6,7}, {2}, {3}, {5}. The color-three contacts are H–{3} and{2}–{5}, so Q_3 is disconnected. Q_1 and Q_2 are connected. This fiber is covered by the one-bad-color theorem but not by the stronger all-connected assumption alone.

## Limits

The contact condition is sufficient, not necessary. Global GT1 concerns one-color moves and remains separate. General rainbow-tree or graceful-fiber width-two connectivity is not asserted. These moves preserve labels, colors and F, but do not preserve the abstract tree type; no implication to the Graceful Tree Conjecture is claimed.

## Exact controls and dependencies

Run `python verify_connected_contacts.py` and `python -O verify_connected_contacts.py`. Outputs agree:1099799 checks and122732 validated moves. The checker verifies both tiny boundaries, every acyclic upper tail leaving three free colors at m=5,6,7 (2,6,22 fibers), and every ordered source/target pair in those fibers. It also replays the improved algorithm over all pairs of the earlier hostile fibers at k=2,3,4. The generic theorem and the all-m three-color conclusion follow from the proof, not from these samples.

The [prior frozen-fiber theorem](../graceful_fiber_obstructions_20261019.md) supplies the original one-color obstruction. The [leaf-buffer theorem](../graceful_two_color_buffer_20261019/GRACEFUL_TWO_COLOR_BUFFER.md) remains valid; this criterion subsumes its contact assumptions and improves the bound. All computational checks in the present packet are self-contained and use only the Python standard library.
