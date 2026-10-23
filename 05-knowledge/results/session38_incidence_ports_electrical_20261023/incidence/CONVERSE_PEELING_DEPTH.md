# Alphabet size versus converse-peeling depth

October 23, 2026. **PROVED:** the elementary tree-depth formula below.
**FINITE-EXACT:** the complete order-eight signature forests for E_43,
E_33, and the ternary sum palette. This depth concerns anonymous signature
MARGINAL equations, not rounds of vertex Weisfeiler--Leman refinement.

In the signature-incidence forest, one parallel round solves every class
edge incident to a degree-one vertex, subtracts its forced mass from the
other endpoint margin, and removes all such edges. Empty components need
zero rounds. A nonempty tree of diameter D needs exactly ceil(D/2) rounds.
For a forest, take the maximum over its components.

Proof: removing all leaves shortens every nontrivial longest path by two
edges, until a single edge or a star is removed in the final round. More
formally, for D>=2 the surviving nonempty tree has diameter D-2; D=1 and
D=2 disappear in one round. Induction gives the formula. The same deletion
is valid for marginal equations because a leaf's sole incident unknown is
its known remaining margin. Parallel edges are not trees and cannot be
silently simplified.

At order eight:

| Signature model | Forest component diameters | Classes solved by round |
|---|---|---|
| E_43, twelve modular symbols | 6,880 diameter-one components | 6,880 |
| E_33, nine modular symbols plus reflection | 6,864 diameter-one; eight diameter-two | 6,880 |
| Ternary sum plus reflection | 6,820 diameter-one; 24 diameter-two; four diameter-three | 6,876, then four |

E_43 already identifies target classes directly, so no separate converse
inference is necessary even though its forest representation has one-edge
components. E_33 uses one layer of reflection marginal elimination. The
ternary sum rule needs two layers on four exceptional classes. Thus the
smaller alphabet trades explicit labels for a small amount of additional
inference; it does not reduce the total class inventory or remove the
source/path/target association.

The four middle class canons are 8393858, 8654866, 25694466, and 34621760.
They form two converse pairs, with Hamiltonian-path counts 191,197,191,197.
For example, the first has forward signature fiber
{8393858,8393861}, while its converse-margin target fiber is
{8397924,8393858}. The other two classes are leaves; after their masses are
read, the middle mass is the remaining marginal. These are fixed class
coordinates, not a new count of four actual ports.

The independent ternary audit directly finds the same parallel peeling
counts (6,876,4), so the depth-two claim does not rest only on component
sizes. A global forest is sufficient, not necessary, for wiring uniqueness;
additional capacities can still resolve cyclic signature systems.
