# Full count-residue coverage mod120 inside a provably sparse ray

**Proved with a finite exact certificate and an all-order pump.**
Every odd residue r modulo120 contains arbitrarily large lower
endpoints H0 of physical pairs (H0,H0+2) on the 24:13 marked ray.
These are two tournaments differing in just the mutual apex arc.
This does not assert all integers in any one of these classes occur.

Let the five suffix atoms, coded as 0,1,2,3,4, be respectively

    P, Q, UZ, UVZ, UWZ.

Their half-gaps are +1,−1,0,−1,+1 by the [carry theorem](ray/RAY.md).
Prefix any suffix by PQ, whose half-gap is zero. The
[wheel receipt](count_wheel.txt) lists one suffix of at most six
atoms for each r; its half-gap is +1 and the lower endpoint is r
modulo120. No enumeration of tournament isomorphism classes is needed.

For completeness, compute modulo120 using the two 12-state
three-run automata from the [apex census](../unit_shields_wheel_followup_20261116/apex/CENSUS.md).
Each automaton counts assignments of spine vertices to three
increasing runs in one proposed apex order. It tracks whether
the first run ends at a legal entry, the third run is already
started, and the middle run is empty or ends at a legal/illegal
entry. The empty middle run uses the mutual arc; a nonempty run
does not, so reading both automata supplies both orientations.
The verifier compiles each atom as a linear map, then independently
rescans each printed literal certificate digit by digit.

Apply the neutral block PQ repeatedly to the initial paired state.
The exact split state at t=1 and t=121 is identical modulo120,
with preperiod 1 and period 120 in this orbit. The state contains
all coordinates that the suffix reads. Therefore replacing the
prefix PQ by (PQ)^(1+120j) leaves **both** output residues
unchanged for every j≥0 and leaves the actual signed gap +2.
The lower endpoints become unbounded by the disjoint-cycle lower
bound in [PORTS](PORTS.md). Each pump adds 2880 sites, so the
spine width also remains fixed modulo120 within each lane.

This is a close Riesel analogue: a finite wheel picks a template
for each target residue, and a repeatable physical pump preserves
that residue. But a Riesel factor forces *every term of its
designated exponent class* to be composite. Here the templates
merely supply sparse subsets of each target class. The exact
entropy bound on the same ray gives O(X^0.071771) distinct
absolute counts through X, no matter how the pump is varied.
Modular completeness of the counts, not merely of the addresses,
still fails to give absolute placement.
