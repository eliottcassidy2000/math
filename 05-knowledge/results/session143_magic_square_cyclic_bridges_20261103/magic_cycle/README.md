# Session143 magic/cycle lane

**Main result:** for `f(k)=Ak^3+Bk^2+Ck+D`, `k=-4,...,4`, a real 3-by-3
magic arrangement exists **if and only if `B=0` and (`A=0` or `C=-15A`)**.
This completely classifies numerical assignments of degree at most three to
the nine centered labels, including repeated values. [Proof](THEOREM.md).

A nonzero third centered moment independently excludes every genuine
quadratic. In particular, no squared nonconstant real root AP works, even
for signed roots and repeated entries.

The `4,8,3,7,2,6,10,5,9,4` cycle is a nine-state `+4 mod9` clock. Its two
ternary digits require a quadratic carry; no relabeling makes it affine over
`F3^2`, and its odd transition cycle cannot embed in a hypercube graph.

**Sharp positive continuation:** `k^3-15k` does admit a magic arrangement.
All 945 center/opposite-pair systems and all 384 odd-cubic signed linear
matching branches are classified exactly. A separate normalized-coefficient
checker verifies the four nonodd candidate rays independently. This is an
ordinary magic-square result, not a magic-square-of-squares construction.

This merges the lanes by completely classifying degree-at-most-three numerical
assignments. **The unrestricted magic square of distinct integer squares
remains OPEN.** Higher-degree numerical assignments and general nonlinear root
assignments are not classified.

## Inheritance and concept board

| Concept | Recovered mechanism | Decisive invariant / boundary |
|---|---|---|
| Magic-square entries | Session57 center and decorated-biclique model | Mean center and nonnegative opposite pairs |
| Square line totals | Session57 quaternion witness | Square totals need not be equal |
| Cyclic versus vector labels | THM-2657/2788 odometer carries | Order9 is impossible for an affine F3-plane action |
| Numerical relabeling | Real polynomials of degree at most three | Cubic centered moment excludes quadratics; exact opposite matchings classify cubics |
| Cube interpretation | Actual transition graph, not point drawing | Odd cycle cannot map to a bipartite graph |

Closest inherited proof mechanisms: center/mean identities and odometer
power-map obstruction. Hostile controls: shifted Lo Shu is magic but its
entrywise square is not; the square-total witness is not magic; suppressing
the carry reduces the order from9 to3. Least-used useful sidecar: the actual
center/opposite-pair involution rather than a square-sum Hamiltonian path.

Seven bounded repository sources (about69KB) are copied with commit/blob and
SHA256 pins. This is targeted inheritance, not a whole-repository duplication
search. No new theorem ID is reserved. The cited carry principle is inherited;
the signed-AP proof is recorded without a literature-originality claim.

Run the commands in [THEOREM.md, reproduction](THEOREM.md#6-reproduction).
All checkers use only Python's standard library and raise explicit errors;
optimized execution does not remove their checks. A portable manifest pins
the whole bundle independently of the source checkout.
