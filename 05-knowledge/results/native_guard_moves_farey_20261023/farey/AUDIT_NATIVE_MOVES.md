# Independent audit: native owner-clock moves

23 October 2026. **VERIFIED:** the final local-move and balanced-rectangle
proofs pass. **FINITE-EXACT:** checker and receipt audit passes for the
four stated order-eight sources, including both selected-port rectangles
and the distinct-port obstruction to amplitude two.

Files: [`../native_moves/NATIVE_GUARDED_MOVES.md`](../native_moves/NATIVE_GUARDED_MOVES.md),
[`../native_moves/check.py`](../native_moves/check.py) and
[`../native_moves/receipt.json`](../native_moves/receipt.json).

## Independent controls

Both final normal and optimized runs pass **70,190 checks**. Separately from the
producer's Kruskal/union-find algorithm, direct threshold breadth-first
search recomputed the connected components of every saved accepted fibre.
This independent path, saved as
[`AUDIT_NATIVE_MOVES_CHECK.py`](AUDIT_NATIVE_MOVES_CHECK.py), passed
**6,308 checks** over **1,576 fibres**:

- Source 1238 has 409 Hamiltonian paths.
- Source 2242 has 407 Hamiltonian paths.
- Source 2246 has 523 Hamiltonian paths.
- Source 2315 has 527 Hamiltonian paths.

The independent census rebuilt each adjacency matrix from its canonical
28-bit word, enumerated all `8!` permutations, and compared the complete
Hamiltonian path set to the receipt. Each fibre was reconstructed by
filtering all those rows by its endpoints, phase, and optional indexed
arc. For each nonnegative cost threshold, an edge was inserted exactly
when the number of unequal owner clocks was at most that threshold. The
first connected threshold and the full partition at one less both matched
the receipt. This validates minimality independently of the producer's
spanning-tree selection.

## Mechanism and scope audit

The two-clock transposition legality criterion is exact. Swapping positions
at distance one reverses their own required arc. Swapping at distance two
would require both orientations through the fixed middle owner. At distance
at least three, the only changed incidences are the predecessor and
successor arcs at the two positions. The original and exchanged paths are
both legal precisely when both moved owners receive from both existing
predecessors and send to both existing successors, with boundary owners
omitted when absent.

At distance exactly three, the four consecutive owners form a strong
four-vertex tournament: the old and new paths force five of the six arcs;
either orientation of the remaining endpoint pair yields scores
`{1,1,2,2}`. The fixed inner arc and actual outer boundary incidences must
still be retained. This is a structural local motif, not an identification
of arbitrary four-vertex graphs with tournament classes.

The actual witness at source 1238 has paths

```text
(7,1,0,4,2,6,5,3)
(7,1,0,5,2,6,4,3)
```

It changes owner clocks 3 and 6, pins endpoints 7 and 3 and the indexed
middle arc `2 -> 6` at clock 4, and changes target owner 881 to 2762. Both
targets have potential residue two while the source has residue one, so
the retained source port phase is seven. The sources' automorphism groups
are trivial; this is not a disguised automorphism-orbit move.

Three boundaries are essential:

1. A move here is an arbitrary **simultaneous** reordering from one accepted
   Hamiltonian presentation to another. Its intermediate permutations need
   not be legal, and no decomposition into single swaps is asserted.
2. With both endpoints pinned there are at most six changed clocks; adding
   a disjoint internal indexed arc leaves at most four. These upper bounds
   are immediate from the number of unpinned positions. The useful finite
   content is the sharp lower obstruction and explicit cross-target
   phase-preserving witness, not mere connectivity when all complete-fibre
   pairs are available.
3. Preserving `x_source-x_target mod8` keeps that individual source port's
   phase. A retargeting does not by itself preserve total target capacity,
   every node's harmonic sum, all other marked characters, or globally
   balanced metagraph wiring. Any such extension needs simultaneous owned
   port accounting. No all-size move budget follows from this census.

## Final balanced-rectangle audit

The selected-source change at `S,R` and the frozen four-node Kirchhoff
calculation are correctly separated. Direct expansion of the removed and
added edge Laplacians gives

```text
Delta L x = w*(x_U-x_V, x_V-x_U, x_S-x_R, x_R-x_S).
```

Thus keeping a source port's phase needs the target-potential congruence;
keeping the whole four-node current additionally needs the source-potential
congruence, with the edge weight retained. Cancelling an even weight modulo
a dyadic modulus would be invalid. The proof does not do so.

Independent inspection of the complete path rows recovers fourteen common
unordered target pairs under endpoint/phase pins and exactly two under
the additional indexed-arc pins: `(881,2762)` and `(2650,2843)`. Direct
reading of the original compressed edge atlas gives the stated pure-black
weights `4,2,4,4` and `4,2,4,6`; no blue weights are present. The exact
potential lists `(97,65,170,226)` and `(97,65,150,86)` give paired defects
`(-112,112,64,-64)` and `(128,-128,64,-64)`. Their greatest retained dyadic
moduli are respectively sixteen and sixty-four, not 256.

The amplitude-two obstruction is genuine under the declared distinct
observer rule. In the entire 409-path source-1238 fibre, exactly one path
targets 2762 and exactly one targets 2843. Two copies of either selected
rectangle would need two different source-1238 paths on that scarce side,
before or after the switch depending on its orientation. They do not
exist, even before additional owner pins are imposed. The weight-level
shadow does allow the chosen double removals, so this is an explicit loss
of actual port ownership under aggregate projection. A single witness
establishes the matching amplitude-one lower bound.

Finally, the target pairs are disjoint. Each rectangle contributes
`+64t,-64t` at its own two target coordinates. Therefore a combined repair
modulo 256 requires each amplitude to be divisible by four independently;
the literal bounds `|t|<=1` leave only the zero combination. This conclusion
is restricted to these two fixed rectangles, their fixed source owners,
and distinct actual observer selections. It neither rules out a larger
block using other target pairs nor constructs altered global native
wiring. The final prose states those boundaries explicitly. No correction
is requested.

Reproduction:

```text
python -B 05-knowledge/results/native_guard_moves_farey_20261023/native_moves/check.py
python -O -B 05-knowledge/results/native_guard_moves_farey_20261023/native_moves/check.py
python -B 05-knowledge/results/native_guard_moves_farey_20261023/farey/AUDIT_NATIVE_MOVES_CHECK.py
python -O -B 05-knowledge/results/native_guard_moves_farey_20261023/farey/AUDIT_NATIVE_MOVES_CHECK.py
```
