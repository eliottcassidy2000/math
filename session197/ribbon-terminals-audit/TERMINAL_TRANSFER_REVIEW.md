# Independent review: distinct-terminal one-seam transfer

**PASS**, no corrections requested, 2026-11-05.

Frozen source in `session197/ribbon-terminals/`:

- PROOF.md: 45383a2eeb947587a6ce35c2367053c81c2db1813da990cc19aa4a4a1e93f91e
- verify.py: 9ada1cc1cfaf33030a89122dcf4ed65103aa23f5a2ad01e5dd5c288cde198a3e
- Validation/normal/optimized output: f99d82d37e7fb4ebbbadc6df98d4de5386deb7ce3a305d06eacfd00e5f6374fd

Both copied-source replays reproduce the frozen receipt exactly. No source file was changed.

## General proof

The first row is the inherited one-reverse-seam path split. The only forbidden avoiding join has source endpoint a and gadget start b, giving -te. In the using case, the source after-run contains a and the gadget before-run contains b, so every other cross join is automatically valid.

For the ending-at-c row, the avoiding case restricts the gadget path to end at c and its forbidden joins to start at b, giving Eh-Fe. In the using case the last gadget run must be nonempty and end at c, so the correct response is the ordered pair of nonempty gadget paths ending at b and c. This is D, not an endpoint product. When b=c the using case is impossible and D=0. F is zero there except for the singleton, where F=1. These conventions correctly recover the same-terminal singleton correction.

For R, adding z that beats only c gives H(T^+)=H(T)-e_c(T)+R_c(T). In the same literal graph, z loses to every source vertex, so regrouping it with G leaves b->a as the sole reverse cross edge. Subtraction and addition of the first two rows yields exactly the displayed third row. This proves three-coordinate closure for all nonempty inputs; no fitted-state assumption or extra deletion coordinate is used.

The width compiler is correctly sufficient, not claimed exhaustive. With the old root last and the entry in the gadget's first w positions, the only reversed cross edge has distance at most w. Internal spans and the returned-last exit are preserved on concatenation.

## Physical clock

The four-vertex gadget has the stated edges, constants and augmented constants. Its entry is the third new vertex and its exit is last, exactly matching the sole cross-boundary reverse edge in the repeated width-three mask word0013. The preserved charge is (-2,1,1), with value-2 on the physical seed. The independently checked first-row identity gives the inhomogeneous second-order recurrence with constant-12. Subtracting1 from H makes it homogeneous; the characteristic polynomial is (x-1)(x^2-18x+5).

This is an observable constant component in a second-order physical family. It is correctly separated from the fixed first-order affine exclusions and makes no inverse-spectrum claim.

## Independent finite evidence

The new `terminal-audit/independent.py` uses literal vertex permutations and a single cut at the distinguished endpoint, rather than any subset DP. It independently computes source profiles, all six gadget constants, augmented-gadget constants, and actual output profiles.

Across every labelled source of order at most three and every labelled gadget of order at most three, with every root/terminal choice, it passes2,349 comparisons:841 coincident-terminal and1,508 distinct-terminal cases. Empty first/last path slots are included by the literal cut boundaries. It also independently recovers the example matrix and the full profiles(5,1,7) and(73,12,132) for the first two actual period graphs.

The source replay separately passes all32,045 comparisons through gadget order four, including8,265 coincident and23,780 distinct terminals. Its actual periodic graphs through order16, matrix identities, span checks and thirty-step recurrence controls also reproduce the frozen certificate.
