# Independent review: variable-offset returning clocks

**PASS. No correction required.** Review date2026-11-05. Frozen source directory `session197-hybrid-constructions/variable-offset-clocks/`:

- `PROOF.md`: `f064c0869f889e1ebfdad06b43b1f901d0cff8af42b5eab0dc3c9ce516e948d8`
- `core.py`: `9f454079dccbff774e3ebd38fd90a871947a54d434d9d3faa9cc5ba36eb7106d`
- `graphs.py`: `3204eeb1e9d2b05cac9bcc21c98f25697d10f925d6ff743393e945ce4e73615c`
- `verify.py`: `14bc8e56f8a7074875eeae96cd1ca73108b163113b562861f490c3c144b8a0a1`
- `VALIDATION.json`: `589682a2455e1e19dc4c8f77f36da4f87a98dd4d194087a7381527f474a880aa`

## Mathematical review

The ordered ports0,1,4,5,6 make L12 lose exactly to old4,5, producing the stated nine-vertex seed. Subsequent R7 and L24 births preserve all five ports and have the described opposite-side commutation. Moving the fixed gate past an R birth also preserves its incidence: the gate beats every root clone in either order. The canonical labeling, including the chronologically indexed but reverse-directed L chain, is consistent.

The root crossing cover uses the fixed old vertices1,4. The L crossing cover uses old4 and the fixed gate B, not a variable module. Both have size two uniformly in the other multiplicity. The positive-multiplicity count therefore lies in the nine-dimensional exponential tensor span. Successive Vandermonde recovery is valid, and the nine physical values determine the displayed formula for all positive m,t.

The t=0 case is checked independently by the complete run vector `(537,540,180)`, giving180*3^m-3. Although this agrees with the formal t=0 value of the final expression, it is not inferred by extending a positive-multiplicity interpolation theorem to an empty block. The repeated initial-gate hostile1483 correctly prevents treating that one-generation tripling identity as a return.

The four-coordinate state is sufficient on the physical family. In terms of X=3^(m+t),Z=3^m,D=2^(t+1), the scalar is

    H=40X+80DZ-20Z+D-5.

Substitution gives exactly the two stated affine updates. Their compositions commute identically, and physical graph commutation establishes the graph-level return independently of this algebra. The R offset10-2^(t+2) changes after L; no fixed ternary-digit interpretation survives that change. The negative compressed coefficient is correctly scoped to the reachable physical states.

Strongness follows from the strong source, the gate's two incidence signs and nonempty transitive substitution into the strong ten-vertex quotient. All boundary cases m=1,t=0 and m=t=1 are included. The minimum value at fixed t is H(1,t), strictly increasing, so the stated inverse search terminates and is sound and complete for this family. Uniqueness and unrestricted spectrum coverage are not claimed.

## Independent controls

`variable-offset-audit/independent.py` imports no producer code. It reconstructs the family directly as a quotient blow-up, independently counts19 graphs through order14 by walk inclusion–exclusion, and solves the rational nine-point tensor system. The coefficient vector in bases `(1,2^m,3^m) tensor (1,2^t,3^t)` is exactly

    (-5,2,0,0,0,0,-20,160,40).

Literal outside-permutation/root-slot enumeration independently recovers complete root vectors at t=0 and t=1. Every one of255 words through length7 has the correct full graph relabeling, strongness and four-coordinate state. Symbolic quotient checks verify both fixed singleton covers. The repeated-B graph independently has H1483. An independent full target inventory through100000 agrees with the inverse test.

Copied-source normal and optimized runs reproduce the frozen receipt byte-for-byte, including the explicit reverse-core gate-label swap. Original files were not modified. The all-order result rests on the uniform covers and exact interpolation, not on these finite controls.
