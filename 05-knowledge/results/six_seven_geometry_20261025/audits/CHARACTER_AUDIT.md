# Independent review of the character/sextic packet

2026-10-25. **Verdict: PASS**, within the stated characteristic-zero
geometric scope and the separate odd-finite-field count scope. Reviewed
`characters/CHARACTER_SEXTICS.md` and its checker without importing its code
into the independent controls. No substantive mathematical discrepancy found.
Concurrent Session60 attribution and the distinction between the two K3
quotients are appropriately explicit. Historical priority is not certified.

1. **Grading and algebra.** The full graded ring over the base polynomial
   ring, with base degree two, is free on all `2^n` squarefree root monomials.
   Its even Veronese is free on the `2^(n-1)` even monomials. On a chart
   `L_j!=0`, the monomial attached to even `S` becomes the squarefree
   monomial supported on `S\{j}` in the normalized roots `x_i/x_j`.
   Every subset of the other indices occurs exactly once. Thus the
   even-Veronese description really sheafifies to the finite algebra in C1;
   it neither doubles the degree nor assumes independent specialized values.

2. **Canonical eigenspaces.** The root cover is finite flat and Gorenstein;
   finite duality gives the stated eigensheaves. The characters have order
   two, so taking duals does not alter their labels. The A1 resolution
   preserves the canonical sections and irregularity. The contribution
   `h^0(O(|S|/2-3))` is zero below six and equals
   `binom(|S|/2-1,2)` otherwise. In particular the n=8 control correctly
   adds three sections at weight eight. The determinant twist in C3
   correctly cancels the global projective sign. The seven coordinate
   sections have complementary-six characters, not singleton characters.

3. **Quotient typing.** A six-root deletion has degree two from X7;
   the six-line double plane has degree 32 from X7 and degree 16 from
   that deletion. At ordinary triple branch points the double plane has
   `w^2=uv(u+v)`, a D4 rational double point, whereas the maximal root
   cover has an A1 point. These models must not share the same node count.
   The normal-model maps and the weaker smooth-model statements in the
   text are correct. The holomorphic-form pullbacks remain well-defined
   after resolving their rational correspondences. No finite morphism
   between minimal K3 surfaces is claimed.

4. **Finite-field fibres.** For a fixed base point and a nonzero scalar
   lambda, the number of affine root vectors is exactly the product of
   `1+chi(lambda L_i)`. The nonzero root vectors modulo Fq* form the
   projective lifts, and that scalar action is free even on zero-coordinate
   strata. Summing lambda kills odd supports, producing C5. Rescaling the
   homogeneous base representative also leaves all even-support terms
   unchanged. Spanning prevents the forbidden all-zero base-value vector.
   Characteristic-three incidence changes do not invalidate this count:
   the manuscript correctly counts the singular equations, not their
   characteristic-zero resolution.

5. **Independent hostile controls.**
   `check_character_support_independent.py` directly enumerates every
   point of P6 over F3, F5, and F7 and tests the four hourglass quadrics.
   It compares that third path with the even-character sum, and on every
   P2 fibre also enumerates every scalar and actual square-root multiplicity.
   This catches projective normalization, zero radicands, and the changed
   small-prime incidence. The producer's full run separately passed its
   canonical n=4..12, 35 chain counts, and all 36 seven-subset controls.

Reproduce the independent receipt with

```text
python -B 05-knowledge/results/six_seven_geometry_20261025/audits/check_character_support_independent.py
python -B -O 05-knowledge/results/six_seven_geometry_20261025/audits/check_character_support_independent.py
```

The checked claims are C1--C5 and their stated consequences. External
Picard/derived-equivalence results are explicitly unused and are not
promoted by this audit.
