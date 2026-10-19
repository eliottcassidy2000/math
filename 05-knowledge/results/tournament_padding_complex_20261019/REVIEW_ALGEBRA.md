# Independent algebra review of the tournament padding complex

19 October 2026. **PASS** for the mathematical statements in `PADDING_COMPLEX.md`. This review independently checks the formulas, including the boundary degree; it does not duplicate the accompanying finite-matrix verifier.

- The operator types are correct: d raises order, h lowers it by adjoining a universal source before evaluation. Both preserve isomorphism-invariant functions.
- Each two-vertex deletion occurs twice, so d squared is zero over F2. In hd, deleting the newly added source contributes exactly the original function; every other deletion is the corresponding dh term. The degree-zero convention C(-1)=0, h0=0 is valid: hd alone is the identity there.
- Hence any element of ker d equals d(hf), giving exactness and a concrete preimage. No assumption that h squared vanishes is needed or asserted.
- For spanning injective embeddings, adding an isolated source vertex is exactly d. Mapping to a target with a new universal source forces the preimage vertex to have indegree zero, giving the stated h formula, including isolated source vertices.
- A tournament automorphism group has no involution, because any transposed pair would reverse its arc. Cauchy's theorem makes the group order odd. Therefore tournament embedding parities are indeed class indicators; the operator identities are on the full function spaces, not an unstated smaller span.
- Rank-nullity and exactness give t_n=r_(n-1)+r_n, and thus the alternating integer rank formula. The stated values through d4 are correct; rank d5=9 is appropriately labeled a theorem from exactness rather than a directly checked six-host matrix.
- For odd n, applying the next d to a putative d_n f=1 gives a contradiction. For even n the contracting identity yields f=c*1+d(hf) when d f=c*1. The transitive-host argument also correctly proves the weaker isolated-vertex Rédei obstruction without invoking exactness.
- The graph-level primitive can be a sum and need not remain a tree. This information loss is explicitly retained. No general tree-classification or integer-count conclusion has been smuggled through the mod-two identities.

The full theorem is an exact and useful operator formulation of the same padding mechanism used in the separately checked eight-tree and nine-tree results.
