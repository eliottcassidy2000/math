# Primary-source and inheritance pins

Checked 2026-10-21 UTC. These are bounded source attributions, not a claim of
full-paper or computational-certificate verification.

## External primary

Jaan Allikvere, *Fourteen and fifteen lonely runners*, arXiv:2609.02604v2.

- Version record: https://arxiv.org/abs/2609.02604v2 . The available record
  gives v2 dated 24 September 2026; v1 was submitted 2 September 2026.
- Exact mathematical location:
  https://arxiv.org/html/2609.02604v2#S3.SS1 . Section 3 defines
  `H=v^perp`, `D=P[-1,1]^n`, `E=P Q^(1/2) B_2^n`, with
  `q_i=1+2u_i-u_i^2` and `u_i=v_i/max_j v_j`. Lemma 3.1 proves inclusion;
  Lemma 3.2 gives the ellipsoidal covolume. The present comparison uses only
  these definitions and local formulas, not the main LRC theorem.
- Normalization: equation (1) and the paragraph following it index `LRC(k)`
  by nonzero speeds, which is `k+1` total runners.
- The paper does not assert that its chosen ellipsoid is maximal. Our
  strict fixed-shape volume comparison is therefore not a correction to
  its stated inclusion lemma. No global constant improvement is asserted.

## Repository inheritance (read-only snapshot)

Repository read at commit `de0463bc8ddf` (startup's cached origin agreed).
No fetch or Git write was performed for this packet. Final publication
requires live freshness reconciliation.

- `01-canon/theorems/THM-4014-lrc14-diagonal-polar-ellipsoid-fastest-coordinate-relation-compression.md`:
  exact pair criterion, compressed covariance, and the ambient-inverse
  hostile. The present extension admits a singular ambient diagonal but
  retains positive definiteness on the projected hyperplane.
- `01-canon/theorems/THM-4009-euclidean-covering-transference-short-relation-compression.md`:
  projected zonotope and speed-dependent inradius; necessary relation
  reductions must not be read as LRC sufficiency.
- `01-canon/theorems/THM-3743-lonely-runner-polyhedron-khinchin-flatness-relation-reduction.md`:
  projected-lattice dual and width context; its older numerical cap is
  explicitly superseded in that file.
- `00-navigation/OPEN-QUESTIONS.md`: the recent preprint's broader proof
  claims have a separately scoped computational audit. This packet does
  not settle or extend that audit.

## Limits

The singular endpoint, uniqueness and translated-ellipsoid proof received
an independent mathematical audit. A separate frozen-artifact review is
still required before publication. This is not a literature novelty or
priority claim. No repository files were modified.
