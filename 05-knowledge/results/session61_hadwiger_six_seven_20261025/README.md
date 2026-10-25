# Hadwiger–Nelson themes inside six/seven geometry

Status: second audited research checkpoint. Classical ingredients are attributed; the new contribution is an explicit comparison with the repository's marked hourglass geometry, not a literature-wide novelty claim.

## Main connection

The hourglass's six outer dual labels are the six units of the Eisenstein lattice. Its two conjugate norm-seven ideals give seven-colour quotients. Retaining the triangular faces before quotienting gives a torus with seven vertices, twenty-one edges and fourteen faces. Each face class is a genuine Fano design.

The two ideal choices and two face classes give **exactly the special four-element orbit of free Campedelli codes** in the inherited hourglass classification. All four select the same complex surface class, whose full automorphism group is `C2 × A4`, by the previously audited Session 60 theorem. This is an explicit map with an incidence and marking sidecar, not a coincidence of the integer seven.

- [Exact lattice and metric bridge](hexagonal/RESULTS.md): exactly two index-seven kernels, `χ(T)=3` but `χ(T²)=7`, and a classical whole-plane seven-colouring proved with strict cell diameter/separation bounds. The safe scale interval is `sqrt(3/7) < h < sqrt(3)/2`.
- [Torus → Fano → Campedelli](torus-fano/TORUS_FANO_CAMPEDELLI.md): explicit face systems, all thirty designs, all dihedral orbits, four selected free kernels, and the inherited geometric consequence.
- [Independent proof and hostile audit](audits/lattice/AUDIT.md): separate standard-library reconstruction, exact metric calculation, and the general rank-two lattice proof. The historical THM-419 repair is now supported by the [full lattice-cap packet](lattice-cap/LATTICE_CAP_REPAIR.md) and [second independent audit](audits/cap/AUDIT.md).
- [External source status](SOURCE_STATUS.md): the September 2026 six-colour lower bound is **CITED**, not independently built or axiom-audited here. None of these bridges improves the plane bound or resolves six versus seven.

## Recovered prior themes

- [Paired palettes, Fano completions and deck correlations](palette-deck/RESULT.md): fifteen six-label pairings match the fifteen genus-two Lagrangians and central Fano pencils. A completion bit and two hidden cohomology bits per cycle remain essential. Specified pairwise sign marginals are blind; triple correlations recover the Fano incidence and inherited surface type.
- [Arithmetic residue/CM boundary](residue-cm/RESULT.md): conjugation-fixed primes are exactly the local-integrality boundary for norm-one reduction. The CM split/inert pattern provides a real arithmetic comparison, while the cyclotomic-165 bound is a bound on the named residue-colouring method, not the target field's actual chromatic number.
- [Repaired lattice cap](lattice-cap/LATTICE_CAP_REPAIR.md): every positive norm layer in every rank-two Euclidean lattice is periodically at most three-colourable; the explicit former residual has odd girth 23. The proof does not rely on the defective historical census.
- [Independent palette/arithmetic audit](audits/palette-arithmetic/AUDIT.md), [independent lattice audit](audits/cap/AUDIT.md), and [correction preimage hashes](CORRECTION_PREIMAGES.json) distinguish proved statements, exact finite controls, inherited geometry and historical failures.

## What the quotient does not preserve

The triangular source graph is three-colourable, but its seven-vertex quotient is `K7`. This graph covering does not preserve chromatic number in both directions: a source colouring need not be periodic under the quotient lattice. Nor is the graph-square clique a unit-distance clique; its edge lengths include `sqrt(3)` and `2`. The Fano incidence is supplied by quotient triangular faces, not by literal Euclidean collinearity of a centred hexagon.

[Concept board](CONCEPT_BOARD.md) · [Publication checks](PUBLICATION_CHECKS.md) · [Exact manifest](MANIFEST.json)

Run `python3 verify_packet.py` from any working directory. Individual mathematical reproduction commands and their precise finite universes are in the linked notes. No integer-square hourglass existence, new K3 Picard rank, or plane six/seven closure is claimed.
