# Independent audit of the primitive bounded-span interface

Result: PASS. Frozen inputs were not changed.

Inputs and SHA-256:

- `session197/ribbon-primitive/PROOF.md`: `37d07f8353bd3ebd9960ef618360e32fc24b5d6c77baea10fe49fc91ce9014b3`
- `session197/ribbon-primitive/verify.py`: `ac507e37a555f8f9ffc73e4f19be904a365e8b8f1020850b0c0dfa1237478d12`
- `session197/ribbon-primitive/VALIDATION.json`: `d0e734c28aa8bf6693ea1350c280bbb17316d07daba9c28fa70be300a442743c`

I checked the pinned distinct-terminal transfer argument, especially the ending row and augmentation/regrouping proof of the cover row. The three primitive-packet matrices, preserved charge, reduced affine maps, determinants, and all three inverse formulas agree algebraically. The narrow image-cone inequalities hold throughout the stated physical domain and are disjoint, so the backward decoder is unambiguous. Its integrality/domain checks and strict height decrease make it both terminating and complete for the pair-state language.

The two observable continuation identities are correct. Their coefficients 48 and 64 leave only a possible power-of-two common divisor; oddness removes it, after which the conserved charge forces observable gcd one. This is not confused with scalar completeness.

The directed-cycle splicing proof preserves strongness for every cell sequence. The exceptional seam distance is at most three, and the explicit words match the physical cross-edge construction. Growth by at least eleven at each cell gives the displayed sublinear power bound on the number of outputs. The residue argument and finite-seed scope are correct.

I reran the full checker normally and under Python optimization. Both outputs match the supplied frozen receipt byte-for-byte: 9,841 word decodings, 62,500 admissible pair tests, all zero-/one-/two-cell graph profiles, and the continuation-gcd controls. No correction or remaining blocker was found. The conclusions concern the three named cells, not all bounded-span tournaments.
