# Core-paper intake: Grothendieck, Banach, SOP, and Crouzeix

> **Status / freshness:** primary arXiv records checked 2026-08-14.  Every
> headline result below is a **CITED VERY RECENT PREPRINT**.  It enters proved
> canon only where a separate self-contained repository theorem is named.

## Saha--Li--Xue et al. — *New Lower and Upper Bounds for the Grothendieck Constant*

- **Primary:** [arXiv:2608.11158v2](https://arxiv.org/abs/2608.11158), revised
  2026-08-12.
- **Imported role:** claims the real bound
  `6pi/11 <= K_G <= 1.7818666069360661`.  The lower mechanism is the sharp
  affine Gaussian-sign strip `b_3>=2b_1-11/6`; the upper mechanism combines an
  odd cubic--quintic correlation map, Gaussian realization, and a certified
  inverse-Wiener majorant.
- **Repo consumers:** [THM-3392](../../01-canon/theorems/THM-3392-bipartite-sign-lift-and-synchronization-loss.md)
  uses only the cited upper number after proving its own optimal factor-two
  synchronization theorem.  [THM-3396](../../01-canon/theorems/THM-3396-four-bit-pairwise-independent-fourier-cone.md)
  is a self-contained finite Walsh analogue; it does not depend on the paper.
- **Does not prove:** a complex Grothendieck bound, tournament equivalence, or
  any LRC/FC/JC statement.  A stronger number advertised in a public code
  README without a populated proof folder is not imported.

## Lu--Yang — *A solution to Banach's isometric conjecture*

- **Primary:** [arXiv:2608.13536v1](https://arxiv.org/abs/2608.13536).
  **CLAIMED SOLUTION, v1.**
- **Imported role:** claims the real finite-dimensional Banach isometric
  conjecture.  Finite moment tensors detect the stabilizer; exact local
  sections form a principal bundle; a positive-degree sphere pullback kills
  the relevant torsion class; signed degree produces adjacent polynomial
  powers; unique factorization forces the squared norm to be quadratic.
- **Does not prove:** the remaining complex version, JC(2), or that killing
  finite monodromy makes a rational coordinate polynomial.  THM-3383 is the
  exact repo hostile: its cover decodes rationally while one target still
  fails two boundary divisibilities.  No canonical theorem depends on the
  preprint.

## Chernikov — *SOP2 = SOP3*

- **Primary:** [arXiv:2608.13291v1](https://arxiv.org/abs/2608.13291), submitted
  2026-08-13.
- **Imported role:** proves a first-order model-theoretic equivalence by a
  mixed-partial-type dichotomy; the witness construction retains
  witness--parameter pairs rather than quotienting to bare vertices.
- **Repo consumers:** [THM-3395](../../01-canon/theorems/THM-3395-small-sheet-typed-cover-star-cochain.md)
  cites that typing move as motivation only; its theorem is instead an
  elementary coset/cochain/CRT/Helly argument.  The `RESERVED`
  [THM-3456 candidate](../../01-canon/theorems/THM-3456-left-permutive-trace-bijection-and-rule30-seed-boundary.md)
  gives a self-contained enriched trace-incidence special case, while the
  [typed-incidence reflection](../../07-reflections/sop2-sop3-typed-incidence-compiler-and-rule30-trace-boundary-codex-20260815.md)
  conditionally repackages the paper's compiler as a set-system lemma.
- **Does not prove:** that an LRC sheet clutter has SOP, that speeds are
  faithful vertices, any named-seed Rule 30 property, or any LRC(14)
  decrement.  It is context rather than a dependency of THM-3456.

## Lorist--Schwenninger — a claimed proof of Crouzeix's conjecture

- **Primary:** [arXiv:2608.03841v1](https://arxiv.org/abs/2608.03841).
  **CLAIMED SOLUTION, v1.**
- **Imported role:** an all-iterate commuting completion of adjoint powers is
  used to obtain scalar Crouzeix constant two.  The repository independently
  extracts and proves the mass-`c` finite-dimensional completion lemma in
  [THM-3390](../../01-canon/theorems/THM-3390-all-iterate-commuting-completion-norm-bound.md).
- **Does not prove for the repo:** complete or matrix-valued Crouzeix, a better
  constant than two, or a raw-compression theorem.  THM-3390's first-commutator
  hostile shows ordinary compression collapses the gate to normality; the
  external headline remains cited, not canonized.

## Shared transfer boundary

The papers expose different forms of one three-gate mechanism:

```text
local compatibility  ->  global closure  ->  target-category effectivity.
```

Gaussian correlations, bundle sections, partial types, and commuting
all-iterate defects solve different first or second gates.  None permits an
untyped transfer to positive LRC covers, polynomial JC coordinates, or
Hadamard completion.  The corresponding exact repo sidecars are respectively
ternary support/higher chaos (THM-3392/3396), boundary divisibility
(THM-3383), typed coset stars (THM-3395), and a distinct positive completion
layer rather than raw compression (THM-3390).
