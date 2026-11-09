# Session 230: creator fibers and structural generation

**Status: ACTIVE RESEARCH CHECKPOINT. No new global completeness claim.**

The task is to follow the first structures realizing an odd Hamiltonian-path
count, and ask which controlled operations generate all tournament classes.
Chan–Pak completeness remains OPEN. First appearances are defined only for
attained counts; 7 and 21 are globally excluded by the cited prior proofs.

## Inheritance and distinct questions

The baseline is main commit `f870c2abf10f08ac8648923e313ca417e5e797f3`.

- Session224 proves that a birth fiber need not be one converse orbit: H=9
  has two self-converse classes at order five. Retain the entire fiber.
- THM-4922 / `first_creators_genealogy_20261107` supplies the two H=63 orphan
  creators and the class-level creator-deletion census through order ten.
- THM-4942 / `first_creator_generation_20261108` supplies a five-operation
  grammar through order eleven. Its numerical census and its proof status
  are being audited here, not silently promoted to an all-order theorem.
- THM-4169 supplies the cited Schmerl–Trotter prime-deletion route.
- THM-4114 and THM-4118 already supply directed-cut/submodular attachment
  responses and show that a gcd-two response lattice need not be a solid
  interval. Reusing these mechanisms is not a novelty claim.

A birth library is minimal in **vertex order for a prescribed count**. It is
not automatically an irredundant algebraic generating set. An all-class
construction is also strictly stronger than realizing each scalar count.

## Portfolio and decisive probes

1. **Anchor: structural grammar.** Audit the finite grammar, investigate
   substitution-prime obstructions, and decompose the H=63 orphan pair.
2. **Niche: sharp connectivity threshold.** Generalize the degree/cut argument
   and test the sharp boundary rather than extrapolating order eleven.
3. **Wildcard: proof objects as operations.** Refine the odd-cycle packing to
   an explicit attachment-cut kernel, crediting prior cut calculus, then
   determine what it can certify about missing values and minimum growth.

The concept board is birth fibers, modular primes, SCC words, cut capacity,
and odd-cycle packings. A result must retain the lost coordinate (class,
attachment, or cycle incidence), not replace it by scalar H.

## Preliminary verified distinctions

The two-vertex transitive tournament TT2 and a directed triangle C3 generate
C3[TT2,1,1], the unique strong four-vertex class with H=5. Thus even this first
creator is redundant once ordinal sum and substitution are allowed. This
small example disproves the claim that the complete first-creator library
is a minimal basis for that grammar. It does not disprove generation.

A fresh independent check of the 97 supplied order-eleven witnesses confirms
their path counts and connectivity and identifies eleven isomorphism classes,
including seven non-first-creator candidates with nontrivial modules. This
is witness soundness, not a fresh replay of the full large census or of all
membership labels in S(10). The detailed receipt will follow in a checked
checkpoint.

## Remaining obligations

Finish independent review of the new structural and kernel statements;
publish self-contained witnesses, optimization-safe checks and receipts;
separate recovered results, proved reductions, finite evidence and OPEN
claims; then verify every checkpoint on remote main. No arbitrary vertex
insertion or arc rewriting is presented as a solution to the research goal.
