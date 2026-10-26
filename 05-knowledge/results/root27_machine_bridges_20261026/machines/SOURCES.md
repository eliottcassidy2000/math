# Primary sources and inspection scope

Checked 2026-10-26 UTC. All links were opened during this session unless
explicitly labelled as a source-linked reference. No external Lean or Coq
build was performed. Web descriptions and repository claims remain CITED.

| Source | Date / identity | Used for |
|---|---|---|
| [BusyBeaverWiki nonhalting atlas](https://wiki.bbchallenge.org/wiki/Non-halting_Turing_machine) | Current page inspected 2026-10-26 | Exact spelling and transition tables for spaghetti sandwich and irregular bell; not a permanent chaos theorem |
| [Bouncer](https://wiki.bbchallenge.org/wiki/Bouncer) | Page states last edited 2026-07-15; oldid 8027 | Two-state native table and transcript |
| [Bell](https://wiki.bbchallenge.org/wiki/Bell) | Page states last edited 2026-03-28; oldid 6789 | Cubic/exponential example tables; displayed exponential transcript and recurrence rejected by our native proof/control |
| [Tony Guilfoyle's decider announcement](https://discuss.bbchallenge.org/t/decider-bouncers/126) | 2022-11-23 | Distinguishes a verified bouncer decision from merely probable bell classification |
| [Bouncer Coq verifier](https://github.com/meithecatte/busycoq/blob/master/verify/Bouncers.v) | Source-linked reference from the Bouncer page | Formal-verifier lead only; not rebuilt or separately audited here |
| [Leng Goldbach25 repository](https://github.com/lengyijun/goldbach_tm) | Current source inspected 2026-10-26 | Author's claimed formal verification and original-machine lineage |
| [Native Goldbach table](https://raw.githubusercontent.com/lengyijun/goldbach_tm/main/GoldbachTm/Tm25/TuringMachine25.lean) | SHA256 `9fa72dad7c4b903747ab81debe6cc1a068e9cf548fe1a10ab0b33375e791e9e0` | Read, parse and offline transcription; definitions of nth_cfg and cfg45 |
| [Goldbach semantic theorem](https://raw.githubusercontent.com/lengyijun/goldbach_tm/main/GoldbachTm/Tm25/Content.lean) | Theorems halt_lemma', halt_lemma_rev | Both directions of the blank-halting/counterexample equivalence; read, not built |
| [Cryptids semantics](https://github.com/google-deepmind/cryptids/blob/main/Cryptids/Semantics/README.md) | Current source inspected 2026-10-26 | External 21-state reduction and provenance; no “25-state current record” claim |
| [Lagarias, Elementary Problem Equivalent to RH](https://websites.umich.edu/~lagarias/doc/elementaryrh.pdf) | Author's July 29, 2001 version; journal publication 2002 | Universal divisor-sum inequality; our interval-witness quantifier deduction is separate |
| [Aaronson, Busy Beaver Frontier](https://www.scottaaronson.com/papers/bb.pdf) | 2020 survey | Native clock convention and halting-oracle equivalence |
| [Cook, Complexity of Theorem-Proving Procedures](https://www.cs.cmu.edu/~15455/resources/Cook1971-complx-thm-proof.pdf) | 1971; retyped primary paper | Uniform polynomial reductions versus bare simulation |
| [Buss–Rosu hidden-algebra paper](https://mathweb.ucsd.edu/~sbuss/ResearchWeb/HiddenAlgebra/incomp.pdf) | Author-hosted primary PDF | Configuration semantics and the Pi^0_2 totality family; standard result, not a classification of fixed Collatz |
| [Kurtz–Simon, generalized Collatz](https://link.springer.com/chapter/10.1007/978-3-540-72504-6_49) | TAMC 2007, LNCS4484 pp542–553; publisher abstract read | CITED Pi^0_2 completeness of their uniform generalization; author PDF fetch failed, so no new full proof audit claimed |

The webpage's machine *tables* and our native interpreter are distinct from
its unverified dynamical nomenclature. A table with no missing instruction
proves nonhalting directly; neither an image nor the word “irregular” is
needed. Conversely no amount of bounded native execution verifies the
Goldbach theorem on all even integers; that direction is supplied only by
the cited mathematical/formal source.
