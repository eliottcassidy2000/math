# The procedural approach foundry: a typed generator of approaches to Collatz-type problems, and where a new mechanism is required

**Status: METHOD / TYPOLOGY (this session's classification; not a theorem).
The exceptional-set probes are FINITE-EXACT and inherited from the
[choice ladder](collatz_procgen_20260922_choice_ladder.md). Barrier
assignments are modelling judgments, each with a printed rationale
(`--explain`). A "live" card only means no barrier in this model blocks
it. Collatz OPEN.**

Script: `04-computation/experiments/collatz_procgen_20260922_foundry.py`.
Output: `collatz_procgen_20260922_foundry.out` (`--probes --explain`).

## 1. Grammar

An approach card is a tuple

```text
PROBLEM x SUBTARGET x MECHANISM x LENS  (+ optional OPERATION)
```

drawn from seven problems, seventeen mechanisms and thirteen lenses. The
problems are Collatz, E-SCC, the `3n-1` sheet, `5n+1` divergence,
Lagarias's rational periodicity conjecture, Mahler's `Z`-numbers and
Erdős's ternary digits of `2^n`. The lenses include 2-adic, 3-adic,
archimedean, adelic, affine-monoid, staircase carry, Beatty clock, base-6
CA, tropical, function field, exceptional dimension, choice relaxation
and sheet. The generator currently emits `473` cards.

Each mechanism carries:

* a set of **control barriers** it is blind to. It cannot separate the
  target from a control where the target is false:
  * SHEET: `3n-1` has three positive cycles, and all odd `b` are
    2-adically conjugate.
  * DRIFT: `5n+1`, `7n+1`.
  * DEFECT: planted density-zero modifications.
  * INTEGRAL: every parity word has a rational 2-adic cycle.
  * UNIFORM: generalized Collatz maps are undecidable.
* **structural requirements**. The one used here is THIN: residue
  certificates, escape induction, choice strategies, Krasikov--Lagarias
  inequalities and finite-state arguments close only when the exceptional
  set is thin.

Each problem splits into sub-targets with the controls each one needs, and
carries its exceptional-set dimension. A card is blocked when the
mechanism is blind to a needed control or fails a structural requirement.
A **hybrid plan** assigns one unblocked mechanism to each sub-target; the
glue between the parts is the remaining proof obligation.

## 2. Output

| problem | sub-target | unblocked real mechanisms |
|---|---|---|
| Collatz (dim `0.95`) | no-divergence | **none** (only the placeholder "some potential function") |
| Collatz | unique-cycle | linear forms in logs, staircase anti-concentration, functional equations, order patterns |
| Collatz | finite check | verification |
| E-SCC (thin at finite levels; dimension open) | Q1, Q2 | escape induction, Krasikov--Lagarias LP, finite-state, transversality, automata (25 hybrid plans) |
| `3n-1` sheet | no-divergence | **none** |
| rational periodicity | no-divergence for all `3x+d` | **none** |
| `5n+1` divergence | one divergent orbit | automata, order patterns, transversality |
| Mahler `Z` (dim `0.585`) | no safe integer code | order patterns |
| Erdős ternary (dim `0.631`) | transversal avoidance | order patterns |

Probe row (FINITE-EXACT, `2^20`): Collatz `27,328` exceptional classes;
minus sheet `27,328` (identical, the SHEET barrier made visible); E-SCC
`664`; `E_S` with `S={6 mod 8}` `2,165`; `5n+1` `232,912`, and with
choice still `133,600` (the DRIFT barrier made visible).

## 3. What the typology says

1. **The cycle half of Collatz has live mechanisms; the divergence half
   has none.** This matches the literature. There are unconditional
   theorems excluding cycles with few blocks (Steiner, Simons--de Weger,
   Hercher: CITED via the portrait lane). Every divergence result is a
   density or log-density statement.
2. **A new idea for no-divergence must do four things at once:**
   * see the drift: use `log 3<2 log 2` quantitatively, or it proves
     nothing about `5n+1`;
   * be defect-aware: act on the exact arithmetic of each orbit, not on
     densities;
   * use a non-uniform arithmetic input, the role Baker's theorem plays for
     cycles;
   * work against an exceptional set of dimension `h(log_3 2)=0.95`, so no
     finite list of escape lemmas can suffice.

   The only listed candidate is "transversality plus a Diophantine input".
   It would say that the parity sequence of a positive integer cannot
   persistently stay above density `log_3 2`, and it is exactly the
   unproved "every-`n`" form of Tao's theorem. The foundry does not supply
   it; it locates it.
3. **Relaxations with thin exceptional sets are the natural training
   ground.** E-SCC passes every barrier its sub-targets need (SHEET is not
   needed, since the relaxed problem is sheet-symmetric and true on both
   sheets). The choice ladder shows that the thinness comes from the
   `6 mod 8` rising-run excursion. HYP-9120 is the concrete next target.
5. **Corrections from the literature lane**
   ([barrier atlas](collatz_procgen_20260922_barrier_atlas.md), primary
   sources):
   * **UNIFORM.** Kurtz--Simon exclude only uniform *complete* methods, not
     sound ones.
   * **THIN/DIMENSION.** This is a heuristic, not a theorem, and "finitely
     many escape lemmas" should read "finitely many parametrized escape
     families". Applegate--Lagarias needed an infinite multiplier family
     for their one class.
   * **SHEET.** Residue statements must be transported by `x->-x`, not
     compared class by class.
   * **Tao.** Tao's theorem is now PROVED blind to SHEET: GGM 2025 covers
     `3N+r` for every odd `r`.
   * **Verification.** The bound is `2^71` (Barina 2025).

   The atlas also finds that no published all-orbits mechanism overcomes
   both SHEET and DRIFT at unbounded complexity. Removing DIMENSION by
   choice (Applegate--Lagarias, `E`) exposes SHEET: the hostile class `-1`
   is the minus sheet's fixed point. This matches section 3 item 2 here.
6. **Caveat on UNIFORM (original wording).** Undecidability blocks only criteria claimed to be
   necessary and sufficient. A decidable *sufficient* condition, such as a
   rewriting-termination interpretation found by search, is not blocked.
   So automata-type searches are marked UNIFORM-blind only as full
   criteria. The table above follows the stricter reading; under the looser
   reading automata-rewriting becomes a live candidate for the Collatz
   no-divergence half. Yolcu--Aaronson--Heule tried this and it has not
   succeeded (CITED, not re-verified here).

## 4. How to extend

Add a problem (sub-targets with needed controls, exceptional dimension), a
mechanism (blind controls with rationale, structural flags, optional
probe), or a lens; rerun. The intended use is as a pre-registration filter:
before investing in an approach, read its card, check which control it must
fail on, and run the attached probe on that control first.
