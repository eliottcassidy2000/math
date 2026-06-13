# Lean 4 Formalization Skeleton: THM-343

`H(T) ≠ 7` for every tournament T.

## Files
- `Thm343.lean` — Lean 4 skeleton with type definitions, lemma statements, and `sorry` placeholders.

## Status
**Skeleton only.** No proofs filled in yet. Use the file as the starting point for incremental formalization.

## Build
Requires Lean 4 + Mathlib4. To validate the skeleton parses:

```sh
# In a Mathlib-aware project:
lake build
```

## Order of attack (suggested)

1. `Tournament` type + decidability (~50 LoC).
2. `induced` (sub-tournament) + decidability (~50 LoC).
3. `DirCycle.isSC` (cycle is SC) — trivial (~30 LoC).
4. `SC_union` — case-analysis on transit through shared vertex (~80 LoC).
5. `cycleCount` definition — straightforward (~30 LoC).
6. `moon_cycle_count` for k=3 only (Moon-Moser bound) — induction (~300 LoC).
7. `camion` — corollary of (6) at k=n (~30 LoC).
8. `OCF` — axiomatize for now (line 1 above).
9. Main theorem `H_ne_seven` — case analysis using above (~150 LoC).

## See also
- `03-artifacts/drafts/lean-formalization-plan-thm343.md` — detailed plan.
- `01-canon/theorems/THM-343-H7-impossible.md` — informal proof, with S6 audit.
- `04-computation/audit_thm343_s6.py` — computational verification of every assumption at small n.
- `05-knowledge/results/audit_thm343_s6.out` — audit results (all OK).
- `05-knowledge/results/audit_n7_exhaustive_s6.out` — exhaustive n=7 H-spectrum (H=7 has 0 occurrences in all 2,097,152 tournaments; also H=21 and H=63).
