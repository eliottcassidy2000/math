# Paper Rebuild Notes

When the time comes to write the new, cleaner paper from scratch, use this document to track ideas about what to include, what to cut, and what the structure should be.

---

## What Definitely Goes In

- **THM-001 (Rédei):** All four proof routes (A–D), clearly labeled as conditional/unconditional
- **THM-002 (OCF formula):** H(T) = I(Ω(T), 2) — the main new result, with full inductive proof conditional on Claim A
- **THM-003 (Claim B):** Proved algebraic companion
- **THM-004 (F1):** Algebraic identity (inshat−1)/2 = #Type-II — should be a lemma
- **THM-005 (F2):** Bijection between Type-II and 3-cycles — should be a lemma
- **THM-007 (F4):** Distribution formula C(L-2, 2k-1) — clean new result, add to paper
- **CONJ-001 (Claim A):** State clearly as the central open problem, with verification table
- **MISTAKE-002 note:** Explicitly state that the exact formula H(T) = B_v + S_v + R_v is FALSE (with evidence)
- **Failure analysis:** Clean explanation of WHY n=6 breaks the per-path identity (using THM-005)

## What to Cut or Restructure

- The LaTeX source is very long and complex with heavy formatting (tcolorbox environments, color definitions, etc.). The new paper should be more streamlined.
- The proof-status table format is good — keep it but simplify the surrounding text
- Remove the "Route dependency summary" box or integrate it as a simple list

## Open Questions to List in the Paper

See `00-navigation/OPEN-QUESTIONS.md` for the current live list. At paper time, the unresolved ones become the open problems section.

## Format Decisions

- Use AMS LaTeX (amsart or similar)
- Markdown for working drafts; convert to LaTeX only at final stage
- All figures should be TikZ (reproducible, no external image dependencies)
- The pin grid figures (IMG-001, IMG-002 in 03-artifacts/images/) should be reproduced in cleaner TikZ

## Title Ideas

- *Parity in Tournaments and the Odd-Cycle Collection Formula*
- *The Odd-Cycle Collection Formula: A New Identity for Hamiltonian Path Counts*
- Keep something about "tiling model" in the title or subtitle

---

## Images to Recreate

From `03-artifacts/images/INDEX.md`: the pin grid diagrams. These are the most visually important and will be needed early.
