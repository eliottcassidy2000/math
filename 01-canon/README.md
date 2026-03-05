# Canon: Admission Rules & Certainty Scale

The `01-canon/` directory contains only mathematics we are collectively most certain of.

---

## Certainty Scale

| Level | Label | Meaning |
|-------|-------|---------|
| 5 | **PROVED** | Complete proof exists in this repo or in a peer-reviewed citation |
| 4 | **PROOF SKETCH** | Proof sketch here is convincing; full proof expected to be routine |
| 3 | **VERIFIED** | Exhaustively or extensively computationally verified; no proof |
| 2 | **CONJECTURE** | Plausible but evidence is limited |
| 1 | **SPECULATION** | Worth tracking but not yet supported |

Only levels 4–5 enter `theorems/` as theorems. Level 3 enters as CONJECTURE with a VERIFIED badge. Levels 1–2 belong in `00-navigation/OPEN-QUESTIONS.md` or `00-navigation/TANGENTS.md`.

---

## Theorem File Template

Every file in `theorems/` follows this format:

```markdown
# [ID]: [Short Name]

**Type:** Theorem | Lemma | Corollary | Conjecture
**Certainty:** [1–5] — [label]
**Status:** PROVED | VERIFIED | OPEN | REFUTED
**Last reviewed:** [INSTANCE-ID] — [DATE]
**Disputes:** [links to court cases, or "none"]
**Tags:** #[tag1] #[tag2]

---

## Statement

[Formal statement. Use standard notation from definitions.md.]

## Proof / Proof Sketch

[Proof or proof sketch. Mark gaps explicitly as **[GAP]**.]

## Verification Record

[Computational verification details if applicable.]

## Notes & History

[Key context, how this was discovered, what it implies.]
```

---

## How to Update a Theorem

- If you improve a proof, edit the file and update **Last reviewed**
- If you find an error, open a court case in `02-court/active/` BEFORE editing
- If a claim is refuted, change Status to REFUTED and explain in Notes
- All edits must include your instance ID
