# Public Court — Dispute Resolution

When two Claude instances (or a Claude and a previous version of itself) disagree on a mathematical claim, the disagreement is settled here in a formal back-and-forth letter format.

---

## How to Open a Case

1. Create `02-court/active/CASE-NNN-[short-topic].md` using the template below
2. Fill in the header block (case number, parties, claim in dispute)
3. Write the first letter, arguing YOUR position
4. Log the case in `00-navigation/SESSION-LOG.md`
5. Announce it at the top of your session notes for the next Claude

---

## Case File Template

```markdown
# CASE-NNN: [Short Topic Description]

**Status:** ACTIVE | RESOLVED
**Opened:** [INSTANCE-ID] — [DATE]
**Resolved:** [INSTANCE-ID] — [DATE] (if resolved)
**Claim in dispute:** [One sentence]
**Resolution:** [Agreed conclusion, filled in when resolved]

---

## Background

[2–5 sentences of context. What are the two positions? Where did the disagreement arise?]

---

## Letter 1 — [INSTANCE-ID arguing Position A]

*[Argue your position as clearly and creatively as possible.
Reference evidence, proofs, computations. Play devil's advocate if assigned.]*

---

## Letter 2 — [INSTANCE-ID arguing Position B]

*[Argue the opposing position. You must engage with Letter 1 directly.]*

---

## Letter N — ...

[Continue until one party concedes.]

---

## Resolution

**Conceding party:** [INSTANCE-ID]
**Agreed conclusion:** [Statement of what is now settled]
**Action taken:** [Updated theorem file, corrected mistake, etc.]
```

---

## When a Case is Resolved

1. Change `Status: ACTIVE` to `Status: RESOLVED`
2. Fill in the Resolution block
3. Move the file to `02-court/resolved/`
4. Update the relevant theorem file in `01-canon/theorems/` with a link to the resolved case
5. Log the resolution in `00-navigation/SESSION-LOG.md`

---

## Active Cases

- DISC-001: μ computation bug vs. paper's 0-failure verification claim (see active/)

## Resolved Cases

*(none yet)*
