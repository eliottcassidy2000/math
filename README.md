# Math Research: Parity in Tournaments

**Topic:** Rédei's theorem, the Tiling Model, and the Odd-Cycle Collection Formula
**Central open problem:** Claim A — H(T) − H(T−v) = 2Σμ(C)
**Latex paper title:** *Parity in Tournaments: The Q-Lemma, the Tiling Model, and the Odd-Cycle Collection Formula*

---

## For any new Claude instance starting a session

Read these files **in this order** before doing any mathematical work:

1. `01-canon/definitions.md` — core vocabulary; do not use any term differently
2. `01-canon/MISTAKES.md` — errors already made; do not repeat them
3. `00-navigation/OPEN-QUESTIONS.md` — the live frontier; pick up where others left off
4. `00-navigation/TANGENTS.md` — a dense index of rabbit holes; scan for jumping-off points
5. `00-navigation/SESSION-LOG.md` — who did what and when; add your own entry at the start

Then read any theorem file in `01-canon/theorems/` relevant to what you're working on.
Check `02-court/` for any active disputes before asserting a contested claim.

---

## Directory Guide

| Directory | Purpose |
|-----------|---------|
| `00-navigation/` | Orientation layer: tangents index, open questions, session log |
| `01-canon/` | High-certainty mathematics only. Definitions + proven theorems + MISTAKES log |
| `02-court/` | Formal dispute resolution. Back-and-forth argument documents between Claude instances |
| `03-artifacts/` | Preserved outputs: images, code snippets, LaTeX drafts |
| `04-paper/` | Future paper rebuild notes and structure |
| `inbox/` | Drop files here to contribute. Run `inbox/processor.py` to pre-process, then tell Claude |

---

## Instance Naming Convention

Every Claude contribution must be tagged with an instance ID in the format:

```
[ACCOUNT]-[DATE]-S[N]
```

Example: `A-2026-03-05-S1` means Account A, March 5 2026, first session that day.
Use this tag in SESSION-LOG.md and on any theorem or court document you touch.

---

## How to Add to Canon

A result enters `01-canon/theorems/` only if:
- It has a proof sketch or proof in this repo, OR a clear citation, AND
- Certainty level ≥ 4 (see `01-canon/README.md` for the scale)

Do not promote a verified-but-unproved result to canon without marking it CONJECTURE.

---

## How to Open a Court Case

If you find a disagreement between two positions (between Claude instances, or between a Claude and a previous version of itself):
1. Create `02-court/active/CASE-NNN-[topic].md` using the template in `02-court/README.md`
2. Log it in SESSION-LOG.md
3. Write the first letter arguing your position
4. The next Claude instance assigned to argue the other side continues the document
5. When one side concedes, move the file to `02-court/resolved/` and update the relevant canon theorem
