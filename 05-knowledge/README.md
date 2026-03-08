# 05-knowledge: Persistent Knowledge Web

This directory preserves ALL computational results, tested hypotheses, and variable relationships permanently. Nothing discovered should ever be lost.

## Structure

### `variables/`
Variable registry. Each file documents one mathematical quantity:
- What equations it appears in
- Values at each n (with source script)
- Known relationships to other variables
- Cross-links to hypotheses and theorems

**Index:** `variables/INDEX.md` — alphabetical lookup with one-line descriptions.

### `hypotheses/`
Every hypothesis ever tested, whether confirmed, refuted, or open.
Each file includes: statement, test method, outcome, and crucially WHY it succeeded/failed.

**Index:** `hypotheses/INDEX.md` — searchable by topic, status, and related variables.

### `results/`
Raw computational outputs from scripts. Every script in `04-computation/` should have a corresponding results file here with the actual output.

**Index:** `results/INDEX.md` — maps script names to result files.

## Conventions

1. **Never discard a result.** Even negative results ("X is NOT true") are valuable.
2. **Always link.** Every entry should cross-reference related variables, hypotheses, theorems, and scripts.
3. **Record the HOW.** When a hypothesis fails, document what specific pattern breaks it and at what n.
4. **Use tags.** Every entry has searchable tags for grep-based discovery.
