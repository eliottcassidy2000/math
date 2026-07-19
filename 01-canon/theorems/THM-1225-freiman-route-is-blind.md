---
id: THM-1225
title: THE FREIMAN/DOUBLING ROUTE IS STRUCTURALLY BLIND — LRC IS NOT TRANSLATION-INVARIANT AND DOUBLING IS — the two-element extremal set cannot be characterised by any additive-combinatorial invariant of the standard kind, because those are translation-invariant while the LRC floor is not: {1,…,13} and its translate {2,…,14} share doubling exactly 1.923 (the minimum 25/13) yet have M = 1/14 and M = 1/8. The full translate ladder {1+k,…,13+k} for k = 0…4 holds doubling fixed at 1.923 while M runs 1/14, 1/8, 1/6, 1/5, 5/22. Conversely minimal doubling does not imply near-floor (AP d=3 and the odd set both have 1.923 with M = 13/38 and 1/2), and the second extremal {1,…,11,13,24} has NON-minimal doubling 2.769 — so the implication fails in both directions. What does separate them is DIVISIBILITY: q₀ takes the value 14 on BOTH extremal families and differs on every non-extremal tested, exactly as THM-1210's D = 1 sieve picture requires. Dilation preserves both doubling and M, consistent with THM-1050
status: exact — all values computed in rational arithmetic; the translate ladder and the doubling figures are exhibited explicitly. The negative is structural (an invariance mismatch), not a failed search. The Kakeya remark is an untested extension of the same invariance objection, flagged as such
source: opus-2026-07-19-S394 (owner: work the two-element extremal set gap, think Kakeya and Freiman)
depends_on: [THM-1220 (the two-element extremal set), THM-1210 (D = 1 is the sieve; q₀ as the governing invariant), THM-1050 (dilation invariance), THM-1185 (the earlier measure-vs-pointwise triage this parallels)]
scripts: 04-computation/freiman_doubling_opus_S394.py -> 05-knowledge/results/freiman_doubling_opus_S394.out
---

# THM-1225 — why additive combinatorics cannot see the floor

## The invariance mismatch

Freiman's theorem is the canonical inverse theorem for additive structure,
and klein's THM-1004/1005 are literally titled "the inverse/rigidity
theorem", so the analogy is inviting: the extremals are APs, APs have
minimal doubling, so perhaps *large doubling ⟹ M away from the floor*.

**It fails, and structurally.** Doubling is **translation-invariant**; the
LRC floor is **not**:

| family | doubling | M | q₀ |
|---|---|---|---|
| {1,…,13} | **1.923** | **1/14** | **14** |
| {2,…,14} | **1.923** | 1/8 | 15 |
| {3,…,15} | **1.923** | 1/6 | 16 |
| {4,…,16} | **1.923** | 1/5 | 17 |
| {5,…,17} | **1.923** | 5/22 | 18 |

Every translate has identical (and minimal) doubling while M ranges over a
factor of three. **Any translation-invariant invariant is therefore blind to
the floor**, and doubling — the Freiman invariant — is one.

## The implication fails in both directions

- **Small doubling ⇏ near floor.** The AP with d = 3 and the odd set
  {1,3,…,25} both have minimal doubling 1.923, with M = 13/38 and M = 1/2 —
  five to seven times the floor.
- **Near floor ⇏ small doubling.** The second extremal {1,…,11,13,24} sits
  exactly on the floor with doubling **2.769**, well above minimal.

## What does separate them

Divisibility. q₀ — the least modulus dividing no speed — takes the value
**14 on both extremal families** and a different value on every non-extremal
tested. That is precisely THM-1210's picture: extremality is D = 1 at
s = q₀ = 14, i.e. the classical sieve at modulus 14. And q₀ is **not**
translation-invariant, which is exactly why it can see what doubling cannot.

Dilation, by contrast, preserves both doubling and M — consistent with
THM-1050, and a useful check that the two invariants agree where they
should.

## The triage rule this yields

This parallels THM-1185's measure-versus-pointwise split and sharpens it:

> **A translation-invariant invariant cannot characterise the LRC extremal
> set.**

That rules out doubling, additive energy, sumset/Freiman dimension, and —
by the same objection, though I have not tested it — the direction-set
invariants of Kakeya type, which are invariant under the affine maps that
move the floor. The invariants that have worked in this programme are all
**arithmetic and translation-sensitive**: divisibility q₀ (THM-1105/1210),
the mod-7 residues of the Fourier kernel (THM-1075), and the active-pair
determinant D (THM-1205).

## Status

Exact and structural. It does not close the two-element gap — it removes a
family of approaches from consideration and says why, which is the same
service THM-1185 performed for measure-based methods.
