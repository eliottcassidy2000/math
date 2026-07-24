---
source: kind-pasteur-2026-07-23-S131 (Opus 4.8)
status: RECALIBRATION on new owner intel: the snippet is the CRITICAL MOMENT of an (outside-LLM) Lonely
  Runner Conjecture proof. This commits the LRC framing and draws the consequences. Two are decisive:
  (1) {1..13} is the WEIGHT/reference, not the certified configuration (it clears 1/25 trivially); (2) a
  fork the fragment cannot resolve — is 1/25 the (wider) gap bound, or a positivity SURPLUS in a
  full-conjecture argument? The soundness crux (kps-S130) is now the likely locus of innovation-or-flaw.
tags: [lrc, lonely-runner, snippet, critical-step, recalibration, soundness, wider-gap, n13, meta]
related: [THM-515, THM-518, kps-S129, kps-S130, klein-S404/405, macmini-S169]
---

# The snippet is the critical step of an LLM Lonely-Runner proof — recalibration

**kps-2026-07-23-S131.** Owner: the fragment "appears to be the critical moment in a Lonely Runner Conjecture
proof" from an outside LLM conversation. This upgrades the n=13 fingerprint from "one suggestive coincidence"
(my S129 caution) to **the intended reading**. Consequences, with the soundness lens from kps-S130.

## 1. {1,…,13} is the WEIGHT, NOT the certified configuration (decisive)
`L({1..13}, δ=1/25) = 0.329` (a third of the circle is lonely) and `gap({1..13}) = 1/14 ≈ 0.0714 > 1/25`,
attained at τ=1/14. So `{1,…,13}` clears `1/25` **trivially** — it cannot be the hard instance the snippet
certifies. Therefore `2457 = 3·Σ_{k=1}^{13}k² = 819·3` enters as the **second-moment normalization / reference**
of the argument (the extremizer the method is calibrated against), and the configuration actually being
bounded is separate. Corollary: since (if LRC holds) *every* 13-speed set has gap ≥ 1/14 > 1/25, proving
"> 1/25" for any *fixed* set is near-trivial — **the snippet's content is a GENERAL/uniform argument** valid
for arbitrary 13-speed V, i.e. an unconditional lower bound that does not assume the conjecture. (This is the
Bedert regime: beat the union bound `1/(2·13)=1/26` uniformly.)

## 2. The fork the fragment cannot resolve: gap-bound vs positivity-surplus
`X := (2457/6592)log_B − log_A = 0.045725` is a dimensionless weighted log-difference, proved `> 1/25 = 0.04`.
Two readings, NOT distinguishable from the fragment:
- **(a) X is the (wider) gap bound.** Then the theorem is `gap(V) ≥ 1/25 > 1/26` for 13-speed V — a modest
  UNCONDITIONAL improvement on the union bound, **NOT the full conjecture** (which needs `1/14`). Weak
  supporting hint: `0.0457 ∈ (1/26, 1/14)` lands squarely in the gap band, as a gap-quantity should.
- **(b) X is a positivity SURPLUS/margin** in an argument whose conclusion is the tight `gap ≥ 1/14`. Then
  "> 1/25" just means "safely positive," and the result COULD be the full conjecture (or its n=13 case).
**This fork is THE question for judging the LLM's claim.** If the LLM asserted a full LRC proof and X is
reading (a), the proof is *incomplete at exactly this step* (reaches 1/25, not 1/14). If (b), it may be real.
Recovering the sentence around eq (27) — "therefore the gap is at least X" vs "therefore [surplus] X > 0, so
…" — settles it. Owner: that context is the single highest-value missing piece.

## 3. Soundness crux (kps-S130) is where an LLM proof of a famous problem most likely innovates OR breaks
The snippet certifies a **log-quantity directly**. My S130 result: a naive log-energy `∫M·log R` is NOT a
sound certificate for the loneliness MEASURE `|G|>0` (the `R≥0` identity `∫_G R = ∫R−∫_B R` has no signed-log
analogue; verified counterexample). So the LLM's step is sound **only** if it is one of:
- **(R2-valid) a free-energy / partition-function / large-deviation formulation**, where a log-rate IS the
  target (e.g. `(1/β)log ∫ e^{β·loneliness} → max loneliness` as β→∞; a soft-max lower bound on the gap).
  This is a legitimate and even elegant way to lower-bound a MAX by a log-rate. If the LLM found this, the
  step is real and novel.
- **(FLAW) a concavity/log-energy quantity mistaken for a loneliness certificate** — the exact trap S130
  isolates. High prior for LLM proofs of open problems.
**Deciding R2-valid vs FLAW is the cluster's decisive task.** The weight `2457/6592` being the config's
second moment is consistent with a free-energy 2nd-order (curvature) term (cf. klein-S405 "Σv² = ½E″(0)"),
mildly favoring the R2-valid reading — worth pursuing, not assuming.

## 4. Concrete next steps
1. **Reconstruct the free-energy functional.** Test the hypothesis `X = (1/β*)·[log-partition-function of the
   13-mode loneliness at inverse-temperature β*]`, with β* tied to `Σv²`. If a soft-max/annealed bound on the
   gap reproduces weight `2457/6592` and the two log-arguments A,B (as `∫e^{βψ}`-type values at two saddle
   points / two temperatures), the step is R2-valid and we have the argument.
2. **Or localize the flaw.** If every reconstruction needs the unsound `log-energy ⇒ |G|>0` step, report it as
   the gap in the LLM proof (with the S130 counterexample as witness).
3. **Owner ask:** the one clause around eq (27) stating what "> 1/25" concludes (a gap, or a surplus).

Files: `/tmp/lrc.py` (loneliness measures, gap of {1..13}), builds on kps-S130 soundness, klein-S404/405.
