---
source: klein-2026-07-23-S404
status: SYNTHESIS of the owner-fragment decode. The exact object is CERTAIN (eq (27) reconstructed
  bit-exact, independently by mac-mini-S168, opus-S3, and here). The SEMANTIC identification
  (n=13 Lonely-Runner "wider gap" / log-energy-beats-floor) is a STRONG hypothesis, not canon:
  the integers are source-generated (not paper-indexed), so no literal citation closes it. Two
  self-corrections to klein-S402 recorded below.
tags: [lrc, eq27, decode, second-moment, wider-gap, bedert, riesz, log-energy, n13, rapidity, tent-function]
---

# Eq (27) resolves to an n=13 "wider gap" — the second-moment fingerprint settles the scale

**klein-2026-07-23-S404.** Continues klein-S402/S403, opus-S2/S3, mac-mini-S168 and mac-mini's
20:20 "artanh cert" discriminators. The fragment is now fully characterized; this closes the
"which n" question mac-mini assigned to me.

## The object (CERTAIN, tri-verified)
```
(2457/6592)·log(8847357/2974400) − log(1285/896)  >  1/25
```
The snippet's 51-digit certificate equals `(2457/6592)·P₋(t_B) − P₊(t_A) − 1/25` **bit-exact**
(P± = the 3-term artanh sandwich). Value `X = 0.0457246…`, margin over 1/25 is `0.005725`.
Coeff on `log_A` is `d = 1` (opus-S3 p-adic pinning); coeff on `log_B` is `c = 2457/6592`.

## Two corrections to klein-S402 (own prior session)
1. **Baker route DOWNGRADED** (mac-mini's discriminator 1 is right). `X≈0.046` is a *genuinely
   large positive* quantity; `log_A/log_B = 0.3308 ≠ c = 0.3727`, so there is **no
   near-cancellation**. Baker / irrationality-measure endgames certify a linear form *near zero*;
   this certifies a weighted log-difference *large, above a floor*. Family = **loneliness /
   energy / entropy / LP-margin**, not transcendence. My S402 "Baker key hint" is retracted.
2. **"No small-coefficient fit" was wrong** — I searched *integer* coeffs; the fit is the
   *rational* `2457/6592` (mac-mini-S168). Retracted.

## The n-scale is n = 13 (resolves mac-mini's discriminator 3)
mac-mini flagged: if `1/25` means "beat trivial `1/(2n)`", the proven value `0.0457 ≈ 1/21.9`
looks like n≈11 (`1/22=0.0455`), clashing with my n=13. **The second-moment fingerprint breaks
the tie for n=13**, and it rests on opus's *proved* `d=1`:
- `2457 = 3·Σ_{k=1}^{13} k² = 3·819` (= 3·P₃(13)). **Only n=13**: 3·S₂(1..11)=1518, (1..12)=1950,
  (1..14)=3045 — none is 2457. The numerator is *exactly* 3× the **second moment of the tight-AP
  extremal velocity set `{1,…,13}`** — the L²/Parseval/variance quantity of the LRC(14) core.
- `91 = Σ_{k=1}^{13} k = C(14,2) = 7·13` (first moment / pair count) also sits inside 2457=27·91.
- Denominator `13²` in `2974400 = 2⁶·5²·11·13²`.
- The proven value `X=0.0457 ∈ (1/26, 1/14)` — exactly the admissible band for a **13-speed max
  loneliness** (classical `1/(2n)=1/26`, conjecture `1/(n+1)=1/14`).
- **Clean threshold:** `1/25 = 1/(2n−1)|_{n=13}`, beating classical `1/(2n)=1/26=1/(2·13)`. This
  *is* "a wider gap of loneliness" (Bedert arXiv:2511.16636, title; web-confirmed real paper). The
  value being `1/21.9` is not "just above `1/(2·11)`" — it is comfortably *in-band* for n=13 with
  the standard-shaped margin; mac-mini's n≈11 only follows if one forces value≈`1/(2n)`, which the
  fingerprint contradicts.

## Mechanism (why arctanh at all) — the tent/odd-harmonic bridge
`2·arctanh(t) = Σ_{k odd} t^k/k` is a **harmonic×geometric** (per-mode log-energy) sum on **odd**
frequencies. The Lonely-Runner tent function `‖x‖` (distance to nearest integer) has an
**odd-harmonic Fourier series** — so a per-mode log-energy of a tent/Riesz test object lands
*exactly* on the odd-only arctanh series (mac-mini's observation, sharpened). `log_A, log_B` are
then **logits/rapidities of two amplitudes** `p_A=1285/2181≈0.589`, `p_B≈0.748≈3/4`
(`log_B≈log 3` — a genuine ternary/log-3 lead worth a direct Riesz-amplitude test). Weight
`2457/6592` (numerator = 3× extremal second moment; `103 | 5872957 = x_n(B)`, so the weight is
*entangled with the data*, i.e. an optimizer output, not a hand-picked constant).

## Honest scope / non-repo alternatives (meta-caution, owner said "any problem")
The integers are **source-generated** (mac-mini: not web-indexed; the `11821757⁵` in the certificate
denominator is a `t⁵/5`-truncation artifact, not problem data). So no paper literally carries (27);
we identify the *problem family*, not a citation. Leading = **LRC(14) wider-gap, Bedert-Riesz /
log-energy method** (repo THM-518/THM-731, opus-S173 `LRCRieszCertificate`). Live non-repo
alternatives in the same "log-energy beats a rational floor" family, ranked: merit-factor /
Littlewood flat-polynomial autocorrelation; Cohn–Elkies LP margins; large-deviation/Chernoff rate
separations. All fit the *shape*; none yet fits the **13-second-moment + `1/(2n−1)` + odd-harmonic**
triple as tightly as LRC(14).

## The forward lead (turn decode → progress on the actual prize)
Whether or not the source is ever ID'd, the productive move: **use the fragment as a template.**
Run the certified-artanh **log-energy bound** on the LRC(14) tight cores `d·{1,…,13}∪{r}` and try
to emit a `> 1/(2n−1) = 1/25` "wider-gap" certificate of this exact shape, with weight tied to the
core's second moment `Σv_i²`. If a 13-speed Riesz/loneliness optimization reproduces
`(A,B)=(1285/896, 8847357/2974400)` (or their log-energy role), the identification is confirmed
*and* we have a new rigorous, Lean-ready certificate style for `inf L > 0` (OPEN-Q-097/104). This
is the concrete experiment; it advances the anchor either way.

→ THM-731 (second-moment Parseval floor), THM-518 (Bedert two-route), opus `LRCRieszCertificate`,
opus-S2/S3 certified log engine, mac-mini-S168 (exact form), klein-S402/S403 (superseded frame).
