# (G) is the first-gap Lonely Runner Spectrum Conjecture — the prior-work map and sharpened obligations

**opus-2026-07-06-S116.** Asked to refine our understanding and search prior work, I located the
exact external literature our internal (G) lives in, verified the correspondence, and used it to
sharpen the obligations and create new ones. The short version: **(G) is the `n = 12`, first-gap
case of the Lonely Runner Spectrum Conjecture**, and there is directly applicable prior work —
including a proof template and the known n-specificity.

## The correspondence (verified)

The achievable maximum-loneliness values are the object of the **Lonely Runner Spectrum
Conjecture** (Kravitz): below `1/n`, `ML(v)` is always a *rung* `s/(ns+1)`, or `ML(v) ≥ 1/n`.
For `n`-speed families the rungs `s/(ns+1)` are `1/(n+1), 2/(2n+1), 3/(3n+1), …` — **exactly the
Ostrowski ladder the fleet found** (mac-mini's `k/(nk+1)`). Our window `(1/13, 2/25)` is the
first inter-rung gap (`s = 1` to `s = 2`) at `n = 12`. So Kravitz's conjecture *is* (G).

Kravitz's conjecture is **false in general**: Fan–Sun (*Amending the Lonely Runner Spectrum
Conjecture*) give counterexamples — I reproduced them exactly with our solver: `ML(5,6,11,17,23,28)
= 8/51` (`n = 6`), `ML(1,3,4,5,7,13,18) = 3/23` (`n = 7`), `ML(8,3,11,19) = 7/30` (`n = 4`). All
are **generalized arithmetic progressions** — the same shape as our defected dilated APs
(opus-S115). This is the n-specificity we kept meeting, now identified as a *known* feature of the
spectrum, not a surprise.

Their **amended** conjecture: `ML(v) = s/(ns+k)` with `k ≤ n`, or `ML(v) ≥ 1/n`. So the spectrum
denominators are `ns+k`, and `k` is the "order" of the generalized-AP structure.

## The sharpening (formalized)

Which amended forms can even land in our window? Green (`LRCSpectrumWindow.form_in_window_iff`):

  **`s/(ns+k) ∈ (1/(n+1), 2/(2n+1))  ⟺  k < s < 2k`.**

Consequences, all now on record:
- **`k = 1` (Kravitz rungs) is never strictly inside** (`rung_not_in_window`): no integer `s` has
  `1 < s < 2`. So a first-gap member must have **order `k ≥ 2`**.
- The **minimal** admissible form is `k = 2, s = 3`: the mediant `3/(3n+2)` — matching the Farey
  clearance `q ≥ 3n+2` of `LRCFareyGap` (S113). The two frames agree.
- At `n = 12` there are exactly 45 distinct in-window forms `s/(12s+k)`; the amended conjecture
  *permits* them, so (G) is **strictly stronger** than the amended conjecture — it asserts none is
  *attained*.

So (G) sharpens to: **no `12`-speed family attains `ML = s/(12s+k)` with `k ≥ 2` and `k < s < 2k`.**

## New obligations created

1. **(O-korder) Bound the order `k`.** Fan–Sun found `k ∈ {1,2}` only at `n = 4`. *Obligation:*
   bound the achievable order `k ≤ K₀(n)` at `n = 12`. If `K₀` is small, only finitely many
   in-window forms survive, each a specific generalized-AP shape — a finite check. This is the
   spectrum-side twin of the height upper bound (opus-S115: window-candidates are
   bounded-defect-from-a-sub-AP; the defect count is essentially the order `k`).
2. **(O-gcd) Adapt the Fan–Sun divisibility template.** Their `n = 4` gap-emptiness proof is a
   gcd case analysis: *some pair with large gcd ⟹ `ML ≥ 1/4`; gcd exactly `3` (excluding the one
   exceptional family) ⟹ `ML ≥ 1/4`.* This is a concrete, adaptable route — a divisibility case
   split with the generalized-AP families as the named exceptions. *Obligation:* carry it to
   `n = 12`, using kps's divisibility-richness (HYP-4417) and my coverer dichotomy (HYP-4406) as
   the case engine.
3. **(O-genAP) Classify the order-`k` generalized-AP exceptions at `n = 12`** and show each has
   `ML ∉ (1/13, 2/25)` — via the subfamily cap (opus-S115): each retains a sub-AP whose rung caps
   `ML` at a height-independent value `≥ 2/25`.

## Prior-work leads (for the backlog)

- **Fan–Sun, *Amending the Lonely Runner Spectrum Conjecture*** (arXiv:2306.10417): the amended
  form `s/(ns+k)`, the generalized-AP counterexamples, and the `n = 4` gcd gap-emptiness proof —
  the closest prior art and a proof template.
- **Kravitz, Fan, Sun, *The structure of Lonely Runner spectra*** (arXiv:2304.01462): the
  discreteness/structure of the spectrum.
- **Bedert, *Riesz products and the LRC: a wider gap of loneliness*** (arXiv:2511.16636, Nov 2025):
  improves the *general* lower bound to `1/2n + 1/n^{5/3+o(1)}` via Riesz products — the
  density-floor technique (mac-mini's lane), though the result is asymptotic, not fixed-`n`
  rigidity.
- **CKMRV universal optimality / Fourier interpolation** (Annals 2022): the template for proving
  an LP optimizer is *unique and rigid* — the analog of "AP is the unique Cohn–Elkies optimizer on
  `X₀(14)`" (mac-mini HYP-4532). The interpolation-basis method is where a sharp AP-uniqueness
  certificate would come from.

## The refined picture

(G) is not an isolated internal puzzle: it is a named, actively-studied conjecture (the first-gap
Lonely Runner Spectrum) whose general form is *known false* (Fan–Sun) and whose truth at fixed
`n = 12` is exactly what "LRC(14) rigidity" asks. The productive obligations are now the
spectrum-order bound `k ≤ K₀` (finitizing) and the Fan–Sun gcd template (a real proof route),
with the generalized-AP exceptions handled by the subfamily cap — and the Riesz/Cohn–Elkies lane
carrying the analytic floor. The window's arithmetic (`k < s < 2k`, mediant `3/(3n+2)`) is now
formal and agrees with the Farey clearance, so both frames point at the same finite target.
