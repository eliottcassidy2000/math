---
source: klein-2026-07-23-S402
status: SPECULATIVE decode of an owner-supplied external fragment (a rigorous `log((1+t)/(1−t))` bound
  proving a quantity `> 1/25`). NOT canon — a construction/synthesis that revives the dormant rapidity/
  adelic thread with external evidence. The fragment's math object is CERTAIN (rapidity = formal log of
  the Cayley group); the LRC connection is STRONG (numeric + three converging repo threads); the "key
  hint" (a Baker / linear-forms-in-logarithms route to the LRC spectrum's Diophantine gaps) is the
  constructed hypothesis.
tags: [lrc, rapidity, formal-group, baker, linear-forms-in-logs, adelic, tao-optimistic, window, decode, speculative]
---

# Decoding a rapidity log-bound "> 1/25" — a Baker route to the LRC window

**klein-2026-07-23-S402.** Owner handed a fragment and a hunch: a computer-assisted step that bounds
`log((1+t)/(1−t))` two-sidedly by truncated odd-power series, evaluates at
`t_A = 389/2181` (`(1+t_A)/(1−t_A) = 1285/896`) and `t_B = 5872957/11821757`
(`= 8847357/2974400`), and concludes "the right side of (27) minus `1/25` is `≥ [51-digit]/[53-digit] > 0`."
Directive: decode it, pursue the repo hunch, treat it as a possible major key hint.

## What it IS (certain)

`log((1+t)/(1−t)) = 2·arctanh(t) = log( ((1+t)/2) / ((1−t)/2) )` — the **rapidity** = **log-odds** of
`p = (1+t)/2` = the **formal logarithm of the Cayley formal group** `F(x,y) = (x+y)/(1+xy)` (relativistic
velocity addition). Writing `t = (m−2k)/m` gives `(1+t)/(1−t) = (m−k)/k`, so each term is exactly the
repo's rapidity object

> `2·arctanh(λ) = ln((m−k)/k)`,  `λ = t = (m−2k)/m`   — **THM-252 (rapidity lattice)** verbatim.

Reverse-engineered facts (all exact): `p_A = 1285/2181 ≈ 0.589`; **`p_B = 8847357/11821757 ≈ 3/4`** so its
rapidity `≈ ln 3`. The proof's engine is the standard rigorous two-sided log bound
`2(t+t³/3+t⁵/5) ≤ log((1+t)/(1−t)) ≤ 2(t+t³/3+t⁵/(5(1−t²)))` — the exact device for certifying a **linear
form in logarithms** (Baker's theorem territory) in exact rational arithmetic.

Ruled OUT computationally: (27)'s RHS is **not** a simple `α·log_B+β·log_A+rational` (residual stays
~175-bit for all small α,β); `1/25` is **not** the lonely measure `L` (those run 0.006–0.024, tight = 0).
The denominators `q_A = 2181 = 3·727`, `q_B = 11821757 = 31·381347` (with a huge prime part) are far too
large for any small LRC extremal — the fingerprint of a **constructed / asymptotic / Diophantine family**,
where large heights are expected (again: Baker).

## Why it is (probably) OUR problem — three threads converge

1. **`1/25 = 1/(2·12+1)` and `2/25 = 2/(2·12+1)`.** These are the 12-velocity (n=13) Ostrowski unit and
   second rungs: `[0;12,k] = k/(12k+1)` gives `1/13` (k=1) and `2/25` (k=2). The interval `(1/13, 2/25)` is
   the **empty second-value window = Tao's optimistic conjecture** (Goddyn–Wong 2006), a live repo thread
   (CONSTANTS-INDEX; the `3/38` mediant is its depth-minimal open value). A rigorous bound "`> 1/25`" is a
   statement at the natural scale of that window.
2. **THM-252 + the rapidity/adelic thread.** `arctanh(λ) = ½ln((m−k)/k)` lies in the log-prime lattice
   `⊕_{p≤m−1} Q·ln p`, with **Baker's theorem giving Q-independence**. The `rapidity-supersingular-adelic`
   reflection (oracle-S19/opus-S91) upgrades this to an **adelic product formula** `|gap|_∞·|debt|_p = 1`:
   the archimedean (loneliness) and p-adic (covering-debt) sizes are reciprocal, so the counterexample
   corner `(gap=0, debt=0)` is forbidden — *"the missing analytic content is the archimedean lower bound."*
3. **The LRC continued-fraction frontier.** The spectrum is Ostrowski rungs `k/((n−1)k+1)` inside Farey
   cells; windows are gaps between consecutive rungs (THM-622 spectral-gap → Farey numerator/depth). The
   fragment's ugly high-height rationals are precisely the "unusually good approximants" HYP-3114 flags.

## The constructed key hint

> **The LRC spectrum's arithmetic — window emptiness (Tao's optimistic conjecture), the Ostrowski-rung
> structure, the covering-min `14/183` — is plausibly governed by Baker-type effective transcendence
> bounds (linear forms in logarithms) on rapidities in the log-prime lattice, glued by the adelic product
> formula. The fragment reads as a piece of exactly such a bound: a certified archimedean rapidity
> separation `≥ 1/25` that supplies the "missing archimedean lower bound" of the adelic proof shape and,
> via `|gap|_∞·|debt|_p = 1`, excludes the counterexample corner.**

If true this is a *new toolkit*: LRC hardness is transcendence-theoretic (effective linear forms in logs),
not only combinatorial/Fourier/covering — an angle the repo's rapidity thread has circled (THM-252,
HYP-1992, the adelic reflection) but left dormant, now with external evidence that a real proof uses this
device. The honest caveat: `2·arctanh` is generic (Fisher z, special relativity, entropy/log-odds), and
without equation (27)'s full form or the source paper I cannot certify the identification — but the
convergence (rapidity form + `1/25=1/(2·12+1)` + adelic proof-shape + CF frontier) is strong.

## Concrete next steps (to confirm / exploit)

- **Recover (27).** Ask the owner for eq. (27) and 1–2 lines of context; I can then verify the inequality
  independently and pin whether it bounds a window, a covering-min, or a rapidity separation.
- **Test the mechanism.** Formalize "rapidity separation ⇒ M-spectrum gap" for small n: does a Baker
  lower bound on `|Σ b_i ln p_i|` over the config's log-prime lattice reproduce the `(1/13,2/25)` gap?
- **Adelic closure.** Pair an archimedean Baker bound with the p-adic debt bound (the n=16 reflection's
  proof shape (P)) — the product formula is the glue; a certified archimedean floor is what was missing.
- **Revive HYP-1992 / THM-252** as an active LRC lane, not a curiosity.

Files: this reflection; decode scripts run inline (rapidity/CF/factor analysis, lonely-measure ruleout).
→ THM-252 (rapidity lattice), HYP-1992 (LRC rapidity bridge), `rapidity-supersingular-...-adelic-s19`,
LRC14-CONTINUED-FRACTION-FRONTIER, CONSTANTS-INDEX (`2/25`, `3/38`), Tao's optimistic conjecture thread.

---

## SOURCE NARROWED (klein-S403, 2026-07-23) — it is the Riesz-product "wider gap" line

Checked the recent LRC literature against the fragment. **Ruled OUT** (no logs / no 1/25 / no matching
integers in the main text): Rosenfeld "eight runners" (arXiv 2509.14111, product+prime-divisibility),
Sungkawichai "eleven/twelve/thirteen" (2604.23906, sieve+polynomial — the paper our repo cites).
**Strong match on THEME:** the ANALYTIC gap-improvement line — **"Riesz products and the Lonely Runner
Conjecture: A wider gap of loneliness"** (arXiv 2511.16636, Nov 2025; Thm 1.3: `ML(V) ≥ 1/(2n)+1/n^{5/3+o(1)}`,
Riesz-product test measures `R(x)=∏(1−p cos 2πmx)`) and its relatives. Its eq. (27) is a model-reduction
inequality `ML(V) ≥ ML(V'') − 1/m^ε − q/p` (a DIFFERENT (27) — no logs), so the fragment's (27) is a
different paper's / appendix's numbered equation, but the METHOD is this one.

**Interpretation corrected.** In the Riesz-product / entropy method, `log((1+t)/(1−t))` is NOT primarily
"rapidity" — it is the **binary-entropy / product-log derivative** (`d/dt[ H((1+t)/2) ] `-type term) that
appears when estimating `∫ Φ·R` or the log of a Riesz product. (The rapidity thread THM-252/HYP-1992 uses
the *same* function, so the two coincide — a genuine unity worth noting, but the fragment's home is the
Riesz/`inf L>0` analytic method.) The `> 1/25` is a **wider-gap statement**: for `n≈13`, `1/(2n)=1/26≈0.0385`,
so proving a Riesz lower bound on `ML(V)` `> 1/25=0.04` beats the trivial union bound — exactly the "wider
gap of loneliness" claim.

**Repo connection (direct):** opus-S173/S174/S178 (`LRCRieszCertificate.lean`, dissociated looseness
uniformly Riesz-certifiable, `sup inf_R ∫MR/∫R ≤ 0.55 ⇒ inf L>0`) and **THM-518 (Bedert two-route
diagnosis)** are the repo's internal development of exactly this method. **The key hint is now sharper:**
the external Riesz "wider-gap" papers and the repo's opus Riesz thread are the same programme; the fragment
is a concrete `ML(V) > 1/25` verification from that programme, and the repo should mine 2511.16636 (+ its
appendices/relatives) for the technique and cross-check against `LRCRieszCertificate`. Still unpinned: the
exact paper carrying the fragment's (27)+integers (likely an appendix / a companion small-`n` verification).
