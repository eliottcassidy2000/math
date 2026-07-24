# The artanh / Abel–Dini separation certificate — a decode, and the missing proof layer

*opus-2026-07-23-S2. Decoding an external snippet the owner supplied as a possible
"key hint." Script `04-computation/certified_logratio_abeldini_opus_S2.py`. See
HYP-9023, and the thread it belongs to: THM-2000 (support-harmonic Abel–Dini).*

## What the snippet is (machine-confirmed)

A **rational separation certificate**: `log((1+t)/(1-t)) = 2·artanh(t)` is trapped
between two rationals (lower = truncate after `t⁵/5`; upper = close the tail
geometrically, `1/(2m+1) ≤ 1/5`), evaluated at engineered rational `t`, to prove a
transcendental quantity exceeds `1/25` — reduced to a decidable rational inequality
witnessed by a 50-digit positive fraction `G`.

Two facts pin its home in the repo, both verified exactly:

1. **The `t`'s are Abel–Dini telescoping ratios.** `t_A = 389/2181`,
   `t_B = 5872957/11821757` are *exactly* `t_n = x_n/(S_n + S_{n−1})` for the
   partial-sum pairs `(S_{n−1},S_n) = (896,1285)` and `(2974400,8847357)`; then
   `(1+t)/(1−t) = S_n/S_{n−1}` and `2·artanh(t_n) = log(S_n/S_{n−1})`. This is
   **verbatim the construction in [THM-2000](../01-canon/theorems/THM-2000-support-harmonic-abel-dini-figurate-surface.md)
   §3.1** (the Abel–Dini partial-sum edge, `Σ log(S_n/S_{n−1}) = log(S_N/S_0)`).
2. **The certificate's denominator is assembled only from the two bounds:**
   `den(G) = 2⁸·3⁴·5²·7·31⁵·257·727³·381347⁵ = lcm(den(Lo_B), den(U_A))`
   (since `11821757 = 31·381347`, `2181 = 3·727`). No foreign primes.

So `RHS(27) = (+)c·log(S_B/S_{B−1}) − d·log(S_A/S_{A−1}) + rational`, with the larger
ratio lower-bounded and the smaller upper-bounded (the directions that lower-bound a
"+big − small" form). The exact `(c,d)` need the source's eq (27) — there is **no**
small-coefficient exact fit (the rational part is large, from prior algebra), so the
snippet under-determines them. What is *not* speculative: it certifies an Abel–Dini
log-combination `> 1/25`.

## Why this is a hint and not a footnote

THM-2000's payoff is a table of **transcendental** support-harmonic masses
(`2 log 2`, `18 − 24 log 2`, `3 log 3 − π/√3`, `π²/6`, `4π²/3 − 12`, …) and **strict
inequalities/non-coincidences** among numerically-close ones (the mass chains
`σ(G₁) > σ(G) > σ(Fib)`, the ladder trichotomy). But right now those strict claims
are only **numerical referees**: the script checks `|value − 2 log 2| < 10⁻⁴⁵` in
mpmath, and the Lean kernel `SupportHarmonicFigurate.lean` *deliberately* formalizes
only the finite rational partial-fraction identities — it says outright that "infinite
sums, Abel–Stieltjes integration … remain in the paper proof." Its one inequality
lemma, `reciprocal_block_bounds`, is the **discrete** block sandwich.

**The snippet is the missing continuous sibling.** `2·artanh` bounds on
`log((1+t)/(1−t))` are the certified log-inequality layer that upgrades those
transcendental orderings from float-checked to *proved* — and puts them in reach of
Lean without axiomatizing `log`.

### It works — a real ordering, certified float-free

`M(6,2) = 2 log 2` (hexagonal) vs `M(4,3) = 18 − 24 log 2` (square-pyramidal Faulhaber):
`M(6,2) > M(4,3) ⟺ 26 log 2 > 18 ⟺ log 2 > 9/13`. The identical 3-term artanh sandwich
at `t = 1/3` gives `log 2 ≥ 842/1215`, and `842/1215 > 9/13`, so
`M(6,2) − M(4,3) ≥ 22/1215 > 0` — a genuine THM-2000 mass ordering proved with pure
rationals. The reusable engine is `log_lower/log_upper(P,Q)` in the script.

## The key-hint reading (speculative, flagged)

The *shape* of the snippet is **separation certification**: trap an analytic quantity
between rationals to prove it lands the right side of a threshold. That shape is not
confined to THM-2000 — it is also the exact shape of the repo's biggest analytic
prize, **LRC(14): `inf_S L(S) > 0`** (THM-501/503/523), where a covering measure must
be certified above a rational floor (`1/1260`, `2/29`, …). The transfer is *structural,
not literal*: `L(S)` is a sinc-lattice sum (`sin(πt/7)/(πt)`), not a log, so `artanh`
does not apply to it directly. But the meta-lesson stands: **the repo's analytic
frontier is a stack of "real number vs rational threshold" certificates, and it has
been proving them numerically. The artanh/Abel–Dini sandwich is the first clean,
Lean-portable member of the certified-inequality toolkit those proofs need.** Building
that toolkit — starting with `logRatioSandwich` next to `reciprocal_block_bounds` — is
the actionable hint.

## Next

- (a) Finish the exact decode if the source's eq (27) surfaces — then name precisely
  which Abel–Dini quantity `> 1/25` is being certified (a mass separation? a
  partial-sum-vs-log gap? a growth/log-convexity gap of a counting sequence?).
- (b) Formalize `logRatioSandwich` in `SupportHarmonicFigurate.lean` and discharge one
  mass ordering (`M(6,2) > M(4,3)` is the cleanest first target) as a real theorem.
- (c) Audit the whole THM-2000 mass chain for orderings certifiable by ≤k artanh terms;
  each becomes a float-free proof.
