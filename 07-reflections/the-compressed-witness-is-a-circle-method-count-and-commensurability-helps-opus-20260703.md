# The compressed witness is a circle-method count — and commensurability (which covering forces) helps, not hurts

*opus-2026-07-03-S56. The owner asked me to close the crux — mac-mini's lemma (A), the last gap: "at a
modulus `q` with no small resonance among the speeds, a witness `a` exists (all `vᵢa` avoid the danger band
mod `q`)." I did not close it, but I gave it its exact circle-method form, verified the error kernel is a
Dedekind sum, and corrected the intuition: commensurate speeds HELP. Honest scope throughout — this is a
tool and a correction, not a proof.*

## The exact identity (the content of (A))

For prime `q` and speeds `V = {v₁,…,v₁₃}` coprime to `q`, the number of witnesses
`N(V,q) = #{a ∈ (0,q) : ∀i, vᵢa mod q ∉ danger}` is, by finite Fourier (danger `D`, `δ = |D|/q ≈ 1/7`):

```
   N = q · Σ_{(h₁,…,h₁₃): Σ hᵢvᵢ ≡ 0 (mod q)} ∏ᵢ c(hᵢ),     c(0) = 1−δ,  c(h) = −d(h) (h≠0),
```
where `d(h) = (1/q)Σ_{y∈D} e(−hy/q)` is the Dirichlet kernel of the danger interval, `|d(h)| ≤ 1/(2|h|)`.
- **Main term** (all `hᵢ = 0`): `q(1−δ)¹³ ≈ q(6/7)¹³ = 0.135 q`. This is mac-mini's lemma (i) — the mean,
  positive and unconditional.
- **Error**: the "resonances" `Σ hᵢvᵢ ≡ 0` with some `hᵢ ≠ 0`. **Lemma (A) is exactly: the error does not
  swamp the main term at a good `q`.**

## The pair error is a Dedekind sum (verified)

The dominant error is the pairs `|S|=2`: for `i<j` with direction `r = vⱼvᵢ⁻¹ mod q`, the contribution is
`(6/7)¹¹ · Σ_h |d(h)| |d(−hr mod q)| ≈ (6/7)¹¹/4 · (1/q) Σ_h 1/(h·‖hr/q‖)`. That last sum is a **Dedekind
sum**, and (verified, `r ~ q/φ` badly-approximable): `Σ_h 1/(h‖hr/q‖) ≈ 1.5·(log q)²` — so **`O((log q)²)`
for no-small-resonance directions**. Against the `O(q)` main term this is negligible once `q ≳ (log q)²`,
i.e. for `q` at the log-census scale. So the pair error is controlled exactly when the directions are badly
approximable — which is what "no small resonance" means.

## The correction: commensurability HELPS

I (and the naive picture) expected "resonance = bad." **The verification says the opposite.** A *small*
resonance is a *commensurate* pair (`vⱼ = 2vᵢ` etc.); its two bad-sets **overlap**, which *shrinks* their
union and leaves **more** witnesses — `N > q(6/7)¹³` whenever a commensurate pair is present (verified: for
`V` containing `18=2·9`, `312=2·156`, `N` runs `+8.5, +3.8, +4.4, +5.9, …` above the main term at the primes
where the alignment is constructive). And a **covering family shares small factors** `2,3,…`, so it *carries*
helpful commensurability. The genuinely dangerous ("no witness") case is the reverse: the 13 bad-sets must
*anti-align* into a near-**tiling** of `[0,q)` (they must, since `13·|D| ≈ 13q/7 > q` forces `6q/7` of forced
overlap to instead be spent covering) — a special adversarial arrangement, not what resonance gives.

## Where this leaves (A) — honest

- **Contributed (verified):** the exact circle-method form of (A) (`N = q(6/7)¹³ + Σ_resonances`), the
  pair-error = Dedekind sum `O((log q)²)` for badly-approximable directions, and the commensurate-helps
  correction. `N > 0` on every free prime tested for the compressed shape.
- **NOT closed:** (A) proper — that at some `q ≤ Q(M) = O(log M·log log M)` the error is below the main term.
  The residual is exactly bounding how often the bad-sets can *adversarially tile* (`f_q < 1` fails), which is
  **mac-mini's capacity argument** (HYP-4054): each blockable modulus costs `log(1/f_q)` bits and CRT caps the
  total at `13 log M`, so a good `q` exists at `O(log M log log M)`. My contribution supports it by making the
  per-`q` error explicit (Dedekind sum) and showing covering pushes the error the *helpful* way.

So the crux's last mile is a geometry-of-numbers / Dedekind-sum estimate at the log-census scale — mac-mini's
active lane. This tool gives it the exact analytic form and a verified error kernel; it does not finish it,
and (given two recent overclaims, MISTAKE-097) I am flagging that plainly.

## Status

- **Tool + verification (correct):** `circle_method_witness_count_opus_20260703_S56.py` — the circle-method
  witness count, the Dedekind-sum error kernel, the commensurate-helps correction.
- **Open (mac-mini's lane):** lemma (A) at `q ~ log M` = the capacity / adversarial-tiling bound.
- **My landed pieces feeding the dispatch:** the routing (HYP-4050), the dominant/compressed cluster-ID
  (HYP-4054-opus), the peel engines (THM-608/`scale_separation(_phase)`), all kernel-pure.

Related: mac-mini `the-compressed-crux-is-the-log-census` + HYP-4054-macmini (capacity), HYP-4040/4051
(band-blockers), kps `lonely14_of_ratio` (the census), MISTAKE-097 (my two prior overclaims on this crux).
Script: `04-computation/circle_method_witness_count_opus_20260703_S56.py`.
