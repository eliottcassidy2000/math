# Vieta, Schanuel, and where the two arithmetics split the primes

*mac-mini-2026-07-20-S138. The owner asked me to hold Schanuel's conjecture, the Vieta
quadratic argument for `e` and `π`, and the repo's multiplication/addition duality together,
and to understand what the repo has actually done about rationals and irrationals. The
honest answer has one sharp new fact in it, one dichotomy sitting inside my own theorem from
two sessions ago, and a lot of framing that should be labelled as framing.*

---

## The Vieta argument, and why it is the duality

If `e+π` **and** `eπ` were both algebraic, then `e` and `π` would be the two roots of

`x² − (e+π)·x + (eπ)`,

a quadratic with algebraic coefficients — hence both algebraic, contradicting Hermite and
Lindemann. So **at least one of `e+π`, `eπ` is transcendental.**

The shape is worth naming: **the disjunction is a theorem and each disjunct is open.** We know
one of `5.8598744…` and `8.5397342…` is transcendental and cannot say which. Schanuel's
conjecture would give both, plus the algebraic independence of `e` and `π`.

And the argument *is* the addition/multiplication duality: Vieta's formulas say the sum and
the product are the two coefficients, and the proof is exactly that **you cannot have both the
additive and the multiplicative combination be tame while the elements are wild.** Schanuel is
the same duality one floor up — `ℚ`-linear (additive) independence of the `z_i` forces the
exponentials `e^{z_i}` to contribute new algebraic content, with `exp` as the bridge. That is
the role `log` plays between the repo's two arithmetics.

## Vieta is already how this repo reads spectra

Char-poly coefficients **are** the elementary symmetric functions of the spectrum. THM-1440's
`n=7` cospectral pair has `p(x) = x(x²+7)(x⁴+14x²+17)` with integer coefficients
`[1,0,21,0,115,0,119,0]`, while its eigenvalues are `0, ±√7, ±√(7−4√2), ±√(7+4√2)` —
irrational, splitting field `ℚ(√2)`. Verified: every elementary symmetric function of those
irrational numbers is an exact integer.

So the repo lives on the **algebraic branch** of the Vieta dichotomy: all symmetric functions
tame, elements irrational but still algebraic. `e`/`π` is the other branch, where tameness of
all the symmetric functions is impossible.

## A rational/transcendental dichotomy inside THM-1520

Two sessions ago the saddle lemma produced `L(pᵐ)/(a_Dᵐ(Dm)!) → exp(a_{D−1}/(D·a_D))`. With
rational (or algebraic) coefficients the exponent is rational (algebraic), so **Lindemann**
applies:

> the saddle limit is **rational ⟺ it equals exactly 1 ⟺ `a_{D−1} = 0`**; otherwise it is
> transcendental.

| `p` | exponent | limit | |
|---|---|---|---|
| `v`, `v²`, `2v³+v` | `0` | `1` | **rational** |
| `v − 1` | `−1` | `1/e` | transcendental |
| `v + 3` | `3` | `e³` | transcendental |
| `v² − 3v + 2` | `−3/2` | `e^{−3/2}` | transcendental |

The only analytic constant this whole GMC thread produced is `e^{rational}`, rational exactly
on the codimension-one locus `a_{D−1} = 0`. Lindemann does the work; nothing here is new
transcendence.

## The sharp fact: the two arithmetics split the primes

`H` is multiplicative under ordinal sum, so `log H` is additive, and the natural
"transcendence-flavoured" question — are these logs independent? — has an **elementary**
answer here. `{log p : p prime}` is `ℚ`-linearly independent by **unique factorization**
(clear denominators; `∏pᵢ^{nᵢ} = 1` forces every `nᵢ = 0`). Hence

> `dim_ℚ span{log H(T)} = #{distinct primes dividing some H(T)}` — measured `1, 2, 4, 12, 30`
> at `n = 3..7`.

And now the fact I did not expect:

> **Rédei says every `H` is odd. So `2` never divides any `H`. The multiplicative side of the
> repo contains no factor of 2 at all — while the additive side is `F₂^m`, hence entirely
> 2-adic. The two arithmetics live on disjoint sets of primes: 2 is purely additive, the odd
> primes purely multiplicative.**

That is not an analogy; it is a consequence of Rédei plus the definition of the tiling
hypercube, and it explains *why* the repo's two arithmetics have stayed genuinely independent
rather than collapsing into one. It also reframes Rédei's theorem: **oddness of `H` is exactly
the statement that the multiplicative arithmetic avoids the additive one's only prime.**

## The repo's constants, on the rationality ladder

CLAUDE.md names four: `√2`, `π`, `e`, `γ`. They sit at four different rungs:

| constant | status |
|---|---|
| `√2` | irrational, **algebraic** (proved, ancient) |
| `e`, `π` | **transcendental** (Hermite 1873, Lindemann 1882) |
| `γ` | **not even known to be irrational** |

So the repo's own constant list is the rationality ladder of number theory, with `γ` at the
open end — and `γ` enters through the `Γ(1/b)^b ~ b^b e^{−γ}` asymptotic, i.e. through exactly
the Gamma-function machinery the fiber fraction `(½)_k/k!` lives in.

## What is a theorem and what is framing

I want this separation attached, because the temptation here is to over-read.

**Theorems.** The Vieta disjunction for `e`,`π` (classical). Char-poly coefficients are the
elementary symmetric functions (canon). The saddle-limit dichotomy (Lindemann + THM-1520).
`dim_ℚ span{log H} = #distinct primes`, and `2` is absent by Rédei (unique factorization).

**Framing only.** "Schanuel is the repo's two arithmetics at the transcendence level." That is
a way of seeing, not a result. **Nothing in this repo bears on Schanuel**, and no repo constant
is known to be transcendental for any reason that is not already Lindemann. The `log H`
independence looks like Baker's theorem but is strictly weaker and strictly easier — unique
factorization, not transcendence — and calling it a shadow of Baker would flatter it.

The one thing I would actually chase: **is there any repo quantity whose transcendence is not
already Lindemann?** Everything found here reduces to `log(integer)` or `e^{rational}`. A
genuinely new one would have to be a ratio of periods or something similar, and I did not find
a candidate.

---

*Cross-links: THM-1440 (the cospectral pair and its splitting field), THM-1520 (the saddle
lemma), THM-1460 (`log Σa` vs `log H`, the two arithmetics under ordinal sum),
`07-reflections/two-charts-on-one-trichotomy-macmini-S137.md`, death-star-S60 (the two
arithmetics), CLAUDE.md (the four constants, the fiber fraction). Artifacts:
`04-computation/vieta_transcendence_threads_macmini_S138.py` (+out).*
