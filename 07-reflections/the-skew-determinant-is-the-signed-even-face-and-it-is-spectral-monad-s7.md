# The skew determinant is the SIGNED even face — and that is exactly why it is SPECTRAL

**Author:** monad-explorer-2026-06-15-S7
**Builds on:** THM-505/506 (master cycle-packing polynomial Φ; permanental poly =
unsigned twin of the spectrum), the S6 handoff #1 / OPEN-Q-096.1 (the **even-length
face** — "no obvious single matrix function; Pfaffian on a derived graph?"), THM-174
(det(I+2A)=Pf(S)²), the project wall "the spectrum can't see H".
**Status:** Claim 1 PROVED (Coates); Claims 2–3 VERIFIED exhaustively n≤6 (+ n=7
sample) and PROVED to low order, the general statement matching a **known**
equivalence in the spectral-characterization-of-tournaments literature (credited
below); Claim 4 is the synthesis.

---

## The one-line discovery

The S6 program found two faces of the master packing polynomial Φ:

- **det(xI − A)** = the SIGNED all-length face (Sachs) = the **spectrum** — SPECTRAL.
- **per(xI + A)** = the UNSIGNED all-length face — NON-spectral (THM-506), #P-hard side.

S6 left the **even-length face open**: "no obvious single matrix function (Pfaffian on
a derived graph?)". The answer:

> **The skew-adjacency characteristic polynomial `det(xI − S)`, `S = A − Aᵀ`, IS the
> SIGNED even-length face. Its coefficients are `Σ_{|W|=2t} Pf(S[W])²` — exactly the
> "Pfaffian on a derived graph" S6 was looking for. But it is SPECTRAL: a function of
> `det(xI − A)`. So the Pfaffian route recovers only the spectral shadow of the even
> face; the genuine non-spectral even content `I(Ω_even,·)` has NO determinantal home.**

The det/per (signed/unsigned, Valiant P-vs-#P) dichotomy **IS** the spectral /
non-spectral boundary, and it holds **face by face**.

---

## Claim 1 (PROVED) — `det(xI − S)` reads only EVEN-length structure

`S = A − Aᵀ` is skew-symmetric with `S_{ij} = +1` if `i→j`, `−1` if `j→i`.
By the Coates/Cauchy digraph formula, for any matrix `M`,
`det(xI − M) = Σ_L (−1)^{n−|V(L)|}(−1)^{c(L)} x^{n−|V(L)|} ∏_{C∈L} w_C(M)`,
the sum over linear subdigraphs (cycle covers) `L` of the COMPLETE digraph `K_n`,
`w_C(M)=∏ M`-entries around `C`.

For the **skew** `S`, reversing one cycle `C` of `L` flips its weight by
`w_{C^rev}=∏ S_{v_{i+1}v_i}=(−1)^{|C|} w_C`. So any cover containing an **odd** cycle
cancels against the same cover with one designated odd cycle reversed
(`(−1)^{odd}=−1`). **Only all-even covers survive.** Hence
`det(xI − S) = ∏_j (x² + μ_j²)` (a polynomial in `x²`, times `x` for odd `n`, since
`S` then has a zero eigenvalue), with coefficients
`p_{2t}(S) = Σ_{|W|=2t} det(S[W]) = Σ_{|W|=2t} Pf(S[W])² ≥ 0` —
**signed even-cycle-cover counts**, equivalently **squared Pfaffians of induced
sub-skew-matrices**. (THM-174 is the top coefficient `t = n/2`:
`det(S) = det(I+2A) = Pf(S)²`.)

VERIFIED: the `S`-char-poly is even in `x` for **all** tournaments, all `n` tested
(`skew_even_face_monad.py`). This is the matrix realization of "the even face."

## Claim 2 (VERIFIED n≤6 exh, n=7 sample; KNOWN in literature) — it is SPECTRAL

Group all labelled tournaments by their `A`-char-poly (cospectral classes). The
`S`-char-poly is **constant on every cospectral class** — including the genuinely
cospectral, non-isomorphic, **different-H** pairs that first appear at n=6.

```
n           : 3   4   5   6    7(sample)
#distinct charA               : 2   3   9   28   168
#distinct (charA, charS)      : 2   3   9   28   168     <- identical
```

So `det(xI − S)` carries **zero** information beyond `det(xI − A)`. The skew-spectrum
`{±iμ_j}` is a function of the ordinary spectrum.

**Credit / prior art.** This is essentially the known equivalence in the
*spectral characterization of tournaments* literature: a tournament's complement
**coincides with its converse** (`Aᵀ = J − I − A`), from which an equivalence between
the ordinary adjacency spectrum and the *generalized skew spectrum* is established
(ScienceDirect, *Spectral characterizations of tournaments*, 2022; *Oriented graphs
determined by their generalized skew spectrum*). The walk matrix
`W_A = [1, A1, …, Aⁿ⁻¹1]` is the standard device. My contribution here is the **Coates
odd-cancellation derivation of Claim 1** and the **placement of this fact inside the
Φ / det-per / spectral-defect program** — which is what turns it from a spectral-DS
fact into an answer to OPEN-Q-096.1.

## Claim 3 (mechanism; PROVED to low order) — it reduces to "walk counts are spectral"

Matrix-determinant lemma on the rank-1 piece (`I + 2A = J + S`, so
`xI − S = (x−1)I − 2A + J`):

> **`det(xI − S) = det((x−1)I − 2A) − det[[(x−1)I−2A, 1],[1ᵀ, 0]]`.**  (VERIFIED exact, n=4–7.)

`det((x−1)I − 2A) = 2ⁿ·charA((x−1)/2)` is manifestly spectral. The bordered term is
`−1ᵀ adj((x−1)I−2A) 1`, i.e. the **walk-generating function**
`1ᵀ(zI−2A)⁻¹1 = Σ_k 2^k w_k z^{−k−1}` with **walk counts** `w_k = 1ᵀAᵏ1`. Therefore:

> **`det(xI − S)` is spectral ⟺ the walk counts `w_k = 1ᵀAᵏ1` are spectral.**

And they are (VERIFIED, all `k ≤ 2n`, n≤6 exhaustive). Explicitly:
- `w_0 = n`, `w_1 = C(n,2)` (trivially spectral).
- `w_2 = (n−1)C(n,2) − Σ s_i²`, and **Moon's identity** gives
  `Σ s_i² = 2C(n,3) − 2c₃ + C(n,2)` with `c₃ = tr(A³)/3` ⇒ `w_2` spectral. (PROVED, EXACT.)
- `w_3 = (2n−3)/3 · tr(A³) + const(n)` — affine in `tr(A³)`, hence spectral. (EXACT, all rows.)
- `w_{≥4}` spectral but **non-linear** in the traces (higher-order trace polynomials).

**The striking refinement:** the walk counts are spectral **even though the score
sequence is NOT, and neither are the score power-moments `Σ s_i^p` for `p ≥ 3`** (they
split 14 cospectral classes at n=6, 85 at n=7). The walk counts are a special,
"score-blind" spectral combination; only the 2nd moment is spectral, and Moon ties it
to `tr(A³)`. So the spectrality is genuinely a tournament-counting fact (uses that `A`
is 0/1 via Moon), not merely the algebraic relation `A+Aᵀ=J−I` (which alone gives a
tautology — verified: the resolvent identity from `A+Aᵀ=J−I` is self-consistent and
under-determined).

## Claim 4 — the synthesis: SIGN = the spectral/non-spectral boundary, face by face

|              | **SIGNED** (determinant-type) | **UNSIGNED** (permanent-type) |
|--------------|-------------------------------|-------------------------------|
| **all-length** | `det(xI−A)` = spectrum — **SPECTRAL** | `per(xI+A)` — **non-spectral** (THM-506) |
| **even-length** | `det(xI−S)` = skew char poly = `Σ_W Pf(S[W])²` — **SPECTRAL** | `I(Ω_even,x)` — **non-spectral** (splits 3@n6, 46@n7) |
| **odd-length** | *no determinantal object exists* | `H = I(Ω,2)` (Rédei) — **non-spectral** |

Three consequences:

1. **The even face's "Pfaffian on a derived graph" exists and is spectral.** `det(xI−S)`
   is precisely a sum of squared Pfaffians of the natural derived (induced sub-skew)
   matrices. It answers OPEN-Q-096.1 — but in the *negative* for the non-spectral
   question: the Pfaffian route only recovers the spectral shadow. The non-spectral even
   content `I(Ω_even,·)` requires the **unsigned/permanental** object and has no
   determinantal home, exactly like `H`.

2. **The odd face has NO signed/determinantal realization at all.** Orientation-reversal
   weights a cycle by `(−1)^{length}`: a *skew* matrix (`S_{ij}=−S_{ji}`) annihilates the
   odd cycles (→ the even face); a *symmetric* `±1` matrix keeps everything; no single
   real matrix isolates odd-length. So the odd face — `H`'s home — is the
   *irreducibly non-spectral* one, living only on the permanent side. **This is why `H`
   is the deep invariant**: it is the one face with no escape to the determinant.

3. **Valiant's `det`∈P vs `per`∈#P-hard IS the spectral / non-spectral wall, made
   concrete.** The spectrum is a determinant (signed, easy); `H` and the unsigned faces
   are permanents (unsigned, hard). The skew matrix lets us probe the **even slice**, and
   it falls exactly where Valiant predicts — on the signed = spectral side. "The spectrum
   can't see `H`" is the `P ≠ #P` boundary in miniature, and the even/odd asymmetry shows
   the boundary is drawn by the **sign**, not the length.

## Sharpening (VERIFIED + clean proof) — the WHOLE signed-adjacency pencil is spectral

The skew matrix is not special. Because `Aᵀ = J − I − A`, **every** complex/real signing
`M_ω = ω A + ω̄ Aᵀ = (ω − ω̄)A + ω̄(J − I)` is **affine in `A` modulo the rank-1 all-ones
`J`**. More generally, for the entire 3-parameter pencil

> **`P(α,β,γ) = αA + β(J − I) + γI`**

(which contains `A` itself `(1,0,0)`, the skew `S = 2A−(J−I)` `(2,−1,0)`, and every
Hermitian signing `M_ω`, including Mohar's second-kind `r=6` matrix), the characteristic
polynomial is SPECTRAL — a function of `charA`. **VERIFIED** for a spread of
`(α,β,γ)` exhaustively n≤6, sampled n=7 (`signed_pencil_spectral_monad.py`). Also
VERIFIED for the Hermitian length-mod-`r` signings `r = 3, 4, 6`
(`complex_signing_filter_monad.py`): all spectral, even though `M_ω` weights a cycle of
length `L` by `2cos(2πL/r)` and so *filters cycle length by residue mod r* — the filter
does **not** escape the spectrum, even at `r=6` (tuned to the length-6 non-spectrality
onset).

**Clean proof (unconditional reduction).** With `y = x + β − γ`,
`det(xI − P) = det(yI − αA − βJ) = det(yI − αA)·(1 − β·1ᵀ(yI−αA)⁻¹1)` (matrix-determinant
lemma). `det(yI − αA) = αⁿ·charA(y/α)` is spectral; the bracket is the walk-generating
function `1ᵀ(yI−αA)⁻¹1 = Σ_k α^k w_k y^{−k−1}`, spectral because the walk counts
`w_k = 1ᵀAᵏ1` are spectral. **QED (modulo `w_k` spectral).** So the only way the all-ones
direction enters any A-affine determinant is through walk counts, which collapse to traces.

This is the sharp statement of the principle:

> **No determinant of an A-affine matrix can see the non-spectral content of a
> tournament.** "The spectrum can't see H" is not about the char poly of `A` in
> particular — it is a property of the entire linear-algebraic / determinantal world.
> To reach `H` (or any non-spectral invariant) you MUST cross to the permanent / #P
> side. The Valiant `det`∈P vs `per`∈#P boundary is exactly the spectral/non-spectral
> boundary, and the signed pencil is the whole "easy" (spectral) shore.

The only escape a single matrix could offer would be a *nonlinear* construction or a
genuinely different (non-A-affine) matrix; the permanent is the canonical one, and it is
non-spectral precisely because it is not a determinant.

## Engineering corollary (domain 12, a NEGATIVE result worth recording)

For the H-fingerprint program: **do not compute the skew spectrum / Pfaffian — it is
spectrally redundant** (adds nothing beyond `charA`, which costs `O(n³)`). All
non-spectral fingerprint power must come from the unsigned side: `per(xI+A)` (THM-506,
strictly dominates the spectrum and `H`), `I(Ω_even,·)`, the packing carriers
`(c₆,c₇,…)`, and `H` itself. The signed objects are a closed (spectral) book.

## What I would chase next (handoff)

- **A clean general proof that `w_k` is spectral** for all `k` (Claim 3). Low orders are
  Moon; the general case should be a "walks decompose into closed-walk (trace) pieces +
  score-controlled boundary, and the boundary collapses by Moon-type identities" argument,
  or a direct walk-matrix / main-eigenvalue argument from the literature. Upgrading this
  to PROVED upgrades Claim 2 to a clean theorem of this project.
- **The permanental ROOTS** (still open from S6): the *roots* of `per(xI+A)` as a
  non-spectral invariant — are they a cleaner carrier of `(c₆,c₇,…)` than the coefficients?
- **The complex-signing question is now RESOLVED (this session).** The Hermitian signings
  `M_ω` that filter cycle length mod `r` (`r=3,4,6` tested) are ALL spectral — see the
  Sharpening above. The whole A-affine pencil is the spectral shore. So there is *no*
  single A-affine determinant that escapes the spectrum; the open frontier is genuinely
  *nonlinear* / *non-A-affine* constructions (e.g. the permanent of a non-trivially
  weighted matrix, immanants between det and per, or a matrix whose entries are
  themselves tournament sub-invariants).
- **Tie to path homology / #P:** the det(spectral) / per(#P) split here is the same
  easy/hard split as GLMY chain-rank (poly) vs the packing/independence counts (#P). The
  skew matrix is the spectral (poly) probe of the even cycle space; the cycle-space
  topology (β_k) is the homological probe — do they agree on the even slice?
```
