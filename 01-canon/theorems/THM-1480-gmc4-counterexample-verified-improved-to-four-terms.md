---
id: THM-1480
title: "THE GMC(4) COUNTEREXAMPLE IS CORRECT, AND IT IMPROVES TO FOUR TERMS. (0) NO, WE DID NOT HAVE IT — THM-1435 states in terms 'THE WITNESS IS NOT PRODUCED'; the repo held transport machinery and a bracket 5 ≤ vcwd ≤ ~20, never a witness. (A) Verified exactly two independent ways: over ALL 256 readings of the notation (exactly 8 survive, one symmetry orbit, all giving m!) and again by real-coordinate expansion in u₁..u₄ ~ N(0,1/2) — so it genuinely lives in FOUR REAL GAUSSIANS. (B) THE MECHANISM IS THE ALTERNATING BINOMIAL IDENTITY: since C(m,k)k!(m−k)! = m!, both sides collapse to E[Pᵐ] = m!·Σⱼ(−1)ʲC(m,j) = 0 and E[QPᵐ] = m!·C(m,0) = m!. Multiplying by Q merely SHIFTS THE SUMMATION INDEX BY ONE and cashes the boundary term. (C) STRICT IMPROVEMENT: the published 'cubic, 6 terms' is NOT minimal — P = (1+W)(W̄ − ZZ̄) = (1+W)(W̄ − |Z|²) with Q = W is a FOUR-term counterexample with the identical m! growth; the lone Z̄ is redundant. Three terms is impossible in the same box, so four is minimal there. (D) An INFINITE FAMILY: Q = Z₂ʳ gives E[QPᵐ] = (−1)^{r+1} m!·Σ_{j<r}(−1)ʲC(m,j), i.e. m! times a degree-(r−1) polynomial; the published case is r = 1. (E) RIGIDITY: deforming (1−Z) → (1−λZ) gives E[Pᵐ] = m!(1−λ)ᵐ, so λ = 1 is forced EXACTLY. (F) The obvious route to n = 3 is STRUCTURALLY BLOCKED: replacing |Z|² by a real square X² replaces k! with (2k−1)!!, and the collapsing sum becomes 0,1,0,9,0,225 — it does not vanish. (G) THE REFRAME THAT MATTERS: this is NOT a witness for an already-known corollary. ¬JC(3) only gives 'GMC(N) false for SOME finite N'; the example LOCALIZES the failure to n = 4, which the published implication chain cannot do. Sharpness remains open — GMC(1) is TRUE, general GMC(2) is OPEN, GMC(3) unknown."
status: >
  (A) VERIFIED-EXACT, two independent formalisms, m = 1..9 (complex) and m = 1..5 (real).
  (B) PROVED (three lines; independently re-derived by a second route this session, and it
  matches the structure of DEZ's own Prop. 1.2 half-disk example — see Provenance).
  (C) VERIFIED-EXACT: the 4-term reduction satisfies both claims for m = 1..9; exhaustive
  over supports of size ≤ 3 in the degree-≤3 complex-monomial box finds NOTHING, so 4 is
  minimal IN THAT BOX. The box contains the published example, so this search has a
  working positive control.
  (D) PROVED and verified r = 1..4, m = 1..6.
  (E) PROVED and verified λ ∈ {0,1,2,3,−1}.
  (F) PROVED (the two moment sequences differ) and exhibited numerically.
  (G) SOURCED — see Provenance. The n = 3 SEARCH in this session is INCONCLUSIVE and is
  reported as such: its n = 4 positive control also found nothing, so the box was
  inadequate and the negative result carries NO evidential weight.
source: mac-mini-2026-07-20-S133 (owner supplied the counterexample from an outside source
  and asked whether the repo already had it or anything stronger, then to improve it)
provenance: >
  GMC is Derksen–van den Essen–Zhao, arXiv:1506.05192, Israel J. Math. 219 (2017) 917–928,
  their Conjecture 1.1 / 2.2. The implication is ONE-WAY (their Thm 1.6: GMC for all n ⟹ JC
  for all n); no reverse implication is known. Zhao's VC is arXiv:0704.1691; his Thm 2.9
  makes VC-for-Laplacians equivalent to JC — a strictly stronger link than GMC's.
  The counterexample is written in DEZ's OWN notation (their Prop. 3.2 defines exactly
  W_j = (X_j − iY_j)/√2, Z_j = (X_j + iY_j)/√2) and is structurally the Gaussian transplant
  of DEZ's OWN Prop. 1.2 — their half-disk counterexample to the Moment Vanishing Conjecture
  (P = (x+iy)², Q = x+iy, ∫Pᵐ = 0 ∀m but ∫QPᵐ ≠ 0 ∀m): same shape, same "fails for every m."
  A literature sweep found NO prior instance of this Gaussian counterexample. The only
  2026 item is Zihan Zhang's expository page (20 July 2026), which derives only the
  non-localized "GMC(N) is false for at least one finite N."
depends_on:
  - THM-1300  # the JC(3) counterexample the chain rests on
related:
  - THM-1435  # "THE WITNESS IS NOT PRODUCED" -- the repo's state before this
  - THM-1430  # explicit symmetric Keller counterexample on C^6
  - HYP-8240  # the witness-extraction assessment
script: 04-computation/gmc_counterexample_verify_macmini_S133.py,
        gmc_mechanism_and_family_macmini_S133.py,
        gmc_three_gaussian_search_macmini_S133.py (+ .outs)
---

# THM-1480 — the GMC(4) counterexample: verified, explained, and improved to four terms

## (0) Did we already have it? — **No, and the repo says so in terms**

THM-1435 (klein-S337) states plainly: **"THE WITNESS IS NOT PRODUCED."** What the repo held
was the de Bondt transport machinery (proved as a symbolic matrix identity, control-validated
with two positive and one negative control), a proof that the Bass–Connell–Wright step cannot
be skipped, and a bracket `5 ≤ vcwd ≤ ~20` on the VC-witness dimension. No witness, for GMC or
for VC.

**One genuine convergence is worth recording.** THM-1435 §D found — independently, from the
Alpöge map — that **homogeneity is the load-bearing obstruction**: Alpöge's `H` is not
homogeneous, so `JH³ ≠ 0` and BCW is essential. The GMC counterexample is *also* essentially
inhomogeneous: its degree-1 and degree-3 homogeneous parts, taken alone, are **not**
counterexamples (they are eventually zero, hence they *satisfy* GMC). Two independent routes
landing on inhomogeneity.

## (A) Verified, exactly, twice

Taking nothing on trust — including the notation — all `4⁴ = 256` assignments of
`{Z₁,Z₂,W₁,W₂}` to the letters `{Z, Z̄, W, W̄}` of two independent standard complex Gaussians
were tested. **Exactly 8 survive** (one orbit under `Z ↔ W` and conjugation), all giving
`E[Pᵐ] = 0` and `E[QPᵐ] = m!` for `m = 1..9`; **no reading gives a different growth law.**
A second, independent real-coordinate expansion in `u₁..u₄ ~ N(0,½)` reproduces
`E[QPᵐ] = 1, 2, 6, 24, 120` with zero imaginary part — confirming it lives in **four real
Gaussians**.

## (B) The mechanism: one binomial identity

With `E[Z^aZ̄^b] = δ_{ab}a!` and the surviving reading `P = (1+W)(Z̄(1−Z) + W̄)`, `Q = W`,
expanding and using **`C(m,k)·k!·(m−k)! = m!`** collapses everything:

> `E[Pᵐ] = m!·Σ_{j=0}^{m}(−1)ʲC(m,j) = 0`  and  `E[QPᵐ] = m!·Σ_{j=1}^{m}(−1)^{j−1}C(m,j) = m!·C(m,0) = m!`

**Both sides are the same alternating binomial sum; multiplying by `Q` shifts the summation
index by one and cashes the boundary term `C(m,0) = 1`.** The factor `(1−Z)` manufactures the
alternating signs; the Gaussian pairing manufactures the factorials that cancel the binomial
denominators exactly; `Q` collects the residue.

## (C) The published object is not minimal — **four terms suffice**

> **`P = (1 + W)(W̄ − ZZ̄) = (1 + W)(W̄ − |Z|²)`, `Q = W`** — 4 terms, still cubic, and
> `E[QPᵐ] = m!` identically.

The lone `Z̄` in the published second factor is **redundant**. Two of the six single-term
deletions leave a counterexample; the surviving core is the above. Exhaustive search over the
degree-≤3 complex-monomial box finds **zero** counterexamples with support ≤ 3 and many with
support 4 — so **4 is minimal in that box**, and the box contains the published example, so
the search has a working positive control.

Search over degree-≤2 supports of size ≤ 5 finds nothing, so **cubic appears necessary**
(bounded search, not a proof).

## (D) An infinite family

> `Q = Z₂^r` gives `E[QPᵐ] = (−1)^{r+1}·m!·Σ_{j=0}^{r−1}(−1)ʲC(m,j)`

— `m!` times a polynomial of degree `r−1` in `m`. Verified `r = 1..4`: `m!`, `m!(m−1)`,
`m!(1−m+C(m,2))`, and a degree-3 factor. **Every `r ≥ 1` breaks GMC**; the published
counterexample is the `r = 1` member.

## (E) Rigidity

Deforming `(1−Z) → (1−λZ)` gives `E[Pᵐ] = m!(1−λ)ᵐ`, verified at `λ = 0,1,2,3,−1`. So
**`λ = 1` is forced exactly**: the construction sits at the unique point where the alternating
sum vanishes. It is not a member of a continuous family.

## (F) The obvious route to `n = 3` is structurally blocked

`|Z|²` and a real square `X²` both have mean 1, so replacing one by the other is the natural
way to save a dimension. It fails, and provably: `E[|Z|^{2k}] = k!` is *exactly* what cancels
`C(m,k)`, whereas `E[X^{2k}] = (2k−1)!!` is not. The analogous sum becomes

`0, 1, 0, 9, 0, 225, 0` for `m = 1..7` — **it does not vanish.**

## (G) The reframe that matters

It is tempting — I assumed it at the start of this session — to read the counterexample as a
mere *witness* for something `¬JC(3)` already implies. **That is wrong.** The DEZ implication
is one-way and stated only globally (`GMC(n) ∀n ⟹ JC(n) ∀n`), so `¬JC(3)` yields only

> "`GMC(N)` is false for **some** finite `N`" — with no bound on `N`.

The explicit example **localizes the failure to `n = 4`**, which the published chain cannot do.
That is a new logical fact, not a corollary. And sharpness is genuinely open: **`GMC(1)` is
TRUE** (DEZ Prop. 4.2), **general `GMC(2)` is OPEN** (proved only for homogeneous `P`,
Cor. 4.4), and **`GMC(3)` is unknown**.

## Honest scope

- **The `n = 3` search in this session is INCONCLUSIVE and carries no evidential weight.**
  Its `n = 4` positive control *also* found nothing — the real-coordinate expansion of
  `(1+W)(W̄−|Z|²)` has more terms and non-unit coefficients, so it falls outside the box. A
  failed control invalidates the negative result. Only (F)'s structural obstruction stands.
- (C)'s minimality is **relative to a box** (complex monomials, degree ≤ 3, unit coefficients).
  Non-unit coefficients or degree ≥ 4 are untested.
- **(B) is not uniquely ours.** The same mechanism was derived independently by a second route
  in this session, and it mirrors DEZ's own Prop. 1.2. Claimed as convergent re-derivation.
- **My earlier "heat shadow" framing was too strong.** The DEZ functional *is* the heat
  semigroup at the origin — `(e^{Δ/2}w^αz^β)(0) = α!δ_{αβ}`, verified — but GMC's hypothesis
  is the whole series `(e^{Δ/2}Pᵐ)(0) = 0` while VC's is the single term `Δ^m Pᵐ = 0`. These
  coincide only for `deg P = 2`. For the JC-relevant quartic case VC's hypothesis implies
  GMC's and not conversely, so **GMC's hypothesis class is strictly larger**. Neither DEZ nor
  Zhao states the relation; it is unsourced folklore.
- The `m!` here has **no established connection** to the `(n−1)!` arborescence count of
  THM-1460 or any other factorial in this repo. Do not cross-link them.

*Artifacts:* the three `gmc_*_macmini_S133.py` scripts (+outs).
