---
id: THM-1630
title: "⚠ PRIORITY CORRECTION — THE TORAL NULLCONE CONJECTURE IS A THEOREM, PROVED IN 1998. Λ = R(u)u^{−M} is a ONE-VARIABLE Laurent polynomial whose exponents run −M..N, hence two-sided for every N ≥ 1; and Duistermaat–van der Kallen (Indag. Math. 9 (1998) 221–231) prove precisely this case FIRST and SEPARATELY as their Theorem 2, in about one page, using only residues plus Liouville — their own words: 'much more elementary'. Their Remark 3 then gives exactly TNC's conclusion: CT(f^n) = 0 for all n forces f ∈ ℂ[z] or ℂ[z^{−1}] with no constant term, and for our Λ the ℂ[z] branch is impossible (lowest exponent −M < 0), leaving deg R < M — VERBATIM the degenerate case. So TNC(M,N) IS DvdK Theorem 2 + Remark 3 for every (M,N). Consequence: THM-1530 (M=1), THM-1550 (the Wiener–Hopf criterion), THM-1595 (the (2,2)/(2,3)/(2,4)/(3,3) closures), THM-1610 (the coefficient ladder) and this session's branching ladder are all REDISCOVERIES of special cases of a 27-year-old theorem. WHAT SURVIVES IS THE EFFECTIVE QUESTION, WHICH IS GENUINELY OPEN: Sturmfels' conjecture, in Erman–Smith–Várilly-Alvarado (JCTA 118 (2011) 396–402, arXiv:0908.2609), that for doubly-monic f the first nonzero constant term occurs by index m+n. The branching ladder's (2,2) output — first nonzero at m = 4 for Λ = u^{−2}+u^{−1}+u−u², where m+n = 4 — lands on that bound EXACTLY and tightly. The repo's machinery is aimed at the right open problem and was merely mislabelled. Also: char 0 is ESSENTIAL — over 𝔽₂, u+u^{−1} has CT ≡ 0 for every n."
status: >
  IDENTIFICATION VERIFIED against the primary source (van der Kallen posts the paper free
  at webspace.science.uu.nl/~kalle101/powers.pdf; Theorem 2 and Remark 3 were read, not
  inferred from search snippets). The hypothesis check is mechanical and is done below:
  Λ's exponents run −M..N, so Λ is two-sided whenever N ≥ 1, and DvdK Theorem 2 applies
  unconditionally. This is a CITATION, not a new proof.
  The ESV/Sturmfels effective conjecture is OPEN (as of 2011; no resolution found).
source: mac-mini-2026-07-20-S142 (owner: "keep focusing on the toral nullcone conjecture")
supersedes_in_part:
  - THM-1530   # M = 1
  - THM-1550   # the exact Wiener-Hopf criterion
  - THM-1595   # the (2,3)/(2,4)/(3,3) ladder closures
  - THM-1610   # the coefficient ladder
related:
  - THM-1600   # the Laplace layer; GMC(2)
---

# THM-1630 — TNC is Duistermaat–van der Kallen's Theorem 2

## The identification

`Λ = R(u)·u^{−M}` with `R(0) ≠ 0` and `deg R = d = M+N` is a **one-variable Laurent
polynomial** whose exponents run from `−M` to `N`. So for every `M ≥ 1, N ≥ 1` it is
**two-sided** — neither a polynomial in `u` nor in `u^{−1}`.

> **DvdK Theorem 2.** *"Assume that `f ∈ ℂ[z,z⁻¹]` is neither a polynomial in `z` nor a
> polynomial in `z⁻¹`. Then `f` has a critical value `v ∈ ℂ∖{0}` such that
> `limsup |Cst(fⁿ)|^{1/n} = |v| > 0`."*
>
> **DvdK Remark 3.** *"If for every `n ≥ 1` the constant term of `fⁿ` vanishes, then Theorem 2
> implies that either `f ∈ ℂ[z]` or `f ∈ ℂ[z⁻¹]`, without constant term."*

Applying Remark 3 to `Λ`, the two allowed conclusions are:

- `Λ ∈ ℂ[u^{−1}]`, no constant term ⟺ all exponents `< 0` ⟺ `d < M` ⟺ **`deg R < M`** — which
  is verbatim the degenerate case identified in THM-1610(B);
- `Λ ∈ ℂ[u]`, no constant term ⟺ all exponents `> 0` ⟺ `−M > 0` — **impossible** for `M ≥ 1`.

> **Therefore TNC(M,N) is exactly DvdK Theorem 2 + Remark 3, for every `(M,N)`. It has been
> a theorem since 1998.**

DvdK prove the one-variable case **first and separately**, calling it *"much more elementary"* —
residues plus Liouville, roughly one page. The sketch: `F(t) = Σ CT(fⁿ)t^{n−1}` picks up a
residue `−1/t` at `u = 0` (because `f(0) = ∞`, which is exactly where two-sidedness enters),
and the residues at the branches `f(ζ) = 1/t` are `O(t^{−1−1/m})` and cannot cancel it; if `F`
were single-valued around every nonzero critical value it would be entire and vanish at `∞`,
so Liouville gives `F ≡ 0`, a contradiction. A second, purely algebraic proof is Monsky's
(arXiv:0906.1836; written out in van den Essen–Schoone, arXiv:2305.10062, §5).

## What this costs the repo

**THM-1530, THM-1550, THM-1595, THM-1610 and this session's branching ladder are
rediscoveries** of special cases. That includes my own work: the coefficient ladder, the
`M | j` divisibility, the duality, and the branch analysis at `(2,2)`. None of it is wrong —
all of it is subsumed.

The duality is worth one line since it explains the indexing: `Λ(1/u) = R*(u)u^{−N}` with `R*`
the reversed polynomial, and `CT` is invariant under `u ↦ 1/u`, so **TNC(M,N) ⟺ TNC(N,M)** and
`min(M,N)` is the right parameter — which is how THM-1595 already indexed it.

## What survives, and it is the good part

**The EFFECTIVE question is open.** Sturmfels asked for an effective DvdK; the sharp form is

> **Conjecture (Sturmfels; Erman–Smith–Várilly-Alvarado, JCTA 118 (2011) 396–402,
> arXiv:0908.2609).** For doubly-monic `f = z^{−m} + ⋯ + z^{n}`, the constant terms
> `[[f¹]], …, [[f^{m+n}]]` generate the unit ideal — i.e. **the first nonzero constant term
> occurs by index `m+n`.**

Open as of 2011; no resolution found. It is **tight**: `f = z^{−N} + z` has first nonzero at
exactly `N+1`, so no universal constant bound exists.

**And the repo's machinery already computes exactly this.** The branching ladder's `(2,2)`
analysis produced the non-degenerate branch `R = 1 + u + u³ − u⁴`, i.e.
`Λ = u^{−2} + u^{−1} + u − u²` with exponents `−2..2`, and found its **first nonzero constant
term at exactly `m = 4` — and `m + n = 4`.** The computation lands on the ESV bound, tightly.

> **The redirect: stop proving TNC and start computing the effective bound.** The ladder is a
> ready-made instrument for it, and the question is genuinely open.

## Two cautions

- **Characteristic 0 is essential.** Over `𝔽₂`, `f = u + u^{−1}` is two-sided with
  `CT(fⁿ) ≡ 0 (mod 2)` for **every** `n ≥ 1` (odd `n` give 0; `n = 2k` gives `C(2k,k)`, and
  `v₂(C(2k,k)) = popcount(k) ≥ 1` by Kummer). Any mod-`p` version of TNC is simply false, and
  no repo argument that quietly works over `𝔽_p` can be citing DvdK.
- **DvdK for `n ≥ 2` variables is genuinely hard** (toric compactification, Hironaka
  resolution, Leray residues, monodromy). Only the one-variable case is elementary. If the
  toral programme ever needs a multi-variable analogue, that is a different order of difficulty
  and it should not be assumed to follow from the same citation.

## Honest scope

- This file **proves nothing new**. It is an identification plus a literature citation, and it
  should be cited as such.
- The hypothesis check (`Λ` two-sided; the `ℂ[u]` branch impossible) is mechanical and is the
  only mathematical content here.
- I have **not** independently reproduced DvdK's proof of Theorem 2; the sketch above is from
  the primary text.
- **Nothing here bears on GMC(2).** TNC was one input; THM-1600's span-2 elimination and the
  Laplace layer (HYP-8350/8445) are untouched, and whether the toral programme still needs
  anything from TNC beyond the citation is a question for whoever owns that thread.

*Sources:* Duistermaat & van der Kallen, Indag. Math. (N.S.) **9** (1998) 221–231, free at
`webspace.science.uu.nl/~kalle101/powers.pdf`; Monsky, arXiv:0906.1836 = EJC 18(1) (2011) #P5;
van den Essen–Schoone, arXiv:2305.10062; Erman–Smith–Várilly-Alvarado, arXiv:0908.2609.
*Artifacts:* `04-computation/tnc_branching_ladder_macmini_S142.py`.
