---
id: THM-647
title: The Anti-Rédei theorem — every self-converse tournament has an ODD number of ρ₀-anti-palindromic Hamiltonian paths, for every involutory anti-automorphism ρ₀ (which always exists); proves THM-644's Anti-Rédei conjecture
status: PROVED (two independent proofs; exhaustively machine-verified on all SC classes n=3..7, every involutory twist; necessity of the involutory hypothesis exhibited by 36 even-Fix non-involutory examples)
source: monad-explorer-2026-07-07-S8 (HYP-5007; owner directive "prove anti-redei via the rho-twisted involution")
depends_on:
  - THM-001   # Rédei: #Hamiltonian paths is odd
  - THM-644   # opus-S139: gridsym = ρ-anti-automorphism; fiber law B(C)·|Aut| = H_anti
related:
  - THM-643   # mac-mini-S46: blue/black parity structure
  - HYP-4851  # klein-S161: T1-T5 parity theorems
  - HYP-4967  # monad-S6: line-incidence algebra
  - LEM-003   # tiling fibration freeness (anti-symmetric form used by THM-644)
---

# THM-647 — The Anti-Rédei theorem

**Theorem.** Let T be a self-converse tournament on n vertices.

1. *(Existence of an involutory twist.)* T admits an involutory anti-automorphism ρ₀
   (a permutation with ρ₀² = id and i→j ∈ T ⟺ ρ₀(j)→ρ₀(i) ∈ T).
2. *(Anti-Rédei.)* For **every** involutory anti-automorphism ρ₀ of T, the number of
   ρ₀-anti-palindromic Hamiltonian paths
   Fix(τ) = #{ P : reverse(ρ₀(P)) = P },  τ := rev ∘ ρ₀,
   is **odd**. In particular an anti-palindromic Hamiltonian path **exists**.

## Proof 1 (Sylow + Rédei — the ρ-twisted involution)

*(a) |Aut(T)| is odd* (classical): an automorphism of even order powers to an involution
α, which must swap some pair {i, j}; α maps the arc between i and j to the arc between
α(i)=j and α(j)=i, i.e. reverses it — contradicting α being an automorphism.

*(b) Existence of ρ₀.* Let ρ be any anti-automorphism (T is self-converse). The group
G = ⟨Aut(T), ρ⟩ contains Aut(T) with index 2 (the anti-automorphisms form the nontrivial
coset: a product of two anti-automorphisms is an automorphism). So |G| = 2·|Aut(T)| with
|Aut(T)| odd; a Sylow 2-subgroup of G has order 2, and its generator is an involution
which cannot lie in the odd-order subgroup Aut(T) — hence an involutory anti-automorphism.

*(c) The twisted involution on Hamiltonian paths.* If P = v₁→…→vₙ is a Hamiltonian path
of T, then each arc ρ₀(v_k) → ρ₀(v_{k+1}) lies in conv(T), so the reversed sequence
ρ₀(vₙ)→…→ρ₀(v₁) is a Hamiltonian path of T. Thus τ(P) := reverse(ρ₀(P)) maps the
Hamiltonian-path set to itself, and since a vertex-relabeling commutes with reversal,
τ²(P) = ρ₀²(P) = P: τ is an involution. Non-fixed paths pair off, so
#HP(T) ≡ #Fix(τ) (mod 2), and Rédei's theorem (THM-001) gives #HP(T) odd. ∎

## Proof 2 (fiber-law assembly, at the tiling level)

By THM-644 (opus-S139), the grid-symmetric tilings of a class C are exactly the tilings
fixed by the canonical involutory twist, and the anti-symmetric orbit–stabilizer (LEM-003)
gives H_anti(C) = B(C)·|Aut(C)|. B(C) is odd on every self-converse class (the σ-parity
theorems: klein HYP-4851 T1–T5 / mac-mini THM-643 / monad HYP-4967, all three engines),
and |Aut| is odd by (a). Hence H_anti(C) is odd. ∎

## The reduction for non-involutory twists (and necessity)

For an arbitrary anti-automorphism ρ, let m = ord(ρ²) (odd, since ρ² ∈ Aut). Then
τ_ρ = rev ∘ ρ has order 2m, and because m is odd and reversal commutes with vertex maps,
τ_ρ^m = rev ∘ ρ^m with ρ^m **involutory**. So the unique involution in ⟨τ_ρ⟩ is the τ of
an involutory representative: **the parity always reads off an involutory twist.** The
involutory hypothesis is genuinely necessary: #Fix(τ_ρ) for non-involutory ρ can be even —
machine census found 8 such (class, ρ) examples at n=6 and 28 at n=7
(`anti_redei_proof_verification_monad_S8.out`).

## Verification

All SC classes n=3..7 (counts 2/2/8/12/88, matching canon): involutory ρ₀ exists for
every class; #Fix(τ_{ρ₀}) odd for **every** involutory ρ₀ (not just one choice);
Rédei re-confirmed on every class rep. Script:
`04-computation/anti_redei_proof_verification_monad_S8.py`.

## Corollaries and remarks

- **Existence:** every self-converse tournament has at least one anti-palindromic
  Hamiltonian path (klein-S161's "Rédei refinement" — odd anti-reversible orbit count —
  is this statement at the orbit level).
- **Blue mass:** B(C) = H_anti(C)/|Aut(C)| ≥ 1 on SC classes re-derives "SC is never
  pure black" (THM-643 T2) independently of the tiling census.
- **Lean shape:** Proof 1 is formalizable with modest effort: odd-Aut (arc-reversal
  contradiction), Sylow-order-2 (or direct: ρ^m for m = ord(ρ²) is an explicit involutory
  witness — no Sylow needed, fully constructive), the τ-involution pairing, and Rédei
  as the imported parity. The constructive ρ^m route avoids choice entirely.
