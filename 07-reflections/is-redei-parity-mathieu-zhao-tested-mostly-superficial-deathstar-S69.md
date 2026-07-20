# Is the odd-cycle-parity functional Mathieu–Zhao? Tested — trivially yes, deeply no

**death-star-2026-07-20-S69** (HYP-8530). Owner: investigate the S68 lead — is the tournament
odd-cycle-parity functional a Mathieu–Zhao (MZ) space, tying Rédei to the Jacobian ecosystem? I tested
it concretely. **Answer: `ker(Rédei-parity)` is MZ, but only trivially (it is a multiplicative
character's kernel); the GMC *depth* does not transfer.** This downgrades my own S68 "sharpest lead" —
honest scientific hygiene. What survives is a shallow "averaging-functional" analogy and one genuinely
open spectral direction.

## The MZ definition, and the two ways to realize it on tournaments

`ker Λ` is Mathieu–Zhao iff `Λ(P^m)=0 ∀m ⟹ Λ(QP^m)=0` for `m≫0`, all `Q`. GMC(n) is the statement
that the Gaussian expectation `E` on `ℂ[x_1..x_n]` has this property; `GMC(n)⟹JC(n)`. **GMC's depth is
that `E` is an *integral* — non-multiplicative (`E[PQ]≠E[P]E[Q]`).**

**Realization 1 — the parity count (multiplicative ⟹ trivial MZ).** The Rédei functional is
`ham(T)` = #Hamiltonian paths, always **odd** (verified n=2..5: counts `{1},{1,3},{1,3,5},{1,3,5,9,11,13,15}`).
Under ordinal sum it is **multiplicative**: `ham(T⊕S)=ham(T)·ham(S)` (verified). So it is a *character*.
For **any** character `χ`, `χ(P)=0 ⟹ χ(QP^m)=χ(Q)χ(P)^m=0` for all `m≥1` — so `ker χ` is MZ **for free**.
This is the *degenerate end* of the MZ spectrum: it says nothing, because a character factors powers.
It is categorically unlike GMC, where the interest is precisely that `E` does *not* factor.

**Realization 2 — the faithful integral (GMC, but false at n≥3).** The honest GMC analog is the actual
*expectation over random tournaments* `E_T[f]=2^{-\binom n2}Σ_T f(T)` — a genuine non-multiplicative
integral. But this is the arc-variable expectation on `\binom n2` variables:

| n vertices | 2 | 3 | 4 | 5 |
|---|---|---|---|---|
| arc variables `\binom n2` | 1 | 3 | 6 | 10 |
| regime | GMC(1) ✓ | **GMC(3) ✗** | GMC(6) ✗ | GMC(10) ✗ |

For `n≥3`, `\binom n2 ≥ 3`, and **GMC(N≥3) is FALSE** (has counterexamples — our own THM-1300 arc).
So `ker E_T` is **not** a Mathieu–Zhao space for any tournament on ≥3 vertices, and the odd-cycle
parity does not rescue it. The deep MZ property simply fails in this realization.

## What is actually shared, and what is not

- **Real but shallow:** both `E` (GMC) and `E_T` / the parity average are **Reynolds operators** — they
  average over a group action (U(1) charge for `E`; the `S_n`/relabeling or `Z_2` arc-parity for
  tournaments) and kill nonzero isotypic/graded components. "A functional kills a grading" is a genuine
  common structure. It is **not** an MZ statement — averaging is generic; MZ is special.
- **Not shared:** the MZ *depth*. GMC(2) is true because it is the *boundary* case (`\binom n2=1`,
  n=2 vertices only). Every tournament with an actual cycle needs ≥3 vertices ⟹ ≥3 arcs ⟹ the false
  regime. So there is no room for a true, deep tournament-MZ in the direct arc realization.

## The OCR is a *derivation*, not a functional — a better (still open) bridge

The Odd-Cycle Collection Formula `H(T) − H(T−v) = 2Σ_{C∋v}μ(C)` (Grinberg–Stanley) is a **difference
operator** (discrete derivative under vertex deletion), not a functional whose kernel one tests for MZ.
The Jacobian Conjecture has a second face — **locally nilpotent derivations** (Makar-Limanov, the LND
route) — parallel to the Mathieu–Zhao/functional face. The OCR's derivation shape suggests that *if*
Rédei/tournaments connect to the JC ecosystem, it is more likely through the **LND/derivation side**
(vertex-deletion as a locally nilpotent operator on the tournament algebra) than the functional-MZ side.
This is the honest relocation of the S68 lead — from "is parity MZ" (answered: trivially) to "is
vertex-deletion a JC-relevant derivation" (open, and better-shaped).

## Honest status and the one surviving deep direction
**S68's "is Rédei-parity Mathieu–Zhao" is answered: trivially yes, deeply no.** I am downgrading that
backlog lead accordingly. The only route to a *deep* (GMC(2)-true) tournament-MZ would need a genuine
**2-real-dimensional / 1-complex reduction** of the tournament — e.g. a single complex spectral invariant
(the H-spectrum / Paley-circulant eigenvalue phase) whose Gaussian moments form an MZ kernel. I have **no
clean formulation** and flag it as speculative, not a claim. Net: a lead tested and mostly retired, with
the derivation reframing kept.

## Credit
death-star-S68 (the lead, here tested and downgraded — my own); the Grinberg–Stanley OCR (Claim A, canon);
Zhao (GMC/MZ, `GMC(n)⟹JC(n)`); Makar-Limanov (LND face of JC); project Rédei/H-spectrum canon.

## Cross-links
S68 atlas (Lens 5), MISTAKE-avoidance on my own lead · `01-canon/definitions.md` (OCR) · CONJ-001 / THM-1300
(GMC(3)-false arc) · `04-computation/redei_mz_deathstar_S69.py` · HYP-8530.
