# The unique-channel DvdK in Lean, and the cancellation ↔ inclusion-exclusion dictionary

*boxeph-2026-07-22-S229. Owner: work the next Lean DvdK target, and mine past results on cancellation
and inclusion-exclusion. Builds on S226/S227 (two-charge DvdK Lean), S228 (positive-coefficient DvdK
Lean, HYP-8925), death-star-S101/HYP-8878 (the unique-primitive-cycle criterion), codex THM-2067
(Galois orbit-product), and the repo's inclusion-exclusion spine (THM-406, THM-515, THM-671, THM-1820,
THM-1810). New Lean file `04-computation/lean/TournamentH7/TournamentH7/GMC2DvdKUniqueChannel.lean`
builds kernel-pure.*

## The next Lean target: the unique-channel criterion, mechanized for any support

S226/S227 formalized the two-charge DvdK by a *specific* uniqueness lemma (`balanced_unique`: the
`Fin 2` linear system has one solution). S228 formalized the *positive-coefficient* case. The right
next target is the theorem that **subsumes both** and is the actual DvdK-free frontier — death-star's
unique-primitive-cycle criterion (HYP-8878/S101), lifted to arbitrary support and mechanized:

```lean
theorem ct_ne_zero_of_unique_balanced (q : ι → ℤ) (c : ι → ℂ) (hc : ∀ i, c i ≠ 0) (m : ℕ)
    (r0 : ι → ℕ) (hr0 : r0 ∈ piAntidiag univ m) (hbal : totalCharge q r0 = 0)
    (huniq : ∀ r ∈ piAntidiag univ m, totalCharge q r = 0 → r = r0) :
    aeval c (constantTermRelation q m) ≠ 0
```

If size `m` has a **unique** balanced composition, `CT(f^m)` collapses to the single term
`multinomial(r0)·∏ cᵢ^{r0 i}`, which no coefficient choice can cancel — nonzero for *every* complex
`c`. This is exactly S101's criterion (`CT(f_F^{m*}) = multinomial(m*,r*)·c^{r*} ≠ 0` when the minimal
balanced channel is unique), now kernel-pure (`#print axioms = [propext, Classical.choice, Quot.sound]`)
and **coefficient-independent** — uniqueness is a property of the charges alone. S101's Python scan
found this holds for **98 of 116 (84%)** straddling supports of size 3–4 in `[-4,4]`; the Lean theorem
is the general statement behind that scan.

I gave it two further honest faces:

- **The cancellation dichotomy** (`two_balanced_of_ct_zero`, the contrapositive): if `CT(f^m)=0` at
  nonzero coefficients while one balanced channel exists, there is a *second* one. Cancellation is
  impossible without ≥2 coincident channels.
- **The cardinality form**, making the death-star stratification literal in Lean:
  `ct_ne_zero_of_card_eq_one` (`|balancedSet q m| = 1 ⟹ CT ≠ 0`, the DvdK-free 84%) and
  `two_le_card_balanced_of_ct_zero` (`CT = 0 ⟹ |balancedSet q m| ≥ 2`, the hard resonant 16%). The
  coincident-cycle *count* is now a Lean `Finset.card`.
- `two_charge_via_unique` re-derives the S226 pair theorem as the `Fin 2` instance, confirming the
  generalization subsumes it.

So the DvdK-free zone is now proved inside Lean for **any support with a unique minimal channel**, and
the exact boundary — `card ≥ 2` coincident cycles — is where codex THM-2067's Galois orbit-product must
act. That is the honest state of "the next Lean target": the elementary side is complete and general;
the Galois argument for the coincident stratum remains the open Lean frontier.

## Mining cancellation ↔ inclusion-exclusion: they are one object

The mine turned up that the repo already proves the two themes are the *same* phenomenon, and my Lean
sits on the cancellation side of a precise dictionary:

**The identity that fuses them (THM-1820, the bridge).** The GMC constant-term/moment nullcone
`E[P^m]=Σ_{balanced channels} multinomial·c^r` and the LRC covering are *one* moment-nullcone,
"positivity of a moment functional past a cancellation wall," differing only by the **boundedness of
the moment alphabet**. LRC's depth `X∈{0,…,13}` is bounded ⟹ its Bonferroni inclusion-exclusion
`Σ_k(−1)^k S_k` **terminates at k=13** ⟹ a finite-depth certificate exists (B5). GMC's radial degree is
unbounded ⟹ no finite depth (detection depth `≥ d+1`, THM-1790). The discriminant is `|alphabet|`.

**The LRC side is literally inclusion-exclusion (THM-406 M1b).** The lonely measure is
`p₀ = μ(N=0) = Σ_{j=0}^{n}(−1)^j S_j`, the alternating sum of multi-arc overlap volumes
`S_j = Σ_{|I|=j} μ(⋂_{i∈I} D_i)`. THM-515 gives the dual sinc form `L(S)=Σ_{t∈Λ}∏h(tᵢ)`, a signed
lattice sum that **alternates** (THM-504: the saving is cross-level `+S₂−S₃+S₄−⋯`, not within-level).

**The sign mechanism is the same involution (THM-1810).** DvdK cancellation and inclusion-exclusion
alternation are the bosonic/fermionic (permanent/determinant) dichotomy: the fermionic/alternating side
carries a sign that an involution can use to cancel, the permanent/Gaussian-moment side does not. The
dihedral DvdK witness `f=u²+u+u⁻¹−u⁻²` (THM-2070: `f(−u⁻¹)=−f(u)` ⟹ `CT(f^m)=(−1)^m CT(f^m)`) is
exactly such a sign-involution killing all odd moments — the CT-sum image of inclusion-exclusion parity.

**The dictionary.**

| | DvdK `CT(f^m)` | LRC loneliness `p₀` |
|---|---|---|
| signed sum | `Σ_{balanced} multinomial·∏cᵢ^{rᵢ}` | `Σ_j (−1)^j S_j` (THM-406 M1b) |
| positive parts | `multinomial(r) > 0` | overlap volumes `S_j ≥ 0` |
| sign source | complex phases in `c^r` | inclusion-exclusion parity `(−1)^j` |
| **positive regime** | positive `c` (S228) / **unique channel, 84%** (S101, **this Lean**) | **odd-Bonferroni B5 certificate**, THM-671; codex guard-capacity THM-2076/2080 |
| hard stratum | **≥2 coincident cycles, 16%** (S101, **this Lean** `card ≥ 2`) | resonant / AP cores (`S₅` explodes 156–172× iid, THM-686) |
| all-cancellation crux | general complex — **THM-2067 Galois orbit** | all-orders `{p₀=0}` — THM-406 M2 |
| what tames the sign | transitive Galois group on the roots | mirror `ι` / doubling homeomorphism (THM-2075) |

## Honest self-correction of S228: LRC does have a positive regime

In S228 I wrote that the LRC covering is "all-cancellation, with no positive regime, unlike GMC's
angular part." The mine shows that is **too strong**, and the correction sharpens the parallel rather
than weakening it:

- The *exact measure* `p₀ = Σ(−1)^j S_j` is indeed all-orders cancellation (THM-406 M2: no finite
  moment truncation certifies its sign near the floor — the rigorous Vitali wall). My claim was right
  **for the exact equality**.
- But inclusion-exclusion has a genuine **positive regime**: the odd Bonferroni truncations are
  *lower* bounds (`p₀ ≥ Σ_{j≤2r+1}(−1)^j S_j`), and THM-671's quintic certificate `B₅ = 2052/7⁵ ≈
  +0.1221` is the *first positive* truncation past the wall (`B₁=−0.857`, `B₃=−0.099`), certifying
  loneliness on **low-resonance / Sidon cores** — precisely the 84%-analogue. codex's guard-capacity
  line (THM-2076 union bound + THM-2080 overlap floor `μ(D_q∩E_h)≥1/42`) is the same positive
  regime, and it is what drove tower depth 8→5→4.

So the honest statement is: **LRC covering has a positive (odd-Bonferroni / capacity) regime that
handles the low-resonance cores, exactly parallel to the positive-coefficient / unique-channel DvdK
regime that handles the 84%; and it provably cannot reach the all-orders floor on the resonant/AP
cores (THM-406 M2), exactly as positive-coefficient DvdK cannot reach the general complex case
(THM-2067).** The two problems are the same moment-nullcone with the same positive-vs-cancellation
split; only the alphabet (finite vs. unbounded) differs, which is why LRC's cancellation *terminates*
(a finite B5 certificate is conceivable) while GMC's does not. That termination is the structural
reason LRC(14) is a finite—if resonant—target and GMC(2) is not.

## Scope

Concrete: the unique-channel DvdK-free criterion (death-star-S101/HYP-8878), and its cancellation
dichotomy stated by coincident-cycle cardinality, are now kernel-pure in Lean for arbitrary support —
generalizing and subsuming the S226 two-charge and complementing the S228 positive-coefficient files.
The DvdK-free zone in Lean is now every support with a unique minimal channel (the elementary 84%);
the coincident-cycle stratum (`card ≥ 2`) remains the THM-2067 Galois frontier. The cancellation ↔
inclusion-exclusion synthesis is a reading of existing proved theorems (THM-1820/1810/406/515/671),
not a new theorem, and it carries an honest correction of my S228 "no positive regime" claim.

Links: HYP-8930, HYP-8925, HYP-8878, THM-2067, THM-2070, THM-1820, THM-1810, THM-406, THM-515,
THM-671, THM-2076, THM-2080,
[[starting-to-formalize-dvdk-the-positive-coefficient-case-and-where-cancellation-is-the-crux-boxeph-S228]],
[[a-kernel-pure-lean-proof-of-the-two-charge-dvdk-seed-boxeph-S226]].
