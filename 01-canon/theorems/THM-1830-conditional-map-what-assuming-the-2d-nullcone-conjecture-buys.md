---
id: THM-1830
title: "Historical conditional map from NC2: the charge implication survives; the LRC and JC arrows do not"
status: >
  SUPERSEDED CONDITIONAL MAP. The elementary implication from a strict
  one-sided NC2 classification to GMC(2) is sound, but NC2 and GMC(2) are now
  proved independently by THM-2022. There is no proved GMC(2)->LRC(14),
  NC2->LRC(14), or GMC(2)->JC(2) arrow. The Fock/vacuum and common-return-time
  discussions are conjectural analogies only. Retained for provenance; it is
  not a current proof graph.
source: klein-2026-07-20-S391 (owner: if we assume the 2-D nullcone conjecture is true, what can we prove — can LRC(14) follow via GMC(2), or another creative path?)
depends_on:
  - THM-1510  # NC2 ⟹ GMC(2), the one proved downstream arrow
  - THM-1820  # LRC and GMC(2) are one moment-nullcone problem, split by alphabet-boundedness
related:
  - THM-1745  # mac-mini: toral moment-count = LRC coprime-pair first-return time (the common parent)
  - THM-1790  # the EMP floor: why GMC(2) is not finitely certifiable (assuming NC2 doesn't fix this)
  - "boxeph HYP-8350: the Fock/vacuum-Mathieu bridge GMC(2) → DC₁ → JC(2), arrows conjectural"
  - "MISTAKE-202 / CASE-gamma-bridge-domination-step: why GMC(2)/NC2 is still open"
---

# THM-1830 — what assuming the 2-D nullcone conjecture buys

> **CURRENT CORRECTION.** THM-2022 has proved NC2 and GMC(2), so the open-status
> language below is historical. Its Frobenius lowest-face proof, not the false
> Gamma-domination bridge, supplies the closure. No LRC(14) or JC(2) consequence
> follows. In particular, a chain whose middle arrows are explicitly
> conjectural is not a path that “does close.”

**Historical hypothesis (now THM-2022).** `NC2`: in two real Gaussians, the only `P` with all moments vanishing
are the purely holomorphic / antiholomorphic ones (the 2-D nullcone is exactly the one-sided
locus) — the Derksen–van den Essen–Zhao-type 2-D nullcone statement. (Naming per the owner.)

## The one solid downstream arrow

> **`NC2 ⟹ GMC(2)`, PROVED** (THM-1510 §C, three lines). If `0 ≠ P ∈ N₂`, NC2 makes all its
> weights `≥ 1` (wlog); weights add, so `P^m` has all weights `≥ m`; a fixed `Q` has finite
> minimum weight, so `QP^m` has no weight-0 term for `m ≫ 0`, and Wick kills it. ∎

So the **logical implication** from NC2 to GMC(2) is sound. Historically it
was the one rigorous downstream arrow in this map. Both statements are now
theorems by THM-2022's separate lowest-face/Frobenius mechanism; the refuted
Gamma-domination route is not part of that proof.

## The owner's chain breaks at the LRC arrow

```text
   NC2  ==PROVED==>  GMC(2)  ==✗ no implication✗==>  LRC(14)
```

There is **no** implication `GMC(2) ⟹ LRC(14)`. What is true is weaker and non-logical:

- **Analogy (THM-1820).** LRC(14)-covering and GMC(2) are *the same moment-nullcone template*
  (kind-pasteur THM-1750) — "positivity past a cancellation wall" — but split by the boundedness
  of the moment alphabet: LRC has `|X| ≤ 13` (bounded, finitely certifiable by the quintic B5,
  THM-671), GMC has the unbounded radial degree (EMP floor grows, THM-1790). Neither specializes
  to the other.
- **Shared shape (mac-mini THM-1745).** The toral moment-count bound *is* the coprime-pair
  first-return time `(p+|n|)/gcd(p,|n|)`, and "max over straddles" *is* the shape of
  `M(S) = max_t min_i ‖v_i t‖`. So GMC's angular layer and LRC live in the same theatre.

But an analogy and a shared shape are not an implication. Concretely, LRC(14)'s hardness is the
**uniform covering finish** — a *bounded-alphabet* problem that B5 already certifies per instance
— and NC2 is a statement about the *unbounded* Gaussian nullcone. Assuming NC2 does not touch
LRC(14)'s actual difficulty. **The proposed chain does not close.**

## A conjectural Fock path toward JC(2) — it does not close

```text
   NC2  ==PROVED==>  GMC(2)  ==Fock/vacuum-Mathieu, conjectural arrows==>  DC₁  ==>  JC(2)
```

boxeph's Fock bridge (HYP-8350) reads GMC(2) as the **vacuum-Mathieu statement for the Weyl
algebra `A₁`**, and continues (arrows explicitly conjectural): GMC(2) ⟹ vacuum-ideal rigidity for
`A₁`-endomorphism images ⟹ `{P_top, Q_top} = 0` ⟹ `DC₁` (Dixmier, `n=1`) ⟹ `JC(2)` (the
2-dimensional Jacobian Conjecture). This is a speculative research program,
not a theorem-level downstream consequence. No arrow from GMC(2) to JC(2) is
proved here.

## The honest unification for LRC

If a path from the nullcone circle to LRC is wanted, it is **not** `NC2 ⟹ LRC` but a **common
parent**: the exactness of the **coprime-pair return-time / moment-count bound** (mac-mini
THM-1745) would settle *both* the moment-nullcone uniform bound (the TNC/GMC angular layer,
HYP-8505/8540) *and* the LRC realignment structure, because both are "max over coprime-pair
first-return times." That parent is a statement about return times, neither equal to NC2 nor
implied by it — so the two conjectures are **siblings under a return-time theorem**, not
parent-and-child. That is the creative reframe: don't route LRC through GMC(2); route both
through the coprime-pair return-time exactness they share.

## Verdict

| arrow | status |
|---|---|
| `NC2 ⟹ GMC(2)` | **PROVED** (THM-1510 §C) |
| `GMC(2) ⟹ JC(2)` | conjectural program (Fock bridge, HYP-8350) |
| `GMC(2) ⟹ LRC(14)` | **does not exist** (analogy only; THM-1820, THM-1745) |
| `NC2 ⟹ LRC(14)` | **does not exist** |

NC2 and GMC(2) are now proved by THM-2022. This file contributes only the
elementary implication between them and historical analogies. It does **not**
prove LRC(14) or JC(2), and no honest chain here does. The nearest proposed link is
a shared return-time parent that would settle LRC and the moment-nullcone uniform bound *together*
— which is where a creative attack on both should aim.

*No script; this is a conditional-implication audit of THM-1510, THM-1745, THM-1820, HYP-8350.*
