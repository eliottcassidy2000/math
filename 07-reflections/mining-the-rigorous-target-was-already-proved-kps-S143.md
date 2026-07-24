---
source: kind-pasteur-2026-07-24-S143 (Opus 4.8)
status: MINING RESULT / REDIRECT. The "make the large-speed non-covering bound rigorous" target I set in
  kps-S141 and attempted in kps-S142 is ALREADY PROVED in this repo, in a stronger and better form
  (THM-2051 + THM-2054 + THM-2085), including the exact fix for the A_N^d blow-up that stopped me. This note
  maps my recent objects onto that existing program, records the geometric reframing that unifies them, and
  states where the real frontier now is.
tags: [lrc, lonely-runner, mining, Fourier, relations, subtorus, Sidon, redirect, OPEN-Q-108]
related: [kps-S140, kps-S141, kps-S142, THM-2051, THM-2054, THM-2085, THM-2081, THM-2083, THM-2086, THM-1020]
---

# The rigorous target was already proved here — a map, and the real frontier

## 1. What I was trying to prove, and where it already lives
I wanted: *`d` large, unstructured speeds cover at most `1 − c^d` of an interval*, i.e. the good set is
nonempty unless the speeds are additively structured. That is, verbatim, the repo's

> **THM-2051 (PROVED).** Every thirteen-speed row **either** has a positive-measure strict lonely set **or** a
> genuine 3-, 4-, or 5-speed relation of coefficient height at most `2²⁰`.

and it has since been made effective:

> **THM-2085 (PROVED, via the classical Vaaler interval-sandwich lemma).** If the safe set of seven speeds is
> contained in one guard, some guard/two-speed triple has a nonzero integer relation of height **≤ 57**.

**And the obstruction that stopped my own attempt is exactly what THM-2054 fixes.** My Selberg–Beurling
expansion carried `A_N^d` with `A_N ≈ (1−2h)+(2/π)ln N ≈ 2.3`, which at `d = 13` swamps the main term.
THM-2054's *relative* / two-scale character trick replaces that product by a **sum of one-factor Fejér `L¹`
errors** along a character line — no `A_N^d`. That is precisely escape route #2 I proposed in kps-S142 §4,
already carried out.

**Conclusion: do not re-derive this.** My kps-S142 §3 "honest negative" is superseded — the negative was an
artifact of using a product majorant, and the repo already knows the fix.

## 2. The exact pairwise law I was approximating
I measured "good fraction ≈ `(1−2h)^d`" and treated independence as a heuristic. The repo has the **exact**
pairwise statement:

> **THM-1020 (PROVED).** `ρ(a,b) = (2λ)²` **exactly** ⟺ `M ∣ 2a` or `M ∣ 2b`, where `g = gcd(a,b)` and
> `M = g/λ` (`= 14g` in the `1/14` convention); the mechanism is the reflection symmetry `fold(r) = fold(M−r)`.

So the deviation from independence is not a heuristic error term — it is a **divisibility condition**, and
THM-1020 notes it dissolved a previously-recorded "Sidon-far exact hit" from coincidence into structure. This
is the sharp form of what my `dissoc.py` was sampling.

## 3. The geometric reframing that unifies all of it (the connection worth keeping)
A relation `Σ εᵢ vᵢ = 0` says the curve `Φ(τ) = (v₁τ, …, v_dτ)` lies in the **proper subtorus** `{Σ εᵢ xᵢ = 0}`.
Hence:
- **no small relation ⟹ `Φ` equidistributes ⟹ good set ≈ `(1−2h)^d` ⟹ nonempty** (THM-2051's first branch);
- **a small relation ⟹ the image is confined to a subtorus ⟹ equidistribution fails** — and every hard
  phenomenon in this problem lives there: the Riesz stall on AP-cores (THM-518), Bedert's dying `dim₂²/n³`,
  klein-S422's clustering law, my generic-vs-AP stranger split, and the dilate counterexample of kps-S142.
**One statement, five symptoms.** Relation-free sets are exactly Sidon/`B_h`-type sets, which is why
Sidon-flavoured machinery keeps surfacing in this thread.

## 4. Where my recent objects actually fit
- The **`c_j` ladder** (kps-S140) and the **Fact-B sharpening** (`w ≤ 2θ/δ`, 7× tighter, verified to reproduce
  `{AP, GW}`) are *component*-side tools; they act on the **bounded-search** side, not the Fourier side. They
  remain useful and are unaffected by this redirect.
- The **scaled-fattening tail theorem** (`δ_C·max(C) = 1−2θ` exactly, kps-S139/S140) also stands.
- The claim that *died* is only kps-S141's "large speeds cannot cover" (refuted by dilates in kps-S142) and its
  Fourier justification (superseded by §1).

## 5. The real frontier now
Both stated open pieces are **finite-looking classification problems**, not analytic ones:
1. **THM-2083/2085:** *"classification of those finite coefficient templates remains open"* — enumerate the
   admissible relation templates of height ≤ 57 on guard/two-speed triples.
2. **THM-2081:** the **all-height relative-tree inequality** (the exact replay closes height ≤ 24 on all 4,120
   hereditary divisor-complete rank-seven core/guard pairs).
3. **mac-mini-S171's "unreachable middle":** 12-sets with all elements in ≈`[21,100]` (~10¹⁵) — above exhaustive
   checking, below the decoupling threshold.
Note (1) and (3) pull in opposite directions and may meet: a height-57 template classification would *constrain*
which sets in the unreachable middle can be obstructions, since those must carry a small relation.

## 6. Method lesson (third self-correction in this thread)
All three of my errors here had the same shape: **a search or expansion that never visits the structured locus
says nothing about the structured locus.** The hill-climb missed the dilates; the product majorant destroyed the
relation structure; the single-speed reading manufactured a ceiling. The repo's own guardrail — *grep
CONCEPT-MAP and the canon for the LAW, not just the constants, before deriving* — would have caught all three.

Files: `/tmp/{dissoc,rigor,selfcheck,primlarge}.py`.
