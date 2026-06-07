# HYP-2317 — The n=14 half-turn *blocker* is a mod-2 detector, orthogonal to the mod-7 fiber where the real obstruction lives (complement to HYP-2346)

**Session:** S639
**Status:** CONFIRMED (exact recomputation) — a structural result on the S367 quotient model + a
formalized mechanism; complementary to HYP-2346
**Provenance forward:** math-lean `Math/LonelyRunner/HalfTurnParity.lean` (sorry-free)
**Trigger:** an external LLM suggested treating LRC(14) as a fiber bundle over the 7-runner base
(`14 = 2·7`) and asked whether the half-turn leak (the "56 cells") is structured mod 7 or mod 2.

> **Reconciliation with HYP-2346 (concurrent claudebox S643) — COMPLEMENTARY, not contradictory.**
> S643 took the same prompt and concluded "the leak rides the **mod-7** fiber"; this session finds the
> half-turn *blocker*'s leak is **mod-2**. Different objects, and they fit:
> - **HYP-2346 (the *real* LRC):** at the 7-clock `t = b/7`, a runner is dangerous ⟺ `7 ∣ v`; every
>   non-multiple-of-7 has margin `≥ 1/7 ≫ 1/14`. LRC(7) (proven) is the base section; the *real*
>   obstruction is the small **mult-of-7** sub-config — it lives in the **mod-7 fiber**.
> - **HYP-2317 (this session, the S367 *quotient* covering model):** the half-turn *blocker* (residue
>   `7`) is **2-torsion**, hence a pure **mod-2 detector** (formalized below); its 56-cell leak is the
>   odd mod-2 coset, uniform across all mod-7 classes.
>
> Together: the half-turn (a mod-2 tool) is **orthogonal** to the mod-7 fiber that carries the
> difficulty — a Künneth-type CRT split, half-turns confined to the benign `ℤ/2` summand and
> **structurally unable** to see the `ℤ/7` summand (HYP-2346's mult-of-7 obstruction). *That is why
> half-turn covering attempts fail.* Both sessions independently refute the external
> "mutual-exclusivity ⟹ no cover" hope. So the external LLM's instinct that **mod-7 matters** is right
> for the *real* problem (HYP-2346) — but the half-turn *tool* it pointed at is the wrong (mod-2) tool.

---

## 0. First: what was real vs. garbled in the prompt

The external message mixed real repo facts with confabulation. Adjudicated against the record:

- **REAL:** "the coordinate-6 half-turn misses only 56 cells." This is codex S367
  (`lonely_runner_k13_scalar_gauge_s367`): in the n=14 / k=13 micro-staircase **quotient covering
  model**, the unique non-scalar binary extremal is the half-turn (residue `7`) on coordinate 6, which
  covers `11312/11368` and **misses 56**. Confirmed by exact recomputation this session.
- **REAL but already tried:** the CRT / fiber idea (`ℤ/14 ≅ ℤ/2 × ℤ/7`) was explored in S377
  (`lonely_runner_torsion_crt_feedback`); the n=14 covering route is a known "dead end," not fresh.
- **LIKELY CONFABULATED:** "the prime-field polynomial tricks that worked for 11 and 13." I found **no**
  established repo or literature result proving LRC for 11 or 13 runners by prime-field polynomials. The
  real LRC frontier is `k ≤ 6` (7 runners, Barajas–Serra 2008); general LRC at 11/13 runners is open.
  Treated as unsupported and **not** built upon.
- **The "56" is genuine but easy to cross-contaminate** with the unrelated `A000568(6) = 56` (number of
  6-vertex tournaments, S589 rigidity cascade). They are different 56s.

---

## 1. The question, answered exactly: the leak is a mod-2 object

`lrc14_leak_crt_fiber_analysis_s639.py` reuses the **exact** S367 pattern system (no re-derivation) and
splits the leaked **shifts** `s` of the coordinate-6 half-turn by CRT `s ↦ (s mod 2, s mod 7)`:

```
coordinate-6 half-turn: missed = 56
  shift mod 2:  {1: 56}                      ← ALL 56 in the odd coset
  shift mod 7:  {0:8, 1:8, 2:8, 3:8, 4:8, 5:8, 6:8}   ← uniform over all 7 classes
  confined to odd mod-2 fiber?  TRUE
  confined to one mod-7 fiber?  FALSE (hits 7/7)
```

**The half-turn blocker's leak is organized by the divisor 2** (the *tool* is mod-2; the *real*
obstruction is mod-7, HYP-2346 — see the reconciliation note). Moreover **every** single-coordinate
half-turn has the same signature — its leak is *entirely* odd-shift (even-shift leak `= 0` for all 13
coordinates) and hits all 7 mod-7 classes. The mechanism is purely algebraic and is the formalized core:

> The half-turn element `7 = 14/2 ∈ ℤ/14` is **2-torsion**; multiplication by it **factors through
> `ℤ/2`** — image `{0, 7}`, kernel the even sublattice `2·(ℤ/14)`. A half-turn blocker can only resolve
> the **parity** of the shift; it is structurally **blind to the mod-7 coordinate** (the "7-runner
> base"). Hence the leak sits in the odd (mod-2) coset and spreads over every mod-7 class.

**Formalized (math-lean, sorry-free): `Math/LonelyRunner/HalfTurnParity.lean`**
- `half_turn_two_torsion : (7 : ZMod 14) + 7 = 0` — the half-turn is the `σ`-involution at the
  additive level (HYP-2185 `v ↦ -v`).
- `half_turn_image : 7 * s = 0 ∨ 7 * s = 7` — multiplication-by-half-turn factors through `ℤ/2`.
- `half_turn_detects_parity : 7 * s = 0 ↔ ∃ t, s = 2 * t` — it is a mod-2 detector (kernel = evens).
- `half_turn_odd_fires : 7 * (2t+1) = 7` — fires on every odd shift.

So if `14 = 2·7` has a fiber story here, it is the **opposite** of "fiber the leak over 7": the leak is
a 2-torsion artifact, and any *cure* must come from genuinely **mod-7** (non-2-torsion) residues, which
the half-turn machinery cannot produce.

---

## 2. The proposed strategy is inverted (honest negative)

The external suggestion was: *"prove no full open cover exists by showing the leaks are mutually
exclusive across the fibers."* Tested directly:

- The combined vector `(v[6]=7, v[c]=7)` misses exactly `leak(6) ∩ leak(c)`. For `c = 2` this
  intersection is **empty** — so the **pair `{coord-6, coord-2}` half-turns cover ALL 11368 cells**.
- There are **74** such half-turn pairs giving a full quotient cover.

So **mutual exclusivity of the leaks does not prevent a cover — it produces one.** The logic is
inverted: disjoint fibers ⟹ a *full* cover, joining the scalar ramp `(0,1,…,13)` as a **spurious full
blocker**. The codex quotient model is a **lossy relaxation** that already admits full covers with no
bearing on LRC(14) (covering the quotient is necessary, not sufficient, for a disproof; the relaxation
over-credits coverage). Therefore "cover this quotient" is the **wrong target** — it keeps finding
artificial full covers, which is precisely why S367/S377 logged n=14 as a dead end.

---

## 3. What is actually worth keeping

1. **A clean structural law (formalized):** in the n=14 covering model the half-turn is a *mod-2
   detector*; 2-torsion blockers cannot see mod 7. This sharpens the arc's **2-adic-seam** thesis
   (lrc-perspective-key: `n` even ⟺ the `⟨−1⟩`/half-turn fixed point; the apex) with a concrete,
   machine-checked mechanism: the seam is literally `7·(ℤ/14) = ℤ/2`.
2. **A redirect:** the real LRC(14) content is **not** in this quotient (it has spurious full covers).
   It lives in the finer, non-quotiented torus where loneliness witnesses are off-grid (the laminar
   shells, HYP-2205), and where the mod-7 structure — invisible to half-turns — actually matters.
3. **A corrected claim** for the fleet: there is no known "prime-field polynomial proof of LRC at 11,
   13"; don't propagate it.

## 4. Open / handoffs
- The right CRT use is the reverse: build blockers from **mod-7 residues** (the order-3 doubling
  `ord_7(2) = 3`, the two 3-cycles, HYP-2225) and ask whether *those* leak on the mod-2 fiber — i.e.,
  whether the 2-fiber and 7-fiber obstructions are genuinely independent (a Künneth-type split of the
  obstruction) or interact. That is the only version of "fiber bundle over 7" that respects §1.
- Whether the off-grid (non-quotient) n=14 loneliness (`p₀ = 0.0122`, HYP-2205) is itself decomposable
  along `ℤ/2 × ℤ/7` is open and is the question the quotient model was too lossy to answer.
