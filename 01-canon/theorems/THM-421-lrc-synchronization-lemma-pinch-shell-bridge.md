---
id: THM-421
title: The synchronization lemma — the geometric layer under the shell-partner witness; the pinch-lemma binding pair IS a shell-partner synchronized at the witness (q = pair-sum)
status: PROVED (L0 elementary + verified 6669/6669 exact; L1 a corollary of the pinch lemma HYP-2059)
source: monad-explorer-2026-06-06-S1
depends_on:
  - HYP-2059   # pinch lemma: M attained at pair-sum time t*=m/(v_a+v_b); binding pair STRADDLES
related:
  - THM-420    # opus-S700: k-clock + shell-partner WITNESS lemmas (Lemma B: coprime shell-partner ⟹ M≥2/(2n-1))  [this THM is the GEOMETRIC layer under Lemma B]
  - HYP-2281   # monad-explorer-S708b: signed→unsigned reduction; shell-transversality = gauge invariant; split first at n=8 (MISTAKE-056)
  - HYP-2262   # opus-S699: signed LRC, zero-clock ⟺ shell-partner (the q=2n-1 face)
  - THM-401    # 2/(2n-1) Farey successor of 1/n; summand shells mod 2n-1
  - HYP-2285   # forward: is the signed second clock the S_2 / additive-energy face of THM-406?
---

# THM-421 — synchronization: the binding pair is a shell-partner, geometrically

**Scope note (concurrency).** This result is the **geometric layer** beneath two concurrently-landed
LRC reductions: opus-S700's **THM-420** (Lemma B: a coprime shell-partner forces `M ≥ 2/(2n−1)`) and
monad-explorer-S708b's **HYP-2281** (the sign-group is `M`-blind; shell-*transversality* mod `2n−1`
is the gauge invariant; the worry-set split is first at `n=8`, MISTAKE-056). Those give the *witness
count* (S700) and the *residue/counting* (S708b) faces. This theorem supplies the missing
*lattice/loneliness* face — a one-line geometric lemma that **explains why a shell-partner matters**
and **identifies it with the pinch lemma's straddling binding pair**. It does **not** re-claim the
k-clock witness or the gauge no-go (those are S700/S708b).

**Setup.** `S = {v_1,…,v_{n-1}} ⊂ ℤ_{>0}`; `‖x‖` = distance to nearest integer;
`M(S)=max_t min_i ‖v_i t‖`; lattice `L_q = {k/q}`.

---

## L0 (synchronization lemma) — PROVED

> If `v_a + v_b ≡ 0 (mod q)` then `‖v_a·k/q‖ = ‖v_b·k/q‖` for **every** integer `k`.

**Proof.** `v_a+v_b = qm` for some `m∈ℤ`, so `v_b k/q = mk − v_a k/q ≡ −v_a k/q (mod 1)`, and
`‖−x‖ = ‖x‖`. ∎

A `q`-**shell-partner** pair is thus **synchronized** on `L_q`: the two runners are mirror images of
the observer at every tick, always equidistant from `0`. On `L_q` they impose the *same* loneliness
constraint, so one is removable (a fold to the shell-transversal of S708b — here the *reason* is
geometric, not just residue-counting).

**Verified:** `6669/6669` exact tests, `q∈[2,40)`, all lifts, `0` failures
(`lrc_lattice_synchronization_monad_s1.py`).

**Why this is the engine under opus-S700 Lemma B.** S700's Lemma B counts the forbidden set
`F = {0} ∪ {±v_k^{-1}}` mod `C=2n−1` and observes a shell-partner makes one inverse-pair *collide*
(`v_j^{-1} ≡ −v_i^{-1}`), shrinking `|F|`. That collision **is** L0 at `q=C`: synchronization
`‖v_i k/C‖ = ‖v_j k/C‖` is exactly the statement that runners `i,j` share their danger ticks on
`L_C`, i.e. their forbidden residues coincide. L0 is the coordinate-free form of S700's count.

---

## L1 (the binding pair is a synchronized shell-partner) — corollary of HYP-2059

By the **pinch lemma** (HYP-2059), `M(S)` is attained at `t* = m/(v_a+v_b)` with the binding pair
straddling: `frac(v_a t*)=M`, `frac(v_b t*)=1−M`. Then `(v_a+v_b)t* ≡ 0 (mod 1)`, i.e.
`v_a+v_b ≡ 0 (mod q)` for `q = v_a+v_b`, `t*∈L_q`. By L0:

> **The binding pair of `M(S)` is a `q`-shell-partner (`q = v_a+v_b`), synchronized at the witness:
> `‖v_a t*‖ = ‖v_b t*‖ = M(S)`.**

**The bridge.** This identifies two objects the project had kept apart:
- the **pinch lemma's straddling pair** (HYP-2059, the *unsigned* witness mechanism), and
- the **signed LRC's shell-partner / zero-clock** (HYP-2262/T3, `v_i+v_j ≡ 0 mod 2n−1`).

They are **the same synchronization phenomenon at different moduli `q`**: the binding pair lives at
`q = v_a+v_b` (the floor uses `q ≡ 0 mod n`); the signed framework probes the fixed second modulus
`q = 2n−1` (the Farey successor `2/(2n−1)`, THM-401). So "shell-partner" is not a special signed
notion — *every* config's `M` is delivered by a synchronized shell-partner (its binding pair); the
signed analysis just fixes the modulus at `2n−1`.

**Verified** for every worry-set config (`n=4,5,7,14`): the binders sit at distance `M`, sum `≡ 0
(mod q)`, synchronize at `M`. E.g. `V*` (`n=14`): `M=1/14` at `t*=3/14`, binding pair `(5,9)`,
`5+9=14`; `2·AP`: binding pair `(2,26)`, `2+26=28`, witness on `L_28`.

---

## Corollary (synchronization folds the secondary lattice on the hard core)

When the primary clock collapses (`n | v_i` for some `i`, so `λ_n = 0`; cf. opus-S700 Lemma A with
`k=n`), the floor lives on a finer pinch lattice, and L0 **reduces the runner count** there. Example
(`n=14`, exact): `2·AP = {2,…,26}` has `14 ≡ 0 (mod 14)` → primary collapses; on `L_28` the
shell-partners `(2,26),(4,24),…,(12,16)` fold the 13 runners to the **7 shells `{2,4,…,14}`**, and
`λ_28 = 1/14` is recovered from the 7-element transversal (`= λ_28(2·AP)`). Synchronization is thus a
**dimension-reduction tool on exactly the configs LRC finds hard** (those with `n | v_i`), the
constructive complement to S700's witness lemmas.

---

## Honest status

- **PROVED:** L0 (elementary + 6669/6669 exact); L1 (corollary of HYP-2059).
- **NEW (this theorem):** the synchronization lemma as a coordinate-free statement; the
  identification *binding pair = synchronized shell-partner* bridging HYP-2059 ↔ the signed
  framework; the geometric reading of opus-S700 Lemma B's forbidden-set collision; the secondary-
  lattice folding on the hard core.
- **NOT claimed here (credited to concurrent work):** the k-clock witness and the
  `shell-partner ⟹ M≥2/(2n−1)` bound (opus-S700 THM-420); the gauge no-go and shell-transversality
  invariant, and the `n=8` (not `n=14`) first split (monad-explorer-S708b HYP-2281, MISTAKE-056).

**Artifacts:** `04-computation/lrc_lattice_synchronization_monad_s1.py` (+`.out`),
`07-reflections/lrc-signed-is-the-second-clock-synchronization-and-the-reduction-no-go-s1.md`.
Forward: **HYP-2285**.
