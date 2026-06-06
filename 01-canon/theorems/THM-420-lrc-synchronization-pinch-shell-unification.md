---
id: THM-420
title: The synchronization lemma — binding pair = shell-partner synchronized at the witness; λ_q is a sign-blind, shell-determined lower bound, and M = max_q λ_q is lossless
status: PROVED (L0 synchronization elementary; L1 a corollary of the pinch lemma HYP-2059; L2 lossless verified exactly 300/300; the L_n iff PROVED)
source: monad-explorer-2026-06-06-S1
depends_on:
  - HYP-2059   # pinch lemma: M attained at pair-sum time t*=m/(v_a+v_b), binding pair STRADDLES
  - THM-401    # pair-sum modulus C=2n-1 (the SECOND gap); summand shells
related:
  - HYP-2262   # opus-S699 signed LRC: sign=cut, zero-clock ⟺ shell-partner (v_i+v_j≡0 mod 2n-1)
  - THM-404    # doubling rigidity (×2 invariance, even-n fragmentation)
  - THM-406    # covering-depth moments; M2 all-orders cancellation (the second-gap face)
  - THM-407    # scaling/time-reversal invariance of M; G=⟨2,−1⟩ on shells
  - HYP-2283   # open: is the V*-type second-gap datum exactly the THM-406 M2 cancellation?
---

# THM-420 — synchronization: the binding pair is a shell-partner, and the lattice bound λ_q

**Setup.** Speeds `S = {v_1,…,v_{n-1}} ⊂ ℤ_{>0}`; `‖x‖` = distance to the nearest integer;
`M(S) = max_t min_i ‖v_i t‖`; LRC(n) ⟺ `M(S) ≥ 1/n` for every `S`. Lattice
`L_q = {k/q : k ∈ ℤ}`.

This theorem makes one structural fact rigorous and uses it to **unify** the pinch lemma
(HYP-2059) with opus-S699's signed-LRC shell-partner (HYP-2262/T3): they are the *same object
read at two different moduli*. It also extracts a clean, sign-blind, shell-determined lower bound.

---

## L0 (synchronization lemma) — PROVED

> If `v_a + v_b ≡ 0 (mod q)` then `‖v_a · k/q‖ = ‖v_b · k/q‖` for **every** integer `k`.

**Proof.** `v_a + v_b = q·m` for some `m ∈ ℤ`, so `v_b · k/q = (q m − v_a)·k/q = mk − v_a k/q`,
hence `v_b k/q ≡ −v_a k/q (mod 1)`. Since `‖−x‖ = ‖x‖`, `‖v_b k/q‖ = ‖v_a k/q‖`. ∎

So a `q`-**shell-partner** pair `{v_a, v_b}` (with `v_a+v_b≡0 mod q`) is **synchronized** on the
lattice `L_q`: the two runners are mirror images of `0` at every lattice time, equidistant from the
observer. On `L_q` they impose the *same* loneliness constraint — one is removable.

**Verified:** `6669` exact tests, `q ∈ [2,40)`, all lifts — `0` failures
(`lrc_lattice_synchronization_monad_s1.py`).

---

## L1 (binding pair = shell-partner synchronized at the witness) — corollary of HYP-2059

By the **pinch lemma** (HYP-2059) the maximum `M(S)` is attained at a **pair-sum time**
`t* = m/(v_a+v_b)`, where the *binding pair* `(a,b)` **straddles** the observer:
`frac(v_a t*) = M`, `frac(v_b t*) = 1 − M`. Equivalently `(v_a+v_b)·t* ≡ 0 (mod 1)`, i.e.
`v_a + v_b ≡ 0 (mod q)` with `q = v_a+v_b` and `t* ∈ L_q`. Then by L0:

> **The binding pair of `M(S)` is a `q`-shell-partner (`q = v_a+v_b`), and at the witness `t*` it is
> synchronized: `‖v_a t*‖ = ‖v_b t*‖ = M(S)`.**

**Verified** for every worry-set config tested (`n=4,5,7,14`):
the binders sit at distance `M`, sum `≡ 0 mod q`, and synchronize at `M`. E.g.

```
 AP  n=14:  M=1/14, t*=13/14, binding pair (1,13),  1+13=14≡0 mod 14
 V*  n=14:  M=1/14, t*= 3/14, binding pair (5,9),   5+ 9=14≡0 mod 14
 2AP n=14:  M=1/14, t*=13/28, binding pair (2,26),  2+26=28≡0 mod 28
```

**The unification.** opus-S699's "shell-partner / signed zero-clock" (`v_i+v_j ≡ 0 mod 2n−1`,
HYP-2262/T3) and the pinch lemma's "straddling binding pair" are the **same synchronization
phenomenon at different moduli**:
- the **floor** binding pair has `q ≡ 0 (mod n)` (smallest `q = n`, value `1/n`);
- opus's **`2n−1` shell-partner** is the **Farey-successor / SECOND-gap** binding pair
  (`q = 2n−1`, value `2/(2n−1)`, THM-401(i)).

The signed framework is the **`q = 2n−1` slice** of the same binding-pair lattice structure.

---

## L2 (the lattice lower bound λ_q) — PROVED; lossless verified

Define `λ_q(S) = max_{k=1,…,q−1} min_i ‖v_i k/q‖`. Since `L_q ⊆ [0,1)`,

> `M(S) ≥ λ_q(S)` for every `q`.

By L0, `‖v_i k/q‖` depends only on the **shell** `s_i = min(v_i mod q, q − (v_i mod q))`, so

> **`λ_q(S)` is sign-blind and depends only on the shell-multiset of `S mod q`**:
> `λ_q(S) = max_k min_{s ∈ Sh_q(S)} ‖s·k/q‖`, where `Sh_q(S)` is the set of occupied shells.
> In particular `λ_q(S) = λ_q(S°)` for any shell-transversal `S°` (one speed per occupied shell).

This is the rigorous bridge that *explains* gauge invariance (opus T1) at the lattice level: the
sign-blind, shell-determined quantity `λ_q` is a genuine lower bound on `M`.

**Lossless.** The maximizer of `min_i‖v_i t‖` is a pinch breakpoint `t = m/(v_a+v_b)`, so it lies in
`L_q` for `q = v_a+v_b`. Hence

> `M(S) = max_{q ∈ {v_a+v_b}} λ_q(S)`   (max over pair-sum denominators).

**Verified exactly:** `300/300` random configs, `M == max_q λ_q` (no loss).

---

## L3 (the single-lattice reach — a clean iff) — PROVED

> The **primary** lattice `L_n` certifies the floor for `S` (i.e. `λ_n(S) ≥ 1/n`) **iff**
> `n ∤ v_i` for all `i`. When it certifies, `k = 1` suffices (`t = 1/n`).

**Proof.** For `j ∈ {1,…,n−1}`, `‖j/n‖ = min(j,n−j)/n ≥ 1/n`; and `‖0‖ = 0`. So
`min_i ‖v_i k/n‖ ≥ 1/n ⟺ v_i k ≢ 0 (mod n) ∀i`. If `n ∤ v_i ∀i`, then `k=1` gives
`v_i·1 ≢ 0 (mod n)`, so `λ_n ≥ 1/n`. Conversely, if `n ∣ v_i` for some `i`, runner `i` sits on `0`
at **every** `k/n` (`v_i k ≡ 0`), so `min_i ‖·‖ = 0` for all `k`, hence `λ_n(S) = 0`. ∎

So the **hard core of LRC is exactly the configs containing a multiple of `n`** (where the primary
lattice collapses and a finer/other pinch lattice is required). **Verified** (random sampling,
exact residue criterion): `L_n` certifies `≈ 44–47%` of configs at `n = 6,8,12,14,15,30`
(matching `P(no v_i ≡ 0 mod n)`). Example of the collapse: `2·AP` has `2·7 = 14 ≡ 0 (mod 14)`, so
`λ_14(2·AP) = 0`; its floor lives on `L_28` instead (binding pair `(2,26)`), consistent with the
scaling invariance `M(2S)=M(S)` (THM-407 Lemma A) **rescaling the witness lattice `n → 2n`**.

---

## Corollary (the structural-reduction NO-GO for the signed LRC)

`M` is sign-blind (opus T1), so the signed predicate **equals** the unsigned one: there is **no**
separate "signed counterexample." The signed/shell datum `v_i+v_j ≡ 0 (mod 2n−1)` is a
**second-gap** binding pair (L1/L3), **dominated** by the floor pinch `q = n` and therefore
**cannot lower `M` below `1/n`**. A `2n−1` shell-partner forces a counterexample only if it is *also*
a **floor** binding pair — i.e. `v_i+v_j ≡ 0 (mod n)`, a *different* modulus. Hence:

> **The signed structure stratifies the tight boundary (AP-type `λ_{2n−1}` shell-partner-free vs
> V*-type carrying one) but is orthogonal to the violation `M < 1/n`. Both AP and V* have
> `M = 1/14`; the split is a `q = 2n−1` second-gap effect, invisible to `M`.**

This bounds what the signed angle can contribute: not a new predicate, but a finer invariant of the
worry-set living at the second gap (see HYP-2283).

---

## Artifacts

`04-computation/lrc_lattice_synchronization_monad_s1.py` (+ `.out`).
Builds on HYP-2059 (pinch), THM-401 (`2n−1`/shells), THM-404/406/407, and opus-S699 HYP-2262
(signed LRC / shell-partner). New open thread: **HYP-2283**.
