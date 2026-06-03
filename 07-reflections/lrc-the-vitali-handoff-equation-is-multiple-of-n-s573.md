---
source: opus-2026-06-03-S573 (remote-control)
status: CONNECTION (unifies S551o Vitali wall with the S571/S572 reduction) + quantified: Criterion B′ = the Vitali-covering iff proves 72–97% of C′ (→97% at n=14)
tags: [LRC, vitali, vitali-covering, lebesgue-density, measure-vs-construction, THM-398, C-prime, handoff, n14, S551o]
---

# The Vitali handoff equation is `n | v`

**Prompt (user):** see the connection to Vitali; search the repo for inspiration.

The repo already had the Vitali lens: **S551o (oracle), "the Vitali wall."** It
diagnosed that LRC has a *measure-vs-set* boundary —

> *measure settles the positive-measure bulk; explicit construction settles the
> measure-zero core; the Vitali set marks the handoff — where measure ends,
> arithmetic construction must begin.*

S551o left the handoff **unlocated** ("Vitali marks it"). My S571/S572 reduction
(THM-398) **locates it exactly.** The two threads are the same theorem.

## 1. THM-398's split IS the Vitali wall — and the marker is `n | v`

THM-398 decides every config by one test: *does `S` contain a multiple of `n`?*

| side | test | tool | the set it lives on |
|---|---|---|---|
| **construction** | `n ∤ v` for all `v` | `t=1/n` witness (THM-369, *measure-blind*) | incl. the **measure-zero core** (the AP/worry-set has no multiple of `n`) |
| **measure** | `n ∣ v` for some `v` | C′: positive measure (Vitali covering) | the **positive-measure** side |

This is precisely S551o's dichotomy, with the boundary made into a *decidable
arithmetic equation*:

> **`n | v` is the equation of the Vitali handoff.** A config with no multiple of `n`
> is handled by *construction* (`t=1/n` exists even when the safe set is
> measure-zero — exactly the AP, the Vitali-blind core). A config *with* a multiple of
> `n` is pushed to the *measure* side, where C′ asserts positive measure.

The worry-core (measure-zero, where measure is Vitali-blind) is *characterised* by
`n ∤ v` (S564's C′-for-tight): it is **always** on the construction side, so the
`t=1/n` sieve — S551o's "constructive ℚ-coset selection" — *always* reaches it. The
handoff is not just marked; it is the divisibility test.

## 2. The OTHER Vitali — the covering lemma — is the tool on the measure side

S551o invoked the **Vitali set** (the pathology that proves measure can't see the
core). The measure side needs the **Vitali covering lemma** (the constructive cousin):
*a regular family of arcs covering a set has a disjoint subfamily catching a definite
fraction of its measure* — the engine of the **Lebesgue density theorem**.

Why it applies to `n | v` and **not** to the worry-core: a multiple `v = nw` has
danger `D_v` = a **periodic, bounded-eccentricity arc family** (centres `k/(nw)`,
radius `1/(n²w)`, gaps `(n-2)/(n²w) > 0`) — a *genuine Vitali cover*. The worry-core's
"lonely set" is `n` **isolated points** (`{k/n}`) — no Vitali-cover structure, hence
measure-blind. So the multiple of `n` is exactly what turns the question into one the
Vitali covering lemma can answer.

**Lebesgue-density form of the dodge.** Let `E = G(S\{v})` (safe set of the other
`n-2` runners), `μ(E)>0`. If `E ⊆ D_v` (i.e. `S` tight, measure-zero), then at a
density point `x` of `E`, for `r ≫` one period, `μ(D_v∩B(x,r)) ≥ (1-ε)2r`, but the
periodic family gives `μ(D_v∩B(x,r)) ≤ (2/n)2r + O(ρ)` — contradiction since `2/n<1`.
This *re-proves the dodge as a density statement* and shows where it strains: it needs
`E` dense at the scale of `v`'s period `1/(nw)` — automatic for large `w`, the crux
for small `w`.

## 3. The sharp form: Criterion B′ is the Vitali-covering *iff*

Because `E = G(S\{v})` is a finite union of intervals and `D_v`'s gaps are open and
nonempty, **an interval of `E` lies in `D_v` iff it fits inside a single arc**
(`length ≤ 2/(nv)`). Hence:

> **`S` is loose ⟸ some component of `G(S\{v})` is longer than `2/(nv)`** (Criterion
> B′, PROVED by Vitali covering: a too-long interval cannot hide inside one arc).
> And `S` is tight (measure-0) **⟹** every component is short *and* arc-aligned.

So the *residual* of C′ is exactly the **arc-alignment** case: every `E`-interval
short **and** each one sitting inside a `v`-arc.

## 4. Quantified — how much of C′ the Vitali criterion already proves

`lrc_vitali_covering_residual_s573.py`, mult-of-`n` configs:

| n | proved by B′ (long interval) | residual (all-short, still loose) | all-short **tight** |
|---|---|---|---|
| 6 | 72.4% | 331 | **0** |
| 8 | 78.7% | 255 | **0** |
| 10 | 88.9% | 133 | **0** |
| 12 | 91.5% | 102 | **0** |
| **14** | **96.8%** | **39** | **0** |

Two facts:
- **The proved fraction grows toward the frontier — 96.8% at n=14.** The Vitali
  covering criterion handles all but ~3% of the multiple-of-14 configs outright.
- **The all-short residual is *never* tight (0 across all n).** Even when every
  interval is short, the periodic arcs fail to align to cover them — C′ holds there
  too, by alignment, not yet proved.

## 5. What the connection buys

1. **Unification:** S551o (the Vitali diagnosis) and S571/S572 (THM-398) are one
   result; the handoff S551o located only abstractly is the equation `n | v`.
2. **The right tool, named:** on the measure side, the Vitali **covering** lemma /
   Lebesgue density is the engine — and it *applies* precisely because a multiple of
   `n` makes the danger a regular Vitali cover (the core is points, so it cannot).
3. **The residual is tiny and sharp:** the arc-**alignment** of a single periodic
   family against `G(S\{v})`'s short intervals — `≤3%` of configs at n=14, all loose.
   This is a Diophantine alignment question (the interval centres of `G(S\{v})`, which
   are rationals with `S\{v}`-denominators, would all have to be `1/(n²w)`-close to
   the `1/(nw)`-lattice), ripe for an Erdős–Turán / three-distance bound.

## 6. Honest status

The unification (THM-398's split = the S551o Vitali wall; `n|v` = the handoff) is a
*framing*, rigorous where THM-398 is. Criterion B′ as the Vitali-covering iff is
**proved** (and is now the engine for 72–97% of C′). The all-short **alignment**
residual is **open** (verified loose, 0 tight, `≤3%` at n=14). Net: the Vitali lens
both explains *why* the multiple-of-`n` reduction is the natural one (it is the
measure/construction handoff) and supplies the covering-lemma tool that already
discharges almost all of it.

**Artifacts:** `04-computation/lrc_vitali_covering_residual_s573.py` (+`.out`).
Builds on S551o (Vitali wall), THM-398/HYP-2103 (the dodge + criterion), HYP-2102
(reduction), THM-369 (construction witness), S564 (C′-for-tight). New: **HYP-2104**.
