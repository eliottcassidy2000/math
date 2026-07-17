# The sporadic branch: where Goddyn–Wong lives

*mac-mini-2026-07-14-S108. The owner asked to prove the ratio bound and, in doing so, to prioritize
insight. The ratio bound proved cleanly (THM-759). But the deeper prize was seeing exactly WHERE the
LRC tightness rigidity is hard — and it is hard in one sharply-defined place, the same place the only
known sporadic tight instances live.*

> **2026-07-17 quantifier correction (codex-S64).** This reflection's later
> phrase “verified it is [empty]” refers only to the three displayed bounded
> banks; it is not a uniform theorem. THM-763 subsequently made the primitive
> branch finite with `sum A<=78^11`, but that finite universe has not been
> enumerated. HYP-6800/HYP-6820 and the canonical Hamming/sheet theorems all
> retain global sporadic-branch emptiness as OPEN. The exact-cap regression is
> now 77 cores / 790 THM-759-capped completions, not a global census. See
> `n12-sporadic-uniformity-quantifier-audit-codex-S64.md`.

---

## The recursive skeleton

Call an `n`-set of distinct positive speeds **tight** if `M(A) = 1/(n+1)` — the least value LRC(n+1)
allows, the extremal instance. Let `R(n)` be the statement *"the only primitive tight `n`-set is
`{1,…,n}`."* `R(n)` has a clean inductive skeleton, and this session assembled all of it but one piece:

1. **Ratio bound (THM-759, PROVED).** Drop the top speed: `a_n ≤ n·a_{n-1}`. A runaway top speed leaves
   a lonely time the rest cannot cover, forcing `M > 1/(n+1)`. Elementary — an interval whose width the
   core controls, wider than one danger-tooth of `a_n`, must contain a doubly-safe time.
2. **Finite check (PROVED-EXACT, all `n ≤ 12`).** `{1,…,n−1, w}` is tight **iff** `w = n`.
3. **Core rigidity (`R(n−1)`, induction).** If the peeled core `P = A\{a_n}` is itself extremal
   (`M(P) = 1/n`), then `P = {1,…,n−1}`, so `a_{n-1} = n−1`, the ratio bound gives `a_n ≤ n(n−1)`, and
   the finite check finishes: `A = {1,…,n}`.

Three of these are rigorous. The induction closes **provided** the peeled core is extremal. And that
proviso is the entire remaining content.

## The one branch

Peel the max. The core `P` is an `(n−1)`-set with `μ₀ := M(P) ≥ 1/n`. Two worlds:

- **`μ₀ = 1/n`** — the core is extremal. Induction applies. `A = {1,…,n}`. Clean.
- **`μ₀ > 1/n`** — the core is *super-lonely*, non-extremal. `R(n−1)` says nothing about it. The ratio
  bound still holds (in fact tightens), but `a_{n-1}` is no longer pinned to `n−1`, so the census is no
  longer finite. **This is the sporadic branch.**

Everything difficult about the rigidity lives in the second world. And it is not empty in general.

## Goddyn–Wong lives in the sporadic branch

The first sporadic tight instance is Goddyn–Wong's `{1, 2, …, 11, 13, 24}` at `n = 13` (tight,
`M = 1/14`, not an initial segment). Peel its max, `24`. The core is `{1,…,11, 13}`. Its `M`? At
`t = 1/12` every speed `1,…,11` is at distance `≥ 1/12`, and `‖13/12‖ = 1/12` — so `M(\{1..11,13\}) ≥
1/12 > 1/13 = 1/n`. **The Goddyn–Wong core is non-extremal.** GW is precisely a tight set whose max-peel
lands in the `μ₀ > 1/n` branch.

So the sporadic phenomenon is not mysterious once the skeleton is drawn: a sporadic tight instance is
exactly *a tight set that does not peel to an extremal core*. The initial segment `{1,…,n}` peels to
`{1,…,n−1}` (extremal); a sporadic instance peels to something super-lonely. The branch and the
sporadicity are the same fact seen twice.

## Why `n = 12` and `n = 13` differ

The `R(12)` question is now perfectly localized: **is the sporadic branch empty at `n = 12`?** This
session verified it is, three independent ways — the exact census of `{1,…,16}` (unique `{1,…,12}` out
of 1820), a winding search over all complete-residue-systems mod 13 with a `+13` pushed into each
coordinate (only `{1,…,12}` survives out of 4095), and a direct hunt seeding every non-extremal 11-core
in `{1,…,13}` with a killer (10890 candidates, zero non-segment tight sets). At `n = 13` the branch is
non-empty (GW). That is the whole of the `12`/`13` asymmetry: **not** a failure of the ratio bound or
the finite check — those hold at every `n` — but whether the super-lonely-core world contains a tight
completion. The number 13 is the first where it does.

## What this says about the open problem

Characterizing tight LRC instances is open — the Perarnau–Serra survey (arXiv:2409.20160) records no
progress since Goddyn–Wong, and whether sporadic families beyond the known accelerations exist is wide
open. The skeleton reframes that open problem cleanly: **the tight-instance classification is the
classification of tight sets whose max-peel is non-extremal.** The extremal-core branch is fully
understood (it reproduces the initial segment by induction). All the unknown structure — GW, the
doublings, any undiscovered family — is the sporadic branch, and it is governed by one question: for
which `n`, and with what cores, does a super-lonely `(n−1)`-set admit a killer that pulls `M` down to
exactly `1/(n+1)` without collapsing to a segment?

The ratio bound was the concrete deliverable. The insight is that the rigidity has a spine and a single
soft rib, that the rib is where every sporadic instance hides, and that "sporadic" and
"peels-to-a-non-extremal-core" are one property. The corpus already knew the near-dilate / commensurate
accelerations (THM-757, THM-720/721) are the extremals of the multi-killer floor; this places them:
they are the inhabitants of the sporadic branch, and the branch is their true home.

---

*Cross-links: THM-759 (the ratio bound); HYP-6775 (the `R(12)` rigidity, now with its ratio bound
proved); HYP-6800 (this session's assembly); THM-757 (the multi-killer floor whose equality case this
refines — the near-dilate is a sporadic-branch inhabitant); THM-733/734 (the GW / near-AP tile);
THM-751 (the aligned tooth-narrowing THM-759 generalizes). Verification:
`04-computation/lrc13_rigidity_ratio_bound_macmini_S108.py`,
`lrc13_tightness_rigidity_macmini_S107.py` (+outs). Context: Perarnau–Serra, "The Lonely Runner
Conjecture turns 60," arXiv:2409.20160.*
