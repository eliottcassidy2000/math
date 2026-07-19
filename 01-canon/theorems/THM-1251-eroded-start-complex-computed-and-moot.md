---
id: THM-1251
title: "THE ERODED START COMPLEX E_k(P), computed exactly for the first time -- and MOOT. (0) PRIORITY: supplier (3) of THM-1162 s3 / THM-1203 s12 lives on the clustered r=5 stratum, and codex's OWN next session closed that stratum uniformly (THM-1214, PROVED, Lean-backed, referee verdict 'clustered r=5 is uniformly lonely'). The s12 residue list is STALE, and my THM-1245 repeated it this morning -- corrected here. (I) Computed anyway, since nobody ever had: min |S(P)| = 31807/194040 EXACTLY at P = (1,2,3,5,7,8,9,11), pinning the constant carried as a measured '0.164' since S128c73 and never frozen. (II) |E_w(P)| > 2/21 holds for ALL 495 cores iff k1 >= 204, and FAILS for exactly 19 cores below that -- so the erosion bound is FALSE as stated on part of the clustered range, and the crude |S| - N*w is vacuous in the worst case. (III) But the 456 undecided (core,k1) residue pairs are NOT obstructions: over ~160,000 killer triples the full safe set is never empty, so the COUNTING was lossy, not the geometry"
status: "(0) VERIFIED against canon: THM-1214 is uncontested, cited approvingly four times in LRC14-PROOF-MAP ('THM-1214 closes the entire eight-core/five-killer clustered stratum'), carries a Lean file and a referee output, and appears in no MISTAKES entry. (I) PROVED by exhaustion -- an exact rational computation over all 495 eight-speed cores. (II) PROVED by exhaustion, and independently reproduced by a second agent to the same threshold and the same 19 cores. (III) VERIFIED over ~160,000 killer triples with k4 <= 235 (the r=5 horn bound), NOT proved -- it is a strong negative result about my own criterion, not a theorem. Nothing here advances LRC(14), which remains open; it retires an obligation and freezes two constants"
source: kind-pasteur-2026-07-19-S128c85 (owner: work the eroded start complex E_k(P))
depends_on:
  - THM-1214    # codex-S78: closes the clustered r=5 stratum, making supplier (3) moot
  - THM-1162    # where supplier (3) was named
  - THM-1203    # whose s12 residue list is stale
related: [THM-1245, THM-1123, THM-1168, THM-1213, MISTAKE-183]
script: 04-computation/eroded_start_complex_kps_S128c85.py, eroded_residue_close_kps_S128c85.py, eroded_residue_direct_kps_S128c85.py (+ .out)
---

# THM-1251 — the eroded start complex, computed and retired

## 0. PRIORITY — supplier (3) is moot

THM-1162 §3 and THM-1203 §12 name three "uniform suppliers" for uniform `r=5`.
Supplier (1), the continuum ceiling `bad <= 2/21`, is codex-S77's THM-1203.
Supplier (3) is:

> "an **eroded core-atlas bound** ensuring that a safe phase outside the bad set
> contains a **complete** `k1`-gap, rather than merely a safe point."

That supplier lives entirely on the clustered `r=5` stratum: `P = {p1<...<p8}` a
core in `{1,...,12}`, killers `13*max(P) < k1 < ... < k5`. **codex's very next
session closed that stratum outright.** THM-1214 (codex-S78):

> **Theorem.** Every family (1) has a time at which all thirteen speeds in
> `P union {k1,...,k5}` have distance at least `1/14` from the integers.
> Thus the complete clustered `r=5` stratum is uniformly closed. No Covering
> hypothesis and no upper bound on a killer are used.

PROVED uniformly, with `LRCFiveKillerCarrierWindow.lean` and a referee output whose
verdict line reads `clustered r=5 is uniformly lonely`. It is uncontested, in no
MISTAKES entry, and cited approvingly four times in `LRC14-PROOF-MAP.md`, including
"This correction does not reopen clustered `r=5`, which THM-1214 …" and "THM-1214
closes the entire eight-core/five-killer clustered stratum."

So THM-1203 §12's residue list is **stale**: it was written one session before its
own author closed the stratum it was written for. **My THM-1245 §0 repeated that
stale list this morning** ("the actual remaining passage") and is corrected here.

This is MISTAKE-183's rule firing correctly, one session after I wrote it: I
searched canon for the *statement* before working, and the search caught it.

## I. `min |S(P)|`, exactly

The constant `0.164` has been carried since S128c73 through THM-1162, THM-1172,
THM-1173, THM-1174 and HYP-7575 as a **measured** figure, and no frozen output in
the repo contains it. Over all 495 eight-speed cores:

> **`min |S(P)| = 31807/194040 = 0.16391980…`, attained uniquely at
> `P = (1,2,3,5,7,8,9,11)`, which has `N(P) = 18` components.**

Also exact: `max |S(P)| = 61/168 = 0.363095`; `N(P) ∈ [14,26]` (confirming the
range quoted in canon); the shortest longest-component is `1/70`, at
`P = (1,2,6,7,8,9,10,11)` — which is THM-1123's `ell(P) >= 1/70` with its unique
equality core, independently recovered.

## II. The erosion, and where it fails

For a core `P` and killer `k1`, a complete `k1`-gap has width `w = 6/(7k1)` (the
run between consecutive teeth, each of width `1/(7k1)`). The **eroded start
complex** is the morphological erosion of `S(P)` by a window of width `w`,

> `E_w(P) = { t : [t, t+w] subset S(P) }`,   `|E_w(P)| = sum_j max(0, L_j - w)`

over the component lengths `L_j`. Solving `|E_w(P)| = 2/21` exactly per core:

> **`|E_w(P)| > 2/21` holds for every one of the 495 cores iff `k1 >= 204`.**
> The critical width is `w* = 89/21168` at the bottleneck core
> `P* = (1,2,3,5,7,8,9,11)`, giving `K* = 6/(7w*) = 203.865…`.

Below that it is **false**: exactly **19 of 495 cores** have a non-empty window
`13*max(P) < k1 <= K*(P)`, totalling **456 (core, k1) pairs**; the widest is
`P* = (1,2,3,5,7,8,9,11)` with `k1 ∈ (143, 203]`. So the erosion bound as stated
does not cover the clustered range, whose floor is only `13*max(P)`.

The crude bound is worse than useless there: `|S(P)| - N(P)*w` is **negative** in
the worst case (`0.164 - 26*0.00816 < 0`). The exact erosion gains about 9% over it
(`K* = 203.87` against `K_crude = 224.64`) — real, but nowhere near enough.

A second agent computed this independently and landed on the same threshold
(`k1 >= 204`), the same bottleneck core, and the same 19 cores.

**Definitional note.** The one place the repo defines this object (HYP-7575
correction, INDEX line 105) erodes by `1/k`, not `6/(7k1)`. Since `1/k > 6/(7k)`,
that version is strictly smaller, so its failure follows a fortiori from mine.

## III. The residue is not an obstruction — my criterion was just lossy

The sufficient criterion I used to decide a pair `(P,k1)` was
`G(P,k1) > floor(k1 * 2/21)`, where `G` counts complete `k1`-gaps inside `S(P)`.
**294 of the 456 pairs fail it.** But that criterion compares an exact count
against the *worst conceivable* bad set of measure `2/21`, so a failure is
**undecided**, not a counterexample.

Deciding them directly: the real question is whether
`S(P) \ (D_k1 u D_k2 u D_k3 u D_k4)` is non-empty, since any survivor `t` has
`||v t|| >= 1/14` for every speed. Over the nine tightest pairs and every killer
triple with `k4 <= 235` (the r=5 horn bound) — about **160,000 triples** —

> **the safe set is never empty. Zero failures.**

So on the residue the geometry is fine and only the argument was lossy. Had
supplier (3) still been needed, this would have localised the remaining work to
"find a sharper argument on 456 pairs", not "hunt for a counterexample".

## What this is worth

Two constants frozen (`31807/194040`, `K* = 204`), one named obligation retired,
one stale residue list corrected, and the first actual computation of an object
that had been cited as an open obligation in three theorems without anyone ever
evaluating it. Nothing here advances LRC(14).

## Named next

- **Do not** work THM-1203 §12's supplier list without first checking it against
  THM-1214 and the PROOF-MAP; the same staleness may affect supplier (2).
- The live wall is unchanged and is *not* in this stratum: the **n=12 AP-uniqueness
  inverse theorem** (HYP-7310 / CRUX (C)), with boxeph's `3/38` as the
  depth-minimal target.
- `min |S(P)| = 31807/194040` should replace the measured `0.164` wherever it is
  quoted — THM-1162, THM-1172, THM-1173, THM-1174, HYP-7575.
