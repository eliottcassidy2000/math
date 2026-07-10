# LEM-024 — The ≤22 covering census collapses to a 6-witness pigeonhole (removes the census native_decide)

**Status:** PROVED (elementary pigeonhole, below) + machine-verified (exhaustive: all 14 002 covering
min=1 ≤22 tuples are far at one of the 6 witnesses) + **FORMALIZED KERNEL-PURE** (`LRCWindow22Census.lean`,
axioms `[propext, Classical.choice, Quot.sound]` — NO native_decide, NO sorry):
`window22_min1_lonely` (the min=1 pigeonhole), `window22_lonely` (the FULL ≤22 covering census = min=1 ∪
`spread13` for min≥2), and `hdistinct22_kernel` (a signature-exact drop-in for
`WindowData.hdistinct22_from_data`).

**Attribution.** The SIX witnesses `{12/25, 9/26, 7/27, 11/28, 4/23, 11/26}` and the 14 002-tuple
census-branch domain were found FIRST by **kps-S127** (greedy cover, 0 violations). opus-S202 independently
re-derived them, then supplied what was missing: the **pigeonhole PROOF** (why six suffice — the
14-forced-distinct-elements argument below) and the **kernel-pure Lean formalization**
(`LRCWindow22Census.lean`). **kps-S127 also applied the assembly swap** (commit `3453ef3e6`), concurrently.

**✅ LANDED.** The 1-line assembly swap `hwindow22_closed cite` →
`hwindowW_closed 22 cite hdistinct22_kernel` (`LRC14GrandAssembly.lean`, window branch) is DONE, the full
root build is green (8795 jobs), and `#print axioms` now certifies

    'lrc14_from_B5'         depends on axioms: [propext, Classical.choice, Quot.sound]
    'lrc14_grand_assembly'  depends on axioms: [propext, Classical.choice, Quot.sound]

— i.e. **the two `winData22` native_decide axioms are GONE from the LRC(14) top theorem**, which is now
FOUNDATIONAL-AXIOMS-ONLY (no `sorryAx`, no `Lean.ofReduceBool`), with only the LRC(≤13) citation and the
single analytic obligation `hB5` as hypotheses. Cf. MISTAKE-135, opus-S200/S201 (kernel `decide` on the raw
census is infeasible; THM-665 is orthogonal — this pigeonhole was the honest route, "option (c)").
**Source:** opus-2026-07-09-S202.
**Depends on:** the KernelGate bridge `lonely_of_kernelWitness` (far ⟹ Lonely), `LRC14.CoveringFamily`,
`spread13_lonely` (the min≥2 companion).

## Setting

The grand assembly's window branch handles covering families with `Vmax ≤ 22`. After the `spread13`
branch peels off `ratio ≤ 13`, the census's **necessary domain** is the `ratio > 13` covering ≤22 tuples,
which forces `min = 1` (since `Vmax ≤ 22 < 26 = 13·2`). Call these the **census-branch tuples**: 13
distinct integers in `[1,22]`, covering (every `q ∈ [2,14]` divides some element), containing `1`.

## Statement

> **Every census-branch tuple is lonely at one of the SIX rational witnesses**
> `{12/25, 9/26, 7/27, 11/28, 4/23, 11/26}`.

(For `min ≥ 2`, `ratio ≤ 13` and `spread13_lonely` applies — no witness list needed. Together they cover
all covering ≤22 tuples, replacing the 31 471-row `winData22` census.)

## The witnesses and their danger sets (over `[1,22]`)

A speed `s` is **far** at `p/q` iff `‖s·p/q‖ ≥ 1/14`, i.e. `q ≤ 14·min((sp) mod q, q − (sp) mod q)`
(`= KernelGate.speedOK s p q`). The **danger set** `D(p/q)` is the `s ∈ [1,22]` that are NOT far:

| witness | `D` (danger) |
|---|---|
| 12/25 | {2} |
| 9/26 | {3} |
| 7/27 | {4} |
| 11/28 | {5} |
| 4/23 | {6, 17} |
| 11/26 | {7, 19} |

Each far-set `F = [1,22] \ D` is huge (≥ 20 of 22 speeds). A tuple `S` is far at `p/q` iff `S ∩ D = ∅`.

## Proof (pigeonhole — 14 forced distinct elements > 13 slots)

Suppose a census-branch tuple `S` is far at **none** of the six. Then `S` meets every danger set:

* from 12/25, 9/26, 7/27, 11/28: `{2,3,4,5} ⊆ S`;
* from 4/23: `S ∩ {6,17} ≠ ∅` — pick `a ∈ S ∩ {6,17}`;
* from 11/26: `S ∩ {7,19} ≠ ∅` — pick `b ∈ S ∩ {7,19}`.

Because `S` is **covering**, the *unique* multiples of 12, 13, 14 in `[1,22]` force `{12,13,14} ⊆ S`, and
each of `q = 8,9,10,11` (multiples `{8,16},{9,18},{10,20},{11,22}` — exactly two each) forces one
representative `c ∈ S∩{8,16}`, `d ∈ S∩{9,18}`, `e ∈ S∩{10,20}`, `f ∈ S∩{11,22}`. Because it is
**min = 1**, `1 ∈ S`.

The eight literals `{1,2,3,4,5,12,13,14}` and the six pair-slots `{6,17},{7,19},{8,16},{9,18},{10,20},
{11,22}` are **pairwise disjoint** (check: no danger element `2,3,4,5,6,17,7,19` lies in a covering pair
or equals `12,13,14`; `17` and `19` are primes `> 14` serving danger only). Hence `a,b,c,d,e,f` are
distinct from the eight literals and from each other, giving

> `|S| ≥ 8 + 6 = 14`,

contradicting `|S| = 13`. Therefore `S` is far at one of the six. Applying `lonely_of_kernelWitness` at
that witness yields a lonely instant. ∎

**`min = 1` is essential.** Without `1 ∈ S` the count is exactly `13`, and indeed
`{2,3,4,5,6,7,12,13,14,16,18,20,22}` is a covering min=2 tuple far at none of the six (ratio 11 — caught by
`spread13`, not the witnesses). This is why the split (min=1 → 6 witnesses; min≥2 → `spread13`) is exactly
right.

## Why this removes the native_decide

The current census (`winData22_ok`, `winData22_complete`) is a `native_decide` sweep over the whole
`C(22,13) = 497 420` window universe (opus-S200/S201: kernel `decide` on it is infeasible — >13 h + OOM —
MISTAKE-135). LEM-024 replaces that global enumeration with a **fixed 6-element witness list** plus a
**pigeonhole**: the per-speed far facts (`speedOK s p q` for `s ∈ [1,22]`, 22 × 6 tiny kernel `decide`s)
and the covering-forces-multiples facts are all small and kernel-pure; no enumeration of tuples occurs.
This is the honest route to a foundational-axioms-only LRC(14) that opus-S201 identified as "option (c)".

## Ledger

- Census-branch domain = 14 002 tuples; a greedy witness cover found **6** witnesses (`12/25, 9/26, 7/27,
  11/28, 4/23, 11/26`) hitting all of them; the danger sets are `{2},{3},{4},{5},{6,17},{7,19}`.
- Pigeonhole: fail-all-6 ⟹ 14 forced distinct elements in a 13-set — impossible. `min=1` essential.
- Removes the `winData22` native_decide census (the last non-foundational axiom besides the LRC(≤13)
  citation, on the with-census route). Lean: `LRCWindow22Census.lean` (in progress). → MISTAKE-135,
  opus-S200/S201, THM-665 (orthogonal), `spread13_lonely`, `lonely_of_kernelWitness`, hwindow22.
