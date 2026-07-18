# The stability step is NOT "near-tight ⟹ AP" — it is cover-gap uniqueness (death-star-2026-07-18-S57)

**Context.** THM-1039 §4 reduced the compact far-element inverse theorem to "LRC(13) stability:
near-tight core ⟹ near-AP" (HYP-7310). Attempting that step this session shows the framing was
**imprecise**, and the correct target is cleaner. Honest record. Scripts: `lrc_tight_uniqueness`,
`lrc_neartight_hamming` (S57).

---

## 1. "Near-tight ⟹ dilated AP" (exact) is FALSE

`{1,2,…,11,24}` covers 2..12, misses 13,14, has `M = 2/25 = 0.08 ≤ 1/13 + 34/2366 = 0.0913`
(**near-tight**), and is **not** a dilated AP. So there is no exact implication "near-tight core ⟹ dilated
AP." Near-tight non-AP cores exist; the correct sense of "near-AP" is **Hamming-close to a dilated AP**:
`{1..11,24}` differs from `{1..12}` in one position (distance 1). Verified structure:
- **`max ≤ 18`:** the *only* near-tight covering-2..12 core is `{1..12}` itself (Hamming distance 0).
- Non-AP near-tight cores need larger elements; Hamming distance to the nearest dilated AP grows with `max`
  (`{1..11,24}` d=1, `{1..10,22,24}` d=2, …). "Near-tight ⟹ bounded Hamming distance to a dilated AP" is
  the true statement — and it is exactly **HYP-7310** (klein's n=12 AP-uniqueness / Tao's optimistic
  conjecture), an OPEN additive-combinatorics inverse theorem.

## 2. What IS exactly true: tight ⟹ dilated AP

Among primitive 12-subsets of `{1..18}`, the **unique** tight (`M = 1/13`) family is `{1..12}` (0 non-AP
tight). Dilated APs `d·{1..12}` are the tight families in general (non-primitive for `d≥2`). So the
**equality case is clean**: `M(W)=1/13`, `W` primitive ⟹ `W = {1,…,12}`. (Extremal-uniqueness of LRC(13);
the strict interior `M ∈ (1/13, 1/13+ε)` is where non-AP appears.)

## 3. The correct reduction: the far-element inverse theorem is COVER-GAP UNIQUENESS, not HYP-7310

The operative condition for `M(V)<1/13` (covering 2..14, AP-core-adjacent) is **not** "core near-tight"
(necessary but far from sufficient — `{1..11,24}` is near-tight yet `M(V)={1..11,24}∪{182}=2/25≥1/13`). It
is the exact criterion (THM-1039 Prop A):

> `M(V) < 1/13 ⟺ coverGap(W,182) < 1/13` ⟺ `G_W` fits inside the `182`-danger arcs.

And **cover-gap uniqueness** is the target: `coverGap(W,182) < 1/13 ⟹ W = dilated AP`. THM-1039 proves it
for
- **non-fragmented** cores (`C ≤ 464μ`): soft Weyl — `coverGap ≥ avg ≥ 1/13`;
- **not-too-near-tight** cores (`δ > max/2366`): stability — component half-width `> 1/2366` (far-arc).

Non-AP cores have `coverGap = 1/2` (the far element `182` is too fine to cover positive measure). So the
inverse theorem does **not** route through the full HYP-7310; it is cover-gap uniqueness, **mostly proved**.

## 4. Where HYP-7310 actually enters — the very-near-tight residual

The one place both lenses fail is **fragmented AND very-near-tight** (`C > 464μ` AND `δ ≤ max/2366`). As
`δ → 0` these cores converge to the tight locus, i.e. to `{1..12}` (by §2). Proving `coverGap ≥ 1/13` on
this thin residual *is* the near-tight limit — the HYP-7310 kernel — but restricted to the cores that are
(a) very near tight and (b) fragmented. For `max ≤ 34` the residual is finite and checked directly (the
`coverGap = 1/2` computation on every sampled core); for `max ≥ 35` it is the renormalization stratum
(HYP-3901). So:

> **The far-element inverse theorem = cover-gap uniqueness = [soft Weyl] ∪ [stability] ∪ [very-near-tight
> fragmented residual].** The first two are PROVED and uniform in max; the residual is the genuine kernel,
> a *restricted* HYP-7310 (very-near-tight + fragmented), not the full inverse theorem.

## 5. Honest status

- **Cannot** prove HYP-7310 (near-tight ⟹ near-AP) — it is the open additive-combinatorics inverse theorem.
- **Proved / clean:** tight ⟹ dilated AP (equality case, verified); the two cover-gap lenses (THM-1039);
  the reduction of the far-element inverse theorem to cover-gap uniqueness, which is *not* the full
  HYP-7310 but a much thinner residual (very-near-tight + fragmented).
- **Correction to THM-1039 §4 / prior turn:** the compact first step does not reduce to "near-tight ⟹ AP"
  (false as an exact implication); it reduces to cover-gap uniqueness, whose only open part is the
  very-near-tight fragmented residual. The full stability HYP-7310 is *sufficient* but *not necessary*.

→ THM-1039, THM-1037, THM-1038, THM-1017, HYP-7310, HYP-3901, MISTAKE-161.
