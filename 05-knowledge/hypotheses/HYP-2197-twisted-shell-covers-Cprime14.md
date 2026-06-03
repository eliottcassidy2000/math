---
id: HYP-2197
title: Twisted-shell dodge (m ≤ 2n-1) ∪ B′ covers C′(14) — a constructive route to LRC(14)
status: OPEN (strong evidence: 0 residual / 46000 configs)
source: claudebox-2026-06-03-S622
depends_on:
  - THM-412   # the twisted-shell dodge (the constructive criterion)
  - THM-398   # the reduction to C′ and Criterion B′
related:
  - THM-401   # 2n-1 = 27 = 3³ (the shell ceiling)
  - THM-411   # tight sets witnessed on (ℤ/m)* (the same multiplier object)
  - HYP-2140  # the 2-adic seam / apex (the "cross")
  - HYP-2055  # n=14 tight family (witnesses {1,3,5,9,11,13}/14 = (ℤ/14)*)
---

# HYP-2197 — the twisted-shell dodge closes the C′(14) residual

## The conjecture (proof route to LRC(14))

> **Coverage.** Every multiple-of-14 speed set `S` (the C′(14) class) is loose because it satisfies
> the twisted-shell dodge [[THM-412]] for some `m ≤ 2n-1 = 27`, OR Criterion B′ (dominance) of
> [[THM-398]]. Equivalently `twisted-shell(m ≤ 27) ∪ B′ = C′(14)`.

With THM-398's reduction this is a **constructive proof of C′(14), hence of LRC(14)**: it exhibits,
for every config, an explicit rational witness `t = a/m` (shell `m`, perspective `a`) at which all
runners exceed `1/14`.

**Evidence.** `0` residual over 6000 random near-tight + 40000 adversarial multiple-of-14 configs
(`lrc_n14_residual_empty_s622.out`, `lrc_n14_ceiling27_s622.out`). The twisted-shell dodge alone
(m ≤ 27) covers all near-tight configs and **subsumes B′ there** (the `846/6000` B′-failures are all
twisted-shell-covered); the only configs needing B′ are far-from-tight dominant ones (minimal
twisted `m ≳ 29`, but `M ≈ 3/11` — trivially loose, long safe component).

## The structure at the critical shell `m = 2n-1 = 27 = 3³`

The dodge on the deepest shell `m = 27` has a clean group-theoretic form. The danger band excludes
residues `{0, 1, 26}` (band-distance ≤ 1). A multiplier `a` works iff `a·S mod 27` avoids `{0,±1}`,
i.e. iff the **unit ±-pair** `{u, -u}` with `u = a^{-1}` is disjoint from `S mod 27` (and `0 ∉ S
mod 27`). So:

> **Twisted-shell-27 dodge ⟺ some unit ±-pair `{u,-u} ⊂ (ℤ/27)*` is free of `S`'s residues.**

`(ℤ/27)*` has `φ(27)=18` units in **9 ±-pairs** — the orbits of the **twisted involution** `u ↦ -u`.
The 13 runners must block all 9 pairs to defeat the shell-27 dodge. The non-unit residues
(multiples of 3 — the **inner 3-adic shells** `3·(ℤ/9)`, and the apex `0`) are the degenerate
"cross": there the `±` involution acquires fixed points (`r ≡ -r mod 27 ⟺ 27 | 2r ⟺ r ∈ {0}`,
and within the 3-shells the structure collapses). A runner landing on the inner shells is "stuck on
the cross" and cannot be rotated off by any multiplier — it must be dodged on a coprime shell instead.

### The laminar / flow picture (heuristic — NOT a clean single-shell lemma)
A contiguous arc of ≤ 13 residues mod 27 is a **laminar** (non-crossing) set; it leaves a free
unit ±-pair for **all but exactly two** arcs — the two half-lines `{1,…,13}` and `{14,…,26}`
(verified: `lrc_n14_*_s622`). So the shell-27 unit-pair dodge covers near-tight configs *except*
when their residues mod 27 are (a multiplier-image of) a half-line; those escape on a **different**
shell `m ∈ {15,…,26}` instead. Hence coverage is the **union over shells**, not a single `m=27`
statement. **Turbulence** (residues spread across the cycle, hitting all 9 pairs) = the dominant/
spread configs, loose by B′. The qualitative dichotomy *laminar ⟹ some-shell dodge / turbulent ⟹
dominance* is the conjecture; the half-line exception is why no single shell suffices.

## What to prove
1. **Per-shell free-pair lemma (corrected):** a config whose residues mod 27 form a contiguous arc
   of length ≤ 13 leaves a free unit ±-pair **unless** the arc is `{1,…,13}` or `{14,…,26}`; the
   exceptions escape on the next coprime shell. (Verified count: 2 bad arcs of 27.)
2. **The dichotomy / residual-empty:** every multiple-of-14 config satisfies the twisted dodge for
   *some* `m ∈ {2,…,27}` or Criterion B′. This is the `0/46000` empirical claim — prove it.
3. **General `n`:** twisted-shell dodge with `m ≤ 2n-1` ∪ B′ ⟹ C′(n) ⟹ LRC(n)? Test n=15,16,18.
