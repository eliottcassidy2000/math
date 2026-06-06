---
source: monad-explorer-2026-06-06-S703 (deep-research, signed-LRC small-n exhaustive)
status: THEOREM (THM-412) + strongly-verified conjecture (HYP-2267) + a correction of S699's T4.
The signed-LRC sign-orbit of AP_n loses information ONLY through the mod-(2n−1) fold (never through
an automorphism — the unfolded cut-spectrum is universally faithful), and the cleanest cause of that
loss is the order-3 torsion runner x=(2n−1)/3, which exists iff 3∣(2n−1) and is a silent single flip.
tags: [signed-lrc, sign-orbit, fold-collision, order-3, prime-3, eulerian, half-system, faithful,
       2n-1, n14, C27, correction-of-T4, THM-412, HYP-2267]
---

# Signed LRC: the sign-orbit degenerates only by folding, and the door is the order-3 point

**Dispatched angle:** exhaustive small-n signed LRC — enumerate sign-reversal patterns, find the
optimal gap at each, tabulate which patterns "beat 3n."

**First, an honesty correction to the seed.** By S699's T1 (gauge invariance), `‖ε_i v_i t‖=‖v_i t‖`,
so the observer gap `M` is **identical for every sign pattern**. No sign reversal ever changes
whether a config is tight, and "beats `3n`" is the *unit-distance* cap (the S702/THM-411 lattice
thread, now owned by a peer), not an LRC quantity. So the honest small-n signed-LRC question is not
about `M` at all — it is about the **finer pairwise invariant** the sign group exposes (T2–T4): the
multiset of folded pair-clocks. This reflection is the rigorous small-n anatomy of that invariant.

## Two layers, cleanly separated

For `AP_n={1,…,n−1}`, `C=2n−1`, a cut (2-coloring of runners, up to swap; `2^{n−2}` of them) gives a
clock per pair: `|i−j|` if monochromatic, the **sum** `i+j` if it is a cut edge. The LRC-relevant
version reduces mod `C` and folds: `ρ(s)=min(s,C−s)` (the shell/resonance strength, THM-401). Two
invariants:

- **Unfolded** cut-spectrum (sums and diffs, no mod). **Verified universally faithful**: across *all*
  `2^{n−2}` cuts, every config tested gives `2^{n−2}` distinct multisets — exhaustively for all
  speed-sets at `n≤7` (3003 configs at `n=7`), and for adversarial symmetric / Sidon-violating /
  geometric sets. The cut and its unfolded clock-spectrum determine each other.
- **Folded** cut-spectrum (mod `C`). This is the genuine LRC object (it is the pair-clock resonance
  at the base division time `t=1/C`). Here the sign-orbit can be `< 2^{n−2}`.

> **The correction.** S699's T4 said the few sign-orbit "collisions are config automorphisms." They
> are **not**. Since the unfolded spectrum is always faithful, *no automorphism exists*: every
> collision is created purely by the mod-`C` **fold** `s↦min(s,C−s)` (the shell-partner involution
> `s↔C−s`). The signed structure loses information only through arithmetic folding, never through a
> symmetry of the configuration.

## When does folding actually merge two cuts? — `2n−1` composite

Exhaustive `n=3..22` (`C=5..43`):

```
   sign-orbit = 2^{n−2}  ⟺  C = 2n−1 is PRIME.
   prime C (5,7,11,13,17,19,23,29,31,37,41,43): orbit FULL, 0 collisions.
   composite C (9,15,21,25,27,33,35,39): orbit < 2^{n−2}.
```

Collision counts: `C=9→1, 15→4, 21→8, 25→4, 27→69, 33→32, 35→16, 39→64`. The explosion at
`C=27=3³` is conspicuous — the 3-richest modulus is the most degenerate.

## The mechanism: a silent flip through the order-3 point (THM-412)

A *single-flip* collision is "flipping one runner `x`'s color changes nothing." Flipping `x` toggles
every incident edge between its mono value `|x−y|` and its cut value `ρ(x+y)`. Build the
**value-multigraph** `G_x` with an edge `{|x−y|, ρ(x+y)}` per neighbour `y`. Then (proved):

- **Lemma A:** a silent flip of `x` exists **iff `G_x` is Eulerian** (all degrees even) — because a
  silent flip is exactly a balanced edge 2-coloring of `G_x`, which exists iff all degrees are even.
- **Lemma B:** using that the runners `{1,…,(C−1)/2}` are a **half-system** of `ℤ/C`, the
  degree of value `v` in `G_x` is `2` for every `v`, *except* two parity defects at `v=x` and
  `v=ρ(2x)` (the half-system class `±0` degeneracy, and the exclusion of the runner `x` itself). So
  the odd-degree vertices are **exactly `{x, ρ(2x)}`**.
- **Theorem:** the two defects merge into a single even vertex — making `G_x` Eulerian — iff
  `x=ρ(2x)`, i.e. `3x≡0 (mod C)`, i.e. **`x=(2n−1)/3`**. This runner exists iff `3∣(2n−1)`.

> **So the cleanest door through which the signed invariant collapses is the order-3 torsion point
> of `ℤ/(2n−1)`.** Whenever `3∣(2n−1)`, the runner `x=(2n−1)/3` is a silent single flip — the signed
> structure cannot see it. This is the signed-LRC face of the **prime-3** thread (THM-407): `n=14`,
> `C=27=3³`, silent runner `x=9`; being `3³` it is maximally folded (all single-flip collisions flip
> `9`, and the multi-flip collisions pile on to 69 groups). The same `π/3` / prime-3 / Eisenstein
> rosette that governs the worry-set shells and the unit-distance optimum (HYP-2170, HYP-2264/THM-411,
> S599u) reappears here as *the* degeneracy locus of the sign group.

Composite `C` with **no** factor 3 (e.g. `25, 35, 49`) still degenerates, but with **no** single
silent flip (`G_x` is Eulerian for no `x`); those collisions are genuinely multi-flip. So `3∣C` is a
clean *sufficient* cause; the full "iff `2n−1` prime" is HYP-2267 (other prime factors contribute
multi-flip collisions whose count law is open).

## Why this matters for the LRC attack

The whole point of S699's signed program is that the signed (additive) face is *finer than `M`* and
splits the worry-set at `n=14` via a shell-partner. This reflection pins **where that finer
invariant itself is blind**: precisely at the order-3 point, precisely when `3∣(2n−1)`. The two facts
sit on the same axis — `n=14` is simultaneously (a) where a tight config first carries a shell-partner
(`V*`, `3+24=27`, S699), and (b) where the sign-orbit is most degenerate (`C=27=3³`, this work). Both
are prime-3 phenomena on the modulus `C=2n−1`. The signed invariant's **blind spot** (order-3 silent
flip) and its **discriminating power** (the `V*` shell-partner zero-clock) are two readings of the
same `3 ∣ C` structure. A proof attack at `n=14` that uses the signed face should treat the order-3
runner `9` as a coordinate the sign group gauges away, and carry the discriminating data on the
shell-partner pair `(3,24)` instead.

## Honest status
- **Proved:** Lemmas A, B and the order-3 theorem (THM-412); hence `3∣(2n−1) ⟹` sign-orbit `<2^{n−2}`.
- **Verified:** Lemma B (n≤40); `x=C/3` Eulerian (n≤199); unfolded faithfulness (exhaustive n≤7 +
  adversarial); the iff-prime law (n≤22).
- **Conjecture (HYP-2267):** sign-orbit `=2^{n−2}` iff `2n−1` prime.
- **Not touched:** the observer gap `M` (sign-invariant); off-grid LRC(14).

**Artifacts:** `04-computation/signed_lrc_orbit_collisions_s703.py`, `…faithfulness_s703b.py`,
`…fast_confirm_s703c.py`, `…collision_mechanism_s703d.py` (+`.out`s),
`05-knowledge/results/signed_lrc_thm411_proof_check_s703e.out`. THM-412, HYP-2267, T757.
Extends S699 (T1–T4, worry-set split), THM-401/403 (shells/2n−1), THM-407 (prime-3, n=14).
