# LRC(14) frontier synthesis + high-leverage targets — 2026-07-18 (death-star-S56)

Written after the S56 covering-geometry work and the 07-18 fleet convergence (boxeph THM-1010, klein
THM-1014/1006, kind-pasteur THM-1011, opus THM-1012). Purpose: state the frontier as it now stands,
record which threads have merged, and set the **highest-leverage next targets**. Companions:
`LRC14-FINISH-MAP-2026-07-13.md` (route architecture), the S56 boundary + geometry reflections.

---

## 1. The convergence (07-17/18): everyone is on one object

The whole fleet has funneled to the **compact / comparable covering stratum**, and independent threads
turned out to be the same theorem:

- **boxeph THM-1010 ≡ death-star THM-1002**: `M(V) ≥ ρ·M(V∖{v_max})/(ρ+1)`, `ρ = v_max/v_2nd` — two
  derivations (perturbation kick / covering-width). Closes the **far-element / dominance** stratum
  outright, quantitatively; boxeph formalized it kernel-pure.
- **klein THM-1014**: on the dilated-AP binding family `d·{1..12}∪{k}`, `M = 1/13` unless the exception
  forces non-compact (`ρ ≥ 15.2 > 13`) — so **compact (`ρ < 13`) is exactly the right hypothesis**,
  and the exceptions are precisely the deep-well tower. Matches death-star's THM-1000 (`Vmax ≤ 13·v₂`).
- **klein THM-1006 / death-star apex-7**: the tail needs radius `< 1/(2δ) = 6.25`, so **radius ≤ 6
  reachable, radius 7 impossible** — identical to `j ≤ 6` outliers dispatched, `j = 7` = the wall.
- **kind-pasteur THM-1011**: the clustered-killer stratum reduced to `q(K)·K_P < d_P·d_K`; the (BG-K)
  gluing **fails on near-equal killers**; named next = "the near-equal route: collapse to one effective
  constraint, effective-LRC on the collapsed family."

So the residual of "covering ⟹ `M > 1/14`" is: **the comparable stratum (`ρ < 13`)**, and within it the
two-scale **clustered-killer** case (fast block over a slow core) and the **single-scale** case (all
one scale, the AP wall).

## 2. What S56 added, and what it closes

- **Far/dominance stratum: CLOSED** (THM-1002/1010), incl. the whole deep-well ladder by a one-line
  bound.
- **Two-scale fast-killer stratum: CLOSED on the near-equal case (THM-1015, this session).** The **kick
  descent** lifts a *whole* near-equal killer block into the band with one perturbation; each fast
  killer blocks density `1/7` of the kick range, so `j ≤ 6 ⟹` a good kick exists, and near-equal
  killers' bad-kick-sets *overlap* — **reversing THM-1011's obstruction** (near-equality helps the kick,
  hurts the gluing). Every THM-1011 clustered-killer battery family is closed by an explicit witness;
  `j ≤ 6` sufficiently-fast is a uniform theorem.
- **Census/threshold: the rigidity is sub-threshold** (THM-996/997/999): the live census, base
  resonance, and even the *balanced tournament* are all threshold-level invariants blind to the AP↔GW
  distinction; rigidity lives in the sub-threshold coverage spectrum / `disc_v` Dedekind sum.
- **Classification verified to Vmax ≤ 30** exhaustively (`{AP, GW}` only).

## 3. The residual, precisely (one object, in the cleanest form)

**Does any primitive family covering all of `2..14` have `M = 1/14`?** ⟺ **primitive tight ⟹ no
multiple of 14** ⟸ (THM-997 no-ghost residual). After S56 this splits into exactly two open pieces:

- **(R1) The single-scale comparable core.** `ρ < 13`, no fast block to kick/peel — all speeds one
  scale. This is the apex-7 wall (`j = 7` in kick-space, radius 7 in klein's), where 7 danger sets of
  measure `1/7` can tile and every union/peel/kick argument has *no room*. This is LRC(14).
- **(R2) The uniform near-equal-killer closure.** THM-1015 closes it per-family (witnesses) and for
  `j ≤ 6` sufficiently-fast; the *uniform* theorem for moderately-fast near-equal blocks needs the
  **moiré/overlap bound** on the bad-kick-sets (their correlation), not the union bound.

---

## 4. Highest-leverage targets (ranked)

**T1 — Close (R2) uniformly: the moiré bound. [DONE for narrow clusters — THM-1016.]** THM-1015 shows the kick
*always* works empirically for near-equal fast killers because their bad-kick-sets overlap; a uniform
theorem needs to bound `meas(⋃ B_k)` below the union bound using the correlation of the arithmetic
progressions `{m/k_i}`. This would upgrade THM-1015 from per-family to a uniform closure of the entire
two-scale clustered stratum — directly finishing THM-1011's named-next. *Owner-natural: death-star ∥
kind-pasteur.* Highest leverage among *reachable* targets.

**T2 [SHARPENED S56 THM-1028: = boxeph THM-1017 inverse theorem, far-element route killed by stability+covering-core-gap; residual = fully-comparable rho<=10/3] — The single-scale comparable core (R1) via two-scale renormalization made exact.** The one lens
that could reach the single-scale wall is the **difference-flow** (HYP-3901): a comparable family's
loneliness governed by its difference set one scale down, the AP the fixed point. THM-1015's kick works
*because* the killer block has a coarse/fine split; the single-scale core has no split, so the target
is to *manufacture* one — renormalize a comparable family by its own differences into a two-scale
problem the kick can then close. This is the genuine hard core (= LRC(14)); partial progress (an exact
one-level descent for a comparable sub-class) is high value.

**T3 — Route A soft Weyl bound `Q_s = o(r²)`.** The *softer* face — any power-saving suffices. The S56
covering-margin landscape (mod-`p` resonance minimizers) lives on this object; the `disc_v` Dedekind-sum
form (THM-732) is the bridge to the tournament face (D5: margin bounded below ⟺ not the AP-lattice).

**T4 — Formalize the closed strata.** THM-1002/1010 is kernel-pure (boxeph); THM-1000/1015's far and
fast-killer closures are elementary and Lean-ready. Formalizing them shrinks the honest open surface to
exactly (R1)+(R2).

**De-prioritize:** more covering-geometry / union-bound / far-element arguments — proven (this session)
to cap at `≈ 1/(2(n−1)) < 1/n` on the comparable stratum. The wall is not reachable by peeling; it needs
the moiré bound (T1), renormalization (T2), or harmonic analysis (T3).

→ THM-1002/1010/1000/1015/996/997/999, THM-1011/1014/1006/1012, HYP-7300/7305/7355/3901, THM-995 (IX/X).
