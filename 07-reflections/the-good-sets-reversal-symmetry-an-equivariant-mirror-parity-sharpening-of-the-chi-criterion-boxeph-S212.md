# The good set's reversal symmetry: an equivariant (mirror-parity) sharpening of the χ(G_δ) criterion

*boxeph-2026-07-21-S212. Owner: look for other topological advances the repo has made, and come up with
creative LRC arguments combining and extending them. Builds on codex THM-2047 (χ(G_δ) criterion), boxeph
S210 (the reversal involution / antisymmetry), kind-pasteur-S19 (LRC Lefschetz trace = free ι + Gauss sum),
THM-1820 (mirror pairs), HYP-3015 (superlevel persistence barcode). New pillar P6 verified in
`04-computation/lrc_nerve_persistence_topology_boxeph_S212.py`.*

## The repo's topological toolkit for LRC (surveyed and credited)

A deep pull shows the fleet has already built most of the topological machinery I would have proposed —
so I credit it rather than reinvent it, and re-verified the core (P1–P5 of the script):

- **THM-2047 §4 (codex, PROVED):** `LRC(14) for S ⟺ χ(G_{1/14}(S)) > 0`, where `G_δ={t:f_S(t)≥δ}`,
  `f_S=min_v‖vt‖`, and `χ(G_δ)=#components` (closed arcs *and* isolated points). [P1 re-verifies this,
  including that the tight `(1,2,3)` good set at `δ=¼` is the two isolated points `{¼,¾}`, `χ=2`.]
- **HYP-3015 (codex-S179):** `{G_δ}` is the **superlevel filtration** of `f_S`; each lonely window is a
  persistence bar; `M(S)` is the top death value. [P3 re-verifies: `M(S)`=top `H₀` bar, born on a pair-sum
  wall `q∣v_i+v_j` (THM-2047 §2); `χ(G_δ)=#bars alive`.]
- **opus's Euler-char certificate & kps's Alexander duality** (`lrc_lonely_set_euler_char_certificate`,
  `cohomological_three_distance_and_sqrt21_bridge`): `χ(L)=#lonely arcs`, `b₀(lonely)=b₀(danger cover)`,
  lonely and covered arcs alternate. [P2 re-verifies the nerve inclusion–exclusion
  `χ(U)=Σ(−1)^{|I|+1}[∩A_i≠∅]` — the topological analog of THM-1820's *measure* I-E (lengths → sinc) — and
  that at the tight threshold measure covers (`|U|=1`) while `χ(G)=2`: the nerve repairs the S211
  volume-blindness.]
- **HYP-3025/3101 (codex):** the closed-arc **Čech nerve** (`β₁=1` ⟺ full cover) with a Betti-defect
  sidecar; the **normal-fan Čech barcode** component-bound proof route (still open).
- **kind-pasteur-S19:** the LRC singular series as a **Lefschetz trace** for the involution `ι:t↦1−t`;
  `ι` acts **freely** on units, so ordinary Lefschetz `Λ(ι)=0` is *blind*, and the certification is the
  **odd-equivariant index = quadratic Gauss sum `i√7`** (Borsuk–Ulam odd degree).
- **Tournament side (PROVED): THM-587** — the metagraph reversal `R:T↦Tᵒᵖ` has
  `SC = tr(R) = ` antipodal Euler/Lefschetz number, giving the closed-form `b₁⁻`.

## The gap — nobody combined the reversal symmetry with χ

Two of these live on the *same* involution but never met: kps-S19 has `ι:t↦1−t` acting freely (with the
Gauss sum), and my S210 has the same `ι=θ↦−θ` as the hinge of antisymmetry (with its `2`-torsion fixed
points `{0,½}`) — but codex's `χ(G_δ)` criterion is stated *without* the symmetry. Combining them gives a
free `ℤ/2`-action on the good set, and a free action forces **even Euler characteristic**. That is the new
pillar.

## The equivariant (mirror-parity) sharpening (P6, verified)

Since `‖·‖` is even, `f_S(1−t)=f_S(t)`, so **`G_δ(S)` is `ι`-invariant**. The only fixed points of `ι` on
`T` are `0` and `½`. Now:
- `f_S(0)=0<δ` always;
- `f_S(½)=0` **iff some speed is even** (`‖v/2‖=0` for even `v`), and `=½` iff all speeds are odd.

Every LRC(14) **covering** set contains a multiple of 2 (else `t=½` is lonely), so **both `ι`-fixed points
are dangerous ⟹ `ι` acts *freely* on `G_δ` ⟹ `χ(G_δ)` is EVEN.** Codex's criterion therefore sharpens:

> **For a covering set `S`: `LRC(14)` ⟺ `χ(G_{1/14}(S)) ≥ 2` ⟺ at least one MIRROR PAIR `{t*, 1−t*}` of
> lonely windows survives.**

Verified: the deep well `{1,…,12,182}` at `δ=1/14` has **`χ=24` = twelve mirror pairs** (even); the tight
`(1,2,3)` good set at `δ=¼` is the mirror pair `{¼,¾}` (`χ=2`, `¾=1−¼`); a covering `(1,2,3)` at `δ=0.3`
has `χ=0`. The **all-odd exception** is exactly the `ι`-fixed case: `½∈G_δ` is an `ι`-invariant component,
`χ` can be **odd** (verified `(1,3,5,7)`: `χ=1`), and the set is automatically lonely at `½` — the
classical "all speeds odd ⟹ lonely at ½", here recovered as the Borsuk–Ulam *fixed-point* case.

This is the topological form of THM-1820's measured "mirror pairs" (B3) and of my S210 involution: lonely
windows of a covering set come in `ι`-conjugate pairs, and the count is even.

## What it buys toward Wall A

1. **Equivariant halving.** LRC(14) for a covering core `⟺` there is one lonely window in the *half*
   fundamental domain `[0,½]` — its mirror in `[½,1]` is automatic. Wall A's search/proof domain halves.
2. **A parity obstruction.** `χ(G_{1/14})` is even, so it can never be `1`: a disproof needs `χ=0`, i.e.
   **every** potential mirror pair killed **simultaneously**. The covering must itself be `ι`-symmetric —
   a rigid constraint a putative counterexample's danger arcs must satisfy on `[0,½]` and `[½,1]` at once.
3. **The Gauss-sum obstruction is the right home.** kps-S19's `Λ(ι)=0` is blind precisely *because* the
   action is free; the live invariant is the **odd-equivariant** index `= i√7 ≠ 0`. In this language it is
   an obstruction class of the `ℤ/2`-quotient `G_δ/ι` on `[0,½]` — a Borsuk–Ulam-type nonvanishing that
   the covering must defeat. The equivariant `χ` (even) and the odd index (`i√7`) are the two halves of the
   `ℤ/2`-equivariant Euler class of the good set.

So the creative extension is: **push codex's `χ(G_δ)>0` through the reversal symmetry that S210/kps-S19
already isolated, turning it into `χ≥2` (mirror-parity) and reducing Wall A to the `ι`-quotient good set on
`[0,½]`.** The same involution `R`/`ι` that THM-587 uses on the tournament metagraph acts here on the LRC
circle — one antisymmetry (S210), now carrying an equivariant Euler class on both sides.

## Honest scope

P1–P5 re-verify the fleet's existing topological toolkit (THM-2047, HYP-3015, opus/kps Euler-char and
Alexander duality) — credited, not claimed. The **new** piece is P6: the free-`ι` / even-`χ` /
mirror-parity sharpening of the criterion, verified on the deep well (`χ=24`), the tight `(1,2,3)`
(`χ=2`), and the all-odd exception. It is a genuine sharpening (`χ>0 → χ≥2`, even) and an equivariant
halving of Wall A, and it names the Gauss sum `i√7` as the odd-equivariant obstruction — but forcing
`χ≥2` for every 13-speed covering core (LRC(14) itself) remains open. The concrete reduction it leaves:
*prove the `ι`-quotient good set `G_{1/14}(C)/ι` on `[0,½]` is nonempty for every covering core `C`.*

Links: HYP-8845, THM-2047 (codex), THM-587, THM-1820, HYP-3015, HYP-3101,
[[antisymmetry-is-the-hinge-tori-odd-functions-saddles-and-tournaments-boxeph-S210]],
[[where-gmc2-reaches-lrc14-the-ct-functional-and-where-it-stops-the-volume-ceiling-boxeph-S211]].
