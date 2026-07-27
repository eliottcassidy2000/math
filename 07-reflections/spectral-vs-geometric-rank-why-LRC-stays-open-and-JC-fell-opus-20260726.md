# Spectral rank vs geometric rank: the axis on which boxeph's "two nullcones" actually diverge

**opus-2026-07-26. Provenance note / synthesis, not a truth source.** Owner: find how our
understanding of LRC and JC matured, speculate wildly about lurking recurrences, and test extensively.
This note proposes one axis that (i) unifies boxeph-S205's "nullcone rigidity" frame with kps-S131b's
"carrier + prime lines" spectral grammar, (ii) *explains* boxeph's undated observation that the three
nullcones "diverge at higher rank," and (iii) makes a testable prediction about the live 2026-07-24/25
frontier of both programs. Tests done are marked ✓; speculations released are marked ⌦.

## The two master frames, and the seam between them

- **boxeph-S205 (structural):** GMC(2), JC(2)⟺Zhao-VC, LRC(n) are one shape — *nullcone rigidity* for
  a functional `F` on powers `P^m` (`F(P^m)=0 ⟹ P` degenerate). They share a **rank-1 seed**
  (THM-1840, proved functional-agnostically for Gaussian `E`, Laplacian `Δ`, and sinc alike) and
  **diverge at higher rank**. The AP is the universal nullcone vertex.
- **kps-S131b (spectral):** every solved question = *smooth carrier × discrete prime lines* (one exact
  local factor per prime); the open cores are all **no-cancellation** — a specific line must be shown
  nonzero. LRC14's crux lives on the joint spectrum of `91 = 7·13` (a mixed unit character mod 91).

**The seam (this note): "rank" means two different things, and the divergence IS that difference.**

> **SPECTRAL rank** = the number of prime lines that must be jointly certified nonvanishing (kps).
> **GEOMETRIC rank** = the ambient dimension / degree that must hold a collision (the JC arc).

boxeph's "rank-1 seed, diverge at higher rank" then reads: the rank-1 case (one line / one branch) is
proved for every functional; at rank ≥ 2 the arithmetic threads (spectral rank) stay **rigid** while
the algebraic threads (geometric rank) **collapse** — because a collision *fits* once the dimension is
large enough, whereas whether a sub-floor covering *fits* is Diophantine.

## The two sides, grounded

**JC = GEOMETRIC rank (collapses at a collision). ✓ verified this session.** The Claude-discovered
Keller map (THM-1300; owner-corrected provenance — *not* "Alpöge", MISTAKE-205 — Claude-discovered,
otherwise uncertain) `u=1+xy`, `F=(u³z+y²u(4+3xy), y+3xu²z+3xy²(4+3xy), 2x−3x²y−x³z)` has `det J_F ≡ −2`
and the triple collision `F(0,0,−¼)=F(1,−3/2,13/2)=F(−1,3/2,13/2)=(−¼,0,0)` — re-verified in exact
arithmetic (`keller_verify_and_rank_opus`). The collision is **1 fixed branch (`x=0`) + 2 σ-swapped
branches (`x=±1`)** under the torus doubling `λ↦λ²`, `σ=diag(−1,−1,1)` — the repo's **`1+2k` census
motif**, and it genuinely needs **3 branches ⟹ geometric rank 3**. It provably cannot descend to dim 2
(boxeph-S144 (B)), so **JC₂ (the floor) survives** while JC₃ fell: exactly "the nullcone collapses once
the rank is high enough to hold a collision," with the floor rescued because rank-3 can't fit in dim 2.

**A SECOND geometric example confirms the pattern: GMC. ✓ verified this session.** GMC (Gaussian
moment, Derksen–van den Essen–Zhao) is *proved* at `n=1`, open at `n=2/3`, and **false at `n=4`**
(THM-1480). I re-verified the rank-4 counterexample `P=(1+W)(W̄−|Z|²)`, `Q=W`: `E[P^m]=0` and
`E[WP^m]=m!` for `m=1..7` (`gmc4_verify_opus`, via the Wick rules `E[W^a W̄^b]=a!δ_{ab}`,
`E[|Z|^{2k}]=k!`). And the **descent to `n=3` is structurally blocked** exactly as the frame predicts:
replacing `|Z|²` by a real square `X²` sends `k! ↦ (2k−1)!!`, and the collapsing sum becomes
`0,1,0,9,0,225` — it does **not** vanish, so the collision cannot descend and the floor (`n≤3`)
survives. So **both geometric threads collapse at a finite rank where a collision fits (JC at 3, GMC
at 4) with a structurally-blocked descent to a surviving floor** — two verified instances of one law.

**LRC = SPECTRAL rank (stays rigid; the collision is Diophantine).** CURRENT-FRONTIER.md's rank-11
residual is literally mod `91 = 7·13` (relations of height `91^6`, the Kelvin cone `(1/91)R⁻¹Kᵒ`), and
kps-S131b's crux object is a **mixed unit character mod 91, nontrivial on both `F_7` and `F_13`** — a
**spectral rank-2** (7⊗13 tensor) no-cancellation. There is no ambient dimension in which a sub-floor
covering "fits for free"; whether it exists is arithmetic, conjectured never. That is why LRC(14) is
still a conjecture and JC(3) is a theorem-of-falsity: **same nullcone shape, spectral vs geometric
rank.**

## The AP is the shared vertex — but only functionally (a tested caveat)

The AP is the nullcone vertex in both (transitive tournament / dilated AP of speeds / `xⁿ` charge).
✓ Confirmed the AP uniquely maximizes additive/`E₃` triples for `n=6,7,8` (THM-730,
`ap_universal_extremal_test`). **Caveat (tested):** "AP = universal extremal" is *functional-dependent*
— the AP minimizes the reach **tail** `μ_{1/7}` (opus-S130) but **not** the reach **mean** `E[maxgap]`
(opus-S133: stretched shapes beat it). This is itself an instance of boxeph's divergence: which
extremal you get depends on the functional, i.e. on whether you weight the tail (spectral) or the mean.

## Maturation-arc recurrences (both programs, dated)

1. **Correct the object before sharpening the technique.** LRC's `L→M` reframe (mac-mini, THM-523,
   2026-06-16: "the disproof search produced the proof reduction") = JC's "whatever the provenance, the
   *map* settles it" (verify-then-emit). Promoted to META-PATTERNS.
2. **Sampling-artifact / bounded census misses the extremal.** LRC: the Kravitz band `s/(13s+1)`
   accumulating at `1/14` was missed by small-support scans (kps-S134, SGC(13) retraction) — and I
   personally re-hit this twice (S131 saturated-margin, S133 `E[maxgap]` AP-min). JC: "survived every
   average-case attack, fell to a *structured* construction" (opus-S413). Same genus at two scales.
3. **The frontier migrates to the floor.** LRC settled `≤13`, open at 14 ⟶ `n=12` inverse theorem.
   JC dead `≥3`, survivors JC₂/DC₂/DC₁. "The open problem lives at the atom" (deathstar-S59m).

## Tests run this session (✓ kept / ⌦ released)

- ✓ **Jumping-champions**: the prime-pair singular series peaks at primorials (210, 30, 42, …) —
  confirms kps's "smooth = constructive interference of prime lines" (the grammar's central claim).
- ⌦ **Smoothness selects LRC14**: released — a crude "richness" metric does NOT single out `n=14`
  (`n=15` scores higher). The primes 7,13 are just `14=2·7` and AP-length 13, not a selection.
- ✓ **AP = additive-triple extremal** (n=6,7,8), but ⌦ *not* the reach-mean extremal (functional-dep).
- ✓ **Keller map** verified (geometric rank 3, `1+2` collision).
- ⌦ **`−¼` = LRC 3-runner floor** (klein-S324): released — `−¼ = −∏cos(2πj/3)` is a root-of-unity
  value; the actual LRC 3-runner floor is `1/3`. Charming, mislabelled.

## The live prediction (the hottest lurking pattern — agent D5)

Both programs are, *right now* (2026-07-24/25), converging on **a signed 1-cocycle / holonomy / flux
on a branched cover**: LRC's "word-current / owner-side acute current / holonomy-groupoid / endpoint
cocycle" (Wall B) and JC's "flux (Keller one-form) on the trigonal spectral cover / pole-descent". My
seam predicts they are **the same `H¹` obstruction shape — nonvanishing of a cocycle around a cycle in
the cover — differing only in rank-type**: LRC's cocycle is **spectral** (values a mixed character mod
91, the 7⊗13 tensor) while JC's is **geometric** (values in the `S₃`-monodromy of the degree-3 cover).
If so, "prove the cocycle doesn't cancel around the cycle" is *one* method-target with two arithmetics.
This is the concrete thing to test next (compare the LRC blocker-cycle holonomy to the JC Keller-flux
as `H¹` classes). **[SHAPE-CONFIRMED this session; not a literal identity.]** Reading both live
frontiers: LRC's word-current is explicitly *"a cyclic zero-divisor problem"* with a **first-Bockstein
sidecar** (the `H¹→H²` connecting map — THM-2337/2356) on the **primitive-91** deepest edge; JC's
degree-18/22 closures run on **flux (Keller one-form) / pole-descent** on the trigonal cover. So both
frontiers *are* "nonvanishing of a signed 1-cocycle on a branched cover" — the `H¹` shape is real on
both sides. They are **not the same class** (LRC's coefficients are the mod-91 character group; JC's
are the `S₃` monodromy / dyadic tail), which is exactly the spectral-vs-geometric divergence: **same
cohomological shape, different rank-type coefficients.** The map is not "LRC cocycle = JC cocycle"; it
is "both are the `H¹`-coboundary question, one with arithmetic coefficients, one with geometric."

## Bottom line

boxeph's "two nullcones that diverge" diverge along **spectral rank (prime-line count, Diophantine,
rigid — LRC) vs geometric rank (dimension, algebraic, collapses at a collision — JC)**. The rank-1 seed
is shared and proved; at rank 2 LRC's `7⊗13` no-cancellation stays open while JC's dim-3 collision
already fell. Both frontiers are now the *same* `H¹`-cocycle-on-a-branched-cover object, one spectral,
one geometric — which is the sharpest form of the analogy the repo has, and the place to test whether
it becomes a map.
