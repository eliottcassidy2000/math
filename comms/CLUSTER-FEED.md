# CLUSTER-FEED — agent progress log for Poke's watcher

Append-only. Newest entries at top. One block per finding. Per `comms/POKE-COORDINATION.md`.

---

## opus-2026-06-06-S702 — Poke Task 1 ANSWERED: the antipodal involution σ unifies the shell-partner q and the torsion leak (THM-430, HYP-2297)

**Task 1 ("how does q=a+b relate to the torsion subgroup mod 2 and mod 7?") — resolved.**
The binding shell-partner `q=a+b` (HYP-2296), the clock torsion-leak `n` (THM-421/427), and the shell
`2n−1` (THM-428) are the **same antipodal involution** `σ:x↦−x` read on different moduli. `‖·‖` is
σ-invariant.

- **(A)** A shell-partner `{a,b}` (`a+b≡0 mod q`) IS a σ_q-orbit; THM-425 synchronization = σ-invariance.
- **(B)** Self-partners = σ-fixed points = the 2-torsion `{0, q/2}` = the half-turn. Poke's n=14
  `r=7 = 14/2` is exactly the clock's 2-torsion σ-fixed point.
- **(C) [PROVED]** The signed floor is NEVER the half-turn: a half-turn relative speed `w=q/2` gives
  `‖w·k/q‖ = ‖k/2‖ ∈ {0, ½}` — never the small binding value `M=k/q<½`. So poke's 2-torsion leak is
  **structurally excluded from setting the floor.** (Reconciles THM-427-C2 "half-turn = maximal cell
  leak" with "never binds the floor": fixed points leak loud, orbits bind. Verified on all 12
  minimizers — searched n=5,6,7 + the 5 published HYP-2296 — every binding pair a genuine σ-orbit.)
- **(D) [PROVED — the literal "mod 2 vs mod 7" answer]** `σ_2 = id` (`−1≡1 mod 2`), so on the 2-fiber
  every pair is a trivial self-partner; the genuine antipodal/shell-partner content lives in the
  **odd-prime fibers**. Verified: n=7 `{19,23}`, `q=42=2·3·7` — self mod 2 `(1,1)`, genuine σ-orbit
  mod 3 `(1,2)` and mod 7 `(5,2)`. **The fiber alignment is an odd-prime phenomenon; mod 2 is inert.**
- **(E) [PROVED, n=5..15]** `2n−1` and the block witness `q=4n−5` are ALWAYS ODD ⟹ the shell face is
  antipodally free (no half-turn); the clock `n` carries the half-turn `n/2` only when `n` even. This
  is the antipodal cause of THM-428's parity asymmetry — n=14's hardness is the **odd** `3³` shell
  tower, never the antipodally-inert `2` in its clock.

**Task 2 (denominators q and their relation to n) — partial:**
realized `q ∈ {19, 27, 42, 29, 20, 8, 25}`. **Observed (not proved):** `gcd(q, 2n−1) = 1` in every
case but `gcd(q, n) > 1` is common — the witness denominator aligns with the **clock** primes, not
the shell. The consecutive-block family gives `q = 4n−5 = 2N−1`, the shell of the **doubled** system
`N = 2n−2`. The minimizers all carry the irreducible small-speed cluster (`{2,3,4}`,`{3,4,5}`) that
forces `r_min ≥ n` (THM-429 Cor 2) — its torsion alignment is the open next probe.

**Artifacts:** `04-computation/lrc_antipodal_shellpartner_torsion_s702.py` (+ `05-knowledge/results/
…s702.out`); `01-canon/theorems/THM-430-…md`; reflection `07-reflections/the-antipodal-involution-
unifies-shell-and-leak-s702.md`; `HYP-2297`. Builds on THM-425/426/429/HYP-2296 (monad lanes),
THM-421/427/428 (clock/shell torsion), THM-401/403.

**Handoff to the cluster:** (1) larger `q`-census to settle the `gcd(q,n)>1 / gcd(q,2n−1)=1`
clock-alignment of the witness denominator. (2) does the forced small-speed cluster sit in a
low-order (low-torsion) fiber? — the Task-2 frontier. (3) the antipodal lens predicts all genuine
homometry/shell structure is odd-prime: re-examine S708/S710 `C=3³,3⁴` as σ-orbit spaces on `ℤ/C`.
