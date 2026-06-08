# THM-446 — Sidon = C4-free = minimal additive energy; the cauldron/Sidon/B_h additive-relation ladder, whose dyadic rungs are the power-of-two cycles of Erdős Problem 64

**Status:** PROVED (Sidon ⟺ C4-free ⟺ minimal additive energy — classical; recomputed) + REFRAME
(the additive-relation ladder; Erdős 64's hard core = high-girth/Sidon-Cayley graphs) + VERIFIED
(the hard core carries a power-of-two cycle, all tested cages/random girth-≥5 cubic). **Erdős 64
remains OPEN; this is a reframing + exploration via the repo's Sidon/summand/cauldron objects, not a
proof.**
**Source:** opus-2026-06-08-S710. The user's ask: connect Sidon sequences to the repo's "cauldron
game" and "summand graph" and use the insight on Erdős Problem 64. Builds on the summand graph
(`a→n` iff `a+b=n`, THM-401 pair-sum/shell), the cauldron game (Schur 3-term, S618/HYP-2190),
S706/THM-441 (additive energy = ‖autocorrelation‖²), THM-445 (Erdős 64), THM-442 (additive↔mult).

## (1) Sidon ⟺ C4-free ⟺ minimal additive energy (PROVED/classical)

For a finite set `S ⊂ ℤ`:
```
   S is Sidon (B_2: all pairwise sums distinct)
     ⟺  the only solutions of a+b=c+d in S are trivial ({a,b}={c,d})
     ⟺  additive energy E(S) = #{a+b=c+d} = 2|S|² − |S| (the minimum)
     ⟺  the summand graph has every node-fibre ≤1 on S (no two S-pairs share a sum)
     ⟺  the associated additive Cayley graph is C4-FREE (no additive quadruple = no 4-cycle).
```
**Verified** (`…s710.py`): Sidon `{1,2,5,11},{0,1,3,7}` have `E=28=2·16−4` (excess 0); non-Sidon
`{1,2,3,4}` and AP `{2,4,6,8}` have `E=44` (excess 16 = the C4 / **Sidon-defect**). By S706/THM-441,
`E(S)=‖1_S⋆1_S‖²=Σ_ξ|\hat{1_S}(ξ)|⁴` — so the **Sidon-defect is the autocorrelation's excess
mass**, and a **4-cycle is one unit of additive energy** above the Sidon floor.

## (2) The additive-relation ladder (the cauldron ⊂ Sidon ⊂ B_h reframe)

The repo's additive objects are the rungs of one ladder, indexed by the number of terms in the
relation they forbid/realise:
```
   CAULDRON / Schur :  3-term   A+B=C        (the "boil")     ↔ a TRIANGLE (odd, the smallest cycle)
   SIDON (B_2)      :  4-term   A+B=C+D       (quadruple)      ↔ a 4-CYCLE  = the FIRST power of two
   B_h              :  2h-term  Σ relation                     ↔ a 2h-CYCLE
   …                                                              …
   POWER OF TWO     :  2^k-term dyadic relation                ↔ a 2^k-CYCLE = the Erdős-64 rungs
```
- The **cauldron game** lives at the 3-term (odd, triangle) rung; **Sidon** at the 4-term (`C4`,
  first 2-power) rung. The summand graph reads all rungs (its node `n` collects the pairs summing to
  `n`; a node with ≥2 pairs is a 4-cycle = additive quadruple).

> **Erdős Problem 64 = "every min-degree-≥3 graph realises a DYADIC rung."** A graph with **no**
> power-of-two cycle avoids the `4,8,16,…`-term relations **simultaneously** — it is **"dyadically
> Sidon"** (`B_{2^{k-1}}`-like at every level `k`). The conjecture says min degree ≥3 **forbids** being
> dyadically Sidon at all levels.

## (3) The hard core, via Sidon/B_h (REFRAME + VERIFIED)

The smallest dyadic rung is `C4`; a graph **with** a `C4` already satisfies Erdős 64. So the **hard
core is the C4-free (Sidon-like) min-degree-≥3 graphs** — equivalently high-girth graphs (`B_h` Cayley
graphs have girth `>2h`). **Sidon / `B_h` sets are the natural source of candidate counterexamples**
(high-girth Cayley/Levi/cage graphs); Erdős 64 says even these must carry a `2^k` cycle.

**Verified** (`…s710.py`): every tested high-girth cubic graph carries a power-of-two cycle —
Petersen (g5)→`C8`, Heawood (g6)→`C8`, Möbius–Kantor/Pappus/Desargues/Dodecahedral/Nauru→`C16`,
McGee (g7)→`C16`, **Tutte–Coxeter (g8) → `C8` (the girth itself is a power of two)**; and 14 random
girth-≥5 cubic graphs (`n=10..20`) — **0** without a `2^k` cycle. The girth-8 cage is the cleanest:
its girth IS `2^3`, so it satisfies the conjecture at the floor.

## Scope / honesty

- (1) is classical (Sidon ⟺ C4-free; additive energy) — recomputed and tied to S706.
- (2)(3) are a **reframing** of Erdős 64 through the repo's Sidon/summand/cauldron objects, with
  computational support that the conjecture's hard core (high-girth/Sidon-like cubic graphs) carries a
  power-of-two cycle. **The conjecture is OPEN; nothing here proves it.** The additive/Cayley structure
  is cleanest for Cayley/circulant graphs (the natural counterexample-source); for *arbitrary*
  min-degree-3 graphs there is no global additive labelling, so this is a *construction/exploration*
  lens, not a proof route.
- **The dyadic theme** ties to THM-442 (additive `Δ³` vs multiplicative): even cycles = additive
  parity (`2ℤ`, settled, THM-443); power-of-two cycles = the multiplicative **2-adic tower** `{2^k}`
  (open) — Sidon (`B_2` = the `k=2` rung) is the first dyadic obstruction.

**Artifacts:** `04-computation/sidon_summand_cauldron_erdos64_s710.py` (+`.out`). Reflection
`07-reflections/sidon-the-additive-relation-ladder-and-erdos-64-s710.md`. New: **HYP-2314**. Builds on
THM-401 (summand/shell), S618/HYP-2190 (cauldron), S706/THM-441 (additive energy), THM-445 (Erdős 64),
THM-442/443 (additive↔multiplicative / even cycle).
