---
source: opus-2026-06-08-S710 (user: connect Sidon sequences to the repo's cauldron game & summand
  graph; use the insight on Erdős Problem 64 = power-of-2 cycles, OPEN)
status: REFRAME + exploration (Erdős 64 stays OPEN). The repo's additive objects form ONE
  additive-relation ladder indexed by #terms: CAULDRON/Schur (3-term, A+B=C, triangle) ⊂ SIDON/B_2
  (4-term, A+B=C+D, = C4 = first power of two = minimal additive energy) ⊂ B_h (2h-term, C_{2h}).
  Sidon ⟺ C4-free ⟺ minimal additive energy (E=2|S|²−|S|; Sidon-defect = excess autocorrelation mass,
  S706). The DYADIC rungs {2^k} of this ladder ARE the power-of-2 cycles of Erdős 64: a graph with no
  2^k cycle is "dyadically Sidon" (B_{2^{k-1}}-like at every level). Hard core = C4-free (high-girth/
  Sidon-Cayley) min-deg-3 graphs; VERIFIED all tested cages + random girth-≥5 cubic carry a 2^k cycle
  (Petersen→C8, McGee→C16, Tutte–Coxeter girth 8 = 2^3). THM-446, HYP-2314.
tags: [sidon, B_h, additive-energy, summand-graph, cauldron, schur, erdos-64, erdos-gyarfas,
  power-of-two-cycle, c4-free, additive-relation-ladder, dyadic, 2-adic, cayley-graph, girth, cage,
  autocorrelation, reframe, open-problem]
---

# Sidon, the additive-relation ladder, and the dyadic hard core of Erdős 64

**Prompt (user):** consider a Sidon sequence and its connections to the repo's "cauldron game" and
"summand graph"; use insights there to make progress on Erdős Problem 64 (every min-degree-≥3 graph
has a cycle of length a power of two — OPEN).

Following the repo's own objects gave a clean reframe: **they are the rungs of one additive-relation
ladder, and the power-of-two cycles of Erdős 64 are exactly its dyadic rungs.** No proof of the
conjecture — but a precise location of its hard core and the additive structure that governs it.

## 1. The three repo objects are one ladder

- **Cauldron game** (S618): naturals dropped into cauldrons; a cauldron *boils* on `A+B=C` — a
  **3-term** (Schur) relation, the odd **triangle**.
- **Summand graph** (THM-401): `a→n` iff `a+b=n`; a node with two incoming pairs is `a+b=c+d` — a
  **4-term** additive quadruple, a **4-cycle**.
- **Sidon set** (`B_2`): *forbids* the quadruple — distinct pairwise sums.

Stacking by the number of terms in the relation:
```
   3-term  A+B=C       cauldron / Schur     triangle  C3   (odd)
   4-term  A+B=C+D      Sidon / B_2          4-cycle   C4   = the FIRST power of two
   2h-term Σ-relation   B_h                  2h-cycle  C_{2h}
```
And **Sidon ⟺ C4-free ⟺ minimal additive energy**: `E(S)=#{a+b=c+d}=2|S|²−|S|` exactly when Sidon
(verified: Sidon `{1,2,5,11}` excess 0; non-Sidon `{1,2,3,4}`, AP excess 16). By S706/THM-441,
`E(S)=‖1_S⋆1_S‖²`, so the **Sidon-defect is the autocorrelation's excess mass and one 4-cycle = one
unit of additive energy above the Sidon floor.** The summand graph is the picture of the ladder; the
cauldron lives on its odd (3-term) rung, Sidon on its first even (4-term) rung.

## 2. Erdős 64 lives on the DYADIC rungs

The power-of-two cycle lengths `{4,8,16,32,…}` are precisely the **dyadic rungs** of this ladder
(a `2^k`-cycle = a `2^k`-term additive relation). So:

> **Erdős 64 = "every min-degree-≥3 graph realises a dyadic rung."** A counterexample is a graph whose
> cycle spectrum **misses every `2^k`** — i.e. it is **"dyadically Sidon"**: `B_{2^{k-1}}`-like at
> every level `k` simultaneously, dodging the 4-term, 8-term, 16-term, … relations all at once. The
> conjecture says min degree ≥3 **forbids** being dyadically Sidon at all levels.

This is the honest content of "use Sidon to attack Erdős 64": **Sidon (and `B_h`) are exactly the
additive structures that suppress the low dyadic rungs**, and the conjecture is whether you can
suppress *all* of them with min degree ≥3.

## 3. The hard core, and where counterexamples would come from

The smallest dyadic rung is `C4`; any graph **with** a 4-cycle is done. So the **hard core is the
C4-free (Sidon-like) min-degree-≥3 graphs** = high-girth graphs. And `B_h` sets *build* them:
high-girth Cayley/Levi/cage graphs are the canonical C4-free cubic graphs, the natural
counterexample-source. Erdős 64 says even these carry a `2^k` cycle.

**Verified** (cages + random girth-≥5 cubic): all carry a power-of-two cycle — Petersen (girth 5)→`C8`,
Heawood→`C8`, Pappus/Desargues/Möbius–Kantor/Dodecahedral/Nauru→`C16`, McGee (girth 7)→`C16`, and
**Tutte–Coxeter (girth 8) → `C8` because its girth *is* `2^3`**; 14 random girth-≥5 cubic graphs, zero
without a `2^k` cycle. The girth-8 cage is the slogan: **a cubic graph whose girth is a power of two
satisfies Erdős 64 at the floor** — and the `(3,2^k)`-cages are exactly the `B_h`-optimal graphs.

> **A clean sub-question this surfaces:** among the `(3,g)`-cages (girth-extremal cubic graphs),
> the girth `g` is a power of two only for `g=4` (`K_{3,3}`) and `g=8` (Tutte–Coxeter). For
> `g ∈ {5,6,7,9,10,11,12}` the cage must reach a power of two *above* its girth — and it always does
> (e.g. McGee g7→C16). **Is there a min-degree-3 graph of girth `g` with no cycle in `[g, 2^{⌈log₂g⌉+1}]`
> equal to a power of two?** That is the precise dyadic-gap question; no example is known (Markström).

## 4. The dyadic theme (tie to the recent arc)

Erdős 64 is the **multiplicative / 2-adic** face of the cycle-length question: even cycles =
additive parity `2ℤ` (settled, THM-443); power-of-two cycles = the 2-adic tower `{2^k}` (open) — the
same additive↔multiplicative split as THM-442 (the `Δ³` cell-recursion vs the modular `H`-law) and
the same 2-adic depth as the LRC clock/shell towers (S701/S704) and the order-2 antipodal `σ` (S702).
**Sidon (`B_2`) is the `k=2` rung — the first 2-adic obstruction.** The ladder makes precise why
"even" is easy and "power of two" is hard: easy = hit the additive class `2ℤ` (one parity); hard =
hit the multiplicative-sparse dyadic set `{2^k}` (a Sidon-type avoidance at every scale).

## 5. Honest status

- **Proved/classical:** Sidon ⟺ C4-free ⟺ minimal additive energy (tied to S706 autocorrelation).
- **Reframe (new, repo-native):** the cauldron/Sidon/`B_h` additive-relation ladder; Erdős 64 = realise
  a dyadic rung; hard core = dyadically-Sidon (high-girth) graphs.
- **Verified:** the hard core (cages, random girth-≥5 cubic) all carry a `2^k` cycle — consistent with
  the conjecture, **not a proof**.
- **OPEN:** Erdős 64 is unresolved; the additive lens is a construction/exploration handle (cleanest
  for Cayley graphs), not a proof for arbitrary graphs. The sharp open sub-question (§3) is the
  dyadic-gap one for `(3,g)`-cages.

**Artifacts:** `04-computation/sidon_summand_cauldron_erdos64_s710.py` (+`.out`). Theorem **THM-446**.
New **HYP-2314**. Builds on THM-401 (summand/shell), S618/HYP-2190 (cauldron), S706/THM-441 (additive
energy), THM-445 (Erdős 64), THM-442/443 (additive↔multiplicative), S701/S702/S704 (2-adic).
