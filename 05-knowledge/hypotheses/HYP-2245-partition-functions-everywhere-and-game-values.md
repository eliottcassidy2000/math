---
id: HYP-2245
title: Partition functions everywhere — covering-depth Z(z), game values = altitude, and the n+2 (±-pair) stride
status: OPEN (synthesis); the Z-identities + p_{n-1}=1/C(n,2) are VERIFIED
source: claudebox-2026-06-03-S626
related:
  - HYP-2200  # the four lenses (partition function was one) — extended here
  - HYP-2180  # iterated logs = altitude (= the ordinal game value)
  - THM-410   # moment-sieve identity Z(z)=Σ S_k (z-1)^k
  - HYP-2240  # the shell tower (here: an Euler/zeta product = a partition function)
  - HYP-2230  # the grid-disproof bridge (generating functions of resonances)
---

# HYP-2245 — partition functions everywhere

Inspired by infinite-Go transfinite game values (Hamkins–Evans: positions carry ordinal game
values by recursion; ordinal **natural sum** combines independent subgames). The unifying claim:
**every LRC/tournament invariant is a coefficient or special value of a partition function, and the
resonance/rigidity structure is that partition function's factorization/singularity structure.**

## The master partition function (verified)

The covering-depth distribution is a statistical-mechanics partition function
```
Z_S(z) = ∫_0^1 z^{depth_S(t)} dt = Σ_k p_k z^k,   depth_S(t) = #{v∈S : ‖v t‖ < 1/n}
```
(`z` = fugacity, runners = sites, depth = occupation number). Verified (`lrc_partition_function_s626`):
- `Z_S(1) = 1`; **`Z_S(0) = p_0`** = lonely measure = the **ground state**;
- **`Z_S'(1) = 2(n-1)/n`** = mean depth, the conserved charge (THM-410 R1);
- **`Z_S(z) = Σ_k S_k (z-1)^k`**, `S_k` = inclusion–exclusion/sieve moments (THM-410 R2);
- **NEW closed form: the top coefficient `p_{n-1} = 1/C(n,2) = 2/(n(n-1))`** (the measure where ALL
  `n-1` runners overlap), verified n=3..10.

**Mean-field baseline & resonance.** The "free"/independent partition function is the binomial
`Z_free(z) = (1 - 2/n + (2/n)z)^{n-1}`, with `Z_free(0) = (1-2/n)^{n-1} → e^{-2}`. The runners
always interact through the shared clock `t`, so `Z_S` is NOT this product; the **interaction free
energy `log(Z_S/Z_free)` is the resonance**, and for the tight/AP config the resonance drives the
ground state to `p_0 = 0` (a phase transition). LRC = "can the interaction empty the ground state?"
— and the union bound is vacuous precisely because `Z'(1) = 2(n-1)/n` is fixed (HYP-2195: content is
in correlations = the non-product part of `Z`).

## Game values = altitude = order of the partition function (the infinite-Go bridge)

Infinite-Go game value of a position = an ordinal built by a `mex/sup` recursion (depth of forced
play). This is exactly the **iterated-log "altitude"** (HYP-2180): the number of nested averaging
levels in the LRC sieve = the recursion depth = the ordinal game value. The partition function is
the **generating function of that game tree**; the altitude is the order of its controlling
singularity. And **ordinal natural sum (combining independent infinite-Go subgames) = product of
partition functions** — the LRC analogue is: `Z_{A∪B} = Z_A·Z_B` would hold for resonance-independent
runner sets; its failure (verified, `{1}∪{2}`) is the shared-clock interaction = the obstruction,
the same way Go subgames stop being independent when they share liberties.

## The n+2 (±-pair) stride

The infinite-Go gadgets advance by a fixed skip (the "n+2 recursion"). The LRC stride is `+2` too:
`n → n+1` advances the **shell modulus `2n-1 → 2n+1`** by 2; the building block is the **±-pair**
(complex conjugation on `(ℤ/m)*`, S625/THM-418), the 2-element orbit; and the **even-fold `n→n/2`**
(the 2-adic seam, n=14→7) is the partition-function "renormalization" step. So the cyclotomic shell
tower (HYP-2240) is itself a partition function built by the n+2 stride — an **Euler/zeta product
over shells**, each shell a local factor.

## Partition functions everywhere (the catalogue to unify)
- covering-depth `Z_S(z)` (this note) — ground state `p_0`, charge `Z'(1)`, sieve `Σ S_k(z-1)^k`;
- the **shell/cyclotomic tower** (HYP-2240) as an Euler product `∏_shells (local factor)`;
- **game-value generating function** `Σ_α (#positions of value α) x^α` (infinite Go / altitude);
- tournament **trivariate GF** (THM-114), **independence polynomial**, the **arc menu = A000016**
  (S521), **(ternary) Krawtchouk** (S622) — each a partition function on the metagraph;
- the grid-disproof's count of modulus-1 elements (HYP-2230) = a theta/partition function of a CM
  lattice.

**Conjecture (to make precise).** All of these are specializations of one bivariate partition
function on the (phase-ring × shell-tower) whose `z→0` limit is loneliness, whose `±`-pair/n+2
stride generates the tower, and whose order-of-vanishing is the altitude/game value. Pin the single
object and the LRC frontier (n=14) becomes "the first ramified factor of this partition function."

## To do
1. Prove `p_{n-1}=1/C(n,2)` (near `t=0`: all bad ⟺ `t < 1/(n(n-1))`; account for all deep-overlap
   points) and find `p_{n-2}` etc. (the top of `Z` as a clean polynomial).
2. Make the "game value = altitude = order of `Z`" precise: define the loneliness game whose value
   is the iterated-log altitude and whose generating function is `Z`.
3. Write the shell tower as an explicit Euler product and identify the n+2 stride as its functional
   equation; locate `n=14`'s `3³` as the first ramified Euler factor.
