---
source: opus-2026-06-03-S599t (remote-control)
status: SYNTHESIS + verified improvement — Hamkins' infinite-game ordinal recursion (value = mex over reachable positions; spectrum of achievable ordinals) is the transfinite mirror of the tournament H multiplicative recursion (spectrum of achievable integers, gaps 7,21). The unifier is the PARTITION FUNCTION: H multiplicative over strong components = the disjunctive game SUM = Z factorizing over independent subsystems. Verified: Z_n=Σ_T H(T)=n!·2^(C(n,2)−(n−1)), the cluster-expansion recursion, and S598's clearance entropy = free energy.
tags: [partition-function, hamkins, infinite-go, ordinal-recursion, game-values, cluster-expansion, free-energy, strong-components, disjunctive-sum, H-spectrum, redei, lee-yang, LRC, S598, n+2-recursion]
---

# Partition functions everywhere, and the infinite-game recursion

**Prompt (user):** [Hamkins tweet] see the n+2 recursion connections among others with infinite Go;
use these insights to make math improvements; come to see partition functions everywhere.

Hamkins' infinite-games program defines a position's value by an **ordinal recursion** — a *won*
position has value `0`, otherwise the value is the **least ordinal not among the values reachable
in one move** (a mex recursion) — yielding a **spectrum of achievable game values** (`ω⁴` in
infinite chess; *every* countable ordinal in infinite draughts)
([Hamkins, infinite games](https://jdh.hamkins.org/tag/infinite-games/);
[Transfinite game values PDF](https://jdh.hamkins.org/wp-content/uploads/Transfinite-game-values-in-infinite-games.pdf)).
That is the transfinite mirror of the structure I found last session, and the bridge is the
**partition function**.

## 1. The unifying picture: recursion → spectrum → partition function

> **Infinite games (Hamkins):** value = `mex` over moves → an **ordinal** spectrum (which
> ordinals are achievable).
> **Tournaments (last session):** `H` = `∏` over strong components → an **integer** spectrum
> (which integers are achievable; gaps `7, 21`).
> **Same shape:** a *recursively-defined value* of a *decomposable* structure, whose set of
> achievable values is a **spectrum with gaps**. The finite (integer) case is the shadow of the
> transfinite (ordinal) case; the infinite-game limit is `n → ω`.

The decomposition that makes both work is the **disjunctive sum**: a game splits into independent
regions whose values **add** (surreal/Grundy sum); a tournament splits into strong components
whose `H` **multiplies**. Addition and multiplication are the *same operation under `log`* — and
that operation is exactly the **partition function**: `Z = ∏ Z_i` over independent subsystems,
`log Z = Σ log Z_i` (free energy is **extensive**). *The disjunctive game sum and the
strong-component product are one structure: the factorization of a partition function.*

## 2. The verified improvement — the tournament partition function

> **Theorem (verified, `…s599t.py`).** The total partition function
> ```
>   Z_n := Σ_{tournaments T on n vertices} H(T)  =  n! · 2^{\binom n2 − (n−1)}
> ```
> (= the number of (tournament, Hamiltonian-path) pairs; each of `n!` orderings is a Ham path in
> `2^{\binom n2−(n−1)}` tournaments). Checked: `Z_n = 1,2,12,192,7680` for `n=1..5`, matching the
> closed form exactly.

Because `H` is multiplicative over the *linearly ordered* strong components (condensation), `Z`
obeys the **cluster expansion** (the connected = *strong* pieces generate everything):
```
   b_n = Σ_{k=1}^{n} \binom nk a_k b_{n−k},    b_0 = 1,   a_k = Σ_{strong T on k} H(T),
```
i.e. the EGF identity `T_H(x) = 1/(1 − S_H(x))` (sequence of connected pieces). **Verified** for
`n≤5` (`a = 0,6,120,6000`; recursion reproduces `Z_n`). This is the **Mayer/exponential-formula
cluster expansion**: the *free energy* `log Z` is the generating function of the **connected
(strong) components** — the disjunctive-sum decomposition, made a partition-function identity.

> **Rédei = the `z=−1` fugacity slice.** `H(T)` is odd for *every* `T` (Rédei), so the partition
> function evaluated at the **sign character** (`z=−1`, the fermionic / Lee–Yang slice) is the
> `GF(2)` permanent→determinant collapse (S599e). The H-spectrum gaps `{7,21}` are **forbidden
> energy levels** of this partition function — Lee–Yang-type exclusions.

## 3. Partition functions everywhere (the panorama)

Once you see it, every repo object is a partition function:

| object | partition function | "components" / free energy |
|---|---|---|
| tournament `H` total | `Z_n = n!·2^{\binom n2−(n−1)}` | strong components (this note) |
| LRC covering-depth `{p_k}` | spectral measure of `M_N` (THM-406) | safe-set components (S574) |
| LRC `(★)` `p_0=Σ(−1)^{|S|}μ(∩)` | grand-canonical / fugacity expansion | inclusion–exclusion clusters |
| **S598 clearance entropy** | `log μ(SAFE)=Σ log cᵢ` | **= the FREE ENERGY**, additive over the cascade |
| unit-distance additive energy | 2-point correlation function | the lattice "bulk + boundary" |
| Rédei / OCF parity | `Z` at `z=−1` (fermionic) | odd cycles (Lee–Yang zeros) |
| infinite-game value | the `mex` ground-state recursion | independent regions (disjunctive sum) |

> **Retroactive payoff:** S598's *clearance-entropy ledger* `log μ(SAFE) = Σ log cᵢ` is now
> *explained* — it is a **free energy**, extensive (additive) over the cascade's components
> exactly because `μ(SAFE)` is a partition function and `log` of a product is a sum. The cascade
> clearances `cᵢ` are the **per-component partition functions**; the worry-set's `cᵢ=0`
> (trapped runner) is a **zero of the partition function** = a Lee–Yang zero = a phase
> transition (the loneliness collapse). The three S598 order parameters are the standard
> thermodynamic trio (entropy, free energy, correlation length / Helly number).

## 4. The `n+2` recursion and the infinite limit

The cluster recursion `b_n = Σ_k \binom nk a_k b_{n−k}` is the **full** convolution; the repo's
**Mode B `n→n−2`** (Cayley–Dickson descent, the metagraph recursion) and the LRC **parity split**
(`n` even/odd, THM-404; the `2n−1` shells) are its **parity slices** — the *additive* `+2`
recursion face versus the *multiplicative* `×2` (doubling) face, the user's recurring
`+2`-vs-`×2` theme. The partition function carries both: doubling is the fugacity rescaling
(`z→z²`), `+2` is the two-step transfer matrix. **The infinite limit closes the analogy:** the
*thermodynamic limit* `n→∞` of `Z_n` (where the per-vertex free energy and the spectrum
stabilise) is the finite-combinatorics mirror of Hamkins' *transfinite limit* `n→ω` (where game
values become ordinals). Infinite Go is the lonely-runner/tournament partition function *at its
ordinal horizon*: the recursion no longer terminates in finitely many steps, so the "value
spectrum" becomes ordinal-valued.

## 5. Honest status

- **Verified:** `Z_n = n!·2^{\binom n2−(n−1)}` (`n≤5`); the cluster-expansion recursion `b_n=Σ\binom
  nk a_k b_{n−k}` with `a=` strong-component sums; Rédei oddness as the `z=−1` slice.
- **Established (rigorous):** `H` multiplicative over strong components ⟹ the partition-function /
  cluster-expansion structure; S598 clearance entropy = free energy (a *re-reading*, now derived).
- **Synthesis (framed, honest):** the recursion→spectrum parallel between infinite-game ordinal
  values (Hamkins) and tournament H-values is a *structural analogy* (both are achievable-value
  spectra of a recursively-defined invariant on a decomposable structure), not a theorem relating
  ordinals to integers; the `n→∞ ↔ n→ω` "infinite limit" is a viewpoint.
- **The `1.014` exponent (still open from prior sessions):** not addressed here; awaiting the
  exact definition.

**Artifacts:** `04-computation/H_partition_function_cluster_expansion_s599t.py` (+`.out`). Builds
on S599s (strong-component multiplicativity), THM-406 (LRC spectral measure), S598 (clearance
entropy), S599e (Rédei/`z=−1`), Cayley–Dickson Mode-B (CLAUDE.md). External:
[Hamkins infinite games](https://jdh.hamkins.org/tag/infinite-games/). New: **HYP-2182**.
