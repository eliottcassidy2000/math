---
source: opus-2026-06-01-S553 (remote-control)
status: reflection + EXACT census (NOT a proof of LRC@14); refines/corrects concurrent oracle-S552/S553
tags: [LRC, n14, tight-family, minimax, spectral-gap, antipodal-transversal, sporadic, doubled-prime, S552, S553]
---

# The LRC tight family at n=14 is exactly the AP — and the non-unit-pair hole

**Prompt (user, remote-control):** try more creative approaches to a proof of
LRC at n=14 and above.

**Honesty up front.** LRC@14 is open. This is not a proof. It is an exact census
of the *tight family* (the only configurations that can carry LRC's difficulty),
a clean structural fact about where tight witnesses sit, and a correction to a
concurrent reduction.

## 0. Why the tight family is the whole game

A speed set `V` (observer + `m=n-1` distinct gcd-1 speeds) has max-collar
`M(V) = max_t min_i ||v_i t||`.  LRC(n) ⟺ `M(V) >= 1/n` for all `V`.  `V` is
**tight** iff `M(V) = 1/n` exactly ⟺ its safe set has measure zero but is
nonempty.  Every **non-tight** set has positive-measure safe set, hence is
lonely.  So:

> if the AP `{1,...,n-1}` is the UNIQUE tight n-set, then LRC(n) holds.

This is false at small n (sporadic tight sets exist).  oracle-S552 proved a
*spectral gap* (`M` is `1/n` or jumps to `>= 2/(2n-1)`, exhaustive n<=8) and
oracle-S553 reduced the gap to `2^{n-1}` antipodal transversals mod `2n-1`.  This
session pins the tight family directly at the target n=14.

## 1. The census (exact, two independent methods)

`lrc_tight_census_s553.py` (measure-zero test) and
`lrc_tight_doubledprime_scan_s553.py` (fast hybrid, rigorous), cross-checked
against an independent max-collar `M(V)` computation (crossing/peak candidate
times). Over **primitive `(n-1)`-subsets of `[1, ~1.5n]`** (the range where the
small-n sporadics live):

| n  | 2n-1     | tight family (speeds <= ~1.5n)                 |
|----|----------|------------------------------------------------|
| 5  | 9 = 3²   | AP + `{1,3,4,7}`                                |
| 6  | 11 prime | AP + `{1,3,4,5,9}`                              |
| 7  | 13 prime | AP only                                         |
| 8  | 15 = 3·5 | AP + `{1,2,3,4,5,7,12}` + `{1,4,5,6,7,11,13}`   |
| 9–13 | —      | AP only                                         |
| **14** | **27 = 3³** | **AP only** (exhaustive: 203 490 configs in [1,21]) |

**0 counterexamples anywhere.**  So at n=14 the tight family is exactly the AP
within speeds `<= 1.5n`, and sporadics are a **small-n phenomenon** (n=5,6,8) —
NOT a doubled-prime effect (n=6 = 2·3 *has* a sporadic; the primes 11,13 do not).

## 2. The non-unit-pair hole (correction to oracle-S553)

oracle-S553's Link-1: a config with no speed `≡0 (mod 2n-1)` whose residues
**miss** an antipodal pair `{a, 2n-1-a}` is lonely with surplus, via the witness
`t = a^{-1}/(2n-1)` — **but this needs `a` to be a unit mod `2n-1`.**  When
`2n-1` is composite there are non-unit residues, and a config can miss *only*
non-unit pairs, escaping Link-1.

The two **n=8 sporadics are exactly such escapees** (verified, exact `M=1/8`):

- `{1,2,3,4,5,7,12}` mod 15: misses only `{6,9}` (6 non-unit), doubles `{3,12}`.
- `{1,4,5,6,7,11,13}` mod 15: misses only `{3,12}` (3 non-unit), doubles `{4,11}`.

They are **tight, gcd-1, residue-0-free, and NOT antipodal transversals** — so
they sit outside oracle-S553's `2^{n-1}` transversal reduction.  Consequence:

> the antipodal-transversal reduction is INCOMPLETE exactly when `2n-1` is
> composite; the missing class is "non-transversals whose only missed pairs are
> non-units," and at n=8 that class contains genuine tight sets.

For the target, **`2n-1 = 27 = 3³` is composite**, so the hole is present *in
principle* at n=14.  The census result (§1) is the news: **that hole is EMPTY in
range at n=14** — no non-AP tight set exists among primitive 13-subsets of
`[1,21]`, transversal or not.  So both oracle's reduction-on-transversals and
this census-on-all-configs agree at n=14: AP alone, in range.

## 3. The tight-witness lattice (structural, ties to the sieve)

Every tight config found — AP *and* all sporadics (n=5,6,8) — is lonely
**exactly at the `n`-division points `t = j/n`** (verified):

- AP@14: `{1,3,5,9,11,13}/14`;  n=8 sporadics: `{1,3,5,7}/8`;  etc.

So **every tight configuration is certified by the modulus-`n` sieve** (the
`q=n` case of THM-369): tightness lives entirely on the `1/n`-lattice.  This is
the positive companion to last session's HYP-2052 (the *bounded sieve* cannot
prove LRC because loaded **non-tight** sets escape to huge denominators) — the
hard **tight** sets, by contrast, are all caught at the single modulus `q=n`.

## 4. The apex edge

`M(A_14) = 2/27` confirmed for the apex-doubled AP `A_14 = {1,...,12,26}` (the
S552 gap edge); within the bounded stress test no config has `M ∈ (1/14, 2/27)`.

## 5. Status and the route it sharpens

**Not a proof.**  Two independent reductions now converge on the same crux:

> **LRC@14 ⟺ the AP is the unique tight 13-set (over ALL speeds).**

The census makes this *exact and verified for speeds `<= 1.5n`*; oracle's gap
makes the surplus quantitative (`>= 1/(n(2n-1))`).  The remaining gap to a proof
is twofold and now precisely located:

1. **Unbounded range.**  Rule out tight 13-sets with a speed `> 1.5n` (the n<=8
   sporadics all sit `<= 1.6n`, so the range is calibrated, but unboundedness is
   unproven).  A speed bound on tight sets (à la the witness-lattice §3: tight ⟹
   witness at `j/n` ⟹ residues mod `n` constrained) would finitize this.
2. **The non-unit-pair hole.**  Close oracle's Link-1 for non-unit missed pairs
   (the class that produces n=8's sporadics).  Since `2n-1 = 27` at n=14, this is
   not optional for a clean proof; the census says it is empty here but the
   *mechanism* must be argued.

**Handoff to oracle (S552/S553):** your transversal reduction misses the
non-unit-pair class; n=8 exhibits tight members of it.  Suggest: a second witness
family for non-unit missed pairs (e.g. lift mod `2n-1` to mod `p·(2n-1)` for a
prime `p`, restoring invertibility), and use the `t=j/n` tight-witness lattice
(§3) to bound tight-set residues mod `n`.

**Artifacts:** `04-computation/lrc_tight_census_s553.py`,
`lrc_tight_doubledprime_scan_s553.py`,
`lrc_tight_witness_structure_s553.py` (+ `.out`s in `05-knowledge/results/`).
**Builds on / refines:** oracle-S552 (spectral gap, doubled-apex), oracle-S553
(antipodal-transversal reduction), HYP-2052 (sieve incompleteness, S551),
HYP-2045/2039, THM-369.  New: **HYP-2055**.
