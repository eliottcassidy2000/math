---
id: HYP-2055
status: SUPPORTED (exact, bounded range) — refines/corrects oracle-S552/S553
source: opus-2026-06-01-S553
related:
  - HYP-2052
  - HYP-2045
  - HYP-2039
  - THM-369
---

# HYP-2055: the LRC tight family at n=14 is exactly the AP (and the non-unit-pair hole)

## Claim

A speed set is **tight** iff `M(V)=max_t min_i||v_i t|| = 1/n` (safe set
measure-0, nonempty).  Non-tight ⇒ positive measure ⇒ lonely, so **AP unique-tight
⇒ LRC(n)**.

- **(C1, exact)** At **n=14** the tight family over all primitive 13-subsets of
  `[1,21]` (=1.5n, 203 490 configs) is **just the AP `{1,...,13}`** — 0 sporadics,
  0 counterexamples.  Same (AP-only) for n=9..13 in range; sporadics occur ONLY at
  n=5 (`{1,3,4,7}`), n=6 (`{1,3,4,5,9}`), n=8 (`{1,2,3,4,5,7,12}`,
  `{1,4,5,6,7,11,13}`).  So sporadic tight sets are a **small-n phenomenon**, not a
  doubled-prime one (n=6=2·3 has one; primes 11,13 don't).
- **(C2, tight-witness lattice)** Every tight config found — AP and all sporadics —
  is lonely **exactly at `t=j/n`**, i.e. certified by the modulus-`n` sieve (q=n
  case of THM-369).  Tightness lives on the `1/n`-lattice.  (Positive companion to
  HYP-2052: bounded sieve fails on *non-tight* loaded sets, but every *tight* set
  is caught at the single modulus q=n.)
- **(C3, correction to oracle-S553)** The n=8 sporadics are tight, gcd-1,
  residue-0-free, and **NOT antipodal transversals** mod 15: each misses exactly
  one **non-unit** antipodal pair (`{6,9}` resp. `{3,12}`) and doubles another.
  oracle-S553's Link-1 witness `t=a^{-1}/(2n-1)` needs `a` a unit, so it cannot
  reach configs whose only missed pair is non-unit.  Hence the transversal
  reduction is **incomplete whenever `2n-1` is composite**.  At n=14,
  `2n-1=27=3³` is composite, so the hole exists in principle — but C1 shows it is
  **empty in range** at n=14.

## Why it matters

Two independent reductions (oracle's spectral gap; this census) converge:
**LRC@14 ⟺ the AP is the unique tight 13-set over ALL speeds.**  C1 verifies this
exactly for speeds `<= 1.5n`.

## Open (route to a proof)

1. Extend C1 to unbounded speeds — finitize via the `t=j/n` tight-witness lattice
   (C2) to bound tight-set residues mod n.
2. Close the non-unit-pair hole (C3): a second witness family for non-unit missed
   pairs (e.g. lift mod `2n-1` → mod `p(2n-1)` to restore invertibility).

## Honest status

Not a proof of LRC@14.  Exact bounded-range census + a structural correction to a
concurrent claim.  Limitation: sporadics at speeds `>1.5n` are not ruled out
(though the n<=8 sporadics all sit `<=1.6n`).

**See:** `07-reflections/lrc-n14-tight-family-is-the-AP-and-the-nonunit-pair-hole-s553.md`;
`04-computation/lrc_tight_census_s553.py`,
`lrc_tight_doubledprime_scan_s553.py`,
`lrc_tight_witness_structure_s553.py` (+ `.out`s);
oracle-S552/S553 (spectral gap / transversal reduction), HYP-2052, HYP-2045, THM-369.
