---
id: HYP-2058
status: MIXED — Part II rigorous (proof-lite); Parts III conjectures with strong evidence
source: opus-2026-06-02-S556
related:
  - HYP-2056
  - HYP-2057
  - HYP-2055
  - THM-369
---

# HYP-2058: LRC@14 proof-lite (structured counterexamples impossible) + the multiple-of-14 tension

## Rigorous (PROVED)

- **Even-fold (Lemma 1):** even `v=2u` ⇒ `||vt||=||u·2t||`, so `M(S) ≤ M(fold(S))`.
- **Antipodal (Lemma 2):** odd `v` ⇒ `||v(t+½)||=||vt+½||`; an odd runner is
  dangerous at ≤1 of any antipodal pair `{t,t+½}`. Corollary: all-odd ⇒ lonely at ½.
- **Thm A (q=14 sieve):** `14∤v_i ∀i` ⇒ `t=1/14` is a witness.
- **Cor A′:** the AP with any subset of its even runners doubled is lonely at 1/14
  (doubling `2m→4m`, `m≤6`, never hits a multiple of 14 since `7∤m`); the only
  tight member is `V*={1..11,13,24}`. ⇒ every S555 Part-A probe is provably not a
  counterexample.
- **Thm B (necessary conditions):** a counterexample (1) contains a multiple of 14;
  (2) that multiple is the even member of the mod-7 singleton class (bridge between
  the mod-2 fold and oracle-S552o's mod-7 CRT — they obstruct on the SAME runner);
  (3) has ≥7 speeds of one parity (13 is odd ⇒ minority ≤6).

## Conjectures with evidence

- **C/C′ tight-witness lattice:** every tight n=14 config is witnessed only on
  `(1/14)ℤ`; cleaner form *tight ⟹ no multiple of 14 ⟹ witnessed at 1/14*.
  Evidence: AP & V* (only tight configs within distance 2 of AP) witnessed at
  `{1,3,5,9,11,13}/14`; 2500 hill-climbs forced to contain a multiple of 14
  reached tightness 0 times (robustly loose).
- **THE TENSION (heuristic for LRC@14):** Thm B(1) — a counterexample MUST contain
  a multiple of 14; C′ — a multiple of 14 keeps a config loose. The forced feature
  is the one that prevents non-loneliness. Target: prove `M(S) ≥ 1/14` whenever
  `14 | v` for some `v ∈ S`.
- **D e≤6 fold reduction (127/127):** `LRC(14)[e≤6] ⟸ some even-good s has no
  odd-split` (S554); open.

## Most promising next theorems
1. **C′:** prove tight ⟹ no multiple of 14 (⇒ tight-witness lattice).
2. **The tension lower bound:** `14|v ⇒ M(S) ≥ 1/14` (would handle Thm B's forced
   runner directly).
3. **D:** prove no-odd-split for e≤6 (⇒ LRC(14) even-minority half, cond. on LRC(7)).

**See:** `07-reflections/lrc-n14-proof-lite-impossibility-and-provable-statements-s556.md`,
`04-computation/lrc_n14_proof_lite_and_brainstorm_s556.py` (+.out),
`05-knowledge/results/lrc_n14_tight_lattice_B3_s556.out`,
`lrc_n14_tight_no_mult14_test_s556.out`; HYP-2056 (even-fold), HYP-2057 (mod-7),
HYP-2055 (V*), THM-369.
