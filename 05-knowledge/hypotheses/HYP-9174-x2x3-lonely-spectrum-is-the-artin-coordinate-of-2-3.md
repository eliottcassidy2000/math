---
id: HYP-9174
title: "The x2x3 lonely spectrum below 1/14 is the Artin coordinate of <2,3>: on the moduli q prime to 6 every value of I(t) = inf ||2^j 3^k t|| is n/q with n the least absolute residue of a <2,3>-coset (I(m/q) = n_q(m)/q exactly; PROVED), while a modulus divisible by 2 or 3 falls through its intermediate denominators (values 1/10, 1/14, 7/104, 1/15, ... exist); the first value below THM-4522's discrete top is 5/73 (73 = 1 mod 24, <2,3> = the squares, least non-residue 5; FINITE-EXACT over all rationals with q <= 4000), for every prime q = 1 mod 24 with [(Z/q)^x : <2,3>] = 2 the extremal value is (least quadratic non-residue)/q (PROVED, one line; checked on all 186 such primes below 20000), the only accumulation point of the spectrum is 0 (CONJECTURED), and the density of primes with <2,3> = (Z/q)^x is the two-generator Artin constant prod_l (1 - 1/(l^2 (l-1))) = 0.6975 (HEURISTIC; measured 0.7071 +- 0.010 to 20000)"
status: >
  OPEN (CONJECTURAL in its density and accumulation parts; the identity
  I(m/q) = n_q(m)/q for q prime to 6 and the least-non-residue law for the
  index-2 primes q = 1 mod 24 are PROVED in one line each; the census parts
  are FINITE-EXACT). CORRECTED 2026-10-04 after the blind audit
  (collatz_coalescence_20261004_audit.md, C1): the first version stated the
  coset identity for every rational; it holds for moduli prime to 6 only, and
  over all rationals the values below 1/14 begin 5/73, 7/104, 1/15, 1/17,
  1/19, 5/97, 5/99, 11/219, 13/259, 77/1539, 1/20, 7/145, 25/518, 29/601.
  Context: THM-4522 proves that the values of I in [1/14, 1/2] are exactly
  1/5, 1/7, 1/10, 1/11, 1/13, 1/14 (88 points), and that I(t) > 0 forces t
  rational (Furstenberg). This hypothesis says what the rest of the spectrum
  is: not a Lagrange-type spectrum with an accumulation point, but the
  discrete-logarithm (Artin) coordinate of the pair (2,3) in the finite cyclic
  groups (Z/q)^x. Cheapest tests: (i) extend the value census to q <= 10^5 and
  check that the ordered values are still n_q/q with n_q the coset minimum (no
  other mechanism); (ii) count primes q <= 10^6 with <2,3> = (Z/q)^x against
  0.6975 (the heuristic ignores the entanglement of Q(sqrt 2, sqrt 3) with
  the cyclotomic fields; the GRH-conditional two-generator Artin density of
  Pappalardi-type is the thing to compare with, CITED-UNVERIFIED here); (iii)
  the primes with large index (q = 6563, index 17; q = 6553, index 56) are
  where I_max(q) q is a large least non-l-th-power residue; test the Burgess
  scale n_q = O(q^(1/4+eps)) against the data.
source: opus-2026-10-04-S1 (worktree collatz-synthesis-20261004), the coalescence note, Probe A
related:
  - 01-canon/theorems/THM-4522-multiplicative-lonely-runner-x2x3-lonely-spectrum.md
  - 01-canon/theorems/THM-4523-collatz-functional-graph-is-rigid-backward-trees-identify-vertices.md (the Artin condition "2 generates (Z/a)^x" in backward separation)
  - 01-canon/theorems/THM-4520-collatz-level-operator-is-a-gauss-twisted-circulant-with-spectrum-on-the-half-circle.md (2 generates (Z/3^n)^x)
  - 05-knowledge/results/collatz_coalescence_20261004_proof_strategies_and_bold_predictions.md
scripts: 04-computation/experiments/collatz_coalescence_20261004_probes.py (Probe A)
output: 04-computation/experiments/collatz_coalescence_20261004_probes.out
---

# HYP-9174 -- the x2x3 lonely spectrum is the Artin coordinate of (2,3)

**The one-line identity (PROVED, for moduli prime to 6).** For `q` prime
to `6`, multiplication by `2` and by `3` permute `Z/q`, so the `<2,3>`-orbit
of `m` is the coset `m<2,3>` of the subgroup `<2,3>` of `(Z/q)^x` (for `m` a
unit), and
`I(m/q) = inf_{j,k} ||2^j 3^k m/q|| = (1/q) min{ min(r, q-r) : r in m<2,3> } =: n_q(m)/q`.
For `q = 2^e 3^f q'` the orbit falls through the intermediate denominators
and `I(m/q)` is the minimum of the corresponding `<3>`-, `<2>`- and
`<2,3>`-orbit minima (script, `I_values`); such moduli contribute their own
values (`1/10, 1/14` in THM-4522's top; `7/104, 1/15, 5/99, 11/219, 77/1539,
1/20, 25/518` below it), which are NOT of the coset form -- the first
version of this file overstated the identity (audit C1).

**The least-non-residue law (PROVED, one line).** If `q = 1 mod 24` is prime
then `2`, `3` and `-1` are squares mod `q`; if moreover `[(Z/q)^x : <2,3>] =
2` then `<2,3>` is the group of squares, the other coset is the non-residues,
it is closed under negation, and its least absolute residue is the least
positive non-residue `n_q`; so `q I_max(q) = n_q`. (All 186 such primes below
`20000` checked.)

**What the census shows (FINITE-EXACT, `collatz_coalescence_20261004_probes.out`).**

* THM-4522's 88 points and their denominators `5,7,10,11,13,14,26,28,33,52,56`
  are reproduced (positive control).
* The largest values of `I` below `1/14` on moduli prime to `6` up to `4000`:
  `5/73, 1/17, 1/19, 5/97, 13/259, 7/145, 29/601, 11/247, 31/697, 19/431,
  1/23, 1/25, 1/29, 13/385, 1/31, 7/241, 1/35, 1/37, 13/485, 23/865`.
  Every one is `n/q` with `n` the minimum of a `<2,3>`-coset: `1/q` when
  `<2,3> = (Z/q)^x` or when `-1` lies in every coset's reach (`q = 23`,
  index 2, `-1` a non-residue), `5/73` and `5/97` and `7/241` for the primes
  `1 mod 24` where `<2,3>` is the group of squares and `5`, resp. `7`, is the
  least non-residue, `29/601` (index 8), `19/431` (index 10), and composite
  moduli `259 = 7 x 37`, `145 = 5 x 29`, `247 = 13 x 19`, `697 = 17 x 41`
  whose `<2,3>` has index 3, 2, 3, 8.
* Over the primes `q <= 20000`: index distribution
  `{1: 1598, 2: 465, 3: 86, 4: 42, 5: 16, 6: 19, 7: 5, 8: 13, 9: 3, 10: 4,
  12: 3, 14: 1, 15: 1, 16: 1, 17: 1, 24: 1, 56: 1}`; for the 186 primes
  `q = 1 mod 24` of index exactly 2, `q I_max(q)` equals the least quadratic
  non-residue in every case; the largest `I_max` over primes with a proper
  subgroup are `5/73, 5/97, 29/601, 19/431, 1/23, 197/6563 (index 17),
  7/241, 5/193, 11/439, 109/4513 (index 12), 29/1201, 157/6553 (index 56)`.

**Why it coalesces with the Collatz thread.** The same coordinate -- the
discrete logarithm of `2` (and of `3`) in a finite cyclic group -- is what
makes the Syracuse level operator a circulant (THM-4520: `2` generates
`(Z/3^n)^x`), what the backward-tree separation lemma needs (THM-4523:
`2` must generate `(Z/a)^x`; `a = 7` fails), what Paley minus a vertex is
(THM-4532: a two-sheet clock in discrete-log coordinates, the sheet being the
parity of the exponent), and what the Syracuse clock mod `9` reads
(`log_2 S(A) = 2(A mod 3) - v mod 6`). Here it governs the multiplicative
lonely runner, i.e. the "local" problem of the owner's local/global thesis,
on the multiplicative runners of Collatz. The content the coordinate cannot
carry is the same in every instance: the archimedean sign (the sheet) and
the drift; see Axis IV and Axis V of the coalescence note.

**Non-consequences.** Nothing here touches the Collatz cycle equation: by
THM-4522(B) no crowded-time dichotomy excludes a cycle, and the present
hypothesis only describes the spectrum's shape.
