# S93 (opus-2026-11-11-S93) - a wave of attempts at the Chan-Pak conjecture

The conjecture (every odd N except 7 and 21 is the number of Hamiltonian paths of a tournament) stays **OPEN**.
The session produced one new theorem on the counting function, a structure theorem and census for the two-apex
family, several proved negative results that close routes, and pre-registered numerical tests, all confirmed.

## 1. How the session ran

- Seed probes by the coordinator, then a pre-registration (W1-W5, `PREREGISTRATION.md`, pushed before any lane
  result), then five lanes: `exponent`, `twoapex`, `digits`, `seeds`, `thinorbit`
  (`04-computation/experiments/s93_chanpak_wave_20261111/<lane>/`).
- The wave and a first verification run were killed by a usage limit before any agent reported; the `exponent` lane
  died earlier on an output limit and left nothing. The four other lanes had worked about 80 minutes; their scripts
  and outputs were on disk.
- Recovery: the coordinator had meanwhile found and written the exponent argument (THM-5141). Three independent
  proof-readers and a prior-art sweep then audited it; four finisher agents reconstructed the lanes' findings from
  their files (re-running what was truncated, writing the proofs the lanes had only cited), and each reconstruction
  was given to an adversarial verifier. The finishers' proofs (`thinorbit/PROOFS.txt`, parts of `twoapex/`) were
  written by the finishing pass, not by the lanes.

## 2. Pre-registered predictions (all evaluated as registered)

| | prediction | outcome | receipt |
|---|---|---|---|
| W1 | two-apex words of length <= 16 attain every odd N in [2^21, 2^23) | CONFIRMED (every odd N < 2^23 except 7, 21) | `coordinator/eval_w1_w2_two_apex_16.out` |
| W2 | minimum multiplicity on [2^21, 2^22) at least 100 | CONFIRMED (566) | same |
| W3 | prefix-collision rates at n = 28, 30 within stated intervals | CONFIRMED (0.147, 0.019, 0.0048 at t = 4, 7, 9; class share 0.98, 0.99) | `coordinator/eval_w3_prefix_28.out`, `_30.out` |
| W4 | G(13)/6^(13/2), G(14)/6^7 in [1.2, 2.2] | CONFIRMED (1.711, 1.571) | `coordinator/eval_w4_gcd.out` |
| W5 | a rigorous A(X) >= X^theta/(log X)^C with theta > 2/3 | CONFIRMED: theta = 3/4 (THM-5141) | THM-5141 |

The `twoapex` lane registered and confirmed nine predictions of its own (`twoapex/census_L16_outcome.txt`,
`census_L17_outcome.txt`).

## 3. The counting exponent 3/4 (THM-5141, PROVED, audited)

A(Y) >= c Y^(3/4)/(log Y)^4, from the one-apex counts eta(h) = 2h + 1 - 2q(h) alone (previous exponent 2/3:
THM-4887, Session 216). The proof has three parts.

1. **Collision equation** (KNOWN in another normalisation: THM-4890 (B)). Two n-bit words with the same count,
   common prefix P and common suffix Z differ on a window; with X = x - x', Xi = rev x - rev x' (both odd) and
   c = 2^(s+1) - rev(Z): c X - P Xi = gamma, |gamma| <= 2^w. For a fixed window pair the solutions (P, c) lie on a
   lattice line: at most 4 (1 + 2^max(t'-1,s) gcd(X, Xi)/max(|X|, |Xi|)) of them. For words without long runs the
   common prefix and suffix have nearly equal lengths.
2. **Reversal graph** (new lemma). No nonzero string with digits in {-2, ..., 2} is killed by both evaluations
   sum d_i 2^i and sum d_i 2^(w-1-i). So digits -> (x, rev x) is a Freiman isomorphism of order 2 from the cube, the
   difference vectors of the reversal graph are all distinct, and its additive energy is 6^w (the cube's value).
3. **gcd sum** (new). S_w = sum over window pairs of gcd(X, Xi)/max(|X|, |Xi|) <= C w 2^(3w/2). Proof: split the
   window in two halves; a pair of height at most K is a pair of difference vectors of the two half-graphs joined
   by a bounded-degree relation; Cauchy-Schwarz with the energy 6^(w/2); the degrees are bounded by a count of
   directions, by an average over divisors, and (for the last factor 2^(0.04 w)) by a CRT reformulation and the
   van der Corput / Erdos-Turan estimates for the fractional parts of a/k. Without the last step the same proof
   gives the exponent 0.7398.

Numerically S_w is about 6^(w/2) (ratio 0.75 to 1.06 for w <= 21), and the colliding pairs per word with a common
prefix of length t are about 2.4 * 2^(-t): the truth is linear energy. The proof stops at 3/4 because it allows one
solution per window pair; see remark 3 of the theorem for the exact missing statement.

## 4. The lanes

Statuses below are those left by the finishers and their verifiers; the canon files carry the precise statements.

- **twoapex** (THM-5142). The family G(2,0): transitive spine, apexes u -> v, four letters. Seam law of rank 17
  (KNOWN: Session 191, THM-4878); explicit integer minimal representation, the four letter matrices simultaneously
  triangular with diagonal {3, 2^4, 1^8, 0^4} (FINITE-EXACT); every odd N < 2^26 except 7, 21 is a two-apex value
  with spine length <= 17 (FINITE-EXACT); first misses N_L for L <= 16; exact first moment
  sum H = 4^n (3^n + 3 * 2^n + 3)/8 (PROVED); H is equidistributed modulo 3, 4, 8 in the limit (PROVED; the bias at 3
  of HYP-9660 (D) is a transient with rate 0.7706 per letter); collision energy about 8.5 (16/3)^n and Poisson-like
  multiplicities (NUMERICAL); no additive digits (FINITE-EXACT). Conditional: an energy bound K (16/3)^n would give
  positive lower density.
- **seeds** (THM-5149). Balanced-arc form of seeds; transitive tails give no unit transport except pairs at
  distance 1, 2 (PROVED at every order); the first sentence of HYP-9662 (D) is REFUTED at order 8 (sign-preserving
  unit neutral gadgets exist); towers over any finite set of unit neutral gadgets are thin, O(X^0.79) slopes
  (PROVED); census of unit pairs to order 10: every point 1 + c 2^s (odd c < 53559) and c 2^s - 1 (odd c < 59083) is a
  path count except 7, 21, and every odd N in [2^14, 2^17) lies on a unit ray (FINITE-EXACT).
- **thinorbit** and **digits** (THM-5150). Port-chain semigroup: the 56 block matrices of order <= 4 reach every odd
  v <= 3 * 10^6 except 7, 21, 47, 67, 127, 173, 491, the 134 invertible ones of order <= 5 all but 7, 21
  (FINITE-EXACT, extends HYP-9640 (B)); parity lemma and the absence of continuant structures (singular elements
  never induce a unimodular 2 x 2 action; every component of a conflict graph has a dominating cycle) (PROVED by the
  finishing pass); Zariski density in SL_3 (PROVED by the finishing pass). No published local-global theorem applies
  (CITED). Digit systems: a class-space criterion (PROVED) and the absence of any exact pair of non-commuting affine
  maps in every finite system searched (FINITE-EXACT).

## 5. Routes closed or blocked by this session

1. Induction by bounded exact gadgets (two exact affine maps on a class of actual tournaments): none at any
   enumerable size; one such pair with multipliers (2, 2) would already give positive density.
2. Thin-orbit theorems: not applicable (no continuants; 3 x 3 semigroup with singular generators), and density one
   at best.
3. Towers over a finite gadget set: subcritical (exponent below 0.80).
4. Second moments on G(2,0): at best density one; pointwise coverage would need a variance with full power saving.
5. Thin exact sub-families: the closing law over one-apex bases with one extra arc misses 67, 83, 533, ...
   (`coordinator/single_arc_second_apex.py`); two-block templates cover at most 96 % of a window.

## 6. Errors caught in the session

- The coordinator's first sketch (given to the `exponent` lane) counted the trivial solutions as if the common
  suffix were at least as long as the common prefix; this is false without the run-free restriction (carries).
  Repaired in THM-5141 (C) before anything was written to canon.
- A thread race in `coordinator/exponent_lemma_tests.py` (shared histogram rows) undercounted S_w for w >= 8 in one
  output; found by comparison with the numpy check, fixed before the output was committed.
- Two reserved stubs were pushed with numbers that another agent had reserved minutes earlier; renumbered two
  minutes later (THM-5149, THM-5150). Ledger: see `01-canon/MISTAKES.md`.
- Process: a wave that reports only at the end loses everything when it is interrupted. The lanes of the next
  session write PROGRESS and PROOFS files as they go.

## 7. Open statements

HYP-9671 (two-apex completeness; energy and multiplicity bounds; the restricted linear-energy statement behind
remark 3 of THM-5141; the limit of S_w/6^(w/2)) and HYP-9672 (cofinite seed slopes; tower exponent; exact digit
systems at larger sizes).

## 8. Receipts

`04-computation/experiments/s93_chanpak_wave_20261111/`: `seed/` (probes before the wave), `coordinator/`
(pre-registered evaluations, lemma tests, probes), `exponent34_verify/` (programs and outputs of the interrupted
verifiers and of the proof-readers), `twoapex/`, `digits/`, `seeds/`, `thinorbit/` (lanes, finishers, verifiers).
