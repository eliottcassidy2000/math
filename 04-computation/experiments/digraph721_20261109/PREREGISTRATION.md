# Pre-registration: mac-mini-2026-11-09-digraph721

Written and committed before the computations that decide them. Outcomes are appended at the end; the predictions are
not edited.

Notation: E = {d(G) : G an even graph (von Brömssen)}, d = number of dominating sets. THM-4996: E contains every odd
number below 512 except 3, 7 and 13. HYP-9575 (B) conjectures that this holds for all odd numbers.

- **R1.** An explicit family of even constructions realises every odd k in [15, 2^24]. The constructions, each with
  a proof of evenness and of its value formula, are:
  - disjoint unions;
  - padded cones;
  - multi-joins of E_r's with one padded child;
  - multi-hub joins x + 2^h (prod (1 + 2^l_i) - 1);
  - partial hubs with tracked pendant structure.

  Probability 0.50.
- **R2.** A complete proof of HYP-9575 (B) (all odd k) is found this session. Probability 0.25.
- **R3.** The coverage in R1, if complete, uses partial hubs essentially: without them some odd k <= 2^24 stays
  unreached. Probability 0.70.

## Outcomes

(appended after the runs; the predictions above are unedited)

- **R1: CONFIRMED.** It holds in a stronger form, for every odd k. The constructions finally used are disjoint
  unions (only in the base), padded cones, E_2-joins, and partial hubs with tracked pendant structure: the hub move
  PH(G, u, l), its cases R and PH(3, 2), and the seeds A and B.
  - `evenclos 26` gives E(n) = F(n) for n <= 26, so every odd k <= M(26) other than 3, 7, 13 and X(n) is realised.
    That range includes [15, 2^24].
  - THM-5055 proves it for all n. Multi-joins of E_r's and multi-hub joins were not needed.
- **R2: CONFIRMED.** THM-5055 is a complete proof of HYP-9575 (B) in the sharp form. opus-2026-11-09-S89 found an
  independent proof concurrently (THM-5040); the two are cross-cited.
- **R3: CONFIRMED.** Without hub moves (`evenclos ... noR`) the closure first fails at 501 = X(9) (order 10) and
  misses 2,096,928 values of F(n) for n <= 26. Without PH(3, 2), R alone keeps the closure exact to n = 26 but must
  fail at even n >= 84 (PROOF.md remark 7.1).
