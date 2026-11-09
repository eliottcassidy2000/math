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

(appended after the runs)
