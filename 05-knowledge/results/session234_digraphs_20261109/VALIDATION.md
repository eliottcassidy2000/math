# Validation and reproduction

**PASS.** The packet separates complete finite enumeration, all-order
structural proofs, inherited classical inputs and open spectrum questions.

From this directory:

```bash
python3 replay.py
python3 replay.py --full
```

The first checks saved independent censuses and reconstructs the core
families. The second reruns the complete enumerations and compares their
outputs with the frozen artifacts. Both work in a temporary copy and
leave the evidence unchanged. Requirements: Python 3 and a C++17 compiler
named `g++`; no third-party Python packages are required.

## Exact universes and controls

1. **Semicomplete matrices:** all `3^(n choose 2)` loopless labelled
   matrices for each `1<=n<=6`. Subset DP is checked against direct
   permutations on every matrix through five and every target through
   six. A separate permutation-scatter implementation independently
   reproduces every histogram entry and all order-six target classes.
2. **Tournament floors:** all `2^(n choose 2)` labelled tournaments for
   `3<=n<=7`. The independent scan gives strong minima `3,5,9,15,25`.
   Moon's cited induced-subtournament theorem and insertion give the
   all-order monotonicity; a finite floor table alone is not extrapolated.
3. **Unrestricted matrices:** all `4^(n choose 2)` loopless labelled
   digraphs through order five. Independent implementations agree on
   every target stratum: order, count, absent pairs, reciprocal pairs,
   strongness and labelled frequency. The scatter implementation saves
   complete unrestricted histograms, but the DP implementation exports
   target strata only: full base-four histogram agreement is **not**
   claimed.
4. **Incidence totals:** for semicomplete inputs the sum of H over all
   matrices is `n! 2^(n-1) 3^(binom(n,2)-n+1)`, and for unrestricted
   inputs it is `n! 2^(n-1) 4^(binom(n,2)-n+1)`. For tournaments it is
   `n! 2^(binom(n,2)-n+1)`. The exact graph and incidence totals pass.
5. **Isomorphisms:** every finite core is canonically labelled by trying
   all vertex permutations, not by a heuristic score signature. Explicit
   automorphism groups validate labelled orbit weights and loop-subset
   orbit counts.
6. **Generators:** all loopless representatives through order eight and
   all loop-decorated representatives through order seven are generated,
   recounted and checked for distinct canonical ordered-SCC signatures.
   These are complete outputs of the proved grammar, not a brute-force
   scan of all input matrices at those larger orders.
7. **Hostile controls:** the doubled triangle has exactly three paths;
   its extra arc is unused. A one-reciprocal-pair order-six graph has
   twenty-one paths. The wrong-target replay gate rejects the seven-path
   matrix when deliberately labelled twenty-one. Diagonal changes leave
   every explicit seed count unchanged.
8. **Unrestricted support theorem:** all 120 order-four seven-path
   targets pass the path-union/hypergraph criterion. Every optional-arc
   extension above each extracted support is independently recounted,
   checking both directions of the criterion, not only original examples.

The mathematical proof is in [PROOFS.md](PROOFS.md), the general theorem
in [GENERAL_GRAMMAR.md](GENERAL_GRAMMAR.md), and review records in
[audit/AUDIT.md](audit/AUDIT.md) and
[final_review/REPORT.md](final_review/REPORT.md).

## Source and portability checks

The final review opened Moon's original Cambridge PDF and the Chan–Pak
arXiv text, confirming the exact cited statements and attribution.
Relative report/data paths are portable. The replay output records the
commands and exact-byte or semantic-JSON comparisons; the manifest binds
the final scientific packet. No temporary executables are published.
