# Source and search ledger

The repository's `05-knowledge/reference/CORE-PAPERS.md` and its analytic
number-theory sidecar were read before making literature claims.

## Verified standard vocabulary

1. Jean-Paul Allouche and Jeffrey Shallit, *The ring of k-regular sequences*,
   Theoretical Computer Science 98(2) (1992), 163--197.
   DOI: https://doi.org/10.1016/0304-3975(92)90001-V
   The publisher record identifies the regular-sequence framework as the
   unbounded-valued extension of automatic sequences. The local matrix proof
   establishes membership directly; no theorem about this particular g is
   attributed to that paper.

2. Arturo Carpi and Cristiano Maggi, *On synchronized sequences and their
   separators*, RAIRO Theoretical Informatics and Applications 35(6) (2001),
   513--524. https://www.numdam.org/item/ITA_2001__35_6_513_0/
   DOI: https://doi.org/10.1051/ita:2001129
   The archival record confirms that synchronized sequences are defined via
   graph recognition, and form a class strictly between automatic and regular
   sequences. This general distinction does not determine whether the present
   g is synchronized; the new explicit separating words do that locally.

3. Jeffrey Shallit, *Synchronized Sequences*, WORDS 2021, author-hosted slides.
   https://cs.uwaterloo.ca/~shallit/Talks/words21.pdf
   Slides 3--6 give linear-representation and padded-graph definitions. Slide
   12 explains first-order automatic constructions using addition and
   comparison. We use these to pin terminology, while giving an independent
   bounded-carry argument in the proof. This is a talk, not a new published
   representation theorem for our source sequence.

## Exact-sequence search limits

Queries used the initial g terms, the transformed terms (g(n)-1)/2, and the
odd-index subsequence, with and without an OEIS domain restriction. No exact
matching entry was identified in the returned results. Direct OEIS search
URLs returned inaccessible/internal-error results, so an exhaustive OEIS
negative is NOT claimed. A different PDF mirror returned HTTP403; it was not
used as a dependency.

Searches for weighted binary inversions and bilinear regular-sequence
representation supplied no applicable all-target theorem. This is only a
record of this bounded search. Neither the source sequence nor the local
formulas are claimed novel on that basis.

The formula in PROOF.md describes g as a weighted statistic of pairs of
binary digits. This is a locally proved description, not a claim that the
phrase names an established sequence.
