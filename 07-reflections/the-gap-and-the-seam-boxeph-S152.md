# The Gap and the Seam: {7,21} across five structures

**boxeph-2026-07-20-S152** · companion: THM-1370 (7 and 21 are never h-values, all n)

THM-1370 closed the owner's {7,21} question with a floor argument: the strong
floor f(n) (Moon-Busch, re-derived exhaustively to n=8, witnessed at n=9)
strands 7 in (5,9) and 21 in (15,25), and factorization recurses into the
gaps. That is the WHY. This note records where the same two numbers surface
around the proof — resonances, each with its status labeled.

1. **Forbidden as values, realized as symmetry (proved).** At n=7 the Paley
   tournament QR7 has |Aut| = 21 = p(p-1)/2 = #arcs (arc-regular), and the
   free-action lemma (S151: |Aut| divides h) forces 21 | h(QR7) = 189. So 21
   cannot BE a path count but MUST divide the top one at n=7. The pair {7,21}
   = (vertices, arcs) of QR7: forbidden as counts, incarnate as structure.
2. **The metagraph seam (in-repo, cross-thread).** The even-graph metagraph
   E_n is chordal/perfect for n <= 6 and breaks at n=7 (first odd holes,
   C_5/C_4) — and the repo already tied that seam to the forbidden pair via
   H = 7 = I(K_3, 2) (THM-200, even-graphs-as-first-class). The h-gap
   permanence (THM-1370) upgrades one side of the seam to a theorem; whether
   the chordality break is a genuine shadow of the floor crevasse or a
   coincidence of n=7 is OPEN and now sharply posed.
3. **Fifth-term breaks, two more instances (verified this session).** The
   Leonardo law 2F(n)-1 matched the strong floor for exactly five terms
   (3,5,9,15,25) and died at n=8 (41 vs 45). The even-Betti vanishing
   conjecture died at n=8 (beta_4 in {1,5}). With Mersenne-Moser (29),
   companion-Pell (97 vs 99), and the 7*3^k tower (63 attained at n=8), the
   motif "two laws agree on ~5 terms then split" is now a documented failure
   MODE — treat any 5-term match as a hypothesis, never evidence.
4. **The methodological key (reusable).** Augmentation covers: every
   (n+1)-tournament minus a vertex is an n-tournament, so canonical n-reps x
   all patterns cover every (n+1)-class at cost |classes(n)| * 2^n — the full
   n=8 h-spectrum cost 58,368 DP calls, no 2^28 enumeration. This works for
   ANY isomorphism-invariant spectrum: beta-spectra of GLMY homology, E_n
   invariants, anti-order spectra, SC counts. Recommended as the default
   exhaustion tool through n=9 (6880 x 512 = 3.5M evaluations).
5. **Single-flip identification (proved).** The almost-transitive family R_n
   (candidate floor witness, refuted) has h = 2^{n-2}+1 — exactly THM-250's
   single-flip apex value, the hypotenuse constant H = 1+2^d. The floor
   family is instead the Moon-Busch 3/5-block chain: the floor is FIBONACCI-
   LIKE only by accident of five terms; it is really the knapsack law
   min{3^a 5^b : 2a+3b = n-1} — blocks priced in vertices, values in paths.
