#!/usr/bin/env python3
"""unexpected_applications_s115.py — Applications nobody has thought of"""
from math import sqrt, log, factorial

print("UNEXPECTED APPLICATIONS OF THE CAYLEY-DELANNOY THEORY")
print("="*60)

# ============================================================
print("\n1. HIRING: IS YOUR INTERVIEW PANEL CONSISTENT?")
print("-"*40)
print("""
5 candidates interviewed by a panel. Each interviewer ranks pairs.
The TOURNAMENT of aggregated preferences has H(T) Hamiltonian paths.

If H ~ E[H] = 7.5: the panel's rankings are essentially random.
  -> Your hiring process is NOT discriminating between candidates.
  -> Either the candidates are truly equivalent, or your evaluation
     criteria are inconsistent.

If H = 1: perfect agreement on a single ordering.
  -> Either genuine skill differences, or groupthink.

ACTION: compute H(T) after each round of interviews.
If H is not decreasing toward 1 as you add data, your process
is not converging on a ranking.

Cost: ZERO (just count pairwise preferences you already have).
""")

# ============================================================
print("2. TASTE TESTING: DOES YOUR WINE RANKING MEAN ANYTHING?")
print("-"*40)
print("""
n wines tasted blind, compared pairwise. Tournament T.
H(T)/E[H] measures whether the tasting detected real quality differences.

For n=6 wines: E[H]=22.5, CV=0.577.
A panel of 10 tasters aggregated: if H < 5, the ranking is REAL.
If H ~ 22, the tasters are guessing.

BONUS: the Fourier decomposition tells you WHERE the disagreement is:
  k=1: pairwise preferences (do tasters agree on A vs B?)
  k=2: transitivity (if A>B and B>C, does A>C hold?)
  k=3: 3-way cycles (Condorcet paradoxes in taste)

A wine ranking with strong k=1 but weak k=2 means: tasters have
consistent PAIRWISE preferences but they DON'T compose transitively.
This is the signature of a multi-dimensional quality space
(e.g., one wine is best on aroma, another on taste, another on finish).
""")

# ============================================================
print("3. COMPILER OPTIMIZATION: WHICH PASS ORDER IS BEST?")
print("-"*40)
print("""
n compiler optimization passes (dead code elimination, inlining,
constant folding, loop unrolling, ...). Order matters!

Run all n! orderings on a benchmark. For each pair of passes (i,j):
  T[i][j] = 1 if pass i before j tends to produce faster code.

H(T) = number of orderings consistent with all pairwise preferences.
If H = 1: there is ONE optimal ordering (do it!).
If H ~ E[H]: pass ordering doesn't matter much (parallelize freely).
If H is small but > 1: a few orderings dominate (search among them).

The Z-score tells you: is optimizing pass order WORTH THE EFFORT?
For n=10 passes: if Z > 2, YES. If Z ~ 0, pass order is noise.
""")

# ============================================================
print("4. SPORTS DRAFT: IS THE DRAFT ORDER DEFENSIBLE?")
print("-"*40)
print("""
After a sports season, teams are ranked for the draft.
The ranking is based on win/loss record, but within tiebreakers
there are pairwise head-to-head results.

Compute H(T) for the head-to-head tournament.
If H is low: the ranking is well-determined by the data.
If H is high: the tiebreakers are essentially arbitrary.

LEGAL IMPLICATIONS: A team challenging their draft position could
use the Z-score to argue that the ranking is not statistically
supported by the head-to-head results.
""")

# ============================================================
print("5. MEDICAL DIAGNOSIS: RANKING TREATMENT OPTIONS")
print("-"*40)
print("""
n treatments compared pairwise in clinical trials.
Treatment A is better than B if it has higher efficacy in head-to-head.

H(T) measures the consistency of the treatment hierarchy.
A Z-score < -2 means: the treatments have a CLEAR ranking.
A Z-score ~ 0 means: the treatment comparisons are NOISY.

This could complement meta-analysis by providing a single
number quantifying 'how rankable are these treatments?'

For n=8 treatments: CV = 0.50.
If the clinical data gives H < 50 (out of E[H] = 315):
  the treatment ranking is statistically supported.
""")

# ============================================================
print("6. LLM EVALUATION: IS CHATGPT REALLY BETTER THAN CLAUDE?")
print("-"*40)
print("""
n language models evaluated by pairwise human preference.
This IS a tournament! (Chatbot Arena, LMSYS, etc.)

Current practice: Elo ratings from pairwise comparisons.
OUR ADDITION: compute H(T) and the Z-score.

If Z << 0: the models have a CLEAR hierarchy.
  Elo ratings are meaningful. Rankings are real.
If Z ~ 0: the models are essentially EQUIVALENT.
  Elo differences are noise. Don't trust the ranking.

For n=10 models on the Chatbot Arena:
  E[H] ~ 7088, CV = 0.447.
  If H < 1000: the ranking is VERY significant (Z < -4).

This is a NEW TOOL for the LLM evaluation community.
Nobody is currently reporting H(T) for chatbot tournaments.
""")

# ============================================================
print("7. DEMOCRACY: QUANTIFYING ELECTION QUALITY")
print("-"*40)
print("""
n candidates in a ranked-choice election.
The pairwise Condorcet matrix is a tournament.

H(T) = number of linear orderings consistent with pairwise margins.
Z-score = significance of the electoral ranking.

A WELL-FUNCTIONING DEMOCRACY should have low H (clear preferences).
A DYSFUNCTIONAL ELECTION has H ~ E[H] (preferences are noise).

For n=5 candidates:
  E[H] = 7.5, CV = 0.632.
  If the Condorcet tournament gives H = 1: unambiguous winner.
  If H = 15: the electorate is maximally split.

This gives a QUANTITATIVE MEASURE of electoral clarity.
Not just 'who won?' but 'how clear was the mandate?'
""")

# ============================================================
print("8. SUPPLY CHAIN: VENDOR RANKING SIGNIFICANCE")
print("-"*40)
print("""
n suppliers evaluated on multiple criteria.
For each pair: which supplier is better overall?

Procurement teams often rank suppliers subjectively.
The Cayley-Delannoy test gives an OBJECTIVE measure:
  'Is this ranking statistically distinguishable from random?'

If NO (Z ~ 0): the suppliers are interchangeable.
  -> Choose the cheapest or most convenient.
If YES (Z < -2): there are real quality differences.
  -> The ranking matters. Pay attention to it.

Cost: ZERO (just use the pairwise comparisons you already made).
""")

print("="*60)
print("SUMMARY: THE CAYLEY-DELANNOY TEST IS USEFUL ANYWHERE")
print("YOU HAVE PAIRWISE COMPARISONS AND WANT TO KNOW IF THE")
print("RESULTING RANKING IS REAL OR NOISE.")
print()
print("Formula: Z = (H - n!/2^{n-1}) / (n!/2^{n-1} * sqrt(2/n))")
print("That's it. One number. Instant answer.")
