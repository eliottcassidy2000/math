#!/usr/bin/env python3
"""partial_bit_applications_s116g.py — Practical applications of partial bit accumulation.

The principle: structure from accumulating sub-bit distinctions.
The mechanism: telescoping products, constrained digit patterns, forbidden totals.
Now: what can we BUILD with this?
"""
from math import sqrt, log, exp, atanh, tanh, log2, factorial, comb
from fractions import Fraction
from itertools import combinations

print()
print("  PRACTICAL APPLICATIONS OF PARTIAL BIT ACCUMULATION")
print()
print("="*70)
print()

# ============================================================
print("  1. RANKING CONFIDENCE FROM PARTIAL COMPARISONS")
print("  " + "-"*50)
print()
print("  PROBLEM: You have n items and k < C(n,2) pairwise comparisons.")
print("  Not all pairs have been compared. How confident is the ranking?")
print()
print("  INSIGHT: Each comparison is a PARTIAL BIT of ranking information.")
print("  A consistent comparison (A>B and B>C implies A>C) adds ~1 bit.")
print("  An inconsistent comparison (A>B>C>A cycle) adds 0 net bits")
print("  but reveals STRUCTURE (the cycle).")
print()
print("  TOOL: Partial Bit Ranking Confidence Score")
print()

def ranking_confidence(comparisons, n):
    """Compute ranking confidence from partial pairwise comparisons.
    comparisons: list of (winner, loser) tuples.
    n: number of items.
    Returns: confidence score in [0, 1] and diagnostic info.
    """
    # Build adjacency
    adj = [[0]*n for _ in range(n)]
    compared = set()
    for w, l in comparisons:
        adj[w][l] = 1
        compared.add((min(w,l), max(w,l)))

    total_pairs = n*(n-1)//2
    coverage = len(compared) / total_pairs

    # Count 3-cycles (contradictions)
    cycles = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (adj[i][j] and adj[j][k] and adj[k][i]) or \
                   (adj[i][k] and adj[k][j] and adj[j][i]):
                    cycles += 1

    # Count transitive triples (confirmations)
    trans = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (min(i,j), max(i,j)) in compared and \
                   (min(j,k), max(j,k)) in compared and \
                   (min(i,k), max(i,k)) in compared:
                    # Check if transitive
                    triple = [i, j, k]
                    wins = [sum(adj[x][y] for y in triple if y != x) for x in triple]
                    if sorted(wins) == [0, 1, 2]:
                        trans += 1

    # Each transitive triple contributes ~1 partial bit
    # Each cycle contributes 0 partial bits (but costs information)
    total_triples = comb(n, 3)
    compared_triples = trans + cycles
    if compared_triples > 0:
        consistency = trans / compared_triples
    else:
        consistency = 1.0

    # Confidence = coverage * consistency
    confidence = coverage * consistency

    return {
        'coverage': coverage,
        'compared': len(compared),
        'total_pairs': total_pairs,
        'cycles': cycles,
        'transitive': trans,
        'consistency': consistency,
        'confidence': confidence,
        'partial_bits': trans * 1.0 + cycles * 0.0,
        'max_bits': log2(factorial(n)) if n <= 20 else n * log2(n),
    }

# Demo with 5 items
print("  DEMO: 5 items, partial comparisons")
print()

# Scenario 1: Complete and consistent
comps_perfect = [(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]
r = ranking_confidence(comps_perfect, 5)
print(f"  Scenario A: Perfect hierarchy (0>1>2>3>4, all compared)")
print(f"    Coverage: {r['coverage']:.0%}, Cycles: {r['cycles']}, "
      f"Consistency: {r['consistency']:.0%}, Confidence: {r['confidence']:.0%}")
print(f"    Partial bits: {r['partial_bits']:.1f} / {r['max_bits']:.1f} max")
print()

# Scenario 2: Incomplete but consistent
comps_partial = [(0,1),(1,2),(2,3),(3,4),(0,2),(1,3)]
r = ranking_confidence(comps_partial, 5)
print(f"  Scenario B: Partial, consistent (chain + shortcuts)")
print(f"    Coverage: {r['coverage']:.0%}, Cycles: {r['cycles']}, "
      f"Consistency: {r['consistency']:.0%}, Confidence: {r['confidence']:.0%}")
print(f"    Partial bits: {r['partial_bits']:.1f} / {r['max_bits']:.1f} max")
print()

# Scenario 3: Complete with one cycle
comps_cycle = [(0,1),(1,2),(2,0),(0,3),(0,4),(1,3),(1,4),(2,3),(2,4),(3,4)]
r = ranking_confidence(comps_cycle, 5)
print(f"  Scenario C: Complete with one 3-cycle (0>1>2>0)")
print(f"    Coverage: {r['coverage']:.0%}, Cycles: {r['cycles']}, "
      f"Consistency: {r['consistency']:.0%}, Confidence: {r['confidence']:.0%}")
print(f"    Partial bits: {r['partial_bits']:.1f} / {r['max_bits']:.1f} max")
print()

# Scenario 4: Minimal comparisons
comps_minimal = [(0,1),(1,2),(2,3),(3,4)]
r = ranking_confidence(comps_minimal, 5)
print(f"  Scenario D: Minimal chain (4 comparisons only)")
print(f"    Coverage: {r['coverage']:.0%}, Cycles: {r['cycles']}, "
      f"Consistency: {r['consistency']:.0%}, Confidence: {r['confidence']:.0%}")
print(f"    Partial bits: {r['partial_bits']:.1f} / {r['max_bits']:.1f} max")
print()

print("  The tool tells you:")
print("  - HOW MUCH of the ranking is determined (coverage)")
print("  - HOW CONSISTENT the data is (no cycles = high consistency)")
print("  - HOW CONFIDENT the ranking is (coverage * consistency)")
print("  - HOW MANY 'partial bits' of ranking information you have")
print()

# ============================================================
print()
print("  2. ADAPTIVE COMPARISON SCHEDULER")
print("  " + "-"*50)
print()
print("  PROBLEM: You need to rank n items with minimum comparisons.")
print("  Which pair should you compare NEXT to gain the most information?")
print()
print("  INSIGHT: Each comparison contributes a partial bit.")
print("  The partial bit is LARGEST when the comparison is MAXIMALLY")
print("  UNCERTAIN — when we don't know who will win.")
print()
print("  ALGORITHM:")
print("  1. Maintain current ranking estimate (from existing comparisons)")
print("  2. For each uncompared pair (i,j):")
print("     - Estimate P(i>j) from transitive closure")
print("     - Compute expected information gain = H(P) = binary entropy")
print("     - This is the 'partial bit value' of the comparison")
print("  3. Compare the pair with HIGHEST expected partial bit value")
print("  4. Update ranking, repeat")
print()
print("  The information gain per comparison is maximized when")
print("  P(i>j) = 1/2 (fair coin), which gives 1 full bit.")
print("  But transitive closure often makes some pairs predictable")
print("  (P near 0 or 1), giving < 1 bit. The scheduler avoids these.")
print()

def binary_entropy(p):
    if p <= 0 or p >= 1: return 0
    return -p*log2(p) - (1-p)*log2(1-p)

print("  Expected information gain per comparison:")
print("  P(i>j)   H(P) bits   Rapidity    Interpretation")
print("  " + "-"*55)
for p in [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]:
    h = binary_entropy(p)
    rap = abs(atanh(2*p-1))
    interp = ""
    if p == 0.5: interp = "MAXIMUM INFO (fair coin)"
    elif rap < 0.2: interp = "near-uncertain, high info"
    elif rap < 1.0: interp = "moderate certainty"
    else: interp = "near-certain, low info"
    print(f"  {p:.2f}       {h:.4f}       {rap:.4f}      {interp}")

print()
print("  Comparisons near the rapidity origin (P~0.5) are most valuable.")
print("  Comparisons far from origin (P~0 or ~1) are nearly redundant.")
print("  The scheduler PRIORITIZES pairs near the center of rapidity space.")
print()

# ============================================================
print()
print("  3. CONTRADICTION DETECTOR FOR SURVEY DATA")
print("  " + "-"*50)
print()
print("  PROBLEM: A survey asks respondents to compare items pairwise.")
print("  Some responses contradict each other (A>B>C>A).")
print("  Which contradictions are 'real' (population-level cycles)")
print("  vs 'noise' (individual inconsistency)?")
print()
print("  INSIGHT: The OCF decomposes H into independent cycle contributions.")
print("  alpha_1 = total cycles, alpha_2 = independent pairs, etc.")
print("  If alpha_1 is high but alpha_2 is low, the cycles are NOT independent")
print("  — they share vertices — suggesting a single confused region.")
print("  If alpha_2 is also high, there are MULTIPLE independent")
print("  contradictions — suggesting genuine population-level disagreement.")
print()
print("  TOOL OUTPUT:")
print("  'Your survey has {alpha_1} pairwise contradictions.'")
print("  '{alpha_2} of these are independent (non-overlapping).'")
print("  'The contradictions cluster around items {X, Y, Z}.'")
print("  'Confidence in overall ranking: {confidence}%.'")
print()

# ============================================================
print()
print("  4. OPTIMAL EXPERIMENT DESIGN VIA PARTIAL BITS")
print("  " + "-"*50)
print()
print("  PROBLEM: Design a clinical trial comparing n treatments pairwise.")
print("  Budget allows k comparisons. Which k pairs to test?")
print()
print("  INSIGHT: The total information from k comparisons is at most k bits.")
print("  But some comparison sets are MORE EFFICIENT because their")
print("  partial bits DON'T OVERLAP (independent information).")
print()
print("  STRATEGY: Choose comparisons whose transitive closures are DISJOINT.")
print("  This maximizes the number of independent partial bits.")
print()
print("  For n items, the minimum comparisons for a complete ranking = n-1")
print("  (a spanning path). But this gives NO redundancy (no cycle detection).")
print("  Adding 1 more comparison creates exactly 1 cycle-detection opportunity.")
print()
print("  The OPTIMAL k for n items:")
print("  k = n-1: minimum for ranking, no contradiction detection")
print("  k = n: one redundant comparison, detects 1 cycle")
print("  k = C(n,2): all pairs, maximum information, maximum cost")
print()
print("  n    n-1    C(n,2)   bits needed    ratio")
print("  " + "-"*50)
for n in [3, 4, 5, 6, 7, 8, 10, 15, 20]:
    pairs = n*(n-1)//2
    bits = log2(factorial(n)) if n <= 20 else n*log2(n) - n/log(2)
    print(f"  {n:3d}   {n-1:3d}    {pairs:5d}      {bits:7.1f}         {(n-1)/bits:.2f}")

print()
print("  The ratio (n-1)/bits_needed -> 0 as n grows.")
print("  You need MANY more comparisons than the minimum to get")
print("  reliable information. The partial bits from a spanning tree")
print("  are necessary but not sufficient.")
print()

# ============================================================
print()
print("  5. RECOMMENDATION SYSTEM: PREFERENCE TOPOLOGY")
print("  " + "-"*50)
print()
print("  PROBLEM: A user rates items pairwise (A vs B, B vs C, etc.).")
print("  Build a recommendation based on the TOPOLOGY of preferences.")
print()
print("  INSIGHT: The user's preference tournament T has H(T) paths.")
print("  High H: the user has WEAK preferences (many valid orderings).")
print("  Low H: the user has STRONG preferences (few valid orderings).")
print()
print("  The OCF gives alpha_1 (how many items form cycles).")
print("  Items in cycles are ones the user is AMBIVALENT about.")
print("  Items NOT in any cycle are DEFINITELY preferred in their ranking.")
print()
print("  RECOMMENDATION:")
print("  - Recommend items that are DEFINITELY above unrated items")
print("    (transitive closure from the known comparisons)")
print("  - Flag items in CYCLES as 'you might also like'")
print("    (the user is ambivalent, so both could work)")
print("  - Compute confidence per item-pair: high for transitive,")
print("    low for cycle-adjacent, unknown for uncompared")
print()

# ============================================================
print()
print("  6. SPORTS: SCHEDULE OPTIMIZER")
print("  " + "-"*50)
print()
print("  PROBLEM: A league has n teams. Design a schedule that")
print("  determines the TRUE ranking in minimum games.")
print()
print("  INSIGHT: Each game is a partial bit. Some games are more")
print("  informative than others (evenly matched teams give ~1 bit,")
print("  mismatched teams give ~0 bits).")
print()
print("  ALGORITHM:")
print("  1. Start with estimated skill ratings (Elo, preseason polls)")
print("  2. Schedule games between teams with CLOSEST ratings")
print("     (maximum partial bit per game)")
print("  3. After each round, update ratings and reschedule")
print("  4. Stop when confidence exceeds threshold")
print()
print("  This is ADAPTIVE SCHEDULING: the schedule changes based on results.")
print("  Traditional round-robin wastes games on mismatched pairs.")
print("  Adaptive scheduling focuses on the INFORMATIVE matchups.")
print()
print("  Expected savings:")
print("  n    Round-robin   Adaptive (est.)   Savings")
print("  " + "-"*50)
for n in [4, 8, 10, 16, 20, 32]:
    rr = n*(n-1)//2
    # Adaptive needs ~n*log2(n) games (information-theoretic minimum)
    adaptive = int(n * log2(n) * 1.5)  # with 50% overhead for noise
    savings = 1 - adaptive/rr
    print(f"  {n:3d}    {rr:5d}          {adaptive:5d}           {savings:.0%}")

print()
print("  For a 32-team league: 496 round-robin games vs ~240 adaptive.")
print("  More than 50% savings while determining the same ranking.")
print()

# ============================================================
print()
print("  7. ML MODEL SELECTION VIA PARTIAL BITS")
print("  " + "-"*50)
print()
print("  PROBLEM: Compare n ML models pairwise on k datasets.")
print("  Which model is best? How confident are we?")
print()
print("  INSIGHT: Each dataset comparison is a partial bit.")
print("  Model A beats Model B on dataset d: one bit.")
print("  If A beats B on MOST datasets: high rapidity.")
print("  If results are mixed: low rapidity (uncertain).")
print()
print("  TOOL: Model Comparison Rapidity")
print("  For each pair (A, B): count wins w, losses l.")
print("  Rapidity = arctanh((w-l)/(w+l)) = arctanh(2*(w/(w+l)) - 1).")
print("  = half the log-odds of A beating B.")
print()
print("  w/total   rapidity    confidence")
print("  " + "-"*40)
for w_frac in [0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0]:
    if w_frac < 1.0:
        rap = atanh(2*w_frac - 1)
    else:
        rap = float('inf')
    conf = "uncertain" if abs(2*w_frac-1) < 0.2 else "moderate" if abs(2*w_frac-1) < 0.6 else "strong" if abs(2*w_frac-1) < 0.9 else "decisive"
    print(f"  {w_frac:.2f}       {rap if rap < 100 else 'inf':>8}    {conf}")

print()

# ============================================================
print()
print("  SUMMARY: FIVE BUILDABLE TOOLS")
print("  " + "-"*50)
print()
print("  1. RANKING CONFIDENCE SCORE")
print("     Input: partial pairwise comparisons")
print("     Output: confidence %, partial bits accumulated, cycle list")
print("     Use: product managers, hiring, any ranking task")
print()
print("  2. ADAPTIVE COMPARISON SCHEDULER")
print("     Input: items to rank, comparison budget")
print("     Output: which pairs to compare next (max info per comparison)")
print("     Use: A/B testing, clinical trials, tournament scheduling")
print()
print("  3. CONTRADICTION DETECTOR")
print("     Input: survey pairwise preference data")
print("     Output: cycle identification, independence analysis, confidence")
print("     Use: market research, UX research, political polling")
print()
print("  4. SPORTS SCHEDULE OPTIMIZER")
print("     Input: teams + initial ratings + game budget")
print("     Output: adaptive schedule maximizing ranking information")
print("     Use: league organizers, esports, debate tournaments")
print()
print("  5. ML MODEL COMPARATOR")
print("     Input: model performance on datasets (pairwise wins/losses)")
print("     Output: rapidity-based confidence ranking of models")
print("     Use: ML practitioners, AutoML pipelines, model registries")
print()
print("  ALL FIVE use the same math: partial bit accumulation")
print("  from pairwise comparisons, measured in rapidity.")
