#!/usr/bin/env python3
"""crystallization_tools_s116n.py — Practical crystallization engines.

Crystallization: iterate a simple rule until structure emerges.
Not computing an answer. Not approximating an answer. GROWING the answer.
The iteration IS the computation. The fixed point IS the result.
"""
from math import log, sqrt, gcd
from collections import Counter, defaultdict

print()
print("  CRYSTALLIZATION ENGINES")
print()
print("="*70)
print()

# ============================================================
print("  ENGINE 1: PREFERENCE CRYSTALLIZER")
print("  " + "-"*40)
print()
print("  Problem: Given noisy, inconsistent pairwise preferences,")
print("  find the STABLE CORE — the preferences that survive iteration.")
print()
print("  Method: Kaprekar-like iteration on tournaments.")
print("  Step: for each triple (i,j,k), if it forms a 3-cycle,")
print("  FLIP the weakest arc (the one with fewest supporting votes).")
print("  Repeat until no more 3-cycles exist or a fixed point is reached.")
print()

def crystallize_preferences(votes, items):
    """Crystallize noisy preferences into a stable ranking.
    votes: list of (winner, loser) tuples (may have repeats)
    items: list of item names
    Returns: (ranking, H_value, iterations, cycles_removed)
    """
    n = len(items)
    idx = {item: i for i, item in enumerate(items)}

    # Count wins
    wins = [[0]*n for _ in range(n)]
    for w, l in votes:
        wins[idx[w]][idx[l]] += 1

    # Build tournament from majority
    T = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if wins[i][j] > wins[j][i]:
                T[i][j] = 1
            elif wins[j][i] > wins[i][j]:
                T[j][i] = 1
            else:
                T[i][j] = 1  # tiebreak by index

    # Iterate: find 3-cycles, flip weakest arc
    iterations = 0
    cycles_removed = 0
    max_iter = 100

    while iterations < max_iter:
        # Find all 3-cycles
        found_cycle = False
        weakest = None
        weakest_margin = float('inf')

        for i in range(n):
            for j in range(n):
                if j == i or T[i][j] == 0: continue
                for k in range(n):
                    if k == i or k == j: continue
                    if T[j][k] == 1 and T[k][i] == 1:
                        # 3-cycle: i->j->k->i
                        # Find weakest arc
                        arcs = [(i,j), (j,k), (k,i)]
                        for a, b in arcs:
                            margin = wins[a][b] - wins[b][a]
                            if margin < weakest_margin:
                                weakest_margin = margin
                                weakest = (a, b)
                        found_cycle = True

        if not found_cycle or weakest is None:
            break

        # Flip the weakest arc
        a, b = weakest
        T[a][b] = 0
        T[b][a] = 1
        cycles_removed += 1
        iterations += 1

    # Compute ranking from final tournament
    scores = [sum(T[i][j] for j in range(n)) for i in range(n)]
    ranking = sorted(range(n), key=lambda i: -scores[i])

    # Count remaining H (simplified: just check if transitive)
    remaining_cycles = 0
    for i in range(n):
        for j in range(n):
            if j == i or T[i][j] == 0: continue
            for k in range(n):
                if k == i or k == j: continue
                if T[j][k] == 1 and T[k][i] == 1:
                    remaining_cycles += 1
    remaining_cycles //= 3  # each cycle counted 3 times (once per vertex)

    return [items[i] for i in ranking], remaining_cycles, iterations, cycles_removed


# Demo
print("  DEMO: 5 restaurants rated by 10 people (noisy)")
items = ["Sushi", "Pizza", "Tacos", "Pasta", "Burger"]
# Simulate noisy votes
import random
random.seed(42)
# True ordering: Sushi > Pizza > Tacos > Pasta > Burger
# But with noise: ~30% of votes are flipped
true_order = {"Sushi": 5, "Pizza": 4, "Tacos": 3, "Pasta": 2, "Burger": 1}
votes = []
for _ in range(10):
    for i, a in enumerate(items):
        for j, b in enumerate(items):
            if i < j:
                if random.random() < 0.7:
                    # Correct vote
                    if true_order[a] > true_order[b]:
                        votes.append((a, b))
                    else:
                        votes.append((b, a))
                else:
                    # Flipped vote (noise)
                    if true_order[a] > true_order[b]:
                        votes.append((b, a))
                    else:
                        votes.append((a, b))

ranking, remaining, iters, removed = crystallize_preferences(votes, items)
print(f"  After {iters} iterations, {removed} weak arcs flipped:")
print(f"  Ranking: {' > '.join(ranking)}")
print(f"  Remaining cycles: {remaining}")
print(f"  True order: Sushi > Pizza > Tacos > Pasta > Burger")
print(f"  Match: {ranking == ['Sushi', 'Pizza', 'Tacos', 'Pasta', 'Burger']}")
print()

# ============================================================
print()
print("  ENGINE 2: CONSENSUS CRYSTALLIZER")
print("  " + "-"*40)
print()
print("  Problem: Multiple rankers provide orderings.")
print("  Find the consensus that CRYSTALLIZES from repeated merging.")
print()
print("  Method: Build a tournament from majority votes across all rankers.")
print("  Then crystallize (Engine 1). The fixed point IS the consensus.")
print()

def consensus_crystallize(rankings, items):
    """Crystallize multiple rankings into consensus.
    rankings: list of lists (each is a ranking of items, best first)
    """
    votes = []
    for ranking in rankings:
        for i, a in enumerate(ranking):
            for b in ranking[i+1:]:
                votes.append((a, b))

    return crystallize_preferences(votes, items)

# Demo: 5 experts rank 7 programming languages
items7 = ["Python", "Rust", "Go", "Java", "C++", "JS", "TS"]
rankings = [
    ["Python", "Rust", "Go", "TS", "JS", "C++", "Java"],
    ["Rust", "Python", "C++", "Go", "TS", "Java", "JS"],
    ["Python", "TS", "Rust", "JS", "Go", "Java", "C++"],
    ["Rust", "Go", "Python", "C++", "TS", "JS", "Java"],
    ["Python", "Rust", "TS", "Go", "JS", "C++", "Java"],
]

ranking, remaining, iters, removed = consensus_crystallize(rankings, items7)
print(f"  5 experts rank 7 languages.")
print(f"  Crystallization: {iters} iterations, {removed} flips.")
print(f"  Consensus: {' > '.join(ranking)}")
print(f"  Remaining disagreements: {remaining} cycles")
print()

# ============================================================
print()
print("  ENGINE 3: FEATURE IMPORTANCE CRYSTALLIZER")
print("  " + "-"*40)
print()
print("  Problem: Given pairwise A/B test results on features,")
print("  crystallize the importance ranking.")
print()
print("  The crystallization naturally handles:")
print("  - Feature A beats B on metric 1 but loses on metric 2")
print("  - These form genuine cycles (tradeoffs, not noise)")
print("  - The crystallizer identifies which arcs to flip (noise)")
print("    vs which cycles to preserve (genuine tradeoffs)")
print()
print("  Method: Run Engine 1 but STOP when all remaining cycles")
print("  have high margin (= genuine tradeoffs, not weak preferences).")
print()

def crystallize_with_tradeoffs(votes, items, margin_threshold=3):
    """Like Engine 1 but preserves strong cycles as genuine tradeoffs."""
    n = len(items)
    idx = {item: i for i, item in enumerate(items)}

    wins = [[0]*n for _ in range(n)]
    for w, l in votes:
        wins[idx[w]][idx[l]] += 1

    T = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if wins[i][j] >= wins[j][i]:
                T[i][j] = 1
            else:
                T[j][i] = 1

    iterations = 0
    noise_removed = 0
    tradeoffs = []

    while iterations < 100:
        weakest = None
        weakest_margin = float('inf')

        for i in range(n):
            for j in range(n):
                if j == i or T[i][j] == 0: continue
                for k in range(n):
                    if k == i or k == j: continue
                    if T[j][k] == 1 and T[k][i] == 1:
                        arcs = [(i,j), (j,k), (k,i)]
                        cycle_min_margin = min(wins[a][b] - wins[b][a] for a,b in arcs)
                        if cycle_min_margin >= margin_threshold:
                            # Genuine tradeoff — preserve it
                            cycle_items = (items[i], items[j], items[k])
                            if tuple(sorted(cycle_items)) not in [tuple(sorted(t)) for t in tradeoffs]:
                                tradeoffs.append(cycle_items)
                        elif cycle_min_margin < weakest_margin:
                            weakest_margin = cycle_min_margin
                            for a, b in arcs:
                                if wins[a][b] - wins[b][a] == cycle_min_margin:
                                    weakest = (a, b)
                                    break

        if weakest is None:
            break

        a, b = weakest
        T[a][b] = 0
        T[b][a] = 1
        noise_removed += 1
        iterations += 1

    scores = [sum(T[i][j] for j in range(n)) for i in range(n)]
    ranking = sorted(range(n), key=lambda i: -scores[i])

    return {
        'ranking': [items[i] for i in ranking],
        'noise_removed': noise_removed,
        'tradeoffs': tradeoffs,
        'iterations': iterations,
    }

# Demo
print("  DEMO: A/B test results on 5 UI features")
features = ["Dark Mode", "Search", "Filters", "Speed", "Layout"]
# Simulate: Speed and Layout form a genuine tradeoff
ab_votes = []
# Speed vs Layout: users are split (genuine tradeoff)
for _ in range(6): ab_votes.append(("Speed", "Layout"))
for _ in range(5): ab_votes.append(("Layout", "Speed"))
# Layout vs Dark Mode: users are split
for _ in range(5): ab_votes.append(("Layout", "Dark Mode"))
for _ in range(6): ab_votes.append(("Dark Mode", "Layout"))
# Dark Mode vs Speed: users prefer Dark Mode
for _ in range(4): ab_votes.append(("Dark Mode", "Speed"))
for _ in range(7): ab_votes.append(("Speed", "Dark Mode"))
# Other comparisons: clear winners
for _ in range(9): ab_votes.append(("Search", "Filters"))
for _ in range(8): ab_votes.append(("Search", "Dark Mode"))
for _ in range(7): ab_votes.append(("Search", "Speed"))
for _ in range(8): ab_votes.append(("Search", "Layout"))
for _ in range(7): ab_votes.append(("Speed", "Filters"))
for _ in range(6): ab_votes.append(("Dark Mode", "Filters"))
for _ in range(6): ab_votes.append(("Layout", "Filters"))

result = crystallize_with_tradeoffs(ab_votes, features, margin_threshold=2)
print(f"  Ranking: {' > '.join(result['ranking'])}")
print(f"  Noise removed: {result['noise_removed']} weak arcs flipped")
print(f"  Genuine tradeoffs found: {len(result['tradeoffs'])}")
for t in result['tradeoffs']:
    print(f"    {t[0]} <-> {t[1]} <-> {t[2]} (cycle preserved as real tradeoff)")
print()

# ============================================================
print()
print("  ENGINE 4: TEXT SUMMARIZATION BY CRYSTALLIZATION")
print("  " + "-"*40)
print()
print("  Problem: Given multiple summaries of the same text,")
print("  find the CRYSTALLIZED summary that all agree on.")
print()
print("  Method: Treat each sentence as an 'item'.")
print("  Each summary RANKS sentences by importance (inclusion order).")
print("  Crystallize the rankings to find the consensus importance.")
print()
print("  This is Engine 2 applied to sentences instead of items.")
print("  The crystallized ranking = the sentences in consensus order.")
print("  The top k sentences = the crystallized summary.")
print()
print("  Not implemented here (needs NLP), but the ENGINE is Engine 2.")
print("  The crystallization principle: iterate majority-merge until stable.")
print()

# ============================================================
print()
print("  ENGINE 5: TEAM FORMATION CRYSTALLIZER")
print("  " + "-"*40)
print()
print("  Problem: N people, form teams of 7.")
print("  Each person has pairwise compatibility scores.")
print("  Find teams that CRYSTALLIZE (stable — no one wants to switch).")
print()
print("  Method:")
print("  1. Random initial teams of 7.")
print("  2. Within each team: run Engine 1 to find the team's 'chemistry'")
print("     (H-value = how many compatible orderings exist).")
print("  3. If a team has high H (many cycles = many conflicts):")
print("     find the person involved in the most cycles.")
print("     SWAP them with someone from a low-H team.")
print("  4. Repeat until all teams have acceptable H.")
print()
print("  The crystallization: teams SELF-ORGANIZE through iterated swaps.")
print("  The fixed point: teams where no swap improves any team's H.")
print("  This is a STABLE MATCHING — but via crystallization, not algorithm.")
print()

# ============================================================
print()
print("  ENGINE 6: ANOMALY DETECTOR BY CRYSTALLIZATION")
print("  " + "-"*40)
print()
print("  Problem: Detect fraudulent or anomalous data in comparisons.")
print()
print("  Method: Crystallize the data (Engine 1).")
print("  The arcs that get FLIPPED during crystallization are the anomalies.")
print("  Arcs flipped early (low margin) = noise.")
print("  Arcs flipped late (high margin) = genuine anomalies.")
print("  Arcs NEVER flipped = consistent data.")
print()
print("  The crystallization process SORTS the data by reliability.")
print("  The first flips = least reliable comparisons.")
print("  The last flips = most surprising genuine disagreements.")
print("  If an arc survives crystallization AND contradicts the final ranking:")
print("  it's a GENUINE ANOMALY — structurally preserved, not noise.")
print()

# Demo
print("  DEMO: Detecting a rigged comparison")
items_demo = ["A", "B", "C", "D", "E"]
votes_demo = []
# Normal data: A > B > C > D > E with some noise
for _ in range(10): votes_demo.append(("A", "B"))
for _ in range(8): votes_demo.append(("B", "C"))
for _ in range(9): votes_demo.append(("C", "D"))
for _ in range(10): votes_demo.append(("D", "E"))
for _ in range(7): votes_demo.append(("A", "C"))
for _ in range(8): votes_demo.append(("B", "D"))
for _ in range(9): votes_demo.append(("C", "E"))
for _ in range(6): votes_demo.append(("A", "D"))
for _ in range(7): votes_demo.append(("B", "E"))
for _ in range(8): votes_demo.append(("A", "E"))
# RIGGED: someone injected E > A votes
for _ in range(12): votes_demo.append(("E", "A"))

ranking, remaining, iters, removed = crystallize_preferences(votes_demo, items_demo)
print(f"  Data includes 12 rigged 'E beats A' votes.")
print(f"  Crystallization: {iters} iterations, {removed} flips.")
print(f"  Final ranking: {' > '.join(ranking)}")
print(f"  Remaining cycles: {remaining}")
if remaining > 0:
    print(f"  The surviving cycle likely involves the rigged arc E > A.")
    print(f"  Anomaly detected: E > A has high margin but contradicts everything else.")
print()

# ============================================================
print()
print("  THE PRINCIPLE")
print("  " + "-"*40)
print()
print("  Every engine uses the same principle:")
print()
print("  1. START with raw, noisy, inconsistent data.")
print("  2. ITERATE a simple local rule (flip weakest contradictions).")
print("  3. STOP when no more weak contradictions exist.")
print("  4. The FIXED POINT is the answer.")
print()
print("  The fixed point is not computed from a formula.")
print("  It is not approximated by optimization.")
print("  It GROWS from the iteration.")
print("  It IS the crystal.")
print()
print("  The crystal preserves GENUINE structure (strong cycles = real tradeoffs)")
print("  and removes NOISE (weak arcs that create spurious cycles).")
print("  This separation of signal from noise happens AUTOMATICALLY")
print("  through the iteration. No threshold needed. No training needed.")
print("  The crystallization process IS the signal-noise separator.")
