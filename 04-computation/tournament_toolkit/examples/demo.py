#!/usr/bin/env python3
"""
Complete demo of the tournament-toolkit package.
Shows all four tools working together on practical scenarios.
"""

from tournament_toolkit import FormalRank, CycleDetector, CartanProbe, SpectralAnalyzer
import numpy as np

print("=" * 60)
print("  TOURNAMENT TOOLKIT — COMPLETE DEMO")
print("=" * 60)
print()

# ================================================================
# DEMO 1: LLM Ranking (FormalRank)
# ================================================================

print("DEMO 1: LLM CHATBOT ARENA RANKING")
print("-" * 50)
print()

ranker = FormalRank()

# Simulated arena data
battles = [
    ("Claude-3.5", "GPT-4o", 156, 134),
    ("Claude-3.5", "Gemini-1.5", 178, 112),
    ("Claude-3.5", "Llama-3.1", 201, 89),
    ("GPT-4o", "Gemini-1.5", 162, 128),
    ("GPT-4o", "Llama-3.1", 195, 95),
    ("Gemini-1.5", "Llama-3.1", 167, 123),
    # The twist: Llama beats Claude in code!
    ("Llama-3.1", "Claude-3.5", 45, 35),
]

for a, b, w, l in battles:
    ranker.add(a, b, wins=w, losses=l)

print(ranker.summary())
print()

# Check ambiguity
H = ranker.ambiguity()
print(f"  Tournament ambiguity H = {H}")
if H and H % 2 == 1:
    print(f"  (Redei check: H is odd ✓)")
print()

# Most informative next comparison
pair = ranker.most_informative_comparison()
print(f"  Most informative next matchup: {pair[0]} vs {pair[1]}")
print()

# ================================================================
# DEMO 2: Fraud Detection (CycleDetector)
# ================================================================

print("DEMO 2: WASH TRADING DETECTION")
print("-" * 50)
print()

detector = CycleDetector()

# Normal trading
for _ in range(20):
    detector.observe("Hedge_Fund_A", "Market_Maker_1")
    detector.observe("Market_Maker_1", "Retail_Investor")
    detector.observe("Retail_Investor", "Hedge_Fund_B")

# Wash trading ring (circular)
for _ in range(15):
    found = detector.observe("Shell_Corp_X", "Shell_Corp_Y")
    found = detector.observe("Shell_Corp_Y", "Shell_Corp_Z")
    found = detector.observe("Shell_Corp_Z", "Shell_Corp_X")  # Cycle!

print(detector.summary())
print()

# ================================================================
# DEMO 3: AI Confidence (CartanProbe)
# ================================================================

print("DEMO 3: AI ATTENTION CONFIDENCE DIAGNOSTIC")
print("-" * 50)
print()

probe = CartanProbe()

# Confident attention (near-transitive)
M_conf = np.array([
    [0, 0.8, 0.9, 0.95],
    [0.2, 0, 0.7, 0.85],
    [0.1, 0.3, 0, 0.75],
    [0.05, 0.15, 0.25, 0],
])

# Uncertain attention (cyclic)
M_unc = np.array([
    [0, 0.6, 0.3, 0.7],
    [0.4, 0, 0.6, 0.3],
    [0.7, 0.4, 0, 0.6],
    [0.3, 0.7, 0.4, 0],
])

for name, M in [("Confident head", M_conf), ("Uncertain head", M_unc)]:
    r = probe.analyze(M)
    print(f"  {name}:")
    print(f"    Verdict: {r['verdict']}")
    print(f"    Tournament fraction: {r['tournament_fraction']:.3f}")
    print(f"    Cooperation fraction: {r['cooperation_fraction']:.3f}")
    print(f"    Entanglement: {r['entanglement']:.3f}")
    print(f"    Combined risk: {r['combined_risk']:.3f}")
    print()

# Compare how the Cartan balance shifts
shift = probe.compare(M_conf, M_unc)
print(f"  Shift from confident → uncertain:")
print(f"    Tournament: {shift['tournament_shift']:+.3f}")
print(f"    Cooperation: {shift['cooperation_shift']:+.3f}")
print(f"    Verdict: {shift['verdict_change']}")
print()

# ================================================================
# DEMO 4: Tournament Quality (SpectralAnalyzer)
# ================================================================

print("DEMO 4: TOURNAMENT QUALITY ANALYSIS")
print("-" * 50)
print()

analyzer = SpectralAnalyzer()

# Paley tournament T_7
qr = {1, 2, 4}
adj_paley = [[0]*7 for _ in range(7)]
for i in range(7):
    for j in range(7):
        if i != j and (j - i) % 7 in qr:
            adj_paley[i][j] = 1

# Transitive tournament
adj_trans = [[0]*7 for _ in range(7)]
for i in range(7):
    for j in range(i+1, 7):
        adj_trans[i][j] = 1

for name, adj in [("Paley T_7 (optimal)", adj_paley),
                   ("Transitive (minimal)", adj_trans)]:
    print(f"  {name}:")
    print(analyzer.quality_report(adj))
    print()

# ================================================================
# DEMO 5: Practical A/B Test
# ================================================================

print("DEMO 5: A/B TEST INTRANSITIVITY CHECK")
print("-" * 50)
print()

ranker2 = FormalRank()
detector2 = CycleDetector()

# Three variants being tested
results = [
    ("Version_A", "Version_B", 52, 48),   # A slightly beats B
    ("Version_B", "Version_C", 55, 45),   # B clearly beats C
    ("Version_C", "Version_A", 51, 49),   # C slightly beats A — cycle!
]

for a, b, w, l in results:
    ranker2.add(a, b, wins=w, losses=l)
    for _ in range(w):
        detector2.observe(a, b)
    for _ in range(l):
        detector2.observe(b, a)

print(ranker2.summary())
print()
print(f"  Cycle detection risk: {detector2.risk():.1%}")
if detector2.n_cycles > 0:
    print(f"  ⚠ INTRANSITIVITY DETECTED: {detector2.n_cycles} cycles found")
    print(f"  The test results are CONTRADICTORY — A>B>C>A")
    print(f"  Recommendation: Run more comparisons or test for confounding variables")
else:
    print(f"  ✓ No intransitivity detected")

print()
print("=" * 60)
print("  All tools demonstrated successfully!")
print("=" * 60)
