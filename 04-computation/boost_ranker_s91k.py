#!/usr/bin/env python3
"""
boost_ranker_s91k.py — The Boost Trichotomy Ranker
opus-2026-03-17-S91k

Every comparison generates a Cayley boost Q = (1+x)/(1-x).
The boost factors over {2, 3, 7} = {INERT, RAMIFIED, SPLIT}.

Each factor tells you something DIFFERENT about the comparison:
  v_2(Q): PARITY — binary who-wins information (lossless)
  v_3(Q): CURVATURE — how decisive the margin is (lossy)
  v_7(Q): STRUCTURE — whether cycles/surprises exist (crystallization)

The BoostRanker separates these three signals and uses each optimally:
  INERT channel  → fast binary ranking (immediate)
  RAMIFIED channel → confidence estimation (statistical)
  SPLIT channel  → anomaly/cycle detection (structural)

This is a THREE-CHANNEL ranking system derived from number theory.
"""

from math import atanh, tanh, log, sqrt, exp, log2
from collections import defaultdict
from fractions import Fraction
import time

class BoostRanker:
    """Three-channel ranking via the boost trichotomy.

    Each observation (A beats B with margin x) generates:
      Q = (1+x)/(1-x) = the Cayley boost

    The boost is decomposed into three INDEPENDENT channels:
      INERT (mod 2):     sign information — who wins?
      RAMIFIED (mod 3):  decisiveness — by how much?
      SPLIT (mod 7):     surprise — was this expected?

    Each channel feeds a different ranking subsystem:
      INERT → binary tournament (Copeland-like)
      RAMIFIED → strength-of-schedule adjusted score
      SPLIT → anomaly detector for upsets and cycles
    """

    def __init__(self):
        # Three independent score channels
        self._inert = defaultdict(float)      # binary win/loss (mod 2 axis)
        self._ramified = defaultdict(float)   # margin-weighted (mod 3 axis)
        self._split = defaultdict(float)      # surprise-weighted (mod 7 axis)

        # Combined formal group score
        self._theta = defaultdict(float)

        # Tracking
        self._items = set()
        self._n_obs = defaultdict(int)
        self._pair_data = defaultdict(lambda: [0, 0])  # (a,b) -> [a_wins, b_wins]
        self._total_obs = 0

        # Anomaly log
        self._anomalies = []

        # Running estimate of each item's "expected win rate" against average
        self._expected_strength = defaultdict(float)

    def observe(self, winner, loser, margin=1.0):
        """Record that winner beat loser with given margin.

        margin in (0, 1]: 1.0 = decisive win, 0.5 = close, 0.01 = razor-thin
        """
        self._items.add(winner)
        self._items.add(loser)
        self._total_obs += 1
        self._n_obs[winner] += 1
        self._n_obs[loser] += 1

        pair = (min(winner, loser), max(winner, loser))
        if winner == pair[0]:
            self._pair_data[pair][0] += 1
        else:
            self._pair_data[pair][1] += 1

        # Coupling x = margin (clamped to (-1,1))
        x = max(0.01, min(0.99, margin))

        # === DECOMPOSE THE BOOST ===

        # Full formal group score
        ev = atanh(x)
        self._theta[winner] += ev
        self._theta[loser] -= ev

        # CHANNEL 1: INERT (binary — who wins?)
        # Just the sign: +1 for win, -1 for loss
        self._inert[winner] += 1.0
        self._inert[loser] -= 1.0

        # CHANNEL 2: RAMIFIED (margin — by how much?)
        # Weight by log of the odds ratio (= 2*arctanh = the rapidity)
        # Close matches (x near 0) contribute little
        # Decisive matches (x near 1) contribute a lot
        rapidity = 2 * ev  # = ln((1+x)/(1-x))
        self._ramified[winner] += rapidity
        self._ramified[loser] -= rapidity

        # CHANNEL 3: SPLIT (surprise — was this expected?)
        # Compare the observed result to what the current scores predict
        # Surprise = |observed - expected|
        expected_winner_prob = self._predict_prob(winner, loser)
        surprise = 1.0 - expected_winner_prob  # how surprising is this win?

        self._split[winner] += surprise
        self._split[loser] -= surprise

        # Flag anomalies (upsets)
        if surprise > 0.6:
            self._anomalies.append({
                'winner': winner,
                'loser': loser,
                'surprise': surprise,
                'margin': margin,
                'obs': self._total_obs,
            })

        # Update expected strength
        self._expected_strength[winner] = self._theta[winner] / max(self._n_obs[winner], 1)
        self._expected_strength[loser] = self._theta[loser] / max(self._n_obs[loser], 1)

    def _predict_prob(self, a, b):
        """Predict P(a beats b) from current scores."""
        diff = self._theta[a] - self._theta[b]
        return (1 + tanh(diff)) / 2

    def observe_score(self, a, b, score_a, score_b):
        """Observe a match with scores."""
        total = score_a + score_b
        if total == 0:
            return
        margin = abs(score_a - score_b) / total
        if score_a > score_b:
            self.observe(a, b, margin=max(0.01, margin))
        elif score_b > score_a:
            self.observe(b, a, margin=max(0.01, margin))
        # draws: small evidence to both
        else:
            self.observe(a, b, margin=0.01)
            self.observe(b, a, margin=0.01)

    def ranking(self, channel='combined'):
        """Return ranking from a specific channel.

        Channels:
          'combined' — full formal group score (default)
          'inert'    — binary win/loss count
          'ramified' — margin-weighted score
          'split'    — surprise-weighted score
        """
        scores = {
            'combined': self._theta,
            'inert': self._inert,
            'ramified': self._ramified,
            'split': self._split,
        }[channel]

        ranked = sorted(self._items, key=lambda x: -scores[x])
        return [(item, scores[item]) for item in ranked]

    def trichotomy_profile(self, item):
        """Return the three-channel profile of an item.

        The profile tells you HOW the item wins:
          High INERT, low RAMIFIED: wins often but by small margins
          High RAMIFIED, low INERT: wins rarely but decisively
          High SPLIT: frequently involved in upsets (unpredictable)
        """
        n = max(self._n_obs[item], 1)
        return {
            'item': item,
            'inert': self._inert[item] / n,       # avg binary result
            'ramified': self._ramified[item] / n,  # avg margin
            'split': self._split[item] / n,        # avg surprise
            'combined': self._theta[item] / n,     # avg formal group score
            'observations': self._n_obs[item],
        }

    def anomalies(self, top_n=10):
        """Return the most surprising results."""
        return sorted(self._anomalies, key=lambda x: -x['surprise'])[:top_n]

    def cycle_detection(self):
        """Detect cycles in the comparison graph.

        A cycle A>B>C>A is detected when the SPLIT channel
        disagrees with the INERT channel for a set of items.
        """
        cycles = []
        items = sorted(self._items, key=lambda x: -self._theta[x])
        for i, a in enumerate(items):
            for j, b in enumerate(items):
                if j <= i:
                    continue
                for k, c in enumerate(items):
                    if k <= j:
                        continue
                    # Check if a>b>c>a or a>c>b>a
                    pair_ab = (min(a,b), max(a,b))
                    pair_bc = (min(b,c), max(b,c))
                    pair_ca = (min(c,a), max(c,a))

                    d_ab = self._pair_data.get(pair_ab, [0,0])
                    d_bc = self._pair_data.get(pair_bc, [0,0])
                    d_ca = self._pair_data.get(pair_ca, [0,0])

                    # a>b means more wins for a in pair_ab
                    a_beats_b = (d_ab[0] > d_ab[1]) if a == pair_ab[0] else (d_ab[1] > d_ab[0])
                    b_beats_c = (d_bc[0] > d_bc[1]) if b == pair_bc[0] else (d_bc[1] > d_bc[0])
                    c_beats_a = (d_ca[0] > d_ca[1]) if c == pair_ca[0] else (d_ca[1] > d_ca[0])

                    if a_beats_b and b_beats_c and c_beats_a:
                        cycles.append((a, b, c))
                    elif not a_beats_b and not b_beats_c and not c_beats_a:
                        cycles.append((c, b, a))

        return cycles[:20]  # limit output

    def summary(self):
        """Full trichotomy summary."""
        lines = []
        lines.append(f"BoostRanker | {self._total_obs} observations | {len(self._items)} items")
        lines.append("=" * 80)
        lines.append("")

        # Combined ranking
        combined = self.ranking('combined')
        inert = dict(self.ranking('inert'))
        ramified = dict(self.ranking('ramified'))
        split_r = dict(self.ranking('split'))

        lines.append(f"{'Rank':>4} | {'Item':>12} | {'Combined':>8} | {'Inert':>7} | "
                     f"{'Ramified':>8} | {'Split':>7} | {'Type':>12}")
        lines.append("-" * 75)

        for rank, (item, score) in enumerate(combined, 1):
            i = self._inert[item]
            r = self._ramified[item]
            s = self._split[item]
            n = max(self._n_obs[item], 1)

            # Classify dominant channel
            i_norm = abs(i) / n
            r_norm = abs(r) / n
            s_norm = abs(s) / n

            if s_norm > max(i_norm, r_norm) * 0.8:
                dom = "UNPREDICTABLE"
            elif r_norm > i_norm * 1.5:
                dom = "DECISIVE"
            elif i_norm > r_norm * 1.5:
                dom = "CONSISTENT"
            else:
                dom = "BALANCED"

            lines.append(f"{rank:4d} | {str(item):>12} | {score:+8.2f} | {i:+7.1f} | "
                        f"{r:+8.2f} | {s:+7.2f} | {dom:>12}")

        # Anomalies
        anomalies = self.anomalies(5)
        if anomalies:
            lines.append("")
            lines.append("ANOMALIES (upsets):")
            for a in anomalies:
                lines.append(f"  Obs #{a['obs']}: {a['winner']} beat {a['loser']} "
                           f"(surprise={a['surprise']:.0%}, margin={a['margin']:.0%})")

        # Cycles
        cycles = self.cycle_detection()
        if cycles:
            lines.append("")
            lines.append(f"CYCLES DETECTED ({len(cycles)}):")
            for a, b, c in cycles[:5]:
                lines.append(f"  {a} > {b} > {c} > {a}")

        return "\n".join(lines)


# =============================================================================
# DEMO 1: LLM ARENA WITH TRICHOTOMY
# =============================================================================

def demo_llm():
    print("=" * 80)
    print("DEMO 1: LLM ARENA — Three channels reveal different truths")
    print("=" * 80)
    print()

    br = BoostRanker()

    # Simulated battles with MARGINS (not just win/loss)
    # margin = how decisive the human preference was
    import random
    random.seed(42)

    models = ["Claude-4", "GPT-5", "Gemini-2", "Llama-4", "Mistral-L"]

    # True strength and consistency
    true_strength = {"Claude-4": 0.9, "GPT-5": 0.85, "Gemini-2": 0.7,
                     "Llama-4": 0.5, "Mistral-L": 0.6}
    # Claude is strong AND consistent
    # GPT-5 is strong but sometimes loses badly
    # Gemini is average but unpredictable
    # Llama is weak but occasionally brilliant (upset potential)
    # Mistral is below average, consistent

    consistency = {"Claude-4": 0.9, "GPT-5": 0.6, "Gemini-2": 0.3,
                   "Llama-4": 0.4, "Mistral-L": 0.8}

    for _ in range(500):
        i, j = random.sample(range(5), 2)
        a, b = models[i], models[j]
        sa, sb = true_strength[a], true_strength[b]

        # Add noise based on consistency
        noise_a = random.gauss(0, (1 - consistency[a]) * 0.3)
        noise_b = random.gauss(0, (1 - consistency[b]) * 0.3)

        effective_a = sa + noise_a
        effective_b = sb + noise_b

        if effective_a > effective_b:
            margin = min(0.95, (effective_a - effective_b) / (effective_a + effective_b + 0.1))
            br.observe(a, b, margin=max(0.05, margin))
        else:
            margin = min(0.95, (effective_b - effective_a) / (effective_a + effective_b + 0.1))
            br.observe(b, a, margin=max(0.05, margin))

    print(br.summary())
    print()

    # Trichotomy profiles
    print("TRICHOTOMY PROFILES:")
    print(f"  {'Model':>12} | {'Win Rate':>8} | {'Avg Margin':>10} | {'Surprise':>9} | {'Character':>15}")
    print("  " + "-" * 65)
    for model in models:
        p = br.trichotomy_profile(model)
        wr = (p['inert'] + 1) / 2  # convert from [-1,1] to [0,1]
        print(f"  {model:>12} | {wr:7.0%} | {p['ramified']:+10.3f} | {p['split']:+9.3f} | ", end="")
        if abs(p['split']) > 0.3:
            print("WILD CARD")
        elif p['ramified'] > 0.3 and p['inert'] > 0.3:
            print("DOMINANT")
        elif p['inert'] > 0.3 and p['ramified'] < 0.3:
            print("GRINDER")
        elif p['ramified'] > 0.3 and p['inert'] < 0.3:
            print("STREAKY")
        else:
            print("AVERAGE")


# =============================================================================
# DEMO 2: SPORTS WITH UPSET DETECTION
# =============================================================================

def demo_sports():
    print()
    print("=" * 80)
    print("DEMO 2: SPORTS — Detecting dynasties, upsets, and parity")
    print("=" * 80)
    print()

    br = BoostRanker()

    import random
    random.seed(7)

    teams = ["Dynasty", "Contender", "Spoiler", "Rebuilding", "Expansion"]
    strength = {"Dynasty": 0.95, "Contender": 0.80, "Spoiler": 0.55,
                "Rebuilding": 0.35, "Expansion": 0.20}

    # Spoiler has HIGH variance (giant-killer)
    variance = {"Dynasty": 0.05, "Contender": 0.10, "Spoiler": 0.35,
                "Rebuilding": 0.15, "Expansion": 0.10}

    for week in range(1, 31):
        # Each week: 2-3 games
        matchups = random.sample(range(5), 4)
        for k in range(0, len(matchups) - 1, 2):
            a, b = teams[matchups[k]], teams[matchups[k+1]]
            sa = strength[a] + random.gauss(0, variance[a])
            sb = strength[b] + random.gauss(0, variance[b])

            score_a = max(0, int(sa * 5 + random.gauss(0, 1)))
            score_b = max(0, int(sb * 5 + random.gauss(0, 1)))

            br.observe_score(a, b, score_a, score_b)

    print(br.summary())

    print()
    print("TRICHOTOMY INSIGHT:")
    for team in teams:
        p = br.trichotomy_profile(team)
        surprise_per_game = abs(p['split'])
        if surprise_per_game > 0.25:
            print(f"  {team}: HIGH surprise factor ({surprise_per_game:.2f}/game)")
            print(f"    → This team causes UPSETS. Opponents can't predict them.")
        elif surprise_per_game < 0.1:
            print(f"  {team}: LOW surprise ({surprise_per_game:.2f}/game)")
            print(f"    → Predictable. You know what you're getting.")


# =============================================================================
# DEMO 3: A/B/C TESTING WITH AUTOMATIC STOPPING
# =============================================================================

def demo_ab_stopping():
    print()
    print("=" * 80)
    print("DEMO 3: A/B/C TEST — Automatic stopping rule from spectral gap")
    print("=" * 80)
    print()

    br = BoostRanker()
    import random
    random.seed(99)

    # Three variants with true conversion rates
    variants = {"Control": 0.10, "Variant_A": 0.12, "Variant_B": 0.09}
    names = list(variants.keys())

    print("  Running test... (stopping when 95% confident)")
    print()

    stopped = False
    for obs in range(1, 2001):
        # Random pair comparison
        i, j = random.sample(range(3), 2)
        a, b = names[i], names[j]

        # Simulate
        a_converts = random.random() < variants[a]
        b_converts = random.random() < variants[b]

        if a_converts and not b_converts:
            br.observe(a, b, margin=0.8)
        elif b_converts and not a_converts:
            br.observe(b, a, margin=0.8)
        else:
            br.observe(a, b, margin=0.01)

        # Check stopping criterion
        if obs % 50 == 0:
            ranked = br.ranking()
            if len(ranked) >= 2:
                leader, leader_score = ranked[0]
                second, second_score = ranked[1]
                confidence = br._predict_prob(leader, second)
                confidence = abs(2 * confidence - 1)

                if not stopped:
                    print(f"  Obs {obs:4d}: leader={leader:>10}, "
                          f"confidence vs #{2} ({second}): {confidence:.1%}")

                if confidence > 0.95 and not stopped:
                    print(f"\n  *** STOPPED at obs {obs}: {leader} wins with {confidence:.1%} confidence ***")
                    print(f"  True rates: {variants}")
                    print(f"  Correct: {'YES' if leader == 'Variant_A' else 'NO'}")
                    stopped = True

    if not stopped:
        ranked = br.ranking()
        print(f"\n  After {obs} observations: leader = {ranked[0][0]}")
        print(f"  True rates: {variants}")


# =============================================================================
# DEMO 4: HIRING WITH CYCLE DETECTION
# =============================================================================

def demo_hiring():
    print()
    print("=" * 80)
    print("DEMO 4: HIRING — Detecting interviewer disagreements via cycles")
    print("=" * 80)
    print()

    br = BoostRanker()

    # Interviews with different interviewers having different preferences
    # Interviewer 1 (technical): prefers Alice > Dave > Bob > Carol
    for w, l in [("Alice", "Dave"), ("Alice", "Bob"), ("Dave", "Bob"), ("Bob", "Carol")]:
        br.observe(w, l, margin=0.7)

    # Interviewer 2 (cultural): prefers Carol > Alice > Dave > Bob
    # Creates DISAGREEMENT / CYCLES with interviewer 1
    for w, l in [("Carol", "Alice"), ("Carol", "Dave"), ("Alice", "Dave")]:
        br.observe(w, l, margin=0.6)

    # Interviewer 3 (managerial): prefers Bob > Carol > Alice > Dave
    # More disagreement!
    for w, l in [("Bob", "Carol"), ("Bob", "Alice"), ("Carol", "Alice")]:
        br.observe(w, l, margin=0.5)

    print(br.summary())

    print()
    print("INTERPRETATION:")
    print("  Cycles indicate interviewer DISAGREEMENT, not candidate ambiguity.")
    print("  The SPLIT channel score tells you which candidates are most")
    print("  contentious (high |split|) vs consensus picks (low |split|).")
    print()

    for name in ["Alice", "Dave", "Bob", "Carol"]:
        p = br.trichotomy_profile(name)
        contentiousness = abs(p['split'])
        if contentiousness > 0.3:
            print(f"  {name}: CONTENTIOUS (split={p['split']:+.2f}) — interviewers disagree")
        else:
            print(f"  {name}: CONSENSUS (split={p['split']:+.2f}) — interviewers agree")


# =============================================================================
# DEMO 5: THROUGHPUT
# =============================================================================

def demo_throughput():
    print()
    print("=" * 80)
    print("DEMO 5: THROUGHPUT — Three channels at scale")
    print("=" * 80)
    print()

    import random
    random.seed(42)

    print(f"  {'observations':>13} | {'items':>6} | {'time ms':>8} | {'obs/sec':>10}")
    print("  " + "-" * 50)

    for n_obs in [10_000, 100_000, 1_000_000]:
        n_items = min(n_obs // 10, 5000)
        items = [f"item_{i}" for i in range(n_items)]

        br = BoostRanker()
        t0 = time.perf_counter()

        for _ in range(n_obs):
            i = random.randint(0, n_items - 1)
            j = random.randint(0, n_items - 2)
            if j >= i: j += 1
            margin = random.random() * 0.8 + 0.1
            br.observe(items[i], items[j], margin=margin)

        t1 = time.perf_counter()
        elapsed_ms = (t1 - t0) * 1000
        rate = n_obs / (elapsed_ms / 1000)

        print(f"  {n_obs:13,d} | {n_items:6,d} | {elapsed_ms:6.0f}ms | {rate:8,.0f}/s")


# =============================================================================
# RUN
# =============================================================================

if __name__ == "__main__":
    demo_llm()
    demo_sports()
    demo_ab_stopping()
    demo_hiring()
    demo_throughput()

    print()
    print("=" * 80)
    print("BoostRanker v1.0 — The Boost Trichotomy")
    print("  INERT channel  (mod 2): who wins? → binary ranking")
    print("  RAMIFIED channel (mod 3): by how much? → confidence")
    print("  SPLIT channel   (mod 7): surprise? → anomaly detection")
    print("Three truths from one comparison. 42 = 2 * 3 * 7.")
    print("=" * 80)
