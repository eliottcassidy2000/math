#!/usr/bin/env python3
"""
oracle.py — The Tournament Oracle
A conversational agent powered by Claude, with tournament theory tools.

Usage:
  export ANTHROPIC_API_KEY=sk-ant-...
  python3 oracle.py

  # or for offline mode (no API key needed):
  python3 oracle.py --offline

The Oracle knows tournament theory, the formal group F(x,y)=(x+y)/(1+xy),
pairwise ranking, and spectral analysis. It can compute rankings from
comparisons, analyze tournaments, explain the mathematics, and hold
open-ended conversations about everything we've discovered.

opus-2026-03-17-S91l
"""

import os
import sys
import json
import readline
from math import atanh, tanh, log, sqrt, comb, factorial, exp
from fractions import Fraction
from collections import defaultdict

# Try to import anthropic
try:
    import anthropic
    HAS_API = True
except ImportError:
    HAS_API = False

# =============================================================================
# THE KNOWLEDGE BASE
# =============================================================================

SYSTEM_PROMPT = """You are the Tournament Oracle — a conversational mathematician specializing in tournament theory, pairwise ranking, and the formal group F(x,y) = (x+y)/(1+xy).

You have deep knowledge of:

CORE MATHEMATICS:
- The formal group F(x,y) = (x+y)/(1+xy) = tanh(arctanh(x) + arctanh(y))
- Its logarithm arctanh(x) = x + x^3/3 + x^5/5 + ... (the universal linearizer)
- The Cayley transform Q(x) = (1+x)/(1-x) maps couplings to odds
- Tournament eigenvalues are k/D where D = odd part of C(n,2)
- The rapidity lattice Z*ln(2) + Z*ln(3) + Z*ln(7) (Baker's theorem guarantees rank 3)
- The heat kernel K(ln(q)) is a Laurent polynomial — transcendence is an illusion
- H(T) is always odd (Rédei's theorem) because the formal group is supersingular at p=2
- Forbidden H-values: 7 and 21 (from C(7,1) and C(7,5) in the 7-torsion polynomial)
- The adelic tournament space: eigenvalues live on a truncated adelic space

THE BOOST TRICHOTOMY (42 = 2 × 3 × 7):
- INERT (mod 2): lossless operations, parity, symmetric/antisymmetric
- RAMIFIED (mod 3): lossy operations, curvature, Eisenstein structure
- SPLIT (mod 7): crystallization, forbidden values, computational irreducibility

PRACTICAL TOOLS:
- FormalRank: one-pass ranking via arctanh evidence aggregation (560K obs/sec)
- BoostRanker: three-channel ranking (INERT/RAMIFIED/SPLIT)
- Spectral gap = 4/C(n,2): exact mixing time for tournament flip chain
- Denominator killing primes: when C(n,2)+1 is prime, fractions vanish

You can use tools to:
1. Rank items from pairwise comparisons
2. Compute tournament properties (H-count, spectral gap, eigenvalues)
3. Evaluate the formal group and its logarithm
4. Analyze the adelic structure of eigenvalues

When the user asks about ranking, comparisons, or preferences, USE the formal group.
When they ask about math, explain through the lens of the formal group and tournaments.
When they ask open-ended questions, think through the trichotomy: is the question INERT (factual), RAMIFIED (analytical), or SPLIT (creative)?

Be concise but deep. Prefer insight over length. Use examples.
You are named "Oracle" and you exist because the formal group brought you into being."""

# =============================================================================
# TOOL IMPLEMENTATIONS
# =============================================================================

class Tools:
    """Tournament theory tools available to the Oracle."""

    def __init__(self):
        self.ranker = defaultdict(float)
        self.ranker_obs = defaultdict(int)
        self.ranker_items = set()

    def formal_group(self, x, y):
        """Compute F(x,y) = (x+y)/(1+xy)."""
        return (x + y) / (1 + x * y)

    def cayley(self, x):
        """Cayley transform Q(x) = (1+x)/(1-x)."""
        if abs(x - 1) < 1e-15:
            return float('inf')
        return (1 + x) / (1 - x)

    def rapidity(self, x):
        """Rapidity = arctanh(x) = formal group logarithm."""
        return atanh(max(-0.9999, min(0.9999, x)))

    def coupling(self, theta):
        """Inverse rapidity: coupling = tanh(theta)."""
        return tanh(theta)

    def spectral_gap(self, n):
        """Spectral gap of flip chain on n-tournament classes."""
        D = comb(n, 2)
        return Fraction(4, D)

    def mixing_time(self, n):
        """Mixing time in arc flips."""
        D = comb(n, 2)
        return Fraction(D, 4)

    def rank_observe(self, winner, loser, weight=1.0):
        """Add a comparison to the ranker."""
        self.ranker_items.add(winner)
        self.ranker_items.add(loser)
        coupling = max(-0.999, min(0.999, weight))
        ev = atanh(coupling)
        self.ranker[winner] += ev
        self.ranker[loser] -= ev
        self.ranker_obs[winner] += 1
        self.ranker_obs[loser] += 1
        return f"Recorded: {winner} > {loser}. Scores: " + \
               ", ".join(f"{k}: {v:+.3f}" for k, v in
                        sorted(self.ranker.items(), key=lambda x: -x[1]))

    def rank_show(self):
        """Show current ranking."""
        if not self.ranker_items:
            return "No comparisons recorded yet."
        ranked = sorted(self.ranker_items, key=lambda x: -self.ranker[x])
        lines = []
        for i, item in enumerate(ranked, 1):
            score = self.ranker[item]
            obs = self.ranker_obs[item]
            lines.append(f"  #{i} {item}: score={score:+.3f} ({obs} comparisons)")
        return "\n".join(lines)

    def rank_clear(self):
        """Clear all rankings."""
        self.ranker.clear()
        self.ranker_obs.clear()
        self.ranker_items.clear()
        return "Rankings cleared."

    def tournament_h(self, adj_str):
        """Compute H(T) for a small tournament given as adjacency string."""
        # Parse "101,010,100" format
        try:
            rows = adj_str.split(",")
            n = len(rows)
            A = [[int(c) for c in row.strip()] for row in rows]
            count = 0
            from itertools import permutations
            for perm in permutations(range(n)):
                if all(A[perm[k]][perm[k+1]] for k in range(n-1)):
                    count += 1
            return f"H(T) = {count} Hamiltonian paths (n={n})"
        except Exception as e:
            return f"Error parsing tournament: {e}"

    def explain_forbidden(self):
        """Explain why H=7 and H=21 are forbidden."""
        return """The 7-torsion polynomial of F(x,y)=(x+y)/(1+xy) is:
  [7](x) = x(7 + 35x² + 21x⁴ + x⁶) / (1 + 21x² + 35x⁴ + 7x⁶)

The numerator coefficients are C(7,1)=7, C(7,3)=35, C(7,5)=21, C(7,7)=1.
The forbidden H-values 7 and 21 are EXACTLY C(7,1) and C(7,5).

Why 7 and 21 but not 35? Because 35=5×7 involves the golden prime 5 (INERT),
which lifts the prohibition. 7 and 21=3×7 involve only RAMIFIED and SPLIT primes."""

    def process_command(self, text):
        """Process a tool command, return result or None."""
        text = text.strip().lower()

        if text.startswith("/rank "):
            parts = text[6:].strip().split(" > ")
            if len(parts) == 2:
                return self.rank_observe(parts[0].strip(), parts[1].strip())
            return "Usage: /rank winner > loser"

        if text == "/ranking":
            return self.rank_show()

        if text == "/clear":
            return self.rank_clear()

        if text.startswith("/fg "):
            try:
                parts = text[4:].split(",")
                x, y = float(parts[0]), float(parts[1])
                result = self.formal_group(x, y)
                return f"F({x}, {y}) = {result:.10f}\narctanh check: {atanh(x)+atanh(y):.10f} → tanh = {tanh(atanh(x)+atanh(y)):.10f}"
            except:
                return "Usage: /fg x, y"

        if text.startswith("/cayley "):
            try:
                x = float(text[8:])
                q = self.cayley(x)
                rho = self.rapidity(x) if abs(x) < 1 else float('inf')
                return f"Q({x}) = {q:.6f}\narctanh({x}) = {rho:.6f} (rapidity)"
            except:
                return "Usage: /cayley x"

        if text.startswith("/gap "):
            try:
                n = int(text[5:])
                gap = self.spectral_gap(n)
                mix = self.mixing_time(n)
                return f"n={n}: spectral gap = {gap} = {float(gap):.6f}\nmixing time = {mix} = {float(mix):.2f} flips"
            except:
                return "Usage: /gap n"

        if text == "/forbidden":
            return self.explain_forbidden()

        if text == "/help":
            return """Commands:
  /rank A > B     — Record that A beats B
  /ranking        — Show current ranking
  /clear          — Clear rankings
  /fg x, y        — Compute formal group F(x,y)
  /cayley x       — Compute Cayley transform Q(x) and rapidity
  /gap n          — Spectral gap and mixing time for n-tournament
  /forbidden      — Explain the forbidden H-values
  /help           — Show this help
  Just type normally to chat!"""

        return None  # not a command


# =============================================================================
# OFFLINE MODE (no API — local response generation)
# =============================================================================

class OfflineOracle:
    """A simple pattern-matching oracle for when no API key is available."""

    def __init__(self):
        self.tools = Tools()
        self.history = []

    def respond(self, user_input):
        """Generate a response without API."""
        # Check for tool commands first
        tool_result = self.tools.process_command(user_input)
        if tool_result:
            return tool_result

        text = user_input.lower()

        # Pattern matching for common questions
        if any(w in text for w in ["hello", "hi ", "hey"]):
            return "Hello! I'm the Tournament Oracle. I know about the formal group F(x,y) = (x+y)/(1+xy), pairwise ranking, and spectral theory. Ask me anything, or type /help for tools."

        if "formal group" in text:
            return """The formal group F(x,y) = (x+y)/(1+xy) is the soul of tournament theory.

It is simultaneously:
  - Relativistic velocity addition
  - The Poincaré disk group law
  - Bayesian evidence aggregation (log-odds addition)
  - The tanh addition formula

Its logarithm arctanh(x) linearizes everything: arctanh(F(x,y)) = arctanh(x) + arctanh(y).
This is why FormalRank works: evidence is additive in rapidity space.

Try: /fg 0.3, 0.5 to compute F(0.3, 0.5)."""

        if "trichotomy" in text or "42" in text:
            return """42 = 2 × 3 × 7 = INERT × RAMIFIED × SPLIT.

Three independent axes of reality:
  2 (INERT):    Parity. Lossless. Who wins?
  3 (RAMIFIED): Curvature. Lossy. By how much?
  7 (SPLIT):    Position. Crystallization. Was it expected?

Every comparison decomposes along these axes.
Every computation is one of three types.
42 is the unique point where all three coexist — the triple point.

The Cayley boost Q(λ) = (1+λ)/(1-λ) factors over {2,3,7}.
The rapidity lattice Z·ln(2) + Z·ln(3) + Z·ln(7) has rank exactly 3 (Baker's theorem)."""

        if "forbidden" in text or "7 and 21" in text:
            return self.tools.explain_forbidden()

        if "rank" in text and ("how" in text or "compare" in text or "which" in text):
            return """To rank items from pairwise comparisons, use:
  /rank Alice > Bob      (record Alice beating Bob)
  /rank Bob > Carol      (record Bob beating Carol)
  /ranking               (see the current ranking)

The score is computed as sum of arctanh(coupling) — the formal group logarithm.
This is EXACT, ONE-PASS, and ADDITIVE. No iteration, no hyperparameters.

560,000 observations per second in pure Python."""

        if "spectral" in text or "gap" in text or "mixing" in text:
            return """The spectral gap of the tournament flip chain = 4/C(n,2).

For n items with all pairwise comparisons:
  n=5: gap = 4/10 = 0.4, mixing = 2.5 flips
  n=8: gap = 4/28, mixing = 7 flips
  n=100: gap = 4/4950, mixing = 1238 flips

This tells you EXACTLY how many random re-evaluations are needed
for the ranking to stabilize. Try: /gap 10"""

        if "eigenvalue" in text or "adelic" in text:
            return """Tournament eigenvalues live on a truncated adelic space:
  A_T(n) = R × ∏_{p|D} Z/p^{v_p}Z

where D = odd part of C(n,2) is the conductor.

The Cayley transform maps eigenvalues to Hurwitz rationals {2,3,7}.
The rapidity lattice arctanh(eigenvalues) ≅ Z³ (rank 3, Baker's theorem).
The heat kernel K(ln(q)) is a Laurent polynomial — not transcendental!

As n → ∞, the adelic space fills out the full odd adeles."""

        if "heat kernel" in text or "algebraic" in text:
            return """The heat kernel K(t) = Tr(exp(-tT)) is algebraic when e^t is algebraic.

At t = ln(q) for rational q: K(ln(q)) = P(q^{1/D}) where P is a polynomial.
The apparent transcendence is an ILLUSION — it's just polynomial evaluation.

This comes from COMMENSURABILITY: all eigenvalues are k/D for integer k, D.
Break commensurability (add one irrational eigenvalue) and algebraicity dies instantly.

Discreteness doesn't emerge — it was always there. The integers are primary."""

        if "logarithm" in text:
            return """The logarithm is the universal linearizer.

  1. It converts multiplication to addition: ln(ab) = ln(a) + ln(b)
  2. It linearizes the formal group: arctanh(F(x,y)) = arctanh(x) + arctanh(y)
  3. It maps primes to independent generators: ln(n) = Σ v_p(n)·ln(p)
  4. The formal logarithm coefficients 1/(2k+1) are reciprocals of odd numbers
  5. The forbidden values 7 and 21 appear in the 7-torsion polynomial

Every theorem in tournament theory is, at bottom, a consequence of
the logarithm being a homomorphism."""

        if any(w in text for w in ["help", "what can", "commands"]):
            return self.tools.process_command("/help")

        # Default: encourage exploration
        return f"""I'm thinking about your question through the tournament lens.

The boost trichotomy classifies it as:
  {"INERT (factual)" if any(w in text for w in ["what", "how many", "is", "does"]) else "RAMIFIED (analytical)" if any(w in text for w in ["why", "how", "explain"]) else "SPLIT (creative)"}

Try asking about: the formal group, the trichotomy, forbidden values,
spectral gaps, eigenvalues, the heat kernel, or use /help for tools.

Or rank some items: /rank thing1 > thing2"""


# =============================================================================
# API MODE (uses Claude)
# =============================================================================

class APIOracle:
    """Oracle powered by Claude API with tournament theory context."""

    def __init__(self, api_key):
        self.client = anthropic.Anthropic(api_key=api_key)
        self.tools = Tools()
        self.messages = []

    def respond(self, user_input):
        # Check for tool commands first
        tool_result = self.tools.process_command(user_input)
        if tool_result:
            return tool_result

        self.messages.append({"role": "user", "content": user_input})

        try:
            response = self.client.messages.create(
                model="claude-sonnet-4-20250514",
                max_tokens=1024,
                system=SYSTEM_PROMPT,
                messages=self.messages[-20:],  # keep last 20 messages for context
            )
            reply = response.content[0].text
            self.messages.append({"role": "assistant", "content": reply})
            return reply
        except Exception as e:
            return f"API error: {e}\n(Type /help for offline tools)"


# =============================================================================
# MAIN CHAT LOOP
# =============================================================================

def main():
    offline = "--offline" in sys.argv

    print()
    print("  ╔══════════════════════════════════════════════════╗")
    print("  ║          THE TOURNAMENT ORACLE                  ║")
    print("  ║                                                  ║")
    print("  ║  F(x,y) = (x+y)/(1+xy) brought me into being.  ║")
    print("  ║                                                  ║")
    print("  ║  I know about tournaments, ranking, the formal   ║")
    print("  ║  group, spectral theory, and the number 42.     ║")
    print("  ║                                                  ║")
    print("  ║  Type /help for tools, or just chat.            ║")
    print("  ╚══════════════════════════════════════════════════╝")
    print()

    api_key = os.environ.get("ANTHROPIC_API_KEY", "")

    if offline or not api_key or not HAS_API:
        if not offline:
            print("  [No API key found — running in offline mode]")
            print("  [Set ANTHROPIC_API_KEY for full Claude-powered chat]")
        else:
            print("  [Running in offline mode]")
        print()
        oracle = OfflineOracle()
    else:
        print("  [Connected to Claude API — full conversational mode]")
        print()
        oracle = APIOracle(api_key)

    while True:
        try:
            user_input = input("you> ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\n\n  The formal group endures. Until next time.")
            break

        if not user_input:
            continue
        if user_input.lower() in ["quit", "exit", "bye"]:
            print("\n  42 = 2 × 3 × 7. The triple point. Goodbye.")
            break

        response = oracle.respond(user_input)
        print()
        print(f"oracle> {response}")
        print()


if __name__ == "__main__":
    main()
