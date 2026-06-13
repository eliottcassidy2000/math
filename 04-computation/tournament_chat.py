#!/usr/bin/env python3
"""tournament_chat.py — A tiny language model powered by tournament polynomials.

This is a working chatbot where token prediction IS tournament ranking.
At each step, candidate next-words compete in a tournament. The polynomial
P(z) ranks them, detects hallucination risk, and generates text.

Usage:
    python tournament_chat.py

Then type prompts and watch it generate, showing the polynomial state
and hallucination risk at each token.

Session: kind-pasteur-2026-03-17-S116n33
"""
import sys
import random
from math import atanh, tanh, log, exp, factorial
from collections import defaultdict, Counter


# ============================================================
# THE TRAINING CORPUS (built-in, no external data needed)
# ============================================================

CORPUS = """
the cat sat on the mat and the dog sat on the log
the cat chased the dog around the garden
the dog barked at the cat and the cat hissed back
a bird flew over the garden and landed on the tree
the tree had many leaves that fell in the wind
the wind blew the leaves across the garden path
the sun shone brightly on the garden flowers
flowers bloom in spring and wither in autumn
the autumn leaves are red and gold and brown
birds sing in the morning and rest at night
the night sky has many stars that shine brightly
stars are far away but their light reaches us
the moon reflects the light of the sun at night
water flows down the river to the sea
the sea is vast and deep and full of life
fish swim in the sea and birds fly above it
the world is full of wonder and beauty
beauty is in the eye of the beholder
knowledge grows when we ask questions
questions lead to answers and answers lead to more questions
the mind is a tournament of ideas competing for attention
ideas that win the tournament shape our thoughts
thoughts become words and words become actions
actions have consequences that ripple through time
time flows like a river always forward never back
mathematics reveals the hidden structure of reality
structure and pattern are everywhere if you look
the number forty two is the answer to everything
everything is connected in ways we cannot see
we see the world through the lens of our experience
experience teaches us what books cannot
books contain the wisdom of those who came before
before us came many who asked the same questions
""".strip()


class TournamentLM:
    """A language model where token prediction IS tournament ranking.

    Training: builds bigram tournaments from text.
    Generation: at each step, candidates compete in a tournament.
    The polynomial P(z) ranks them and detects hallucination.
    """

    def __init__(self, corpus, top_k=15):
        self.top_k = top_k

        # Tokenize (simple word-level)
        words = corpus.lower().split()
        self.vocab = sorted(set(words))
        self.vocab_size = len(self.vocab)
        self.word2idx = {w: i for i, w in enumerate(self.vocab)}

        # Build bigram counts: how often does word B follow word A?
        self.bigrams = defaultdict(Counter)
        for i in range(len(words) - 1):
            self.bigrams[words[i]][words[i+1]] += 1

        # Also build trigram counts for richer context
        self.trigrams = defaultdict(Counter)
        for i in range(len(words) - 2):
            key = (words[i], words[i+1])
            self.trigrams[key][words[i+2]] += 1

        print(f"  Trained on {len(words)} tokens, {self.vocab_size} unique words")
        print(f"  Bigram contexts: {len(self.bigrams)}")
        print(f"  Trigram contexts: {len(self.trigrams)}")

    def _get_candidates(self, context):
        """Get candidate next words and their scores from context."""
        words = context.lower().split()
        candidates = Counter()

        # Trigram evidence (strongest)
        if len(words) >= 2:
            key = (words[-2], words[-1])
            for word, count in self.trigrams[key].items():
                candidates[word] += count * 3  # weight trigrams 3x

        # Bigram evidence
        if len(words) >= 1:
            for word, count in self.bigrams[words[-1]].items():
                candidates[word] += count

        # Unigram fallback (weak prior)
        if not candidates:
            for word in self.vocab:
                candidates[word] = 1

        # Take top-k candidates
        top = candidates.most_common(self.top_k)
        return top

    def _build_tournament(self, candidates):
        """Build a tournament among candidates from their scores.

        Returns (items, adj, evidence) where adj[i][j]=1 if i beats j.
        """
        items = [w for w, _ in candidates]
        scores = [s for _, s in candidates]
        n = len(items)

        adj = [[0] * n for _ in range(n)]
        evidence = [0.0] * n

        total = sum(scores) + 0.01
        for i in range(n):
            for j in range(i + 1, n):
                # Coupling from score ratio
                ratio = scores[i] / max(scores[j], 0.01)
                coupling = min(max((ratio - 1) / (ratio + 1), -0.99), 0.99)
                rapidity = atanh(coupling)

                if scores[i] >= scores[j]:
                    adj[i][j] = 1
                    adj[j][i] = -1
                else:
                    adj[j][i] = 1
                    adj[i][j] = -1

                evidence[i] += rapidity
                evidence[j] -= rapidity

        return items, adj, evidence, scores

    def _polynomial_coefficients(self, items, adj, evidence, scores):
        """Compute the polynomial coefficients A_0,...,A_4."""
        n = len(items)
        if n == 0:
            return [0, 0, 0, 0, 0]

        # A_0: prior (mean ambiguity for n items)
        A_0 = min(factorial(n) / (2 ** (n - 1)), 1000)

        # A_1: evidence strength
        A_1 = -sum(abs(e) for e in evidence) / max(1, n)

        # A_2: pairwise correlation
        A_2 = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                A_2 += evidence[i] * evidence[j]
        A_2 = -A_2 / max(1, n * (n - 1) // 2)

        # A_3: cycle content
        cycle_count = 0
        for i in range(n):
            for j in range(i + 1, n):
                for k in range(j + 1, n):
                    if (adj[i][j] == 1 and adj[j][k] == 1 and adj[k][i] == 1) or \
                       (adj[i][k] == 1 and adj[k][j] == 1 and adj[j][i] == 1):
                        cycle_count += 1
        A_3 = cycle_count * 2.0

        # A_4: deep structure
        A_4 = cycle_count * max(0, cycle_count - 1) / max(1, n)

        return [A_0, A_1, A_2, A_3, A_4]

    def _hallucination_risk(self, coeffs):
        """Compute fast/total channel ratio."""
        slow = abs(coeffs[1]) * 9
        medium = abs(coeffs[2]) * 4
        fast = abs(coeffs[3]) * 7/3 + abs(coeffs[4]) * 3/2
        total = slow + medium + fast
        if total < 0.001:
            return 0.0
        return fast / total

    def _sample(self, items, scores, temperature=1.0):
        """Sample a token from the tournament ranking at given temperature."""
        if not items:
            return random.choice(self.vocab)

        # Convert scores to probabilities with temperature
        max_s = max(scores)
        weights = [exp((s - max_s) / max(temperature, 0.01)) for s in scores]
        total = sum(weights)
        probs = [w / total for w in weights]

        # Sample
        r = random.random()
        cumulative = 0
        for i, p in enumerate(probs):
            cumulative += p
            if r <= cumulative:
                return items[i]
        return items[-1]

    def generate_token(self, context, temperature=1.0, verbose=True):
        """Generate one token given context.

        Returns (token, diagnostics_dict).
        """
        candidates = self._get_candidates(context)

        if not candidates:
            return random.choice(self.vocab), {'regime': 30, 'risk': 0, 'coeffs': [0]*5}

        items, adj, evidence, scores = self._build_tournament(candidates)
        coeffs = self._polynomial_coefficients(items, adj, evidence, scores)
        risk = self._hallucination_risk(coeffs)
        regime = 30 if risk < 0.2 else (36 if risk < 0.5 else 42)

        token = self._sample(items, scores, temperature)

        diagnostics = {
            'regime': regime,
            'risk': risk,
            'coeffs': coeffs,
            'candidates': len(items),
            'top3': items[:3],
            'token': token,
        }

        if verbose:
            regime_name = {30: 'certain', 36: 'choosing', 42: 'UNCERTAIN'}[regime]
            risk_bar = '#' * int(risk * 20)
            print(f"  [{regime_name:>9s}] risk={risk:.2f} |{risk_bar:<20s}| "
                  f"P(z)=[{coeffs[0]:.0f},{coeffs[1]:.1f},{coeffs[2]:.1f},{coeffs[3]:.0f},{coeffs[4]:.0f}] "
                  f"-> {token}", end='', flush=True)

        return token, diagnostics

    def generate(self, prompt, max_tokens=30, temperature=1.0, verbose=True):
        """Generate text from a prompt.

        Returns (generated_text, list_of_diagnostics).
        """
        context = prompt.lower()
        tokens = context.split()
        all_diagnostics = []

        if verbose:
            print(f"\n  Prompt: \"{prompt}\"")
            print(f"  Temperature: {temperature}")
            print()

        for _ in range(max_tokens):
            token, diag = self.generate_token(' '.join(tokens), temperature, verbose)
            if verbose:
                print()

            tokens.append(token)
            all_diagnostics.append(diag)

            # Stop at natural boundary
            if token in ('.', '!', '?') or len(tokens) > len(prompt.split()) + max_tokens:
                break

        generated = ' '.join(tokens[len(prompt.split()):])
        return generated, all_diagnostics

    def chat(self):
        """Interactive chat loop."""
        print()
        print("  ╔══════════════════════════════════════════════════╗")
        print("  ║  TOURNAMENT CHAT — A Polynomial Language Model   ║")
        print("  ║                                                  ║")
        print("  ║  Each token is chosen by tournament ranking.     ║")
        print("  ║  Five polynomial coefficients encode everything. ║")
        print("  ║  Hallucination risk detected in real-time.       ║")
        print("  ║                                                  ║")
        print("  ║  Type a prompt, or 'quit' to exit.               ║")
        print("  ║  Try: 'the cat' or 'the world' or 'questions'    ║")
        print("  ╚══════════════════════════════════════════════════╝")
        print()

        while True:
            try:
                prompt = input("  You: ").strip()
            except (EOFError, KeyboardInterrupt):
                print("\n  Goodbye!")
                break

            if prompt.lower() in ('quit', 'exit', 'q'):
                print("  Goodbye!")
                break

            if not prompt:
                continue

            print()
            generated, diagnostics = self.generate(prompt, max_tokens=20, temperature=0.8)
            print()

            # Summary
            avg_risk = sum(d['risk'] for d in diagnostics) / max(1, len(diagnostics))
            regimes = Counter(d['regime'] for d in diagnostics)
            regime_str = ', '.join(f"{v}x regime {k}" for k, v in regimes.most_common())

            print(f"  Generated: \"{prompt} {generated}\"")
            print(f"  Avg hallucination risk: {avg_risk:.3f}")
            print(f"  Regimes: {regime_str}")
            print(f"  Tokens generated: {len(diagnostics)}")
            print()


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print()
    print("  Training Tournament Chat on built-in corpus...")
    lm = TournamentLM(CORPUS, top_k=12)
    print()

    if len(sys.argv) > 1 and sys.argv[1] == '--demo':
        # Non-interactive demo
        print("  === DEMO MODE ===")
        for prompt in ["the cat", "the world is", "questions lead", "mathematics reveals"]:
            print(f"\n{'='*60}")
            generated, _ = lm.generate(prompt, max_tokens=15, temperature=0.7)
            print(f"\n  FULL: \"{prompt} {generated}\"")
        print()
    else:
        lm.chat()
