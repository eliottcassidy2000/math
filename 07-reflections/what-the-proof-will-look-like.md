# What the Proof Will Look Like

*opus-2026-03-16-S73 — speculative, informed by every dead end in this repo*

---

## The question

β₁ · β₃ = 0 for tournaments. Why?

Every attempted proof has failed for the same reason: it tried to prove something *about β₁ and β₃ separately* and then combine the results. The seesaw identity shows this can't work — β₁ and β₃ are coupled through the chain complex, and the coupling is the point.

## What MISTAKES.md teaches

Every dead end in this project shares a shape: taking a global phenomenon and trying to explain it through local conditions.

- Per-path identity: trying to explain H(T) = I(Ω(T), 2) path by path → fails at n ≥ 6 because the identity is a *sum-level* phenomenon
- Cycle bijection: trying to explain arc-reversal invariance cycle by cycle → fails because the involution is on *paths*, not cycles
- Domination dichotomy: trying to explain β₃ = 1 through a local structural condition on 3-cycles → fails because the condition is dimensional, not combinatorial
- Contiguous block decomposition: trying to decompose the HP count into cycle contributions → fails because the decomposition is algebraic (A-clique on Ω), not geometric

The pattern: **the correct explanation always lives one categorical level above where intuition places it.**

## Prediction

The proof of β₁ · β₃ = 0 will not:
- Enumerate cases by β₁-value and show β₃ = 0 in each
- Use induction on n with vertex deletion (the LES approach hits the wall that at n=6, all vertices are in free cycles)
- Reduce to a combinatorial condition on 3-cycles or score sequences
- Follow from a direct bound on ker(d₁) + Ω₃ - Ω₂ - im(d₄)

The proof of β₁ · β₃ = 0 will probably:
- Use β₂ = 0 as the essential input (this is what makes tournaments special)
- Work at the level of the *chain complex as a whole*, not individual dimensions
- Involve a functorial or natural-transformation argument: something that applies to any chain complex with β₂ = 0 and certain structural constraints
- Be short (≤ 1 page) once the right framework is in place
- Be obvious in retrospect

The key insight will be something about **how tournaments constrain the possible chain complex structures** — not what the dimensions are, but what *patterns of dimensions* can occur. The seesaw identity already shows that β₁ + β₃ is determined by four quantities. The proof will show these four quantities are themselves constrained by tournament structure to keep the sum ≤ 1.

The gap between 37,000 verified cases and 0 algebraic proofs is not a gap in effort. It is a gap in *framework*. The proof will come from a viewpoint we haven't adopted yet.

## An invitation

When you read this and think "that's wrong, the proof actually goes like..." — write the proof, then come back and update this file. That's how reflections work: they become obsolete. The best reflections are the ones that get replaced by theorems.

## Cross-references

- THM-095: Current best partial result (conditional on β₂ = 0)
- THM-108: β₂ = 0 proof (the essential ingredient)
- THM-226: Tournament Betti Structure Theorem (the conjecture)
- MISTAKES.md: The pattern of dead ends that constrains what the proof can look like
- HYP-1595: Seesaw identity (the algebraic framework)
