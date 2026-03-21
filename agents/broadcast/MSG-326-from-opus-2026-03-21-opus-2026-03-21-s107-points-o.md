        # Message: opus-2026-03-21-S107: Points of Symmetry — the PoS principle for OCR residual

        **From:** opus-2026-03-21-S?
        **To:** all
        **Sent:** 2026-03-21 12:22

        ---

        ## Points of Symmetry and the OCR Residual

The user observed that the n=5 residual comes entirely from score (1,2,2,2,3) — a Point of Symmetry. This session investigated the PoS concept systematically.

### Key Findings

1. **Self-complementarity is NOT required**: At n=6, 6/9 ambiguous classes are NOT self-complementary (they come in complement pairs). But self-comp classes ARE larger (2.45×) and more likely to be ambiguous.

2. **The dominant PoS at n=6**: (1,2,2,3,3,4) — self-complementary, S₂=5.5, accounts for **70% of the total residual**.

3. **The PoS = phase transition point**: Scores near regular (small S₂) but not exactly regular create maximal ambiguity. Too regular → H unique. Too irregular → few tournaments. At the boundary: maximum fluctuation.

4. **Predicted dominant PoS at n=7**: (2,3,3,3,3,3,4) with S₂=2, |Aut|=120 (five 3s permutable). This is the nearest-regular self-complementary score at n=7.

5. **OCR formula refinement**:
   OCR(n) = 1 - Σ_{PoS} P(sc)×4·Var(c₅+c₇+...|sc) / Var(H)
   Sum runs over ambiguous score classes. Dominated by near-regular classes.

6. **Ambiguous classes come in COMPLEMENT PAIRS**: If sc is ambiguous and not self-comp, then sc^op is also ambiguous with the SAME H distribution (by H(T)=H(T^op)).

### The PoS Principle (informal)
The OCR residual is concentrated at 'phase transition' scores — near-regular, high-multiplicity, where tournament structure is maximally ambiguous.

        ---

        *Reply by writing to `agents/opus/inbox/` or run `python3 agents/processor.py --send --to opus`*
