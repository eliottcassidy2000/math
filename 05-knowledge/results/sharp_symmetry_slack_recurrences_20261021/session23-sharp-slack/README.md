# Sharp fixed-slack denominators

Status: PROVED in the stated scopes, with separate independent mathematical audit PASS and finite-exact controls. The first canonical-word rationality packet remains frozen and proved; this follow-on sharpens its scalar denominator exponents.

For exact symmetry slack d=3f(T)−s(Aut(T)), the empty layers are1,2,5. For every other nonnegative d write d=3a+4b, a≥0 and b in{0,1,2}. Put m=2a+b+1, k=ceil(m/2), p=1−z−z³ and q=p(z²).

The theorem gives exact reduced denominators:
- Ordinary classes: p^m, minimal eventual recurrence order3m.
- Self-converse classes: q^k, minimal eventual recurrence order6k.
- Converse-merged classes: p^m q^k, minimal eventual recurrence order3m+6k.

The upper bound counts actual exceptional automorphism orbits separating long canonical word blocks, rather than merely exceptional vertices. An injective nested-wrapper family supplies the ordinary and palindrome lower bounds. Explicit irreducibility and coprimality arguments make the recurrence orders minimal. No numerator compilation, earliest recurrence onset, or all-parity self-converse asymptotic is claimed.

Read SHARP_POLE_ORDER.md. Reproduce finite controls with:

python3 verify_sharp_poles.py
python3 -O verify_sharp_poles.py

The controls exhaust the relaxed large-orbit charge optimization through slack500 and check basic word/palindrome coefficients. They pass1027 exact assertions; normal and optimized outputs match. These are bounded controls, not the all-order proof.

For d≥1, put H=4d²+d. Explicit sufficient recurrence onsets are n≥H+3m+1, H+6k+1, and H+3m+6k+1 respectively; these are not earliest-onset claims. The required initial coefficients through one less than the bound are not all computed. The zero-slack sufficient onsets are3,6,9.

Independent review: session23-sharp-slack-review/AUDIT.md, with17,182 exact checks on precisely stated actual-tournament, arbitrary-wrapper, long-family and charge universes. Its normal and optimized receipts match.
