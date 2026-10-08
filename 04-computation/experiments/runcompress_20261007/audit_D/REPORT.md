# Audit D — THM-4604 (Collatz words and the reducible locus) and results note §4

Independent audit, 2026-10-07, for session mac-mini-2026-10-07-twoanchor (continuation). Saved by the main session from the auditor's final message; the harness blocked subagent report writes.
* Script: `thm4604_audit.py` (+ `.out`), about 1.8·10^5 checks in 58 families.

## Verdict

**The algebra is CONFIRMED:**
* the conventions, the cocycle identity and positivity (B_w is also odd and prime to 3);
* cyclic coboundaries;
* the centraliser torus, by distinct eigenvalues;
* the trace formula 2cosh(δ/2);
* the Cayley cubic, exactly;
* the example pair;
* Markov at κ = −2, with integer points 3·Markov;
* the mutation as (u, v) → (u, uv), verified on the depth-10 Christoffel tree;
* faithfulness of w → G_w (injective on all 65,535 words with Σw ≤ 16).

## Corrections (applied)

* **C1. Anti-homomorphism.** G_uv = G_v G_u is an anti-homomorphism for concatenation. B_w is odd and prime to 3. B is non-split on the whole monoid (c_1 = −1 ≠ c_2 = 1).
* **C2. Coboundary and centraliser.** The condition 2^(Σw) ≠ 3^|w| is automatic. The centraliser is the torus fixing c_w and ∞ because the eigenvalues are distinct; the anchored states are its elements with multiplier 3^k. Equal-count group elements can be unipotent.
* **C3. Merges.** Distinct words with equal counts have equal traces but always different carries, since (|w|, Σw, B_w) determines w. A merge is an incidence F_u(n) = F_v(h). Example: 483 and 469 have words (1,7) and (7,1) and merge at 17.
* **C4. Where friezes live.**
  * Markov sits at κ = −2.
  * Conway–Coxeter friezes are quiddities with product −I; since tr[M(a), M(b)] = 2 + (a−b)², they are not at −2.
  * No G_w is ±I. Read as quiddities, however, valuation words can close friezes: (1,1,1) and (1,2,2,2,2,2,2,1,7) are the U-words of 7183 and 2583211, which merge at 24245.
* **C5. Title.** "lie on the reducible locus".
* **C6. Status.** KNOWN in substance: Böhm–Sontacchi 1978; Lagarias 1985 and 1990; Goldman 2009 Prop. 2.3.1; Culler–Shalen 1983; Cohn 1955; Goldman 2003. The minor identity, the carry exchange relation and faithfulness are PROVED.
* **C7. The Reading was false in places.**
  * Fricke traces are defined on the whole character variety.
  * Frieze entries are minors of a configuration and are NOT blind to merges: det(v_i, v_j) = 3^i 2^(A_i) B_(w[i+1..j]), and the Plücker relation becomes the carry exchange relation `B_xy B_yz = B_y B_xyz + 3^|y| 2^(Σy) B_x B_z`.
  * The correct dichotomy is that trace coordinates see only the semisimplification, while minor coordinates see the cocycle.
  * The reason frieze moves are not rewrites: ear moves are relations of a ↦ M(a), and a ↦ G_a has no relations.
* **C8 and C9.** Results note §0 and §4 corrected accordingly.

## Prior art

* Böhm–Sontacchi 1978; Lagarias 1985, 1990.
* Goldman 2003, 2009.
* Cohn 1955; Bowditch 1998.
* Conway–Coxeter 1973.
* Berstel–Lauve–Reutenauer–Saliola 2008.
* Adjacent: Fernández–Ibáñez, arXiv:2607.24844 (Christoffel words and Collatz parity).
