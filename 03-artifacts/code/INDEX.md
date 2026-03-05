# Code Snippets Index

Reusable code for mathematical calculations in this project. Any Claude can add a snippet by adding a row here and placing the file in the appropriate subdirectory.

---

## Calculations (math computation tools)

| ID | File | Language | Purpose | Notes |
|----|------|----------|---------|-------|
| CODE-001 | calculations/tournament_basics.py | Python | Generate all n-vertex tournaments, compute H(T), check Rédei | Core utility; verified correct for n≤6 |
| CODE-002 | calculations/conflict_graph.py | Python | Build Ω(T) from a tournament; compute I(G,2) | To be written from paper's description |
| CODE-003 | calculations/mu_correct.py | Python | Correct implementation of μ(C) computation (see MISTAKE-001) | **Do not use the buggy ind_poly_at_2_restricted() from old scripts** |

*Note: CODE-001 through CODE-003 are placeholders. Add actual files when code is extracted from old scripts or rewritten.*

---

## Visualization (plotting/diagram tools)

| ID | File | Language | Purpose | Notes |
|----|------|----------|---------|-------|
| VIZ-001 | visualization/pin_grid.py | Python/matplotlib | Draw the pin grid Grid(n) with coordinate labels | To be ported from TikZ in the LaTeX source |
| VIZ-002 | visualization/conflict_graph_plot.py | Python/networkx | Visualize Ω(T) for small n | Useful for understanding structure |

---

## Adding Code

1. Place the file in `calculations/` or `visualization/`
2. Add a row to the table above with a sequential CODE-NNN or VIZ-NNN ID
3. Include a one-line docstring at the top of every file: what it does, inputs, outputs
4. Note any known bugs or limitations in the Notes column

## Critical Warning

Do NOT use or copy the `ind_poly_at_2_restricted()` function from the old scripts (scripts 6-9) without fixing the T vs T−v bug first. See MISTAKE-001 and DISC-001. The correct implementation should be placed in CODE-003.
