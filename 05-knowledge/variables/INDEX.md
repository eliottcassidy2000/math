# Variable Registry — Index

Alphabetical index of all mathematical quantities in this project.
Each variable has its own file with full cross-references.

**How to use:** Find a variable below, click to its detail file for equations, values, and links.
**How to add:** Create `variables/VARNAME.md` from the template at bottom, add a row here.

---

## Quick Lookup

| Variable | File | One-line description | Key equations |
|----------|------|---------------------|---------------|
| `A[i][j]` | [tournament-matrix.md](tournament-matrix.md) | Adjacency matrix of tournament T | A[i][j] + A[j][i] = 1 |
| `alpha_k` / `a_k` | [alpha-k.md](alpha-k.md) | Count of disjoint k-tuples of odd cycles | H = 1 + 2*sum(alpha_k) |
| `B[i][j]` | [signed-matrix.md](signed-matrix.md) | Signed adjacency: B = 2A - J | B[i][j] in {-1, +1} |
| `bc33` | [bc33.md](bc33.md) | Disjoint 3-cycle pairs | bc33 = alpha_2 at n=7 |
| `bc35_w` | [bc35-w.md](bc35-w.md) | Weighted (3-cycle, 5-cycle) disjoint pairs | Appears in c_2 at n=9 |
| `c_k` | [fourier-coefficients.md](fourier-coefficients.md) | Coefficients in W-polynomial expansion | W(T,r) = sum c_k * r^k |
| `D_k` | [D-k.md](D-k.md) | Degree-2k Fourier component of W | W = sum D_k; D_0 = c_0 |
| `F_k` | [forward-counts.md](forward-counts.md) | Count of HPs with exactly k forward edges | F(T,x) = sum F_k * x^k |
| `fwd(P)` | [forward-edges.md](forward-edges.md) | Number of forward edges in HP P | fwd(P) + bwd(P) = n-1 |
| `H(T)` | [hamiltonian-paths.md](hamiltonian-paths.md) | Total Hamiltonian path count | H = tr(M); H = I(Omega,2); H odd (Redei) |
| `I(G,x)` | [independence-poly.md](independence-poly.md) | Independence polynomial of graph G | I(Omega(T), 2) = H(T) (OCF) |
| `M[a,b]` | [transfer-matrix.md](transfer-matrix.md) | Transfer matrix entry | M symmetric; tr(M) = H; M[a,b] = sum_S (-1)^|S| E_a B_b |
| `mu(C)` | [mu-C.md](mu-C.md) | Weight of odd cycle C in OCF | mu(C) = I(Omega(T-v) restricted, 2) |
| `Omega(T)` | [omega-graph.md](omega-graph.md) | Conflict graph of odd cycles | Vertices = odd directed cycles; edges = shared vertex |
| `s_e` | [edge-variable.md](edge-variable.md) | Centered edge variable | s_e = A_e - 1/2; s_e in {-1/2, +1/2} |
| `S(T)` | [signed-hp-permanent.md](signed-hp-permanent.md) | Signed HP permanent | S(T) = sum_P prod B[P_i][P_{i+1}] |
| `t_k` | [cycle-counts.md](cycle-counts.md) | Number of directed k-cycles in T | t_3 = C(n,3)/4 iff regular |
| `W(T,r)` | [W-polynomial.md](W-polynomial.md) | W-polynomial / weighted HP count | W(T,r) = sum_P prod(r + s_e) |

---

## Template for new variable files

```markdown
# Variable: NAME

**Symbol:** `name`
**Type:** integer / polynomial / matrix / graph
**Defined in:** `01-canon/definitions.md` line N; THM-XXX

## Definition
[Precise mathematical definition]

## Values at small n
| n | Range | Notes | Source script |
|---|-------|-------|--------------|
| 3 | ... | ... | `04-computation/...` |

## Equations involving this variable
- EQ-1: [equation] — source: THM-XXX
- EQ-2: [equation] — source: HYP-YYY

## Relationships to other variables
- Related to [VAR](link): [how]

## Tags
#tag1 #tag2 #tag3
```
